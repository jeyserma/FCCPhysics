"""
Create the GMAD input files for BDSIM from MAD-X Twiss output files in tfs format.

LCC_v2 -- see README.md for the list of changes with respect to LCC_v1.
"""

import numpy as np
import pandas as pd
import sys
import pybdsim
import pymadx
import ROOT
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot, make_axes_locatable
from pathlib import Path
import re

TWISS_FILE = "GMAD/fcc_ee_z_b1_twiss_end.tfs"

# Element at which the beamline is cut for the two-stage runs. Stage 1 ends
# just before it, stage 2 starts at it, so there is no gap and no overlap.
SPLIT_ELEMENT = "QF1BL"

# Stage-2 beam file, written by make_stage2_input.py. BDSIM 1.7.7 accepts
# pdgid and weight as columns (verified against libbdsim), which is what lets
# a mixed-species dump be replayed. Weights are all 1.0 in the uniform
# replication scheme -- the 1/N is absorbed into nPrimariesTotal instead, so
# every downstream consumer keeps working with a single global weight.
# One file per species, and stage 2 is run once per species.
#
# This is NOT cosmetic. BDSIM 1.7.7 converts the energy column using the
# REFERENCE beam particle's mass, ignoring the per-line pdgid. With
# beam particle e+ that silently drops every photon below 511 keV (with an
# E column) or corrupts its energy (with an Ek column) -- verified by injecting
# a decade scan and reading the sampler back. Splitting by species and setting
# the reference particle to match makes the energies round-trip exactly.
STAGE2_INPUT = {"e+": "stage2_input_ep.dat", "gamma": "stage2_input_gamma.dat"}

# Ek, not E: the column is kinetic energy, which is what makes the photon
# energies survive once the reference particle is right.
STAGE2_FORMAT = "pdgid:x[m]:xp[rad]:y[m]:yp[rad]:Ek[GeV]:t[ns]:weight"



def read_tfs_header(tfs_file=TWISS_FILE, keys=("EX", "EY", "SIGE", "ENERGY")):
    """
    Read the @-parameters of a MAD-X tfs file, e.g.

        @ EX               %le    0.000000000700000

    Returns a dict {key: float}. Used so that the emittances and the energy
    spread of the halo distribution always follow the TWISS file that is
    actually being converted, instead of being hardcoded (LCC_v1 hardcoded
    emitx=0.7e-9 and emity=2.6e-12, which happened to match this lattice).
    """
    out = {k: None for k in keys}
    with open(tfs_file) as f:
        for line in f:
            if line.startswith("*"):   # start of the column definitions
                break
            if not line.startswith("@"):
                continue
            fields = line.split()
            if len(fields) >= 4 and fields[1] in out:
                out[fields[1]] = float(fields[3])
    missing = [k for k, v in out.items() if v is None]
    if missing:
        raise ValueError(f"Could not read {missing} from the TWISS header of {tfs_file}")
    return out


def truncated_exponential(rng, size, scale, lo, hi):
    """
    Draw `size` values from p(J) ~ exp(-J/scale) restricted to [lo, hi], by
    inverting the CDF. Exact, one uniform random number per particle, no
    rejection and no resampling -- every draw is an independent point.
    """
    u = rng.uniform(size=size)
    a, b = np.exp(-lo/scale), np.exp(-hi/scale)
    return -scale*np.log(a - u*(a - b))


def generate_4d_distribution(self, twiss_file, idx_start, tfs_file=TWISS_FILE):
    """
    Generate the beam-halo distribution at the entrance of the tracked beamline
    (the TWISS row before idx_start, i.e. the optics at the start of the first
    tracked element).

    Model -- unchanged with respect to LCC_v1: an exponential tail in the
    Courant-Snyder invariant (the "action" J) of each plane, cut to the region
    between an inner and an outer amplitude,

        p(Jx) ~ exp(-Jx/scale_x)  on [ (3.5 sigma_x)^2 , (xtail sigma_x)^2 ]
        p(Jy) ~ exp(-Jy/scale_y)  on [ (4.0 sigma_y)^2 , (ytail sigma_y)^2 ]

    with scale_x = (3.5^2 emitx)/xweight and scale_y = (4^2 emity)/yweight, so
    that xweight/yweight set the tail slope relative to the inner cut. Smaller
    weight = longer tail. For reference, a Gaussian beam has scale = 2*emit.

    Method -- changed with respect to LCC_v1. v1 threw 250*ngenerate uniform
    points in a 4-D box, kept the ~0.4% inside the annulus, and drew ngenerate
    of them with replacement with probability proportional to exp(...). Because
    the box extends to ytail=513 sigma_y while the vertical weight scale is only
    ~16 sigma_y, the weights were so concentrated that 200000 macroparticles
    contained only ~1069 distinct phase-space points (Kish effective sample size
    ~71) and needed a 3.2 GB peak RSS. Sampling the truncated exponential
    directly by inverting its CDF gives the identical distribution with
    ngenerate distinct particles, one random number each. See README.md.
    """
    rng = np.random.default_rng(self._seed)
    size = self._ngenerate

    row = twiss_file.iloc[idx_start-1]
    alfx, alfy = row['ALFX'], row['ALFY']
    betx, bety = row['BETX'], row['BETY']
    dispx, dispxp = row['DX'], row['DPX']

    header = read_tfs_header(tfs_file)
    emitx, emity = header['EX'], header['EY']
    sigmaE, energy = header['SIGE'], header['ENERGY']

    haloNSigmaXInner, haloNSigmaXOuter = 3.5, self._xtail
    haloNSigmaYInner, haloNSigmaYOuter = 4.0, self._ytail

    # amplitude (action) limits of the sampled annulus, in metre-radian
    JxLo, JxHi = haloNSigmaXInner**2 * emitx, haloNSigmaXOuter**2 * emitx
    JyLo, JyHi = haloNSigmaYInner**2 * emity, haloNSigmaYOuter**2 * emity

    # tail slopes, defined relative to the inner cut (as in LCC_v1)
    scaleX = JxLo / self._xweight
    scaleY = JyLo / self._yweight

    Jx = truncated_exponential(rng, size, scaleX, JxLo, JxHi)
    Jy = truncated_exponential(rng, size, scaleY, JyLo, JyHi)

    # uniform betatron phase, then the Courant-Snyder transform back to (x, x'),
    # which satisfies gamma*x^2 + 2*alfa*x*x' + beta*x'^2 = J by construction
    phix = rng.uniform(0.0, 2.0*np.pi, size=size)
    phiy = rng.uniform(0.0, 2.0*np.pi, size=size)

    sampled_dx  = np.sqrt(Jx*betx)*np.cos(phix)
    sampled_dxp = -np.sqrt(Jx/betx)*(alfx*np.cos(phix) + np.sin(phix))
    sampled_dy  = np.sqrt(Jy*bety)*np.cos(phiy)
    sampled_dyp = -np.sqrt(Jy/bety)*(alfy*np.cos(phiy) + np.sin(phiy))

    # ------------------------------------------------------------------
    # ENERGY SPREAD AND DISPERSION: deliberately left at the LCC_v1 values,
    # so that LCC_v2 differs from LCC_v1 in the sampling method ALONE and the
    # effect of that change can be measured on its own. Both are known to be
    # wrong -- see README.md sections 2 and 3 for the full analysis.
    #
    #   * SIGE in the TWISS header is dE/E, so the absolute spread should be
    #     SIGE * ENERGY = 45.6 MeV. The line below uses 1.0 MeV, i.e.
    #     sigma_delta = 2.2e-5 instead of 1.0e-3, a factor 45.6 too small.
    #
    #   * BDSIM applies "userfile" coordinates literally and never adds the
    #     dispersive orbit (dispx/dispxp in the beam block are used only by
    #     gausstwiss), so x should get + dispx*delta and xp + dispxp*delta.
    #     It does not here.
    #
    # The two belong together: with sigma_delta = 2.2e-5 the dispersive offset
    # is 7.7 um instead of 349 um, i.e. invisible next to the 813 um betatron
    # sigma_x. Fixing either alone changes nothing measurable.
    #
    # To restore the physics, replace the single line below with:
    #
    #     delta = rng.normal(0.0, sigmaE, size=size)
    #     sampled_dx  = sampled_dx  + dispx*delta
    #     sampled_dxp = sampled_dxp + dispxp*delta
    #     sampled_E   = energy*(1.0 + delta)
    #
    # sigmaE, energy, dispx and dispxp are read above purely so that this is a
    # copy-paste change.
    # ------------------------------------------------------------------
    sampled_E = rng.normal(45.6, 1.0e-3, size=size)   # LCC_v1 behaviour, see above

    return sampled_dx, sampled_dxp, sampled_dy, sampled_dyp, sampled_E

def add_before_last_semicolon(file_path, new_element):
    # Read the existing file contents
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Print the file contents before modification
    #print("File Contents Before Modification:")
    #for line in lines:
        #print(line.strip())
    # Make sure the file is not empty
    if len(lines) > 0:
        # Process the last line
        last_line = lines[-1].strip()
        if ';' in last_line:
            # Split the last line at the last semicolon
            parts = last_line.rsplit(';', 1)
            modified_line = parts[0] + ',\n'+ new_element + ';' + parts[1]
            # Update the last line in the list
            lines[-1] = modified_line + '\n'
    # Reopen the file in write mode and write the modified content
    with open(file_path, 'w') as file:
        file.writelines(lines)
    # Print the file contents after modification
    #print("\nFile Contents After Modification:")
    #for line in lines:
        #print(line.strip())

class MDIStudy:
    """
    Run a study, generating the GMAD files with pybdsim and run BDSIM through pybdsim
    REBDSIM can also be ran with pydbsim

    Example:
    >>> import long_z_lattice_core
    >>> s = MDIStudy()
    >>> s.genGMAD()
    >>> s.runStudy()

    Optional parameters can be set afterwards like changing the mask and collimators apertures.
    """

    def __init__(self, **kwargs):
        self._ngenerate = kwargs.get("ngenerate", 5000)          # number of primary particles i.e. positrons
        self._nruns     = kwargs.get("nruns", 1)                 # number of iterations

        self._runKey    = kwargs.get("runKey", '')               # key to name the output file

        self._seed      = kwargs.get("seed", 12)                 # a particular seed, for reproducability
        self._userfile  = kwargs.get("userfile", 1)          # 1:Gaussian beam, 2: Uniform halo, 3: Exponential halo, 4: Injection particles at FFQ

        self._roundmask = kwargs.get("roundmask", 3)            # type of mask closer to the IP, 1: round mask (GDML), 2: elliptical mask (GDML), 3: jaw mask (BDSIM), 4: elliptical mask (BDSIM)
        self._maskA     = kwargs.get("maskA", [0.015, 0.015, 0.007, 0.013])                # masks aperture
        self._extraA	= kwargs.get("extraA", ["11.05e-3", "12.08e-3", "20.28e-3", "7.86e-3", "22.07e-3"])	 # Aperture of BWL Hor, QC0L vert and hor(1 and 2) and QC2L vert an hor collimators

        self._X0        = kwargs.get("X0", 0)	                 # horizontal beam centroid displacwithBement
        self._XP0       = kwargs.get("XP0", 0)               # horizontal beam centroid angle
        self._Y0        = kwargs.get("Y0", 0)                # vertical beam centroid displacement
        self._YP0       = kwargs.get("YP0", 0)               # vertical beam centroid angle
        
        self._deltaS    = kwargs.get("deltaS", 0)             # Longitudinal position of the injected beam centroid wrt IP

        self._xtail     = kwargs.get("xtail", 10)             # horizontal halo width
        self._ytail     = kwargs.get("ytail", 10)             # vertical halo width
        self._xweight   = kwargs.get("xweight", 0)           # horizontal tails slope, i.e. lifetime
        self._yweight   = kwargs.get("yweight", 0)           # vertical tails slope, i.e. lifetime

        self._withSol   = kwargs.get("withSol", True)            # include the (anti-)solenoid field map
        self._withCorr  = kwargs.get("withCorr", False)            # include the (anti-)solenoid field map with orbit correctors
        self._withBC1L  = kwargs.get("withBC1L", 0)             # include the three dipoles before the IP (where the solenoid SR hits the beam pipe)
        self._withDip   = kwargs.get("withDip", True)             # include the three dipoles after the IP (where the solenoid SR hits the beam pipe)

        self._repo      = kwargs.get("repo", "DATA/")              # where the simulation outputs will be saved

        # Two-stage ("dual sampling") support. 0 = single stage, identical to
        # LCC_v2. 1 = track from the start of the tracked region up to the split
        # element and sample everything crossing that plane. 2 = start at the
        # split element, read that dump back as the beam, and track to QD0AL.
        self._stage     = kwargs.get("stage", 0)
        # Which species stage 2 injects: "e+" or "gamma". See STAGE2_INPUT.
        self._stage2_species = kwargs.get("stage2_species", "e+")
        # Strip the tunnel from the stage-2 model (see genGMAD). Stage 2 only.
        self._slimGeometry = kwargs.get("slimGeometry", False)

        self._optics    = kwargs.get('optics', False)
        self._traj      = kwargs.get('traj', False)
        self._bpabs     = kwargs.get('bpabs', False)

    def runOptics(self):
        """
        Generate a rebdsim file containing the optical functions.
        """
        print('Loading optics in BDSIM...\n')
        pybdsim.Run.Bdsim('GMAD/input.gmad', outfile='output', ngenerate=5000)
        print('Generating rebdsim optics file...\n')
        pybdsim.Run.RebdsimOptics('output.root', "optics_{}_{}_{}_{}.root".format(self._X0,self._Y0,self._XP0,self._YP0))

    def plotOptics(self, sOffset):
        """
        Generate plots of the orbit and optical functions.
        """
        d = pybdsim.Data.Load("optics_{}_{}_{}_{}.root".format(self._X0,self._Y0,self._XP0,self._YP0))

        fig = plt.figure(figsize=(5, 4), dpi=200)
        ax = plt.subplot()
        plt.plot(d.optics.S()-sOffset, 1e6*d.optics.Mean_x(), lw=1, label='X')
        plt.plot(d.optics.S()-sOffset, 1e6*d.optics.Mean_y(), lw=1, label='Y')
        plt.legend(); 
        #ax.set_ylim(-max(np.abs(ax.get_ylim())), max(np.abs(ax.get_ylim())))
        plt.grid(ls='--'); ax.set_xlabel('Distance from the IP [m]'); ax.set_ylabel(r'X/Y Orbit [$\rm{\mu}$m]')
        pybdsim.Plot.AddMachineLatticeFromSurveyToFigure(fig, d.model, sOffset=-sOffset)
        plt.xlim(-300, 300)
        plt.savefig('plotOrbit.png', dpi=300, bbox_inches='tight')

        fig = plt.figure(figsize=(5, 4), dpi=200)
        ax = host_subplot(111, figure=fig); ax2 = ax.twinx()
        ax.plot(d.optics.S()-sOffset, d.optics.Beta_x(), lw=1, label=r'$\rm{\beta}_x$')
        ax.plot(d.optics.S()-sOffset, d.optics.Beta_y(), lw=1, label=r'$\rm{\beta}_y$')
        ax2.plot(d.optics.S()-sOffset, 100*d.optics.Disp_x(), lw=1, label=r'D$_x$')
        ax.set_xlabel('Distance from the IP [m]'); ax.set_ylabel(r'Beta [m]')
        ax2.set_ylabel('Dispersion [cm]')
        plt.legend(); ax.grid(ls='--'); 
        #ax2.set_ylim(-25, 85); ax.set_ylim(-5, 105)
        pybdsim.Plot.AddMachineLatticeFromSurveyToFigure(fig, d.model, sOffset=-sOffset)
        plt.xlim(-300, 300)
        plt.savefig('plotOptics.png', dpi=300, bbox_inches='tight')
        # plt.show()

    def genGMAD(self):
        """
        Generate a set of GMAD files to the particular specification of this study as defined
        by the passed parameters when intiating this instance.
        """

        #######################
        # Read TWISS file
        #######################
        print("Read TWISS file\n")
        MADX_TWISS_HEADERS_SKIP_ROWS = 50
        MADX_TWISS_DATA_SKIP_ROWS = 52
        
        headers = pd.read_csv(TWISS_FILE, skiprows=MADX_TWISS_HEADERS_SKIP_ROWS,
                        nrows=0, sep=r"\s+")        
        headers.drop(headers.columns[[0, 1]], inplace=True, axis=1)
        twiss_file = pd.read_csv(TWISS_FILE,
                         header=None,
                         names=headers.columns.values,
                         na_filter=False,
                         skiprows=MADX_TWISS_DATA_SKIP_ROWS,
                         sep=r"\s+")
        twiss_file.index.name = 'NAME'

        #######################
        # Search the indexes
        #######################

        # Find an IP NOT on the edge of the sequence
        idx_IP = twiss_file.reset_index()[twiss_file.reset_index()['NAME']=="IP"].index[0]
        # Find the first dipole BEFORE the ip to start the sequence conversion
        found, idx_start = 0, idx_IP
        while found<1+self._withBC1L:
            if twiss_file.iloc[idx_start]['KEYWORD'] == "RBEND":
                found += 1
            idx_start -=1

        # Find the first dipole AFTER the IP to start the sequence conversion
        # Can be adapted to get any dipole after the IP to have a longer beam line
        # WARNING: idx_stop computed here is overwritten further down with
        # idx(QD0AL)+1, so the beamline always ends 2.4 m before the IP and
        # withDip currently has NO effect. Everything downstream of QD0AL
        # (MASK_QC1L, DRIFT_L0/L1, DRIFT_SOL, DRIFT_R1) is therefore defined in
        # input_components.gmad but never placed in the sequence. See README.md.
        found, idx_stop = 0, idx_IP
        while found<self._withDip:
            if twiss_file.iloc[idx_stop]['KEYWORD'] == "RBEND":
                found += 1
            idx_stop +=1

        drift_pre_mask, drift_pre_mask2, drift_pre_ip, drift_post_ip = twiss_file.iloc[idx_IP-3].name, twiss_file.iloc[idx_IP-11].name, twiss_file.iloc[idx_IP-1].name, twiss_file.iloc[idx_IP+1].name

        IR = twiss_file.iloc[idx_start:idx_stop]
        col_name = IR[IR["KEYWORD"]=="COLLIMATOR"].index

        print("Write collimator settings file\n")
        print(self._extraA)
        f = open("collimatorSettings.dat", "w")
        lines = [
            "# Collimator Settings\n",
            "name\tmaterial\txsize[m]\tysize[m]\n",
            f"{col_name[0]}\tinermet180\t{self._extraA[0]}\t30e-3\n",
            f"{col_name[1]}\tinermet180\t{self._extraA[1]}\t30e-3\n",
            f"{col_name[2]}\tinermet180\t{self._extraA[2]}\t30e-3\n",
            f"{col_name[3]}\tinermet180\t30e-3\t{self._extraA[3]}\n",
            f"{col_name[4]}\tinermet180\t{self._extraA[4]}\t30e-3\n",
            f"{col_name[5]}\tinermet180\t15e-3\t15e-3\n",
            f"{col_name[6]}\tinermet180\t15e-3\t15e-3\n",
        ]
        f.writelines(lines)
        f.close()
        cols = pybdsim.Data.Load("collimatorSettings.dat")

        #######################
        # Aperture information
        #######################

        ap = pymadx.Data.Aperture("GMAD/fcc_ee_z_aperture.tfs")
        aa = pd.DataFrame(ap.data.values(), columns=ap.columns)[['NAME', 'N1',  'APERTYPE', 'APER_1', 'APER_2']].set_index('NAME')
        ap = ap.RemoveBelowValue(5e-3)

        found, idx_start_FF = 0, idx_IP
        while found<5: # FOR TTBAR SHOULD CONSIDER 7 TO INCLUDE THE OTHER 2 QF
            if twiss_file.iloc[idx_start_FF]['KEYWORD'] == "QUADRUPOLE":
                found += 1
            idx_start_FF -=1
        found, idx_stop_FF = 0, idx_IP
        while found<5: # FOR TTBAR SHOULD CONSIDER 7 TO INCLUDE THE OTHER 2 QF
            if twiss_file.iloc[idx_stop_FF]['KEYWORD'] == "QUADRUPOLE":
                found += 1
            idx_stop_FF +=1

        magnet_geometry = {}

        last_drift_idx = None

        for i in range(idx_start_FF, idx_stop_FF):
            row = twiss_file.iloc[i]
            kw = row.KEYWORD
            name = str(row.name)

            if kw == "DRIFT":
                last_drift_idx = i
                continue

            if kw == "QUADRUPOLE":
                # quad aperture
                aper = aa.loc[row.name]["APER_1"].drop_duplicates().values[0]

                # 1) add quad itself
                magnet_geometry[name] = {"apertureType": "circular", "aper1": aper}

                # 2) add the drift immediately BEFORE this quad (if it exists)
                if last_drift_idx is not None:
                    drift_row = twiss_file.iloc[last_drift_idx]
                    drift_name = str(drift_row.name)
                    magnet_geometry[drift_name] = {"apertureType": "circular", "aper1": aper}

                    # optional: reset so only the *immediately* preceding drift is used once
                    last_drift_idx = None     

        print("Start converting the beamline into GMAD files\n")


        idx_QC1L1 = twiss_file.reset_index()[twiss_file.reset_index()['NAME']=="QD0AL"].index[0]
        idx_stop = idx_QC1L1 +1

        # --- two-stage split -------------------------------------------------
        idx_split = twiss_file.reset_index()[
            twiss_file.reset_index()['NAME'] == SPLIT_ELEMENT].index[0]
        if self._stage == 1:
            # stop just before the split element; the last element of the stage-1
            # line is then the drift whose exit is the split plane, and that is
            # where the sampler goes (see self._stage1_sampler below).
            idx_stop = idx_split
            self._stage1_sampler = twiss_file.iloc[idx_split - 1].name
            print(f"STAGE 1: {twiss_file.iloc[idx_start].name} -> {SPLIT_ELEMENT} "
                  f"(exclusive), sampling all species at '{self._stage1_sampler}'\n")
        elif self._stage == 2:
            idx_start = idx_split
            print(f"STAGE 2: {SPLIT_ELEMENT} -> QD0AL, beam read from the stage-1 dump\n")

        if self._userfile == 4:
            print("Using injection particles at FFQ, thus cutting the beamline before the FFQ\n")
            idx_QC2L2 = twiss_file.reset_index()[twiss_file.reset_index()['NAME']=="QF1BL"].index[0]
            idx_start = idx_QC2L2
        if self._userfile == 5:
            print("Using injection particles at B0BL, thus cutting the beamline before the B0BL\n")
            idx_B0BL = twiss_file.reset_index()[twiss_file.reset_index()['NAME']=="B0BL"].index[0]
            idx_start = idx_B0BL

        a, b = pybdsim.Convert.MadxTfs2Gmad(TWISS_FILE,
                                            "GMAD/input",
                                            linear = True,
                                            #aperturedict = ap,
                                            collimatordict = cols,
                                            startname = idx_start,
                                            stopname = idx_stop,
                                            samplers = None,
                                            defaultAperture = 30e-3,
                                            userdict=magnet_geometry,
                                            )

        print("Modifying GMAD sequence file\n")

        with open("GMAD/input_sequence.gmad", "r") as file:
            text= file.readlines()
        i = 0
        while i<len(text):
            if drift_pre_ip in text[i]:
                text[i] =text[i].replace(f"{drift_pre_ip}", "DRIFT_SOL")
            if drift_post_ip in text[i]:
                text[i] =text[i].replace(f"{drift_post_ip}", "DRIFT_R1")
            if col_name[5].replace(".", "") in text[i]:
                text[i] =text[i].replace(col_name[5].replace(".", ""), "MASK_QC2L")
            if drift_pre_mask2 in text[i]:
                text[i] =text[i].replace(f"{drift_pre_mask2}", "DRIFT_L2")
            if drift_pre_mask in text[i]:
                text[i] =text[i].replace(f"{drift_pre_mask}", "DRIFT_L1")
            if col_name[6].replace(".", "") in text[i]:
                if self._roundmask==2:
                    text[i] =text[i].replace(col_name[6].replace(".", ""), "MASK_QC1L")
                else:
                    text[i] =text[i].replace(col_name[6].replace(".", ""), "MASK_QC1L, DRIFT_L0")
            i+=1
        with open("GMAD/input_sequence.gmad", "w") as file:
            file.writelines(text)

        # Create the extraphysics file
        f = open("GMAD/emextraphysics.mac", "w")
        lines = [
        "/cuts/setLowEdge 10 eV\n",
        "/physics_lists/em/SyncRadiation true\n",
        "/physics_lists/em/GammaToMuons true\n",
        "/physics_lists/em/PositronToMuons true\n",
        "/physics_lists/em/PositronToHadrons true\n",
        "/physics_lists/em/MuonNuclear true\n",
        "/physics_lists/em/GammaNuclear true\n",
        ]
        f.writelines(lines)
        f.close()

        print(f'Compute optics :{self._optics}\nCompute trajectory {self._traj}\nFully absorbing beam pipe: {self._bpabs}')

        if self._optics:
            lines = [
                'precisionRegion: cutsregion, prodCutPhotons=1e-6, prodCutElectrons=1e-3, prodCutPositrons=1e-3;\n',
                '! physics options - full physics\n',
                'option, physicsList="",\n',
                'stopSecondaries=1,\n',
                'beampipeMaterial="Cu",\n',
                'beampipeRadius=30e-3,\n',
                'beampipeThickness=3e-3;\n',

                'sample, all;'
                ]
        elif self._traj and not self._bpabs:
            lines = [
                'precisionRegion: cutsregion, prodCutPhotons=1e-6, prodCutElectrons=1e-3, prodCutPositrons=1e-3;\n',
                '! physics options - full physics\n',
                'option, physicsList="G4FTFP_BERT",\n',
                'geant4PhysicsMacroFileName="emextraphysics.mac",\n',
                'minimumKineticEnergy = 7e-3,\n',
                'particlesToExcludeFromCuts ="22",\n',
                'magnetGeometryType="none",\n',
                'includeFringeFields=1,\n',
                'beampipeMaterial="Cu",\n',
                'apertureType="pointsfile:FCC_30mm.dat:mm",\n',
                #'beampipeRadius=35e-3,\n'
                'beampipeThickness=3e-3,\n',
                # Stage 2 covers only the last 6.8 m and samples on-axis at QD0AL.
                # The tunnel and its 2 m of soil are the single most expensive
                # thing to navigate, and every injected particle pays that cost
                # on entry. Dropping them for stage 2 is only safe if nothing
                # scattering off the tunnel returns to the sampler -- validated
                # in README section 9; DO NOT copy this to stage 1 or a
                # single-stage run, where the tunnel matters.
                ('buildTunnel=0,\n' if (self._stage == 2 and self._slimGeometry)
                 else 'buildTunnel=1,\n'),
                'tunnelOffsetX=-0.3,\n',
                'tunnelOffsetY=-0.42,\n',
                'tunnelAper1=2.75,\n',
                'tunnelFloorOffset=1.45,\n',
                'tunnelThickness=0.1,\n',
                'tunnelSoilThickness=2,\n',
                'tunnelVisible=0,\n',
                'tunnelIsInfiniteAbsorber=1,\n',
                ]
        elif self._traj and self._bpabs:
            lines = [
                'precisionRegion: cutsregion, prodCutPhotons=1e-6, prodCutElectrons=1e-3, prodCutPositrons=1e-3;\n',
                '! physics options - full physics\n',
                'option, physicsList="G4FTFP_BERT",\n',
                'geant4PhysicsMacroFileName="emextraphysics.mac",\n',
                'minimumKineticEnergy = 7e-3,\n',
                'particlesToExcludeFromCuts ="22",\n',
                'magnetGeometryType="none",\n',
                'beamPipeIsInfiniteAbsorber=1,\n',
                'storeElossLocal=1,\n',
                'storeElossLinks=1,\n',
                'storeCollimatorHits = 1,\n',
                'storeCollimatorHitsAll = 1,\n',
                'storeCollimatorInfo=1,\n',
                'includeFringeFields=0,\n',
                'beampipeMaterial="Cu",\n',
                'apertureType="pointsfile:FCC_30mm.dat:mm",\n',
                #'beampipeRadius=35e-3,\n'
                'beampipeThickness=3e-3,\n',
                ('buildTunnel=0,\n' if (self._stage == 2 and self._slimGeometry)
                 else 'buildTunnel=1,\n'),
                'tunnelOffsetX=-0.3,\n',
                'tunnelOffsetY=-0.42,\n',
                'tunnelAper1=2.75,\n',
                'tunnelFloorOffset=1.45,\n',
                'tunnelThickness=0.1,\n',
                'tunnelSoilThickness=2,\n',
                'tunnelVisible=0,\n',
                'tunnelIsInfiniteAbsorber=1,\n',
                'storeTrajectory = 1,\n',
                'storeTrajectoryParticleID = "22",\n', # Store only photons
                'storeTrajectoryStepPoints  = 1,\n', # Store the first tracked point
                'storeTrajectoryStepPointLast  = 1,\n', # Store the last tracked point
                'trajectoryFilterLogicAND = 1;\n', # Exclude the primaries
                ]
        else:
            lines = [
                'precisionRegion: cutsregion, prodCutPhotons=1e-6, prodCutElectrons=1e-3, prodCutPositrons=1e-3;\n',
                '! physics options - full physics\n',
                'option, physicsList="G4FTFP_BERT",\n',
                # 'option, physicsList="xray_reflection em synch_rad em_extra",\n',
                # 'xrayAllSurfaceRoughness=100e-9,\n',
                'geant4PhysicsMacroFileName="emextraphysics.mac",\n',
                'useGammaToMuMu = 1, \n',
                'minimumKineticEnergy = 7e-3,\n',
                'particlesToExcludeFromCuts ="22",\n',
                'magnetGeometryType="none",\n',
                'includeFringeFields=1,\n',
                'beampipeMaterial="Cu",\n',
                'apertureType="pointsfile:FCC_30mm.dat:mm",\n',
                #'beampipeRadius=35e-3,\n',
                'beampipeThickness=3e-3,\n',
                ('buildTunnel=0,\n' if (self._stage == 2 and self._slimGeometry)
                 else 'buildTunnel=1,\n'),
                'tunnelOffsetX=-0.3,\n',
                'tunnelOffsetY=-0.42,\n',
                'tunnelAper1=2.75,\n',
                'tunnelFloorOffset=1.45,\n',
                'tunnelThickness=0.1,\n',
                'tunnelSoilThickness=2,\n',
                'tunnelVisible=0,\n',
                'tunnelIsInfiniteAbsorber=1;\n',
                # NB: in LCC_v1 the two lines below were missing the separating
                # commas and were silently concatenated into a single string.
                # The resulting GMAD was still valid, but only by accident.
                # 36 mm matches the QD0AL sampler of LCC_v2, where the r < 18 mm
                # acceptance cut is applied anyway. The STAGE-1 sampler must be
                # much wider: it is a phase-space handoff, not an acceptance
                # plane, and anything it fails to record is silently lost from
                # the sample. At 36 mm the recorded photons pile up against
                # |x| = 17.98 mm and the two-stage yield came out 6.5 % low.
                ('option, samplerDiameter=500*mm;\n' if self._stage == 1
                 else 'option, samplerDiameter=36*mm;\n'),
                # Stage 1 samples EVERY species at the split plane -- the e+ that
                # stage 2 will replay, and the photons and shower products already
                # made upstream, which must be carried across or they are lost.
                (f'sample, range={self._stage1_sampler};\n' if self._stage == 1
                 else 'sample, range=QD0AL, partID={22};\n'),
                ]
        with open("GMAD/input_options.gmad", "w") as file:
            file.writelines(lines)

        print("Modifying GMAD component file\n")
        if self._withSol:
            with open("GMAD/input_components.gmad", "a") as file:
                file.write('DRIFT_L1: drift, l=0.04, aper1=0.015, apertureType="circular", fieldAll="d1field";\n')
                file.write('DRIFT_L0: drift, l=0.24, aper1=0.015, apertureType="circular", fieldAll="d0b";\n')
                file.write('DRIFT_SOL: element, fieldAll="detectorfield", geometryFile="gdml:CC_geometry.gdml", l=2*2.1000000000014825, stripOuterVolume=0;\n')
                file.write('DRIFT_R1: drift, l=0.3, apertureType="circular", aper1=15e-3, fieldAll="d2field";\n')
                file.write(f'MASK_QC2L: ecol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[0]}, ysize={self._maskA[1]};\n')
                file.write(f'DRIFT_L2: drift, l=0.710000000001173, aper1=0.025, apertureType="circular";\n')
            if self._roundmask==2:
                with open("GMAD/input_components.gmad", "a") as file:
                    file.write('MASK_QC1L: element, l=0.06, geometryFile="gdml:elliptical_mask7_9.gdml", stripOuterVolume=0, markAsCollimator=1, fieldAll="maskfield";\n')
            elif self._roundmask==3:
                with open("GMAD/input_components.gmad", "a") as file:
                    file.write(f'MASK_QC1L: jcol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[2]}, ysize=0.015, fieldAll="dmask";\n')
            else:
                with open("GMAD/input_components.gmad", "a") as file:
                    file.write(f'MASK_QC1L: ecol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[2]}, ysize={self._maskA[3]}, fieldAll="dmask";\n')
            if self._withCorr:
                with open("GMAD/input_options.gmad", "a") as file:
                    file.write('detectorfield: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad;\n')
                    file.write('maskfield: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.03;\n')
                    file.write('d1field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.26+.02;\n')
                    file.write('d2field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=-2.1000000000014825-.15;\n')
                    file.write('dmask: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.25;\n')
                    file.write('d0b: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.12;\n')
            else:
                with open("GMAD/input_options.gmad", "a") as file:
                    file.write('detectorfield: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad;\n')
                    file.write('maskfield: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.03;\n')
                    file.write('d1field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.26+.02;\n')
                    file.write('d2field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=-2.1000000000014825-.15;\n')
                    file.write('dmask: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.25;\n')
                    file.write('d0b: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.1000000000014825+.12;\n')
        else:
            with open("GMAD/input_components.gmad", "a") as file:
                file.write('DRIFT_L1: drift, l=0.04, aper1=0.015, apertureType="circular";\n')
                file.write('DRIFT_L0: drift, l=0.24, aper1=0.015, apertureType="circular";\n')
                file.write('DRIFT_SOL: element, geometryFile="gdml:CC_geometry.gdml", l=2*2.1000000000014825, stripOuterVolume=0;\n')
                file.write('DRIFT_R1: drift, l=0.3, apertureType="circular", aper1=15e-3;\n')
                file.write(f'MASK_QC2L: ecol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[0]}, ysize={self._maskA[1]};\n')
                file.write(f'DRIFT_L2: drift, l=0.710000000001173, aper1=0.025, apertureType="circular";\n')
                if self._roundmask==2:
                    with open("GMAD/input_components.gmad", "a") as file:
                        file.write('MASK_QC1L: element, l=0.06, geometryFile="gdml:elliptical_mask7_9.gdml", stripOuterVolume=0, markAsCollimator=1, fieldAll="maskfield";\n')
                elif self._roundmask==3:
                    with open("GMAD/input_components.gmad", "a") as file:
                        file.write(f'MASK_QC1L: jcol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[2]}, ysize=0.015;\n')
                else:
                    with open("GMAD/input_components.gmad", "a") as file:
                        file.write(f'MASK_QC1L: ecol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[2]}, ysize={self._maskA[3]};\n')

        if self._stage == 2:
            # STAGE 2 -- the beam is the stage-1 dump, replayed. Every species
            # crossing the split plane is injected, so the photons and shower
            # products already produced upstream are carried across rather than
            # lost; pdgid carries the species and t the absolute arrival time, so
            # the time reference at QD0AL stays that of the full beamline.
            with open("GMAD/input_beam.gmad", "r") as file:
                text = file.readlines()
            for i, line in enumerate(text):
                if "gausstwiss" in line:
                    text[i] = line.replace("gausstwiss", "userfile")
            with open("GMAD/input_beam.gmad", "w") as file:
                file.writelines(text)
            # Force the reference particle to the species being injected.
            species = self._stage2_species
            for i, line in enumerate(text):
                if "particle" in line and "=" in line:
                    text[i] = re.sub(r'particle\s*=\s*"[^"]*"',
                                     f'particle="{species}"', text[i])
            with open("GMAD/input_beam.gmad", "w") as file:
                file.writelines(text)

            add_before_last_semicolon("GMAD/input_beam.gmad",
                                      f'\tdistrFile = "{STAGE2_INPUT[species]}"')
            add_before_last_semicolon("GMAD/input_beam.gmad",
                                      f'\tdistrFileFormat = "{STAGE2_FORMAT}"')

        elif self._userfile==1: # GAUSSTWISS BEAM
            print("Modifying GMAD beam input file\n")
            with open("GMAD/input_beam.gmad", "r") as file:
                text= file.readlines()
            # Match whatever number pybdsim wrote rather than the literal "0.0":
            # it writes e.g. "Xp0=-0.0", so the LCC_v1 string replace for Xp0
            # never fired and the XP0 argument was silently ignored.
            for i, line in enumerate(text):
                line = re.sub(r'X0=[-+0-9.eE]+\*m', f'X0={self._X0}*m', line)
                line = re.sub(r'Y0=[-+0-9.eE]+\*m', f'Y0={self._Y0}*m', line)
                line = re.sub(r'Xp0=[-+0-9.eE]+', f'Xp0={self._XP0}', line)
                line = re.sub(r'Yp0=[-+0-9.eE]+', f'Yp0={self._YP0}', line)
                text[i] = line
            with open("GMAD/input_beam.gmad", "w") as file:
                file.writelines(text)

        elif self._userfile==2: # HALOSIGMA BEAM
            with open("GMAD/input_beam.gmad", "r") as file:
                text= file.readlines()
            i = 0
            while i<len(text):
                if "gausstwiss" in text[i]:
                    text[i] =text[i].replace("gausstwiss", "halo")
                i+=1
            with open("GMAD/input_beam.gmad", "w") as file:
                file.writelines(text)
            add_before_last_semicolon("GMAD/input_beam.gmad", '\thaloNSigmaXInner = 3.5')
            add_before_last_semicolon("GMAD/input_beam.gmad", '\thaloNSigmaYInner = 4.0')
            add_before_last_semicolon("GMAD/input_beam.gmad", f'\thaloNSigmaXOuter = {self._xtail}')
            add_before_last_semicolon("GMAD/input_beam.gmad", f'\thaloNSigmaYOuter = {self._ytail}')

        elif self._userfile==3: # USER BEAM - EXP HALO
            with open("GMAD/input_beam.gmad", "r") as file:
                text= file.readlines()
            i = 0
            while i<len(text):
                if "gausstwiss" in text[i]:
                    text[i] =text[i].replace("gausstwiss", "userfile")
                i+=1
            with open("GMAD/input_beam.gmad", "w") as file:
                file.writelines(text)
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFile = "inputfile.dat"')
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFileFormat = "x[m]:xp[rad]:y[m]:yp[rad]:E[GeV]"')
            bb = generate_4d_distribution(self, twiss_file, idx_start)
            bb = pd.DataFrame(bb).T
            bb.to_csv('GMAD/inputfile.dat', sep='\t', index=False, header=False)

        elif self._userfile==4: # Injection beam at FFQ integrated over 100 turns
            with open("GMAD/input_beam.gmad", "r") as file:
                text = file.readlines()

            # Replace gausstwiss -> userfile
            for i, line in enumerate(text):
                if "gausstwiss" in line:
                    text[i] = line.replace("gausstwiss", "userfile")

            # MODIFY S0 TO THE FFQ POSITION (make sure it's a scalar)
            # S_QC2L2 = twiss_file.iloc[idx_QC2L2].S
            # S_START = twiss_file.iloc[idx_start].S
            # L_QC2L2 = twiss_file.iloc[idx_QC2L2].L
            # Considering that the beam line won't be as in twiss file
            # S0 = float(S_QC2L2 - S_START - L_QC2L2 + self._deltaS)
            S0 = float(self._deltaS)

            beam_idx = None
            for i, line in enumerate(text):
                # catches: "beam," or "beam,   X0=..." etc.
                if line.lstrip().startswith("beam,"):
                    beam_idx = i
                    break

            if beam_idx is None:
                raise ValueError("No 'beam,' definition found in GMAD/input_beam.gmad")

            # Insert a new line right AFTER the 'beam,' line
            indent = "\t"  # match your file formatting; could also use spaces
            text.insert(beam_idx + 1, f"{indent}S0={S0}*m,\n")

            # 4) Write back
            with open("GMAD/input_beam.gmad", "w") as file:
                file.writelines(text)
  
            # FILE IS ALREADY PREPARED OUTSIDE WITH THE NAME injectionpart_ffq.dat --> see xutil
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFile = "injectionpart_ffq.dat"')
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFileFormat = "x[m]:xp[rad]:y[m]:yp[rad]:E[GeV]"')
        elif self._userfile==5: # Injection beam at B0BL integrated over 100 turns
            with open("GMAD/input_beam.gmad", "r") as file:
                text = file.readlines()

            # Replace gausstwiss -> userfile
            for i, line in enumerate(text):
                if "gausstwiss" in line:
                    text[i] = line.replace("gausstwiss", "userfile")
            # 4) Write back
            with open("GMAD/input_beam.gmad", "w") as file:
                file.writelines(text)

            # tab['s','ipg.1']-4.03734033    
            # FILE IS ALREADY PREPARED OUTSIDE WITH THE NAME injectionpart.dat --> see xutil
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFile = "injectionpart_b0bl.dat"')
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFileFormat = "x[m]:xp[rad]:y[m]:yp[rad]:E[GeV]"')
        else:
            raise ValueError("Userfile must be 1, 2, 3 or 4")  # keep the original gausstwiss beam

        return -(twiss_file.iloc[idx_start-1]['S']-twiss_file.iloc[idx_IP]['S'])

    def runBDSIM(self):
        """
        Runs BDSIM, must call genGMAD() before running (unless files are provided manually).
        If providing manually the main gmad must be in the directory as 'GMAD/input.gmad'.
        """
        # runKey lets stage 2 write one file per species without a seed clash;
        # it is "" everywhere else, so single-stage naming is unchanged.
        outfile = f"output_{self._seed}{self._runKey}"
        runOptions = f"--seed={self._seed}"

        # run bdsim
        pybdsim.Run.Bdsim('GMAD/input.gmad', outfile, ngenerate=self._ngenerate, options=runOptions, batch=True)


