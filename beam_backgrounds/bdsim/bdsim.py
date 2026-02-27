"""
Create the GMAD input files for BDSIM from MAD-X Twiss output files in tfs format.
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
import argparse
import json

def generate_4d_distribution(self, twiss_file, idx_start):
    np.random.seed(self._seed)
    size = self._ngenerate
    alfx=twiss_file.iloc[idx_start-1]['ALFX'] #50.02904476961108
    alfy=twiss_file.iloc[idx_start-1]['ALFY']
    betx=twiss_file.iloc[idx_start-1]['BETX'] #7861.70688987599
    bety=twiss_file.iloc[idx_start-1]['BETY'] #285.5641279561109
    emitx = 0.74e-9
    emity = 2.5e-12
    gammax = (1.0 + alfx * alfx) / betx
    gammay = (1.0 + alfy * alfy) / bety

    haloNSigmaXInner = 3.5
    haloNSigmaXOuter = self._xtail
    haloNSigmaYInner = 4
    haloNSigmaYOuter = self._ytail

    sigmaX     = np.sqrt(emitx * betx)
    sigmaY     = np.sqrt(emity * bety)
    sigmaXp    = np.sqrt(gammax * emitx)
    sigmaYp    = np.sqrt(gammay * emity)

    emitInnerX = haloNSigmaXInner**2 * emitx
    emitInnerY = haloNSigmaYInner**2 * emity
    emitOuterX = haloNSigmaXOuter**2 * emitx
    emitOuterY = haloNSigmaYOuter**2 * emity

    xMax  = haloNSigmaXOuter * sigmaX
    yMax  = haloNSigmaYOuter * sigmaY
    xpMax = haloNSigmaXOuter * sigmaXp
    ypMax = haloNSigmaYOuter * sigmaYp

    # Generate random displacement values for x, y, dxp, and dyp
    init_size = 250*size
    dx = xMax * (1 - 2 * np.random.uniform(size=init_size))
    dy = yMax * (1 - 2 * np.random.uniform(size=init_size))
    dxp = xpMax * (1 - 2 * np.random.uniform(size=init_size))
    dyp = ypMax * (1 - 2 * np.random.uniform(size=init_size))

    # Calculate beam parameters for x and y directions
    emitXSp = gammax * dx**2 + 2. * alfx * dx * dxp + betx * dxp**2
    emitYSp = gammay * dy**2 + 2. * alfy * dy * dyp + bety * dyp**2

    # Calculate exponential weights
    exp_weights = np.exp(- self._xweight * emitXSp / emitInnerX - self._yweight * emitYSp / emitInnerY)

    # Apply the condition to filter based on emitXSp and emitYSp ranges
    condition = (
        (abs(emitXSp) > emitInnerX) &
        (abs(emitYSp) > emitInnerY) &
        (abs(emitXSp) < emitOuterX) &
        (abs(emitYSp) < emitOuterY)
    )

    valid_indices = np.where(condition)[0]

    normalised_exp_weight = exp_weights[valid_indices]/sum(exp_weights[valid_indices])

    # Sample random values based on the normalized exponential weights and the condition
    sampled_indices = []

    while len(sampled_indices) < size:
        sampled_indices.extend(np.random.choice(valid_indices, size=size - len(sampled_indices), replace=True, p=normalised_exp_weight))

    sampled_indices = np.array(sampled_indices)
    sampled_dx = dx[sampled_indices]
    sampled_dxp = dxp[sampled_indices]
    sampled_dy = dy[sampled_indices]
    sampled_dyp = dyp[sampled_indices]
    sampled_E = np.random.normal(45.6, 1.0e-3, size=size)

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
        self._extraA	= kwargs.get("extraA", ["21.8e-3", "6.98e-3", "26.8e-3", "22.8e-3", "8.03e-3", "19.3e-3"])	 # Aperture of BWL Hor, QC0L vert and hor(1 and 2) and QC2L vert an hor collimators

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
        ax2.set_ylim(-20, 20)
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
        
        headers = pd.read_csv('GMAD/fcc_ee_z_b1_twiss_end.tfs', skiprows=MADX_TWISS_HEADERS_SKIP_ROWS,
                        nrows=0, sep=r"\s+")        
        headers.drop(headers.columns[[0, 1]], inplace=True, axis=1)
        twiss_file = pd.read_csv('GMAD/fcc_ee_z_b1_twiss_end.tfs',
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
        idx_IP = twiss_file.reset_index()[twiss_file.reset_index()['NAME']=="IPG.1"].index[0]
        # Find the first dipole BEFORE the ip to start the sequence conversion
        found, idx_start = 0, idx_IP
        while found<1+self._withBC1L:
            if twiss_file.iloc[idx_start]['KEYWORD'] == "RBEND":
                found += 1
            idx_start -=1

        # Find the first dipole AFTER the IP to start the sequence conversion
        # Can be adapted to get any dipole after the IP to have a longer beam line
        found, idx_stop = 0, idx_IP
        while found<self._withDip:
            if twiss_file.iloc[idx_stop]['KEYWORD'] == "RBEND":
                found += 1
            idx_stop +=1


        drift_pre_mask, drift_pre_mask2, drift_pre_ip, drift_post_ip = twiss_file.iloc[idx_IP-3].name, twiss_file.iloc[idx_IP-11].name, twiss_file.iloc[idx_IP-1].name, twiss_file.iloc[idx_IP+2].name

        IR = twiss_file.iloc[idx_start:idx_stop]
        col_name = IR[IR["KEYWORD"]=="COLLIMATOR"].index

        print("Write collimator settings file\n")
        print(self._extraA)
        f = open("collimatorSettings.dat", "w")
        lines = [
            "# Collimator Settings\n",
            "name\tmaterial\txsize[m]\tysize[m]\n",
            f"{col_name[0]}\tinermet180\t{self._extraA[0]}\t30e-3\n",
            f"{col_name[1]}\tinermet180\t30e-3\t{self._extraA[1]}\n",
            f"{col_name[2]}\tinermet180\t{self._extraA[2]}\t30e-3\n",
            f"{col_name[3]}\tinermet180\t{self._extraA[3]}\t30e-3\n",
            f"{col_name[4]}\tinermet180\t30e-3\t{self._extraA[4]}\n",
            f"{col_name[5]}\tinermet180\t{self._extraA[5]}\t30e-3\n",
            f"{col_name[6]}\tinermet180\t15e-3\t15e-3\n",
            f"{col_name[7]}\tinermet180\t15e-3\t15e-3\n",
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
        while found<5: # 4 is the number of quadrupole to form the doublet in the LCCO lattice
            if twiss_file.iloc[idx_start_FF]['KEYWORD'] == "QUADRUPOLE":
                found += 1
            idx_start_FF -=1
        found, idx_stop_FF = 0, idx_IP
        while found<5: # 4 is the number of quadrupole to form the doublet in the LCCO lattice
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


        idx_QC1L1 = twiss_file.reset_index()[twiss_file.reset_index()['NAME']=="QC1L1.2"].index[0]
        idx_stop = idx_QC1L1 +1

        if self._userfile == 4:
            print("Using injection particles at FFQ, thus cutting the beamline before the FFQ\n")
            idx_QC2L2 = twiss_file.reset_index()[twiss_file.reset_index()['NAME']=="QC2L2.2"].index[0]
            idx_start = idx_QC2L2-2

        a, b = pybdsim.Convert.MadxTfs2Gmad("GMAD/fcc_ee_z_b1_twiss_end.tfs",
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
            if col_name[6].replace(".", "") in text[i]:
                text[i] =text[i].replace(col_name[6].replace(".", ""), "MASK_QC2L")
            if drift_pre_mask2 in text[i]:
                text[i] =text[i].replace(f"{drift_pre_mask2}", "DRIFT_L2")
            if drift_pre_mask in text[i]:
                text[i] =text[i].replace(f"{drift_pre_mask}", "DRIFT_L1")
            if col_name[7].replace(".", "") in text[i]:
                if self._roundmask==2:
                    text[i] =text[i].replace(col_name[7].replace(".", ""), "MASK_QC1L")
                else:
                    text[i] =text[i].replace(col_name[7].replace(".", ""), "MASK_QC1L, DRIFT_L0")
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
                'buildTunnel=1,\n',
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
                'buildTunnel=1,\n',
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
                'geant4PhysicsMacroFileName="emextraphysics.mac",\n',
                'minimumKineticEnergy = 7e-3,\n',
                'particlesToExcludeFromCuts ="22",\n',
                'magnetGeometryType="none",\n',
                'includeFringeFields=1,\n',
                'beampipeMaterial="Cu",\n',
                'apertureType="pointsfile:FCC_30mm.dat:mm",\n',
                #'beampipeRadius=35e-3,\n',
                'beampipeThickness=3e-3,\n',
                'buildTunnel=1,\n',
                'tunnelOffsetX=-0.3,\n',
                'tunnelOffsetY=-0.42,\n',
                'tunnelAper1=2.75,\n',
                'tunnelFloorOffset=1.45,\n',
                'tunnelThickness=0.1,\n',
                'tunnelSoilThickness=2,\n',
                'tunnelVisible=0,\n',
                'tunnelIsInfiniteAbsorber=1;\n',
                #'csample, range=QC1L12, partID={22};'
                'sample, range=QC1L12, partID={22};'
                ]
        with open("GMAD/input_options.gmad", "w") as file:
            file.writelines(lines)

        print("Modifying GMAD component file\n")
        if self._withSol:
            with open("GMAD/input_components.gmad", "a") as file:
                file.write('DRIFT_L1: drift, l=0.04, aper1=0.015, apertureType="circular", fieldAll="d1field";\n')
                file.write('DRIFT_L0: drift, l=0.24, aper1=0.015, apertureType="circular", fieldAll="d0b";\n')
                file.write('DRIFT_SOL: element, fieldAll="detectorfield", geometryFile="gdml:CC_geometry.gdml", l=2*2.100247523199509, stripOuterVolume=0;\n')
                file.write('DRIFT_R1: drift, l=0.3, apertureType="circular", aper1=15e-3, fieldAll="d2field";\n')
                file.write(f'MASK_QC2L: ecol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[0]}, ysize={self._maskA[1]};\n')
                file.write(f'DRIFT_L2: drift, l=0.659999999998263, aper1=0.025, apertureType="circular";\n')
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
                    file.write('maskfield: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.03;\n')
                    file.write('d1field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.26+.02;\n')
                    file.write('d2field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=-2.100247523199509-.15;\n')
                    file.write('dmask: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.25;\n')
                    file.write('d0b: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map_corr.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.12;\n')
            else:
                with open("GMAD/input_options.gmad", "a") as file:
                    file.write('detectorfield: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad;\n')
                    file.write('maskfield: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.03;\n')
                    file.write('d1field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.26+.02;\n')
                    file.write('d2field: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=-2.100247523199509-.15;\n')
                    file.write('dmask: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.25;\n')
                    file.write('d0b: field, type="bmap3d", magneticFile = "bdsim3d:3D_field_map.dat", axisAngle=1, axisY=1, angle=-15*mrad, z=2.100247523199509+.12;\n')
        else:
            with open("GMAD/input_components.gmad", "a") as file:
                file.write('DRIFT_L1: drift, l=0.04, aper1=0.015, apertureType="circular";\n')
                file.write('DRIFT_L0: drift, l=0.24, aper1=0.015, apertureType="circular";\n')
                file.write('DRIFT_SOL: element, geometryFile="gdml:CC_geometry.gdml", l=2*2.100247523199509, stripOuterVolume=0;\n')
                file.write('DRIFT_R1: drift, l=0.3, apertureType="circular", aper1=15e-3;\n')
                file.write(f'MASK_QC2L: ecol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[0]}, ysize={self._maskA[1]};\n')
                file.write(f'DRIFT_L2: drift, l=0.659999999998263, aper1=0.025, apertureType="circular";\n')
                if self._roundmask==2:
                    with open("GMAD/input_components.gmad", "a") as file:
                        file.write('MASK_QC1L: element, l=0.06, geometryFile="gdml:elliptical_mask7_9.gdml", stripOuterVolume=0, markAsCollimator=1, fieldAll="maskfield";\n')
                elif self._roundmask==3:
                    with open("GMAD/input_components.gmad", "a") as file:
                        file.write(f'MASK_QC1L: jcol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[2]}, ysize=0.015;\n')
                else:
                    with open("GMAD/input_components.gmad", "a") as file:
                        file.write(f'MASK_QC1L: ecol, horizontalWidth=0.046, l=0.02, material="W", region="precisionRegion", xsize={self._maskA[2]}, ysize={self._maskA[3]};\n')

        if self._userfile==1: # GAUSSTWISS BEAM
            print("Modifying GMAD beam input file\n")
            with open("GMAD/input_beam.gmad", "r") as file:
                text= file.readlines()
            i = 0
            while i<len(text):
                if "X0" in text[i]:
                    text[i] =text[i].replace("X0=0.0*m", f'X0={self._X0}*m')
                if "Y0" in text[i]:
                    text[i] =text[i].replace("Y0=0.0*m", f'Y0={self._Y0}*m')
                if "Xp0" in text[i]:
                    text[i] =text[i].replace("Xp0=0.0", f'Xp0={self._XP0}')
                if "Yp0" in text[i]:
                    text[i] =text[i].replace("Yp0=0.0", f'Yp0={self._YP0}')
                i+=1
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
            S_QC2L2 = twiss_file.iloc[idx_QC2L2].S
            S_START = twiss_file.iloc[idx_start].S
            L_QC2L2 = twiss_file.iloc[idx_QC2L2].L
            # Considering that the beam line won't be as in twiss file
            S0 = float(S_QC2L2 - S_START - L_QC2L2/2 + self._deltaS)

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

            # tab['s','ipg.1']-4.03734033    
            # FILE IS ALREADY PREPARED OUTSIDE WITH THE NAME injectionpart.dat --> see xutil
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFile = "injectionpart.dat"')
            add_before_last_semicolon("GMAD/input_beam.gmad", '\tdistrFileFormat = "x[m]:xp[rad]:y[m]:yp[rad]:E[GeV]"')
        else:
            raise ValueError("Userfile must be 1, 2, 3 or 4")  # keep the original gausstwiss beam

        return -(twiss_file.iloc[idx_start-1]['S']-twiss_file.iloc[idx_IP]['S'])

    def convert_hepevt(self):
        infile = f"output_{self._seed}.root"
        outfile = f"output_{self._seed}.hepevt"

        # Read the dist file
        d = pybdsim.Data.Load(infile)
        part2 = pybdsim.Data.SamplerData(d, -1) # Only one sampler, hence take the last one, 0 is the initial positron disitribution
        part = pd.DataFrame([part2.data['x'], part2.data['xp'], part2.data['y'], part2.data['yp'], part2.data['energy'], part2.data['p'], part2.data['zp'],
        part2.data['T'], part2.data['partID'], part2.data['mass']]).T
        part.columns = ['x', 'xp', 'y', 'yp', 'E', 'p', 'zp', 't', 'partID', 'm']
        part['r'] = np.sqrt(part['x']**2 + part['y']**2)
        part['rp'] = np.sqrt((part['x']+part['xp']*(2.2+6))**2 + (part['y']+part['yp']*(2.2+6))**2)
        # Filtering
        # in energy > 2keV 
        # in radius removing photons beyond R=15mm
        # removing all photons R<9mm 8.2m downstream the sampler equivalent to s=+6.0m
        part = part[(part['partID']==22)&(part['E']>2e-6)&(part['r']<18e-3)&(part['rp']>9e-3)]
        part.reset_index(drop=True, inplace=True)
        result = part
        # Creating file in hepevt format following https://hugonweb.com/hepevt/
        # Loop through the photons and add each to a new vertex
        f = open(outfile, "w")
        lines = [
                # Write total number of events in the file
                f'{len(result)}\n'
            ]
        f.writelines(lines)

        for i in range(len(result)):
            # Conversion from sampler units to hepevt units
            # 'time' is in ns AT THE SAMPLER, converted to mm/c w.r.t. the IP
            # 'p' is in GeV/c
            # 'E' is in GeV
            # 'xp', 'yp', 'zp' are in fractional part of the momentum p, convert to GeV/c
            # 'x', 'y' are in m, convert to mm
            # z is -2.2m, location of the sampler
            partID = result.iloc[i]['partID']
            px, py, pz, E = result.iloc[i]['xp']*result.iloc[i]['p'], result.iloc[i]['yp']*result.iloc[i]['p'], result.iloc[i]['zp']*result.iloc[i]['p'], result.iloc[i]['E']
            x, y, z, t = 1e3*result.iloc[i]['x'], 1e3*result.iloc[i]['y'], -2.2e3, 1e3*(result.iloc[i]['t']*1e-9*2.997924588e8-(d.model['QC1L12']['SEnd']+2.4))

            lines = [
                # Create an active photon particle (PDG ID for photon is 22)
                # <Status> <PDG ID> <1st Mother> <2nd Mother> <1st Daughter> <2nd Daughter> <Px> <Py> <Px> <E> <Mass> <x> <y> <z> <t>
                # where Px/Py/Pz are in GeV/c, E is in GeV, and M is in GeV/c^2. x/y/z are in mm and t is in mm/c
                f'1 {partID} 0 0 0 0 {px} {py} {pz} {E} {0} {x} {y} {z} {t}\n'
            ]
            f.writelines(lines)
        f.close()


    def runBDSIM(self):
        """
        Runs BDSIM, must call genGMAD() before running (unless files are provided manually).
        If providing manually the main gmad must be in the directory as 'GMAD/input.gmad'.
        """
        outfile = f"output_{self._seed}"
        runOptions = f"--seed={self._seed}"

        # run bdsim
        pybdsim.Run.Bdsim('GMAD/input.gmad', outfile, ngenerate=self._ngenerate, options=runOptions, batch=True)





if __name__ == "__main__":
    # Parse command-line arguments for xweight and yweight
    parser = argparse.ArgumentParser(description="Run MDIStudy with varying xweight and yweight.")
    parser.add_argument("--ngenerate", type=int, help="Number of particles", default=5000)
    parser.add_argument("--userfile", type=int, help="1--> beam core simulation; 2--> beam halo simulation; 3--> beam halo with exponential density", default=1)
    parser.add_argument("--roundmask", type=float, help="1-->round mask 7mm; 2-->elliptical mask 7x8.5; else-->jaw mask", default=3)
    parser.add_argument("--xtail", type=int, help="Horizontal halo width in sigma", default=10)
    parser.add_argument("--ytail", type=int, help="Vertical halo width in sigma", default=51)
    parser.add_argument("--xweight", type=float, help="X-weight value", default=1.2)
    parser.add_argument("--yweight", type=float, help="Y-weight value", default=0.06)
    parser.add_argument("--withDip", type=int, help="Add post IP dipoles")
    parser.add_argument("--seed", type=int, help="Seed number", default=123)
    parser.add_argument("--withSol", type=int, help="1-->Including solenoid 0-->Without solenoid", default=0)
    parser.add_argument("--withCorr", type=int, help="1-->Including correction 0-->No correction", default=0)
    parser.add_argument("--cfg", type=Path, help="Optional JSON config file")

    # two-pass parse: first read --cfg if present
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--cfg", type=Path)
    ns_pre, _ = pre.parse_known_args(sys.argv)

    if ns_pre.cfg:
        print("Loading args from JSON file\n")
        cfg = json.loads(ns_pre.cfg.read_text())
        parser.set_defaults(**cfg)
    
    args = parser.parse_args()



    # Create an MDIStudy instance including the specified beam halo xweight and yweight
    mdi_study = MDIStudy(seed=args.seed, withSol=0, withCorr=0, ngenerate=args.ngenerate, withDip=True, maskA=[15e-3, 15e-3, 7e-3, 8.5e-3], userfile=args.userfile, roundmask=args.roundmask, xtail=args.xtail, ytail=args.ytail, xweight=args.xweight, yweight=args.yweight)
    # Generate the BDSIM model
    mdi_study.genGMAD()
    mdi_study.runBDSIM()
    mdi_study.convert_hepevt()
