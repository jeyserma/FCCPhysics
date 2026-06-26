
import sys, os, glob, math
import ROOT
import logging
import time
import functions
import config as gpconfig
import argparse
from pathlib import Path
import json
from scipy.constants import c, micro, nano, pi, milli, micro


logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)

ROOT.gROOT.SetBatch(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

ROOT.EnableImplicitMT(4) # use all cores
ROOT.DisableImplicitMT() # single core

# load libraries
header = Path(__file__).resolve().parent / "functions.h"
ROOT.gInterpreter.Declare(f'#include "{header}"')


layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

bins_theta = (36, 0, 180)
bins_phi = (72, -180, 180)
bins_p = (100000, 0, 10000)
bins_z = (int(max_z/2), -max_z, max_z)
bins_layer = (10, 0, 10)

bins_p_log10 = (600, -2, 4)


def meta_info(reader):
    hists = []
    global str_out, json_dict

    # beam-related quantities
    npart = reader.get_metdata('beam1_pars_particles')
    npart = reader.get_metdata('beam2_pars_particles')
    sigmax = reader.get_metdata('beam1_pars_sigma_x') # default nm
    sigmay = reader.get_metdata('beam1_pars_sigma_y') # default nm
    sigmaz = reader.get_metdata('beam1_pars_sigma_z') # default um
    beta_x = reader.get_metdata('beam1_pars_beta_x') # default um
    beta_y = reader.get_metdata('beam1_pars_beta_y') # default um
    energy = reader.get_metdata('beam1_pars_energy') # default GeV

    str_out += f"BEAM RELATED QUANTITIES\n"
    str_out += f"Energy (GeV)             {energy:.2f}\n"
    str_out += f"Bunch density            {npart:.2e}\n"
    str_out += f"sigmax (nm)              {sigmax:.2f}\n"
    str_out += f"sigmay (nm)              {sigmay:.2f}\n"
    str_out += f"sigmaz (um)              {sigmaz:.2f}\n"
    str_out += f"betax (mm)               {beta_x/1000.:.2f}\n"
    str_out += f"betay (mm)               {beta_y/1000.:.2f}\n"

    json_dict['beams'] = {}
    json_dict['beams']['energy'] = energy
    json_dict['beams']['npart'] = npart
    json_dict['beams']['sigmax'] = sigmax
    json_dict['beams']['sigmay'] = sigmay
    json_dict['beams']['sigmaz'] = sigmaz
    json_dict['beams']['betax'] = beta_x/1000.
    json_dict['beams']['betay'] = beta_y/1000.

    # GP related quantities

    cutx = reader.get_metdata('cut_x') # nm
    cuty = reader.get_metdata('cut_y') # nm
    cutz = reader.get_metdata('cut_z') # um
    minx_ = reader.get_metdata('min_x_') # nm
    maxx_ = reader.get_metdata('max_x_') # nm
    miny_ = reader.get_metdata('min_y_') # nm
    maxy_ = reader.get_metdata('max_y_') # nm
    minz_ = reader.get_metdata('min_z_') # nm
    maxz_ = reader.get_metdata('max_z_') # nm
    nx = reader.get_metdata('n_x')
    ny = reader.get_metdata('n_y')
    nz = reader.get_metdata('n_z')
    nt = reader.get_metdata('n_t')

    minx_calc_, maxx_calc_ = -(nx-2)/(nx)*cutx, (nx-2)/(nx)*cutx
    miny_calc_, maxy_calc_ = -(ny-2)/(ny)*cuty, (ny-2)/(ny)*cuty
    minz_calc_, maxz_calc_ = -cutz*1000, cutz*1000

    infx = cutx / sigmax
    infy = cutx / sigmay
    infz = cutz / sigmaz

    str_out += "\n"
    str_out += f"GRID INPUTS\n"
    str_out += f"nx                       {nx}\n"
    str_out += f"ny                       {ny}\n"
    str_out += f"nz                       {nz}\n"
    str_out += f"cutx (nm)                {cutx:.2f} (={infx:.2f} sigmax)\n"
    str_out += f"cuty (nm)                {cuty:.2f} (={infy:.2f} sigmay)\n"
    str_out += f"cutz (um)                {cutz:.2f} (={infz:.2f} sigmaz)\n"

    str_out += "\n"
    str_out += f"PRIMARY GRID DIMENSIONS\n"
    str_out += f"min/max x (nm)           {minx_:.2f}/{maxx_:.2f}\n"
    str_out += f"min/max x (nm)           {miny_:.2f}/{maxy_:.2f}\n"
    str_out += f"min/max x (um)           {minz_:.2f}/{maxz_:.2f}\n"
    str_out += f"min/max x (nm), calc     {minx_calc_:.2f}/{maxx_calc_:.2f}\n"
    str_out += f"min/max x (nm), calc     {miny_calc_:.2f}/{maxy_calc_:.2f}\n"
    str_out += f"min/max x (um), calc     {minz_calc_:.2f}/{maxz_calc_:.2f}\n"

    json_dict['grid'] = {}
    json_dict['grid']['nx'] = nx
    json_dict['grid']['ny'] = ny
    json_dict['grid']['nz'] = nz
    json_dict['grid']['cutx'] = cutx
    json_dict['grid']['cuty'] = cuty
    json_dict['grid']['cutz'] = cutz
    json_dict['grid']['minx'] = minx_
    json_dict['grid']['miny'] = miny_
    json_dict['grid']['minz'] = minz_


    str_out += "\n"
    str_out += f"LUMINOSITY DIAGNOSTICS\n"
    lumi_ee = reader.get_metdata('lumi_ee')/10000. # 1/cm2
    lumi_ee_high = reader.get_metdata('lumi_ee_high')/10000. # 1/cm2 same as ee?
    lumi_eg = reader.get_metdata('lumi_eg')/10000.# 1/cm2
    lumi_ge = reader.get_metdata('lumi_ge')/10000.# 1/cm2
    lumi_gg = reader.get_metdata('lumi_gg')/10000.# 1/cm2
    lumi_gg_high = reader.get_metdata('lumi_gg_high')/10000. # same as gg?
    upsmax = reader.get_metdata('upsmax')
    npairs = reader.get_metdata('pairs_Pairs/bunch_crossing')
    
    str_out += f"lumi_ee (1/cm2)          {lumi_ee:.2e}\n"
    str_out += f"lumi_eg (1/cm2)          {lumi_eg:.2e}\n"
    str_out += f"lumi_ge (1/cm2)          {lumi_ge:.2e}\n"
    str_out += f"lumi_gg (1/cm2)          {lumi_gg:.2e}\n"
    str_out += f"upsmax (nm)              {upsmax:.2e}\n"

    json_dict['luminosity'] = {}
    json_dict['luminosity']['lumi_ee'] = lumi_ee
    json_dict['luminosity']['lumi_eg'] = lumi_eg
    json_dict['luminosity']['lumi_ge'] = lumi_ge
    json_dict['luminosity']['lumi_gg'] = lumi_gg
    json_dict['luminosity']['upsmax'] = upsmax

    nstep = int(nz)*2
    h_lumi = ROOT.TH1D("h_lumi", "Luminosity profile;Step;Luminosity (1/cm^{2}/step)", nstep, 0, nstep)
    for istep in range(1, nstep):
        lumi_step = reader.get_metdata(f'lumi_total[{istep}]')/10000.# 1/cm2
        lumi_step1 = reader.get_metdata(f'lumi_peak[{istep}]')/10000.# 1/cm2 # same as total?
        h_lumi.SetBinContent(istep, lumi_step)
    hists.append(h_lumi)

    str_out += "\n"
    str_out += f"IPC RESULTS\n"
    # https://gitlab.cern.ch/jaeyserm/guinea-pig/-/blob/master/src/resultsCPP.cc
    npairs = reader.get_metdata('pairs_n_pairs')
    pairs_e_pairs = reader.get_metdata('pairs_e_pairs') # total energy, GeV
    pairs_n_BW = reader.get_metdata('pairs_n_BW')
    pairs_e_BW = reader.get_metdata('pairs_e_BW') # GeV
    pairs_n_BH = reader.get_metdata('pairs_n_BH') # somehow different than counted in the pair files?
    pairs_e_BH = reader.get_metdata('pairs_e_BH') # GeV
    pairs_n_LL = reader.get_metdata('pairs_n_LL')
    pairs_e_LL = reader.get_metdata('pairs_e_LL') # GeV
    pairs_n1 = reader.get_metdata('pairs_n.1')
    pairs_b1 = reader.get_metdata('pairs_b.1')
    pairs_n2 = reader.get_metdata('pairs_n.2')
    pairs_b2 = reader.get_metdata('pairs_b.2')
    str_out += f"Number of pairs          {npairs}\n"
    str_out += f"Total pair energy ()     {pairs_e_pairs}\n"
    str_out += f"Number of BW pairs       {pairs_n_BW}\n"
    str_out += f"Energy of BW pairs       {pairs_e_BW}\n"
    str_out += f"Number of BH pairs       {pairs_n_BH}\n"
    str_out += f"Energy of BH pairs       {pairs_e_BH}\n"
    str_out += f"Number of LL pairs       {pairs_n_LL}\n"
    str_out += f"Energy of LL pairs       {pairs_e_LL}\n"

    
    json_dict['ipcs'] = {}
    json_dict['ipcs']['pairs_n'] = npairs
    json_dict['ipcs']['pairs_e'] = pairs_e_pairs
    json_dict['ipcs']['pairs_n_BW'] = pairs_n_BW
    json_dict['ipcs']['pairs_e_BW'] = pairs_e_BW
    json_dict['ipcs']['pairs_n_BH'] = pairs_n_BH
    json_dict['ipcs']['pairs_e_BH'] = pairs_e_BH
    json_dict['ipcs']['pairs_n_LL'] = pairs_n_LL
    json_dict['ipcs']['pairs_e_LL'] = pairs_e_LL



    return hists


def analysis(df, reader, ptype):

    # get meta information and define spatial binning
    n_x = reader.get_metdata('n_x')
    n_y = reader.get_metdata('n_y')
    n_z = reader.get_metdata('n_z')
    min_x_ = reader.get_metdata('min_x_') # inner grid, nm
    max_x_ = reader.get_metdata('max_x_') # inner grid, nm
    min_y_ = reader.get_metdata('min_y_') # inner grid, nm
    max_y_ = reader.get_metdata('max_y_') # inner grid, nm
    min_z_ = reader.get_metdata('min_z_') # inner grid, nm
    max_z_ = reader.get_metdata('max_z_') # inner grid, nm

    if ptype == 0:
        infl_x = 1.0 # adapt according to production volume wrt inner grid
        infl_y = 0.15
        infl_z = 0.15
        bins_x = (int(n_x)*5, infl_x*min_x_/1e3, infl_x*max_x_/1e3) # um
        bins_y = (int(n_y)*5, infl_y*min_y_/1e0, infl_y*max_y_/1e0) # nm
        bins_z = (int(n_z)*5, infl_z*min_z_/1e6, infl_z*max_z_/1e6) # mm

    else:
        infl_xy = 12. # according to grids8
        infl_z = 1.5

        bins_x = (int(n_x)*5, infl_xy*min_x_/1e6, infl_xy*max_x_/1e6) # mm
        bins_y = bins_x # square grid
        bins_z = (int(n_z)*5, infl_z*min_z_/1e6, infl_z*max_z_/1e6) # mm       

    bins_t = (int(n_z)*2, 0, n_z*2)







    hists = []


    # MC particle kinematics
    df = df.Alias("pairs", "Pairs0" if ptype == 0 else "Pairs")
    df = df.Alias("pairs_type", "Pairs0Process" if ptype == 0 else "PairsProcess")
    df = df.Define("sel_bw", "pairs_type == 0")
    df = df.Define("sel_bh", "pairs_type == 1")
    df = df.Define("sel_ll", "pairs_type == 2")


    df = df.Define("pairs_p4", "makeP4Vector(pairs, 0.015)")
    #df = df.Define("pairs_p4", "makeP4Vector(pairs, 0.0)")
    df = df.Define("pairs_E", "get_E(pairs_p4)*1000") # MeV
    df = df.Define("pairs_p", "get_p(pairs_p4)*1000") # MeV
    df = df.Define("pairs_pt", "get_pt(pairs_p4)*1000") # MeV
    df = df.Define("pairs_theta", "get_theta(pairs_p4, true, false)")
    df = df.Define("pairs_phi", "get_phi(pairs_p4, true)")
    df = df.Define("pairs_n", "size(pairs)/2")
    



    hists.append(df.Histo1D(("pairs_E", "Energy distribution;E (MeV);Entries", *bins_p), "pairs_E"))
    hists.append(df.Histo1D(("pairs_p", "Momentum distribution;p (MeV);Entries", *bins_p), "pairs_p"))
    hists.append(df.Histo1D(("pairs_pt", "Transverse momentum distribution;p_{{T}} (MeV);Entries", *bins_p), "pairs_pt"))
    hists.append(df.Histo1D(("pairs_theta", "Polar angle;#theta (deg);Entries", *bins_theta), "pairs_theta"))
    hists.append(df.Histo1D(("pairs_phi", "Azimuthal angle;#phi (deg);Entries", *bins_phi), "pairs_phi"))
    hists.append(df.Histo1D(("pairs_n", "Number of pairs;Pair multiplicity;Entries", *(10000, 0, 10000)), "pairs_n"))
    hists.append(df.Histo1D(("pairs_type", "IPC type;IPC type;Entries", *(5, 0, 5)), "pairs_type"))
    
    df = df.Define("pairs_E_log10", "log10(pairs_E)")
    df = df.Define("pairs_p_log10", "log10(pairs_p)")
    df = df.Define("pairs_pt_log10", "log10(pairs_pt)")
    hists.append(df.Histo1D(("pairs_E_log10", f"Energy distribution;log_{{10}}(E) (MeV);Entries", *bins_p_log10), "pairs_E_log10"))
    hists.append(df.Histo1D(("pairs_p_log10", f"Momentum distribution;log_{{10}}(p) (MeV);Entries", *bins_p_log10), "pairs_p_log10"))
    hists.append(df.Histo1D(("pairs_pt_log10", f"Transverse momentum distribution;log_{{10}}(p_{{T}}) (MeV);Entries", *bins_p_log10), "pairs_pt_log10"))

    # 2d theta-pT log plot
    df = df.Define("pairs_theta_rad", "get_theta(pairs_p4, false, true)") # from 0 to pi/2
    df = df.Define("pairs_theta_rad_log10", "log10(pairs_theta_rad)")
    
    hists.append(df.Histo2D(("pairs_theta_pt", f"#theta-p_{{T}};log_{{10}}(#theta) (rad);log_{{10}}(p_{{T}}) (MeV)", *((500, -4, 1) + bins_p_log10)), "pairs_theta_rad_log10", "pairs_pt_log10"))

    if ptype == 0:
        df = df.Define("pairs_x", "get_x(pairs)/1e3") # um
        df = df.Define("pairs_y", "get_y(pairs)") # nm
        df = df.Define("pairs_z", "get_z(pairs)/1e6") # mm

        hists.append(df.Histo1D(("pairs_x", "x distribution;x (#mum);Entries", *bins_x), "pairs_x"))
        hists.append(df.Histo1D(("pairs_y", "y distribution;y (nm);Entries", *bins_y), "pairs_y"))
        hists.append(df.Histo1D(("pairs_z", "z distribution;z (mm);Entries", *bins_z), "pairs_z"))

        hists.append(df.Histo2D(("pairs_xy", "xy distribution;x (#mum);y (nm)", *(bins_x + bins_y)), "pairs_x", "pairs_y"))
        hists.append(df.Histo2D(("pairs_xz", "xz distribution;x (#mum);z (mm)", *(bins_x + bins_z)), "pairs_x", "pairs_z"))
        hists.append(df.Histo2D(("pairs_yz", "yz distribution;y (nm);z (mm)", *(bins_y + bins_z)), "pairs_y", "pairs_z"))

    if ptype == 1:
        df = df.Define("pairs_x", "get_x(pairs)/1e6")
        df = df.Define("pairs_y", "get_y(pairs)/1e6")
        df = df.Define("pairs_z", "get_z(pairs)/1e6")

        hists.append(df.Histo1D(("pairs_x", "x distribution;x (mm);Entries", *bins_x), "pairs_x"))
        hists.append(df.Histo1D(("pairs_y", "y distribution;y (mm);Entries", *bins_y), "pairs_y"))
        hists.append(df.Histo1D(("pairs_z", "z distribution;z (mm);Entries", *bins_z), "pairs_z"))

        hists.append(df.Histo2D(("pairs_xy", "xy distribution;x (mm);y (mm)", *(bins_x + bins_y)), "pairs_x", "pairs_y"))
        hists.append(df.Histo2D(("pairs_xz", "xz distribution;x (mm);z (mm)", *(bins_x + bins_z)), "pairs_x", "pairs_z"))
        hists.append(df.Histo2D(("pairs_yz", "yz distribution;y (mm);z (mm)", *(bins_y + bins_z)), "pairs_y", "pairs_z"))

    
    # per production process
    for p in ['bw', 'bh', 'll']:

        df = df.Define(f"pairs_E_{p}", f"pairs_E[sel_{p}]")
        df = df.Define(f"pairs_p_{p}", f"pairs_p[sel_{p}]")
        df = df.Define(f"pairs_pt_{p}", f"pairs_pt[sel_{p}]")
        df = df.Define(f"pairs_theta_{p}", f"pairs_theta[sel_{p}]") 
        df = df.Define(f"pairs_phi_{p}", f"pairs_phi[sel_{p}]")
        df = df.Define(f"pairs_n_{p}", f"size(pairs_E_{p})/2")
        
        hists.append(df.Histo1D((f"pairs_E_{p}", "Energy distribution;E (MeV);Entries", *bins_p), f"pairs_E_{p}"))
        hists.append(df.Histo1D((f"pairs_p_{p}", "Momentum distribution;p (MeV);Entries", *bins_p), f"pairs_p_{p}"))
        hists.append(df.Histo1D((f"pairs_pt_{p}", "Transverse momentum distribution;p_{{T}} (MeV);Entries", *bins_p), f"pairs_pt_{p}"))
        hists.append(df.Histo1D((f"pairs_theta_{p}", "Polar angle;#theta (deg);Entries", *bins_theta), f"pairs_theta_{p}"))
        hists.append(df.Histo1D((f"pairs_phi_{p}", "Azimuthal angle;#phi (deg);Entries", *bins_phi), f"pairs_phi_{p}"))
        hists.append(df.Histo1D((f"pairs_n_{p}", "Number of pairs;Pair multiplicity;Entries", *(10000, 0, 10000)), f"pairs_n_{p}"))

        df = df.Define(f"pairs_E_log10_{p}", f"pairs_E_log10[sel_{p}]")
        df = df.Define(f"pairs_p_log10_{p}", f"pairs_p_log10[sel_{p}]")
        df = df.Define(f"pairs_pt_log10_{p}", f"pairs_pt_log10[sel_{p}]")
        hists.append(df.Histo1D((f"pairs_E_log10_{p}", f"Energy distribution;log_{{10}}(E) (MeV);Entries", *bins_p_log10), f"pairs_E_log10_{p}"))
        hists.append(df.Histo1D((f"pairs_p_log10_{p}", f"Momentum distribution;log_{{10}}(p) (MeV);Entries", *bins_p_log10), f"pairs_p_log10_{p}"))
        hists.append(df.Histo1D((f"pairs_pt_log10_{p}", f"Transverse momentum distribution;log_{{10}}(p_{{T}}) (MeV);Entries", *bins_p_log10), f"pairs_pt_log10_{p}"))

        df = df.Define(f"pairs_theta_rad_log10_{p}", f"pairs_theta_rad_log10[sel_{p}]")
        hists.append(df.Histo2D((f"pairs_theta_pt_{p}", f"#theta-p_{{T}};log_{{10}}(#theta) (rad);log_{{10}}(p_{{T}}) (MeV)", *((500, -4, 1) + bins_p_log10)), f"pairs_theta_rad_log10_{p}", f"pairs_pt_log10_{p}"))

        df = df.Define(f"pairs_x_{p}", f"pairs_x[sel_{p}]") # um
        df = df.Define(f"pairs_y_{p}", f"pairs_y[sel_{p}]") # nm
        df = df.Define(f"pairs_z_{p}", f"pairs_z[sel_{p}]") # mm

        if ptype == 0:
            hists.append(df.Histo1D((f"pairs_x_{p}", "x distribution;x (#mum);Entries", *bins_x), f"pairs_x_{p}"))
            hists.append(df.Histo1D((f"pairs_y_{p}", "y distribution;y (nm);Entries", *bins_y), f"pairs_y_{p}"))
            hists.append(df.Histo1D((f"pairs_z_{p}", "z distribution;z (mm);Entries", *bins_z), f"pairs_z_{p}"))

            hists.append(df.Histo2D((f"pairs_xy_{p}", "xy distribution;x (#mum);y (nm)", *(bins_x + bins_y)), f"pairs_x_{p}", f"pairs_y_{p}"))
            hists.append(df.Histo2D((f"pairs_xz_{p}", "xz distribution;x (#mum);z (mm)", *(bins_x + bins_z)), f"pairs_x_{p}", f"pairs_z_{p}"))
            hists.append(df.Histo2D((f"pairs_yz_{p}", "yz distribution;y (nm);z (mm)", *(bins_y + bins_z)), f"pairs_y_{p}", f"pairs_z_{p}"))

        if ptype == 1:
            hists.append(df.Histo1D((f"pairs_x_{p}", "x distribution;x (mm);Entries", *bins_x), f"pairs_x_{p}"))
            hists.append(df.Histo1D((f"pairs_y_{p}", "y distribution;y (mm);Entries", *bins_y), f"pairs_y_{p}"))
            hists.append(df.Histo1D((f"pairs_z_{p}", "z distribution;z (mm);Entries", *bins_z), f"pairs_z_{p}"))

            hists.append(df.Histo2D((f"pairs_xy_{p}", "xy distribution;x (mm);y (mm)", *(bins_x + bins_y)), f"pairs_x_{p}", f"pairs_y_{p}"))
            hists.append(df.Histo2D((f"pairs_xz_{p}", "xz distribution;x (mm);z (mm)", *(bins_x + bins_z)), f"pairs_x_{p}", f"pairs_z_{p}"))
            hists.append(df.Histo2D((f"pairs_yz_{p}", "yz distribution;y (mm);z (mm)", *(bins_y + bins_z)), f"pairs_y_{p}", f"pairs_z_{p}"))


    # projection on barrel cylinder
    df = df.Define("pairs_x_mm", "get_x(pairs)/1e6") # mm
    df = df.Define("pairs_y_mm", "get_y(pairs)/1e6")
    df = df.Define("pairs_z_mm", "get_z(pairs)/1e6")
    df = df.Define("pairs_px", "get_px(pairs_p4)*1e3") # MeV
    df = df.Define("pairs_py", "get_py(pairs_p4)*1e3")
    df = df.Define("pairs_pz", "get_pz(pairs_p4)*1e3")
    df = df.Define("pairs_q", "get_q(pairs)")
    df = df.Define("pairs_acc", "hitsCylinder(pairs_x_mm, pairs_y_mm, pairs_z_mm, pairs_px, pairs_py, pairs_pz, pairs_q, 2, 13.0, 100)")

    df = df.Define("pairs_pt_log10_acc", "pairs_pt_log10[pairs_acc]")
    df = df.Define("pairs_theta_rad_log10_acc", "pairs_theta_rad_log10[pairs_acc]")
    hists.append(df.Histo2D(("pairs_theta_pt_acc", f"#theta-p_{{T}};log_{{10}}(#theta) (rad);log_{{10}}(p_{{T}}) (MeV)", *((500, -4, 1) + bins_p_log10)), "pairs_theta_rad_log10_acc", "pairs_pt_log10_acc"))


    return hists


def get_hist(hname, hists):
    for h in hists:
        if h.GetName() == hname:
            return h
    return None

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--accelerator", type=str, help="Accelerator config", default="FCCee_Z_GHC_V23")
    parser.add_argument("--parameter_set", type=str, help="Parameter set", default="Z256_2T_grids8")
    parser.add_argument("--campaign", type=str, help="Campaign (as in config.py)", choices=["ipc", "ipc_studies", "ipc_sensitivity_studies"], default="ipc")
    parser.add_argument("--maxfiles", type=int, help="number of files", default=-1)
    args = parser.parse_args()


    #FCCee_Z_GHC_V25p1
    #FCCee_Z_LCC_V105_v2_50ns
    #FCCee_Z_LCC_V105_v2_25ns
    #FCCee_Z_LCC_V105
    #FCCee_Z_GHC_V25p3_4
    #FCCee_TOP_GHC_V25p1
    #FCCee_WW_GHC_V25p1
    #FCCee_ZH_GHC_V25p1
    #FCCee_Z_4IP_GHC_V24p4
    #FCCee_Z_CDR
    #FCCee_Z_GHC_V23


    accelerator = args.accelerator
    parameter_set = args.parameter_set
    maxfiles = args.maxfiles
    campaign = args.campaign
    json_dict = {}
    str_out = ""
    
    output_dir = f"/home/submit/jaeyserm/public_html/fccee/guineapig/validation/{campaign}/{accelerator}/{parameter_set}/"
    os.system(f"mkdir -p {output_dir}")

    cfg = getattr(gpconfig, campaign)
    input_dir = cfg[accelerator][parameter_set]['dir']
    logger.info(f"Run over input directory {input_dir}")
    input_files = glob.glob(f"{input_dir}/*.root")
    logger.info(f"Found {len(input_files)} input files")
    if maxfiles > 0:
        input_files = input_files[:maxfiles]

    


    logger.info(f"Get meta information")
    time0 = time.time()
    reader = functions.GuineaPigReader(input_files)
    hists_meta = meta_info(reader)
    time1 = time.time()
    logger.info(f"Meta information done in {int(time1-time0)} seconds")


    logger.info(f"Start analysis")
    time0 = time.time()
    df = ROOT.RDataFrame("events", input_files)
    
    hists0 = analysis(df, reader, 0)
    hists1 = analysis(df, reader, 1)
    df.Count()
    time1 = time.time()
    logger.info(f"Analysis done in {int(time1-time0)} seconds")


    logger.info(f"Write output")
    output_file = f"{output_dir}/output.root"
    fout = ROOT.TFile(output_file, "RECREATE")
    for h in hists_meta:
        h.Write()

    fout.mkdir("pairs0")
    fout.cd("pairs0")
    for h in hists0:
        h.Scale(1./len(input_files)) # normalize per BX
        h.Write()
    fout.mkdir("pairs")
    fout.cd("pairs")
    for h in hists1:
        h.Scale(1./len(input_files)) # normalize per BX
        h.Write()

    fout.cd()
    p = ROOT.TParameter(int)("nevents", len(input_files))
    p.Write()

    fout.Close()


    ## get some IPC numbers
    json_dict['ipcs']['pairs_n_average'] = get_hist("pairs_n", hists0).GetMean()
    json_dict['ipcs']['pairs_n_average_BW'] = get_hist("pairs_n_bw", hists0).GetMean()
    json_dict['ipcs']['pairs_n_average_BH'] = get_hist("pairs_n_bh", hists0).GetMean()
    json_dict['ipcs']['pairs_n_average_LL'] = get_hist("pairs_n_ll", hists0).GetMean()
    json_dict['ipcs']['pairs_n_95pctquant'] = functions.get_quantile(get_hist("pairs_n", hists0))
    json_dict['ipcs']['pairs_n_95pctquant_BW'] = functions.get_quantile(get_hist("pairs_n_bw", hists0))
    json_dict['ipcs']['pairs_n_95pctquant_BH'] = functions.get_quantile(get_hist("pairs_n_bh", hists0))
    json_dict['ipcs']['pairs_n_95pctquant_LL'] = functions.get_quantile(get_hist("pairs_n_ll", hists0))

    json_dict['nevents'] = len(input_files)
    

    with open(f"{output_dir}/summary.txt", "w") as f:
        f.write(str_out)

    json_str = json.dumps(json_dict, indent=4)
    with open(f"{output_dir}/summary.json", "w") as f:
        f.write(json_str)

    logger.info(f"Output saved to {output_file}")
