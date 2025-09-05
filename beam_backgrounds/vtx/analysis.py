
import sys, os, glob, math
import ROOT
import logging
import time

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

ROOT.EnableImplicitMT(8) # use all cores
#ROOT.DisableImplicitMT() # single core

# load libraries
ROOT.gInterpreter.Declare('#include "functions.h"')


layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

bins_theta = (36, 0, 180)
bins_phi = (72, -180, 180)
bins_p = (10000, 0, 10)
bins_z = (int(max_z/2), -max_z, max_z)
bins_layer = (10, 0, 10)

def analysis(input_files, output_file):

    hists = []
    df = ROOT.RDataFrame("events", input_files)

    # MC particle kinematics
    df = df.Define("mc_pdgid", "get_pdgid(MCParticles, true)")
    df = df.Define("mc_genstatus", "get_generatorStatus(MCParticles)")
    df = df.Define("mc_sel_str", "mc_pdgid == 11 && mc_genstatus == 1")
    df = df.Define("mc_sel", "MCParticles[mc_sel_str]")
    
    df = df.Define("mc_sel_theta_rad_f", "get_theta(mc_sel, false)")
    df = df.Define("mc_sel_theta_rad", "tf_theta(mc_sel_theta_rad_f)") # from 0 to pi/2
    df = df.Define("mc_sel_pt", "get_pt(mc_sel)")
    df = df.Define("mc_sel_theta_rad_log10", "log10(mc_sel_theta_rad)")
    df = df.Define("mc_sel_pt_log10", "log10(mc_sel_pt)")
    hists.append(df.Histo1D(("mc_sel_theta", "", *(500, -5, 5)), "mc_sel_theta_rad"))
    hists.append(df.Histo1D(("mc_sel_pt", "", *bins_p), "mc_sel_pt"))
    hists.append(df.Histo2D(("theta_pt_mc_all", "", *((400, -4, 0) + (400, -4, 0))), "mc_sel_theta_rad_log10", "mc_sel_pt_log10"))


    # Vertex hits (only non-secondary ones)
    df = df.Define("isSecondary", "isProducedBySecondary(VertexBarrelCollection)")
    df = df.Define("isPrimary", "!isSecondary")
    df = df.Define("VertexBarrelCollectionPrimary", "VertexBarrelCollection[isPrimary]")

    df = df.Define("hits_x", "getSimHitPosition_x(VertexBarrelCollectionPrimary)")
    df = df.Define("hits_y", "getSimHitPosition_y(VertexBarrelCollectionPrimary)")
    df = df.Define("hits_z", "getSimHitPosition_z(VertexBarrelCollectionPrimary)")
    df = df.Define("hits_r", "getSimHitPosition_r(VertexBarrelCollectionPrimary)")
    df = df.Define("hits_theta", "getSimHitPosition_theta(VertexBarrelCollectionPrimary, true)")
    df = df.Define("hits_phi", "getSimHitPosition_phi(VertexBarrelCollectionPrimary, true)")
    df = df.Define("hits_layer", "getSimHitLayer(hits_r, {14, 36, 58})")

    df = df.Define("hits_layer0", "hits_layer == 0")
    df = df.Define("hits_z_layer0", "hits_z[hits_layer0]")
    df = df.Define("hits_r_layer0", "hits_r[hits_layer0]")
    df = df.Define("hits_theta_layer0", "hits_theta[hits_layer0]")
    df = df.Define("hits_phi_layer0", "hits_phi[hits_layer0]")


    # MC particles on first layer
    df = df.Define("MCVtxBarrel", "getMCParticle(VertexBarrelCollection, MCParticles, _VertexBarrelCollection_particle.index)")
    df = df.Define("MCVtxBarrelPrimary", "MCVtxBarrel[isPrimary]")
    df = df.Define("MCVtxBarrelPrimaryLayer0", "MCVtxBarrelPrimary[hits_layer0]")
    df = df.Define("MCVtxBarrelPrimaryLayer0_pdgid", "get_pdgid(MCVtxBarrelPrimaryLayer0, true)")
    df = df.Define("MCVtxBarrelPrimaryLayer0_genstatus", "get_generatorStatus(MCVtxBarrelPrimaryLayer0)")
    df = df.Define("MCVtxBarrelPrimaryLayer0_sel_str", "MCVtxBarrelPrimaryLayer0_pdgid == 11 && MCVtxBarrelPrimaryLayer0_genstatus == 1")
    df = df.Define("MCVtxBarrelPrimaryLayer0_sel", "MCVtxBarrelPrimaryLayer0[MCVtxBarrelPrimaryLayer0_sel_str]")

    df = df.Define("MCVtxBarrelPrimaryLayer0_theta", "get_theta(MCVtxBarrelPrimaryLayer0_sel, true)")
    df = df.Define("MCVtxBarrelPrimaryLayer0_theta_rad_f", "get_theta(MCVtxBarrelPrimaryLayer0_sel, false)")
    df = df.Define("MCVtxBarrelPrimaryLayer0_theta_rad", "tf_theta(MCVtxBarrelPrimaryLayer0_theta_rad_f)") # from 0 to pi/2
    df = df.Define("MCVtxBarrelPrimaryLayer0_p", "get_p(MCVtxBarrelPrimaryLayer0_sel)")
    df = df.Define("MCVtxBarrelPrimaryLayer0_pt", "get_pt(MCVtxBarrelPrimaryLayer0_sel)")

    df = df.Define("MCVtxBarrelPrimaryLayer0_theta_rad_log10", "log10(MCVtxBarrelPrimaryLayer0_theta_rad)")
    df = df.Define("MCVtxBarrelPrimaryLayer0_pt_log10", "log10(MCVtxBarrelPrimaryLayer0_pt)")


    # histograms
    hists.append(df.Histo1D(("z", "", *bins_z), "hits_z"))
    hists.append(df.Histo1D(("theta", "", *bins_theta), "hits_theta"))
    hists.append(df.Histo1D(("phi", "", *bins_phi), "hits_phi"))

    hists.append(df.Histo1D(("z_layer0", "", *bins_z), "hits_z_layer0"))
    hists.append(df.Histo1D(("theta_layer0", "", *bins_theta), "hits_theta_layer0"))
    hists.append(df.Histo1D(("phi_layer0", "", *bins_phi), "hits_phi_layer0"))

    hists.append(df.Histo2D(("z_phi", "", *(bins_z + bins_phi)), "hits_z", "hits_phi"))
    hists.append(df.Histo2D(("theta_phi", "", *(bins_theta + bins_phi)), "hits_theta", "hits_phi"))

    hists.append(df.Histo2D(("z_phi_layer0", "", *(bins_z + bins_phi)), "hits_z_layer0", "hits_phi_layer0"))
    hists.append(df.Histo2D(("theta_phi_layer0", "", *(bins_theta + bins_phi)), "hits_theta_layer0", "hits_phi_layer0"))

    hists.append(df.Histo1D(("layer", "", *bins_layer), "hits_layer"))
    hists.append(df.Histo2D(("z_phi_layer0", "", *(bins_z + bins_phi)), "hits_z_layer0", "hits_phi_layer0"))

    hists.append(df.Histo1D(("MCVtxBarrelPrimaryLayer0_theta", "", *bins_theta), "MCVtxBarrelPrimaryLayer0_theta"))
    hists.append(df.Histo1D(("MCVtxBarrelPrimaryLayer0_p", "", *bins_p), "MCVtxBarrelPrimaryLayer0_p"))
    hists.append(df.Histo1D(("MCVtxBarrelPrimaryLayer0_pt", "", *bins_p), "MCVtxBarrelPrimaryLayer0_pt"))
    hists.append(df.Histo2D(("theta_pt_mc_layer0", "", *((400, -4, 0) + (400, -4, 0))), "MCVtxBarrelPrimaryLayer0_theta_rad_log10", "MCVtxBarrelPrimaryLayer0_pt_log10"))

    fout = ROOT.TFile(output_file, "RECREATE")
    for h in hists:
        h.Write()
    fout.Close()



if __name__ == "__main__":

    input_dir, output_file = "/ceph/submit/data/user/k/kudela/beam_backgrounds/aug24_ddsim/CLD_o2_v07_0T/FCCee_Z_4IP_04may23_FCCee_Z256_0T", "output.root"
    input_dir, output_file = "/ceph/submit/data/user/k/kudela/beam_backgrounds/aug24_ddsim/CLD_o2_v07_2T/FCCee_Z_4IP_04may23_FCCee_Z256_0T", "output.root"

    logger.info(f"Start analysis")
    input_files = glob.glob(f"{input_dir}/output_*.root")
    logger.info(f"Found {len(input_files)} input files")

    logger.info(f"Start analysis")
    time0 = time.time()
    analysis(input_files, output_file)
    time1 = time.time()

    logger.info(f"Analysis done in {int(time1-time0)} seconds")
    logger.info(f"Output saved to {output_file}")

