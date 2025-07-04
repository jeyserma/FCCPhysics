
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

ROOT.EnableImplicitMT() # use all cores
#ROOT.DisableImplicitMT() # single core

# load libraries
ROOT.gInterpreter.Declare('#include "functions.h"')
ROOT.gInterpreter.Declare('#include "digitizer_simHit.h"')


layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

bins_theta = (36, 0, 180)
bins_phi = (72, -180, 180)
bins_z = (int(max_z/2), -max_z, max_z)
bins_layer = (10, 0, 10)

bins_occupancy = (100, 0, 100)
bins_occupancy_avg = (1000, 0, 100)



def analysis(input_files, output_file):

    hists = []
    df = ROOT.RDataFrame("events", input_files)


    df = df.Define("hits_x", "getSimHitPosition_x(VertexBarrelCollection)")
    df = df.Define("hits_y", "getSimHitPosition_y(VertexBarrelCollection)")
    df = df.Define("hits_z", "getSimHitPosition_z(VertexBarrelCollection)")
    df = df.Define("hits_r", "getSimHitPosition_r(VertexBarrelCollection)")
    df = df.Define("hits_theta", "getSimHitPosition_theta(VertexBarrelCollection, true)")
    df = df.Define("hits_phi", "getSimHitPosition_phi(VertexBarrelCollection, true)")
    df = df.Define("hits_edep", "getEnergyDeposition(VertexBarrelCollection)")
    df = df.Define("hits_cellID", "getCellID(VertexBarrelCollection)")
    
    df = df.Define("hits_isSecondary", "isProducedBySecondary(VertexBarrelCollection)")
    df = df.Define("hits_layer", "getSimHitLayer(hits_r, {14, 36, 58})")


    df = df.Define("hits_selection", "hits_layer == 0 && !hits_isSecondary") # layer 0, not secondary
    df = df.Define("hits_sel_x", "hits_x[hits_selection]")
    df = df.Define("hits_sel_y", "hits_y[hits_selection]")
    df = df.Define("hits_sel_z", "hits_z[hits_selection]")
    df = df.Define("hits_sel_edep", "hits_edep[hits_selection]")
    df = df.Define("hits_sel_cellID", "hits_cellID[hits_selection]")

    hists.append(df.Histo1D(("z", "", *bins_z), "hits_z"))
    hists.append(df.Histo1D(("theta", "", *bins_theta), "hits_theta"))
    hists.append(df.Histo1D(("phi", "", *bins_phi), "hits_phi"))


    df = df.Define("geo_layer1", "BarrelGeometry(14.f, 109.f, 25.f, 25.f)") # radius (mm), z-extent (mm), pixel pitch x (um), pixel pitch y (um)
    df = df.Define("digis_layer1", "DigitizerSimHitBarrel(geo_layer1, hits_sel_x, hits_sel_y, hits_sel_z, hits_sel_edep, -1)")
    df = df.Define("digis_layer1_occ", "getOccupancyBarrel(geo_layer1, digis_layer1, 4, 4)")
    df = df.Define("digis_occupancy_max", "digis_layer1_occ[0]")
    df = df.Define("digis_occupancy_avg", "digis_layer1_occ[1]")

    hists.append(df.Histo1D(("digis_occupancy_avg", "", *bins_occupancy_avg), "digis_occupancy_avg"))
    hists.append(df.Histo1D(("digis_occupancy_max", "", *bins_occupancy), "digis_occupancy_max"))

    fout = ROOT.TFile(output_file, "RECREATE")
    for h in hists:
        h.Write()
    fout.Close()



if __name__ == "__main__":

    input_dir, output_file = "/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z256", "digitizer.root"
    input_dir, output_file = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23/", "digitizer.root"
    

    logger.info(f"Start analysis")
    input_files = glob.glob(f"{input_dir}/*.root")
    #input_files = [input_files[0]]
    logger.info(f"Found {len(input_files)} input files")

    logger.info(f"Start analysis")
    time0 = time.time()
    analysis(input_files, output_file)
    time1 = time.time()

    logger.info(f"Analysis done in {int(time1-time0)} seconds")
    logger.info(f"Output saved to {output_file}")

