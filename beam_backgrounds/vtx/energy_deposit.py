from podio import root_io
import glob
import pickle
import argparse
import functions
import math
import ROOT
ROOT.gROOT.SetBatch(True)


parser = argparse.ArgumentParser()
parser.add_argument('--calculate', help="Calculate", action='store_true')
parser.add_argument('--plots', help="Plot the energy deposits", action='store_true')
parser.add_argument("--maxFiles", type=int, default=1e99, help="Maximum files to run over")
args = parser.parse_args()

##########################################################################################
# this file is for plotting the number of hits over energy deposited in each layer
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_o2_v05/guineaPig_andrea_June2024_v23/"
files = glob.glob(f"{folder}/*.root")


# layer_radii = [14, 23, 34.5, 141, 316] # IDEA vertex approximate layer radii
layer_radii = [14, 36, 58] # CLD vertex approximate layer radii
nLayers = len(layer_radii)

if args.calculate:

    hist_energy = ROOT.TH2D("energy", "", nLayers, 0, nLayers, 200, 0, 200)
    hist_dedx = ROOT.TH2D("dedx", "", nLayers, 0, nLayers, 100, 0, 1000)
    hist_dedx_noSecondary = ROOT.TH2D("dedx_noSecondary", "", nLayers, 0, nLayers, 100, 0, 1000)

    for iF,filename in enumerate(files):

        print(f"starting {filename} {iF}/{len(files)}")
        podio_reader = root_io.Reader(filename)

        events = podio_reader.get("events")
        for event in events:
            for hit in event.get("VertexBarrelCollection"):
                radius_idx = functions.radius_idx(hit, layer_radii)

                edep = 1000000*hit.getEDep() # convert to keV
                path_length = hit.getPathLength() # mm
                mc = hit.getMCParticle()

                if mc.getGeneratorStatus() != 1:
                    continue # particles not input into geant

                hist_energy.Fill(radius_idx, edep)
                hist_dedx.Fill(radius_idx, edep/path_length)
                if not hit.isProducedBySecondary(): # mc particle not tracked
                    hist_dedx_noSecondary.Fill(radius_idx, edep/path_length)

        if iF > args.maxFiles:
            break


    fOut = ROOT.TFile("energy_deposit.root", "RECREATE")
    hist_energy.Write()
    hist_dedx.Write()
    hist_dedx_noSecondary.Write()
    fOut.Close()

if args.plots:

    outdir = "/home/submit/jaeyserm/public_html/fccee/beam_backgrounds/vtx/CLD_o2_v05/"
    fIn = ROOT.TFile("energy_deposit.root")

    hist_energy = fIn.Get('energy')
    hist_dedx = fIn.Get('dedx')
    hist_dedx_noSecondary = fIn.Get('dedx_noSecondary')

    # plot them for all layers
    for i,r in enumerate(layer_radii):
        l = i+1 # layer
        # we slice the histograms for each layer
        hist_energy_layer = hist_energy.ProjectionY(f"energy_layer{l}", l, l)
        hist_dedx_layer = hist_dedx.ProjectionY(f"dedx_layer{l}", l, l)
        hist_dedx_noSecondary_layer = hist_dedx_noSecondary.ProjectionY(f"dedx_noSecondary_layer{l}", l, l)
    
        functions.plot_hist(hist_energy_layer, f"{outdir}/energy_layer{l}.png", f"Energy layer {l}", xMin=0, xMax=100, xLabel="Energy (keV)")
        functions.plot_hist(hist_dedx_layer, f"{outdir}/dedx_layer{l}.png", f"Energy layer {l}", xMin=0, xMax=1000, xLabel="dE/dx (keV/mm)")
        functions.plot_hist(hist_dedx_noSecondary_layer, f"{outdir}/dedx_noSecondary_layer{l}.png", f"Energy layer {l}", xMin=0, xMax=1000, xLabel="dE/dx (keV/mm)")