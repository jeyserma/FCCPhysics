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
# this file is for plotting the number of hits in a 2D map of phi and z and purely as a
# function of phi and theta
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_o2_v05/guineaPig_andrea_June2024_v23/"
files = glob.glob(f"{folder}/*.root")


# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# max_z = 96 # IDEA first layer

layer_radii = [14, 36, 58] # CLD approximate layer radii
max_z = 110 # CLD first layer

z_step = 2



if args.calculate:

    hist_theta_phi = ROOT.TH2D("theta_phi", "", int(180/5), 0, 180, int(360/5), 0, 360)
    hist_z_phi = ROOT.TH2D("z_phi", "", int(max_z/2), -max_z, max_z, int(360/5), 0, 360)

    nEvents = 0
    for i,filename in enumerate(files):

        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)

        events = podio_reader.get("events")
        for event in events:
            nEvents += 1
            for hit in event.get("VertexBarrelCollection"):
                radius_idx = functions.radius_idx(hit, layer_radii)
                if radius_idx != 0: # consider only hits on the first layer
                    continue

                if hit.isProducedBySecondary(): # remove mc particle not tracked
                    continue

                ph = functions.phi(hit.getPosition().x, hit.getPosition().y) * (180 / math.pi)
                th = functions.theta(hit.getPosition().x, hit.getPosition().y, hit.getPosition().z) * (180 / math.pi)
                z = hit.getPosition().z

                hist_z_phi.Fill(z, ph)
                hist_theta_phi.Fill(th, ph)

        if i > args.maxFiles:
            break

    # normalize the histograms over the number of events, to get the average number of hits per event
    hist_z_phi.Scale(1./nEvents)
    hist_theta_phi.Scale(1./nEvents)

    fOut = ROOT.TFile("hitmaps.root", "RECREATE")
    hist_z_phi.Write()
    hist_theta_phi.Write()
    fOut.Close()


if args.plots:

    outdir = "/home/submit/jaeyserm/public_html/fccee/beam_backgrounds/vtx/CLD_o2_v05/"
    fIn = ROOT.TFile("hitmaps.root")

    hist_z_phi = fIn.Get('z_phi')
    hist_theta_phi = fIn.Get('theta_phi')


    functions.plot_2dhist(hist_z_phi, f"{outdir}/z_phi.png", f"Hit maps first layer", xMin=-max_z, xMax=max_z, xLabel="z (mm)", yLabel="Azimuthal angle (deg)")
    functions.plot_2dhist(hist_theta_phi, f"{outdir}/theta_phi.png", f"Hit maps first layer", xMin=0, xMax=180, xLabel="Theta (deg)", yLabel="Azimuthal angle (deg)")

    # plot the projections on phi and theta
    hist_theta = hist_theta_phi.ProjectionX()
    hist_phi = hist_theta_phi.ProjectionY()
    functions.plot_hist(hist_theta, f"{outdir}/theta.png", f"Hit maps first layer", xMin=0, xMax=180, xLabel="Theta (deg)")
    functions.plot_hist(hist_phi, f"{outdir}/phi.png", f"Hit maps first layer", xMin=0, xMax=360, xLabel="Azimuthal angle (deg)")