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
parser.add_argument("--maxFiles", type=int, default=1e99, help="Maximum files to run over")
args = parser.parse_args()


##########################################################################################
#  this file is for calculating the average and max occupancies in the first layer
##########################################################################################

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z/"
folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23/"

folder = "/ceph/submit/data/group/fcc/ee/detector/VTXStudiesFullSim/CLD_guineaPig_andrea_June2024_v23_vtx000/"
#folder = "/ceph/submit/data/user/k/kudela/beam_backgrounds/CLD_o2_v05/FCCee_Z_4IP_04may23_FCCee_Z256"
files = glob.glob(f"{folder}/*.root")


# layer_radii = [14, 23, 34.5, 141, 316] # IDEA approximate layer radii
# max_z = 96 # IDEA first layer

layer_radii = [14, 36, 58] # CLD approximate layer radii
layer_radii = [[13.0000, 13.2850], [14.2850, 14.5700]]
layer_radii = [[13.0000, 14.2850], [14.2850, 15.57]]
max_z = 110 # CLD first layer

totalHists = 0

if args.calculate:

    cellCounts = {}
    maxHits = []
    nEvents = 0
    for i,filename in enumerate(files):

        print(f"starting {filename} {i}/{len(files)}")
        podio_reader = root_io.Reader(filename)
        cellCounts_per_event = {}
        events = podio_reader.get("events")
        totHits = 0
        for event in events:
            nEvents += 1
            for hit in event.get("VertexBarrelCollection"):
                
                if hit.isProducedBySecondary(): # remove mc particle not tracked
                    continue
                
                edep = 1000000*hit.getEDep() # convert to keV
                cellID = hit.getCellID()
                
                #print(hit.rho(), cellID, edep)
                radius_idx = functions.radius_idx(hit, layer_radii)
                
                if radius_idx != 0: # consider only hits on the first layer
                    continue

                

                
                
                # the cellID corresponds to a cluster of pixels = module
                # the readout is done per module
                
                if not cellID in cellCounts:
                    cellCounts[cellID] = 0
                if not cellID in cellCounts_per_event:
                    cellCounts_per_event[cellID] = 0
                # increment the hit count for this module
                # strictly speaking all the primaries+secondaries within the dR should be treated as 1 hit
                cellCounts[cellID] += 1
                cellCounts_per_event[cellID] += 1
                zzz = hit.getPosition().z
                #if(abs(zzz) > 108):
                #    print(zzz)
                
                totalHists += 1
                totHits += 1
            
            # max cellcount
            if len(cellCounts_per_event) > 0:
                max_ = max([cellCounts_per_event[c] for c in cellCounts_per_event])
            else:
                max_ = 0
            maxHits.append(max_)
            print(max_)

        #print(cellCounts)
        #print(totHits)

        if i > args.maxFiles:
            break

    print(f"CellIDs: {len(cellCounts)}")

    print(f"totalHists: {totalHists}")
    print(f"totalHists (norm): {totalHists/nEvents}")

    # normalize the hits over the number of events
    print(cellCounts)
    cellCounts_norm = [cellCounts[c]/nEvents for c in cellCounts]
    print(cellCounts_norm)
    
    max_hits_per_events = sum(maxHits) / len(maxHits) # has a distribution! --> plot
    max_hits = max(cellCounts_norm)
    avg_hits = sum(cellCounts_norm)/len(cellCounts_norm)

    safety_factor = 3
    cluster_size = 5
    
    print(f"Max hits (per-event): {max_hits_per_events}")
    print(f"Max hits (modules-average): {max_hits}")
    print(f"avg_hits: {avg_hits}")

    # this corresponds to the occupancy per module
    print(f"Maximum occupancy: {max_hits*safety_factor*cluster_size}")
    print(f"Average occupancy: {avg_hits*safety_factor*cluster_size}")

