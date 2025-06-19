

import sys, os, glob, math
import ROOT
import logging

logging.basicConfig(format='%(levelname)s: %(message)s')
logger = logging.getLogger("fcclogger")
logger.setLevel(logging.INFO)


ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

ROOT.EnableImplicitMT() # use all cores
#ROOT.DisableImplicitMT() # single core

# load libraries
ROOT.gSystem.Load("libFCCAnalyses")
fcc_loaded = ROOT.dummyLoader()
ROOT.gInterpreter.Declare("using namespace FCCAnalyses;")
ROOT.gInterpreter.Declare('#include "functions.h"')



bins_p = (250, 0, 250)
bins_res = (1000, -0.1, 0.1)

def analysis(input_files, output_file):

    df = ROOT.RDataFrame("events", input_files)

    df = df.Alias("MCRecoAssociations0", "_MCRecoAssociations_rec.index")
    df = df.Alias("MCRecoAssociations1", "_MCRecoAssociations_sim.index")


    # select charged reco/gen particles
    df = df.Define("reco_q", "FCCAnalyses::ReconstructedParticle::get_charge(ReconstructedParticles)")
    df = df.Define("gen_q", "FCCAnalyses::MCParticle::get_charge(Particle)")
    df = df.Define("reco_q_sel", "reco_q != 0")
    df = df.Define("gen_q_sel", "gen_q != 0")
    df = df.Define("reco_charged", "ReconstructedParticles[reco_q_sel]")
    df = df.Define("gen_charged", "Particle[gen_q_sel]")
    df = df.Define("reco_charged_p", "FCCAnalyses::ReconstructedParticle::get_p(reco_charged)")
    df = df.Define("reco_charged_n", "FCCAnalyses::ReconstructedParticle::get_n(reco_charged)")
    df = df.Define("gen_charged_p", "FCCAnalyses::MCParticle::get_p(gen_charged)")
    df = df.Define("gen_charged_n", "FCCAnalyses::MCParticle::get_n(gen_charged)")

    # hadronic resolution = sum of energy of all particles
    df = df.Define("reco_e", "FCCAnalyses::ReconstructedParticle::get_e(ReconstructedParticles)")
    df = df.Define("reco_e_tot", "Sum(reco_e)")
    df = df.Define("qq_res", "(reco_e_tot - 91.2)/ 91.2")
    qq_res = df.Histo1D(("qq_res", "", *bins_res), "qq_res")
    reco_e_tot = df.Histo1D(("reco_e_tot", "", *bins_p), "reco_e_tot")

    fout = ROOT.TFile(output_file, "RECREATE")
    qq_res.Write()
    reco_e_tot.Write()
    fout.Close()



if __name__ == "__main__":

    input_files, output_file = ["samples/IDEA_2T_Zuu_ecm91p2.root"], "output/IDEA_2T_Zuu_ecm91p2.root"
    #input_files, output_file = ["samples/IDEA_3T_Zuu_ecm91p2.root"], "output/IDEA_3T_Zuu_ecm91p2.root"
    #input_files, output_file = ["samples/CLD_2T_Zuu_ecm240.root"], "output/CLD_2T_Zmumu_ecm240.root"

    logger.info(f"Start analysis")
    analysis(input_files, output_file)
    logger.info(f"Done! Output saved to {output_file}")

