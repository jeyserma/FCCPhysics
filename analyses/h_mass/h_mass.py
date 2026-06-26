
import ROOT
import array
import numpy as np

ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB

ecm = 240 # 240 365

fraction = 1


processListBkg = {

    f'p8_ee_WW_ecm{ecm}':                  {'fraction':fraction},
    f'p8_ee_WW_mumu_ecm{ecm}':             {'fraction':fraction},
    f'p8_ee_WW_ee_ecm{ecm}':               {'fraction':fraction},
    f'p8_ee_ZZ_ecm{ecm}':                  {'fraction':fraction},
    f'wz3p6_ee_uu_ecm{ecm}':               {'fraction':fraction},
    f'wz3p6_ee_dd_ecm{ecm}':               {'fraction':fraction},
    f'wz3p6_ee_cc_ecm{ecm}':               {'fraction':fraction},
    f'wz3p6_ee_ss_ecm{ecm}':               {'fraction':fraction},
    f'wz3p6_ee_bb_ecm{ecm}':               {'fraction':fraction},
    f'wz3p6_ee_tautau_ecm{ecm}':           {'fraction':fraction},
    f'wz3p6_ee_mumu_ecm{ecm}':             {'fraction':fraction},
    f'wz3p6_ee_ee_Mee_30_150_ecm{ecm}':    {'fraction':fraction},
    f'wz3p6_ee_nunu_ecm{ecm}':             {'fraction':fraction},

    f'wzp6_egamma_eZ_Zmumu_ecm{ecm}':      {'fraction':fraction},
    f'wzp6_gammae_eZ_Zmumu_ecm{ecm}':      {'fraction':fraction},
    f'wzp6_gaga_mumu_60_ecm{ecm}':         {'fraction':fraction},

    f'wzp6_egamma_eZ_Zee_ecm{ecm}':        {'fraction':fraction},
    f'wzp6_gammae_eZ_Zee_ecm{ecm}':        {'fraction':fraction},
    f'wzp6_gaga_ee_60_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_gaga_tautau_60_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_nuenueZ_ecm{ecm}':           {'fraction':fraction},
}

processListSignal_240 = {

    f'wzp6_ee_mumuH_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-100MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-50MeV_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-100MeV_ecm{ecm}': {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-lower-1pc_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-higher-1pc_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-lower-6pc_ecm{ecm}':    {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-higher-6pc_ecm{ecm}':   {'fraction':fraction},

    f'wzp6_ee_eeH_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-100MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-50MeV_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-100MeV_ecm{ecm}': {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_eeH_BES-lower-1pc_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_eeH_BES-higher-1pc_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_eeH_BES-lower-6pc_ecm{ecm}':    {'fraction':fraction},
    f'wzp6_ee_eeH_BES-higher-6pc_ecm{ecm}':   {'fraction':fraction},

    f'wz3p6_ee_eeH_BES-higher-1pc_ecm{ecm}_E2':     {'fraction':fraction},
    f'wz3p6_ee_eeH_BES-lower-1pc_ecm{ecm}_E2':      {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-higher-25MeV_ecm{ecm}_E2':    {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-lower-25MeV_ecm{ecm}_E2':     {'fraction':fraction},
    f'wzp6_ee_eeH_ecm{ecm}_E2':                     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-100MeV_ecm{ecm}_E2':    {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-50MeV_ecm{ecm}_E2':     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-100MeV_ecm{ecm}_E2':     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-50MeV_ecm{ecm}_E2':      {'fraction':fraction},

    f'wzp6_ee_mumuH_BES-higher-1pc_ecm{ecm}_3T':     {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-lower-1pc_ecm{ecm}_3T':      {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-higher-25MeV_ecm{ecm}_3T':    {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-lower-25MeV_ecm{ecm}_3T':     {'fraction':fraction},
    f'wzp6_ee_mumuH_ecm{ecm}_3T':                     {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-100MeV_ecm{ecm}_3T':    {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-50MeV_ecm{ecm}_3T':     {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-100MeV_ecm{ecm}_3T':     {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-50MeV_ecm{ecm}_3T':      {'fraction':fraction},
    
    f'wz3p6_ee_eeH_BES-higher-1pc_ecm{ecm}_3T':     {'fraction':fraction},
    f'wz3p6_ee_eeH_BES-lower-1pc_ecm{ecm}_3T':      {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-higher-25MeV_ecm{ecm}_3T':    {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-lower-25MeV_ecm{ecm}_3T':     {'fraction':fraction},
    f'wzp6_ee_eeH_ecm{ecm}_3T':                     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-100MeV_ecm{ecm}_3T':    {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-50MeV_ecm{ecm}_3T':     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-100MeV_ecm{ecm}_3T':     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-50MeV_ecm{ecm}_3T':      {'fraction':fraction},


    f'wzp6_ee_mumuH_BES-higher-1pc_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-lower-1pc_ecm{ecm}_SiTracking':      {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-higher-25MeV_ecm{ecm}_SiTracking':    {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-lower-25MeV_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wzp6_ee_mumuH_ecm{ecm}_SiTracking':                     {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-100MeV_ecm{ecm}_SiTracking':    {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-50MeV_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-100MeV_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-50MeV_ecm{ecm}_SiTracking':      {'fraction':fraction},
    
    f'wz3p6_ee_eeH_BES-higher-1pc_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wz3p6_ee_eeH_BES-lower-1pc_ecm{ecm}_SiTracking':      {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-higher-25MeV_ecm{ecm}_SiTracking':    {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-lower-25MeV_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wzp6_ee_eeH_ecm{ecm}_SiTracking':                     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-100MeV_ecm{ecm}_SiTracking':    {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-50MeV_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-100MeV_ecm{ecm}_SiTracking':     {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-50MeV_ecm{ecm}_SiTracking':      {'fraction':fraction},

}

processListSignal_365 = {

    f'wzp6_ee_mumuH_ecm{ecm}':                  {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-lower-100MeV_ecm{ecm}':  {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-lower-50MeV_ecm{ecm}':   {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-higher-100MeV_ecm{ecm}': {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-lower-1pc_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_mumuH_BES-higher-1pc_ecm{ecm}':   {'fraction':fraction},

    f'wzp6_ee_eeH_ecm{ecm}':                  {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-lower-100MeV_ecm{ecm}':  {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-lower-50MeV_ecm{ecm}':   {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-higher-100MeV_ecm{ecm}': {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-higher-50MeV_ecm{ecm}':  {'fraction':fraction},
    f'wzp6_ee_eeH_BES-lower-1pc_ecm{ecm}':   {'fraction':fraction},
    f'wzp6_ee_eeH_BES-higher-1pc_ecm{ecm}':   {'fraction':fraction},

}


if ecm == 240:
    processList = processListSignal_240 | processListBkg
if ecm == 365:
    processList = processListSignal_365 | processListBkg

processList = {

    f'wz3p6_ee_mumuH_mH-lower-25MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_mumuH_noBES_ecm{ecm}':                        {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-50MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-50MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-100MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-higher-25MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-100MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wz3p6_ee_eeH_mH-lower-25MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wz3p6_ee_mumuH_mH-higher-25MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_eeH_mH-lower-100MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-higher-50MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_mumuH_mH-lower-50MeV_noBES_ecm{ecm}':                  {'fraction':fraction},
    f'wzp6_ee_eeH_mH-higher-100MeV_noBES_ecm{ecm}':                  {'fraction':fraction},

}




inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "../h_zh/utils.h"]


# output directory
outputDir   = f"output/h_mass/histmaker/ecm{ecm}/"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 64

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 if ecm == 240 else 3.12e6


# define histograms
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (200, 0, 200) # 1 GeV bins
bins_recoil = (20000, 0, 200) # 10 MeV bins 
bins_recoil_fine = (2000, 120, 140) # 10 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (400, 0, 4)
bins_phi = (400, -4, 4)
bins_aco = (400, 0, 4)

bins_count = (50, 0, 50)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 5)
bins_dR = (1000, 0, 10)

bins_cat = (10, 0, 10)
bins_resolution = (10000, 0.95, 1.05)

bins_m_fine = (100, 120, 130) # 100 MeV bins


def build_graph_ll(df, hists, dataset, ch):

    df = df.Alias("Lepton0", "Muon#0.index" if ch=="mumu" else "Electron#0.index")

    # all leptons (bare)
    df = df.Define("leps_all", "FCCAnalyses::ReconstructedParticle::get(Lepton0, ReconstructedParticles)")
    df = df.Define("leps_all_p", "FCCAnalyses::ReconstructedParticle::get_p(leps_all)")
    df = df.Define("leps_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(leps_all)")
    df = df.Define("leps_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(leps_all)")
    df = df.Define("leps_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(leps_all)")
    df = df.Define("leps_all_no", "FCCAnalyses::ReconstructedParticle::get_n(leps_all)")
    df = df.Define("leps_all_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(leps_all, ReconstructedParticles)") 
    df = df.Define("leps_all_p_gen", "FCCAnalyses::gen_p_from_reco(leps_all, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle)")

    # cuts on leptons
    df = df.Define("leps", "FCCAnalyses::ReconstructedParticle::sel_p(20)(leps_all)")
    df = df.Define("leps_p", "FCCAnalyses::ReconstructedParticle::get_p(leps)")
    df = df.Define("leps_pt", "FCCAnalyses::ReconstructedParticle::get_pt(leps)")
    df = df.Define("leps_theta", "FCCAnalyses::ReconstructedParticle::get_theta(leps)")
    df = df.Define("leps_phi", "FCCAnalyses::ReconstructedParticle::get_phi(leps)")
    df = df.Define("leps_q", "FCCAnalyses::ReconstructedParticle::get_charge(leps)")
    df = df.Define("leps_no", "FCCAnalyses::ReconstructedParticle::get_n(leps)")
    df = df.Define("leps_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(leps, ReconstructedParticles)")
    df = df.Define("leps_sel_iso", "FCCAnalyses::sel_iso(0.25)(leps, leps_iso)") # 0.25

    # momentum resolution
    df = df.Define("leps_all_reso_p", "FCCAnalyses::leptonResolution(leps_all, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, 0)")
    df = df.Define("leps_reso_p", "FCCAnalyses::leptonResolution(leps, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, 0)")
    df = df.Define("leps_reso_theta", "FCCAnalyses::leptonResolution(leps, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, 1)")
    df = df.Define("leps_reso_phi", "FCCAnalyses::leptonResolution(leps, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, 2)")


    # baseline selections and histograms
    hists.append(df.Histo1D((f"{ch}_leps_all_p_noSel", "", *bins_p_mu), "leps_all_p"))
    hists.append(df.Histo1D((f"{ch}_leps_all_p_gen_noSel", "", *bins_p_mu), "leps_all_p_gen"))
    hists.append(df.Histo1D((f"{ch}_leps_all_theta_noSel", "", *bins_theta), "leps_all_theta"))
    hists.append(df.Histo1D((f"{ch}_leps_all_phi_noSel", "", *bins_phi), "leps_all_phi"))
    hists.append(df.Histo1D((f"{ch}_leps_all_q_noSel", "", *bins_charge), "leps_all_q"))
    hists.append(df.Histo1D((f"{ch}_leps_all_no_noSel", "", *bins_count), "leps_all_no"))
    hists.append(df.Histo1D((f"{ch}_leps_all_iso_noSel", "", *bins_iso), "leps_all_iso"))

    hists.append(df.Histo1D((f"{ch}_leps_p_noSel", "", *bins_p_mu), "leps_p"))
    hists.append(df.Histo1D((f"{ch}_leps_theta_noSel", "", *bins_theta), "leps_theta"))
    hists.append(df.Histo1D((f"{ch}_leps_phi_noSel", "", *bins_phi), "leps_phi"))
    hists.append(df.Histo1D((f"{ch}_leps_q_noSel", "", *bins_charge), "leps_q"))
    hists.append(df.Histo1D((f"{ch}_leps_no_noSel", "", *bins_count), "leps_no"))
    hists.append(df.Histo1D((f"{ch}_leps_iso_noSel", "", *bins_iso), "leps_iso"))



    #########
    ### CUT 0 all events
    #########
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut0"))


    #########
    ### CUT 1: at least a lepton with at least 1 isolated one
    #########
    df = df.Filter("leps_no >= 1 && leps_sel_iso.size() > 0")
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut1"))


    #########
    ### CUT 2 :at least 2 OS leptons, and build the resonance
    #########
    df = df.Filter("leps_no >= 2 && abs(Sum(leps_q)) < leps_q.size()")
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut2"))


    # build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
    df = df.Define("zbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, ecm, false)(leps, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define("zll", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result[0]}") # the Z
    df = df.Define("zll_leps", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result[1],zbuilder_result[2]}") # the leptons
    df = df.Define("zll_m", "FCCAnalyses::ReconstructedParticle::get_mass(zll)[0]")
    df = df.Define("zll_p", "FCCAnalyses::ReconstructedParticle::get_p(zll)[0]")
    df = df.Define("zll_theta", "FCCAnalyses::ReconstructedParticle::get_theta(zll)[0]")
    df = df.Define("zll_recoil", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zll)")
    df = df.Define("zll_recoil_m", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil)[0]")
    df = df.Define("zll_category", "FCCAnalyses::polarAngleCategorization(0.8, 2.34)(zll_leps)")

    df = df.Define("zll_leps_p", "FCCAnalyses::ReconstructedParticle::get_p(zll_leps)")
    df = df.Define("zll_leps_theta", "FCCAnalyses::ReconstructedParticle::get_theta(zll_leps)")
    df = df.Define("leading_p_idx", "(zll_leps_p[0] > zll_leps_p[1]) ? 0 : 1")
    df = df.Define("subleading_p_idx", "(zll_leps_p[0] > zll_leps_p[1]) ? 1 : 0")
    df = df.Define("zll_leading_p", "zll_leps_p[leading_p_idx]")
    df = df.Define("zll_subleading_p", "zll_leps_p[subleading_p_idx]")
    df = df.Define("zll_leading_theta", "zll_leps_theta[leading_p_idx]")
    df = df.Define("zll_subleading_theta", "zll_leps_theta[subleading_p_idx]")

    df = df.Define("missingEnergy", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("cosTheta_miss", "FCCAnalyses::get_cosTheta_miss(missingEnergy)")

    df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(zll_leps)")
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(zll_leps)")




    #########
    ### CUT 3: Z mass window
    #########
    hists.append(df.Histo1D((f"{ch}_zll_m_nOne", "", *bins_m_ll), "zll_m"))
    df = df.Filter("zll_m > 86 && zll_m < 96")
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut3"))


    #########
    ### CUT 4: Z momentum
    #########
    hists.append(df.Histo1D((f"{ch}_zll_p_nOne", "", *bins_p_mu), "zll_p"))
    if ecm == 240:
        df = df.Filter("zll_p > 20 && zll_p < 70")
    if ecm == 365:
        df = df.Filter("zll_p > 50 && zll_p < 150")
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut4"))


    #########
    ### CUT 5: recoil cut
    #########
    hists.append(df.Histo1D((f"{ch}_zll_recoil_nOne", "", *bins_recoil), "zll_recoil_m"))
    df = df.Filter("zll_recoil_m < 140 && zll_recoil_m > 120")
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut5"))


    #########
    ### CUT 6: cosThetaMiss
    #########
    hists.append(df.Histo1D((f"{ch}_cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosTheta_miss"))
    df = df.Filter("cosTheta_miss < 0.98")
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut6"))


    ########################
    # Final histograms
    ########################
    hists.append(df.Histo1D((f"{ch}_leps_p", "", *bins_p_mu), "leps_p"))
    hists.append(df.Histo1D((f"{ch}_leps_pt", "", *bins_p_mu), "leps_pt"))
    hists.append(df.Histo1D((f"{ch}_zll_p", "", *bins_p_mu), "zll_p"))
    hists.append(df.Histo1D((f"{ch}_zll_m", "", *bins_m_ll), "zll_m"))
    hists.append(df.Histo1D((f"{ch}_zll_recoil_m", "", *bins_recoil), "zll_recoil_m"))

    hists.append(df.Histo1D((f"{ch}_cosThetaMiss", "", *bins_cosThetaMiss), "cosTheta_miss"))
    hists.append(df.Histo1D((f"{ch}_acoplanarity", "", *bins_aco), "acoplanarity"))
    hists.append(df.Histo1D((f"{ch}_acolinearity", "", *bins_aco), "acolinearity"))

    hists.append(df.Histo2D((f"{ch}_zll_recoil_m_cat", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m", "zll_category"))
    hists.append(df.Histo2D((f"{ch}_leps_reso_p_cat", "", *(bins_resolution + bins_cat)), "leps_reso_p", "zll_category"))


    ########################
    # Systematics
    ########################

    # muon momentum scale
    df = df.Define("leps_scaleup", "FCCAnalyses::lepton_momentum_scale(1e-5)(leps)")
    df = df.Define("leps_scaledw", "FCCAnalyses::lepton_momentum_scale(-1e-5)(leps)")

    df = df.Define("zbuilder_result_scaleup", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, ecm, false)(leps_scaleup, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define("zbuilder_result_scaledw", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, ecm, false)(leps_scaledw, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define("zll_scaleup", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result_scaleup[0]}")
    df = df.Define("zll_scaledw", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result_scaledw[0]}")
    df = df.Define("zll_recoil_scaleup", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zll_scaleup)")
    df = df.Define("zll_recoil_scaledw", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zll_scaledw)")
    df = df.Define("zll_recoil_m_scaleup", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_scaleup)[0]")
    df = df.Define("zll_recoil_m_scaledw", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_scaledw)[0]")

    hists.append(df.Histo2D((f"{ch}_zll_recoil_m_cat_scaleup", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_scaleup", "zll_category"))
    hists.append(df.Histo2D((f"{ch}_zll_recoil_m_cat_scaledw", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_scaledw", "zll_category"))


    # sqrt uncertainty
    df = df.Define("zll_recoil_sqrtsup", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm + 0.002)(zll)")
    df = df.Define("zll_recoil_sqrtsdw", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm - 0.002)(zll)")
    df = df.Define("zll_recoil_m_sqrtsup", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_sqrtsup)[0]")
    df = df.Define("zll_recoil_m_sqrtsdw", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_sqrtsdw)[0]")

    hists.append(df.Histo2D((f"{ch}_zll_recoil_m_cat_sqrtsup", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_sqrtsup", "zll_category"))
    hists.append(df.Histo2D((f"{ch}_zll_recoil_m_cat_sqrtsdw", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_sqrtsdw", "zll_category"))


    return hists




def build_graph(df, dataset):

    hists = []

    '''
    # full/fast sim helpers
    if 'FullSim' in dataset.name:
        # in FullSim, both the reco and gen particles are produced with a crossing angle
        df = df.Define("ReconstructedParticles", "FCCAnalyses::unBoostCrossingAngle(PandoraPFOs, -0.015)")
        df = df.Define("Particle", "FCCAnalyses::unBoostCrossingAngle(MCParticles, -0.015)")

        df = df.Define("leps_all", f"FCCAnalyses::sel_type({'13' if args.flavor == 'mumu' else '11'}, ReconstructedParticles)")
        df = df.Define("photons", "FCCAnalyses::sel_type(22, ReconstructedParticles)")

        df = df.Alias("Particle0", "_MCParticles_parents.index")
        df = df.Alias("Particle1", "_MCParticles_daughters.index")
        df = df.Alias("MCRecoAssociations0", "_RecoMCTruthLink_rec.index")
        df = df.Alias("MCRecoAssociations1", "_RecoMCTruthLink_sim.index")
    else:
        df = df.Alias("Particle0", "Particle#0.index")
        df = df.Alias("Particle1", "Particle#1.index")
        df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
        df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

        df = df.Alias("Lepton0", "Muon#0.index" if args.flavor == 'mumu' else "Electron#0.index")
        df = df.Define("leps_all", "FCCAnalyses::ReconstructedParticle::get(Lepton0, ReconstructedParticles)")

        df = df.Alias("Photon0", "Photon#0.index")
        df = df.Define("photons", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")
    '''

    df = df.Define("ecm", "240" if ecm == 240 else "365")
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    df = df.Define("cut0", "0")
    df = df.Define("cut1", "1")
    df = df.Define("cut2", "2")
    df = df.Define("cut3", "3")
    df = df.Define("cut4", "4")
    df = df.Define("cut5", "5")
    df = df.Define("cut6", "6")
    df = df.Define("cut7", "7")
    df = df.Define("cut8", "8")
    df = df.Define("cut9", "9")
    df = df.Define("cut10", "10")
    df = df.Define("cut11", "11")
    df = df.Define("cut12", "12")
    df = df.Define("cut13", "13")
    df = df.Define("cut14", "14")

    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

    #if "HZZ" in dataset: # remove H(ZZ) invisible decays from HZZ
    #    df = df.Define("hzz_invisible", "FCCAnalyses::is_hzz_invisible(Particle, Particle1)")
    #    df = df.Filter("!hzz_invisible")
    #if "p8_ee_WW_ecm" in dataset: # remove muons/electrons from inclusive WW
    #    df = df.Define("ww_leptonic", "FCCAnalyses::is_ww_leptonic(Particle, Particle1)")
    #    df = df.Filter("!ww_leptonic")

    build_graph_ll(df, hists, dataset, "mumu")
    build_graph_ll(df, hists, dataset, "ee")
    return hists, weightsum