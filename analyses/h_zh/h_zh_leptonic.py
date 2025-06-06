
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

processListSignal = {

    f'wzp6_ee_qqH_Hbb_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_Hcc_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_Hss_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_Hgg_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_Haa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_HZa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_HWW_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_HZZ_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_qqH_Hmumu_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_qqH_Htautau_ecm{ecm}':       {'fraction':fraction},
    f'wz3p6_ee_qqH_Hinv_ecm{ecm}':         {'fraction':fraction},

    f'wzp6_ee_ssH_Hbb_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_Hcc_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_Hss_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_Hgg_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_Haa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_HZa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_HWW_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_HZZ_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ssH_Hmumu_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_ssH_Htautau_ecm{ecm}':       {'fraction':fraction},
    f'wz3p6_ee_ssH_Hinv_ecm{ecm}':         {'fraction':fraction},

    f'wzp6_ee_ccH_Hbb_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_Hcc_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_Hss_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_Hgg_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_Haa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_HZa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_HWW_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_HZZ_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_ccH_Hmumu_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_ccH_Htautau_ecm{ecm}':       {'fraction':fraction},
    f'wz3p6_ee_ccH_Hinv_ecm{ecm}':         {'fraction':fraction},


    f'wzp6_ee_bbH_Hbb_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_Hcc_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_Hss_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_Hgg_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_Haa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_HZa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_HWW_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_HZZ_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_bbH_Hmumu_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_bbH_Htautau_ecm{ecm}':       {'fraction':fraction},
    f'wz3p6_ee_bbH_Hinv_ecm{ecm}':         {'fraction':fraction},


    f'wzp6_ee_eeH_Hbb_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_Hcc_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_Hss_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_Hgg_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_Haa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_HZa_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_HWW_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_HZZ_ecm{ecm}':           {'fraction':fraction},
    f'wzp6_ee_eeH_Hmumu_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_eeH_Htautau_ecm{ecm}':       {'fraction':fraction},
    f'wz3p6_ee_eeH_Hinv_ecm{ecm}':         {'fraction':fraction},

    f'wzp6_ee_mumuH_Hbb_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_Hcc_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_Hss_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_Hgg_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_Haa_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_HZa_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_HWW_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_HZZ_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_mumuH_Hmumu_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_mumuH_Htautau_ecm{ecm}':     {'fraction':fraction},
    f'wz3p6_ee_mumuH_Hinv_ecm{ecm}':       {'fraction':fraction},

    f'wzp6_ee_tautauH_Hbb_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_Hcc_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_Hss_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_Hgg_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_Haa_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_HZa_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_HWW_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_HZZ_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_tautauH_Hmumu_ecm{ecm}':     {'fraction':fraction},
    f'wzp6_ee_tautauH_Htautau_ecm{ecm}':   {'fraction':fraction},
    f'wz3p6_ee_tautauH_Hinv_ecm{ecm}':     {'fraction':fraction},


    f'wzp6_ee_nunuH_Hbb_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_Hcc_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_Hss_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_Hgg_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_Haa_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_HZa_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_HWW_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_HZZ_ecm{ecm}':         {'fraction':fraction},
    f'wzp6_ee_nunuH_Hmumu_ecm{ecm}':       {'fraction':fraction},
    f'wzp6_ee_nunuH_Htautau_ecm{ecm}':     {'fraction':fraction},
    f'wz3p6_ee_nunuH_Hinv_ecm{ecm}':       {'fraction':fraction},

}




processList = processListSignal | processListBkg
#processList = processListSignal

inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "utils.h"]


# output directory
outputDir   = f"output/h_zh_leptonic/histmaker/ecm{ecm}/"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 64

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 if ecm == 240 else 3e6


# define histograms
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (200, 0, 200) # 1 GeV bins
bins_recoil = (20000, 0, 200) # 10 MeV bins 
bins_recoil_fine = (500, 100, 150) # 100 MeV bins 
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

ROOT.EnableImplicitMT(nCPUS) # hack to deal correctly with TMVAHelperXGB
tmva_helper_mumu = TMVAHelperXGB(f"FCCPhysics/analyses/h_zh/bdt_leptonic/xgb_bdt_mumu_{ecm}_converted.root", "ZH_Recoil_BDT")
#tmva_helper_mumu = TMVAHelperXGB("FCCPhysics/analyses/h_zh/trainings/xgb_bdt_mumu_new_v0_converted.root", "ZH_Recoil_BDT")
tmva_helper_ee = TMVAHelperXGB(f"FCCPhysics/analyses/h_zh/bdt_leptonic/xgb_bdt_ee_{ecm}_converted.root", "ZH_Recoil_BDT")

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
    df = df.Define("leps_theta", "FCCAnalyses::ReconstructedParticle::get_theta(leps)")
    df = df.Define("leps_phi", "FCCAnalyses::ReconstructedParticle::get_phi(leps)")
    df = df.Define("leps_q", "FCCAnalyses::ReconstructedParticle::get_charge(leps)")
    df = df.Define("leps_no", "FCCAnalyses::ReconstructedParticle::get_n(leps)")
    df = df.Define("leps_iso", "FCCAnalyses::coneIsolation(0.01, 0.5)(leps, ReconstructedParticles)")
    df = df.Define("leps_sel_iso", "FCCAnalyses::sel_iso(0.25)(leps, leps_iso)") # 0.25


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

    # remove H->mumu/ee candidate leptons
    df = df.Define("zbuilder_result_Hll", "FCCAnalyses::resonanceBuilder_mass_recoil(125, 91.2, 0.4, ecm, false)(leps, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define("zll_Hll", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result_Hll[0]}") # the Z
    df = df.Define("zll_Hll_m", "FCCAnalyses::ReconstructedParticle::get_mass(zll_Hll)[0]")
    df = df.Define("zll_leps_Hll", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{zbuilder_result_Hll[1],zbuilder_result_Hll[2]}") # the leptons
    df = df.Define("zll_leps_dummy", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{}") # the leptons
    df = df.Define("leps_to_remove", "return (zll_Hll_m > (125-3) && zll_Hll_m < (125+3)) ? zll_leps_Hll : zll_leps_dummy")
    df = df.Define("leps_good", "FCCAnalyses::ReconstructedParticle::remove(leps, leps_to_remove)")


    # build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
    # technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
    df = df.Define("zbuilder_result", "FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, ecm, false)(leps_good, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
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
    #df = df.Filter("zll_m > 73 && zll_m < 120") # ILC
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut3"))


    #########
    ### CUT 4: Z momentum
    #########
    hists.append(df.Histo1D((f"{ch}_zll_p_nOne", "", *bins_p_mu), "zll_p"))
    if ecm == 240:
        df = df.Filter("zll_p > 20 && zll_p < 70")
    if ecm == 365:
        df = df.Filter("zll_p > 50 && zll_p < 150")
    #df = df.Filter("zll_p > 20 && zll_p < 55")
    #df = df.Filter("zll_p > 10 && zll_p < 70") # ILC
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut4"))


    #########
    ### CUT 5: recoil cut
    #########
    hists.append(df.Histo1D((f"{ch}_zll_recoil_nOne", "", *bins_recoil), "zll_recoil_m"))
    #df = df.Filter("zll_recoil_m < 140 && zll_recoil_m > 120")
    df = df.Filter("zll_recoil_m < 150 && zll_recoil_m > 100") # 120
    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut5"))


    #########
    ### CUT 6: cosThetaMiss, for mass analysis
    #########
    hists.append(df.Histo1D((f"{ch}_cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosTheta_miss"))
    #df = df.Filter("cosTheta_miss < 0.98")
    #df = df.Filter("cosTheta_miss < 0.99975")  # ILC-equivalent cut

    hists.append(df.Histo1D((f"{ch}_cutFlow", "", *bins_count), "cut6"))

    if ch == "mumu":
        df = tmva_helper_mumu.run_inference(df, col_name="mva_score")
    elif ch == "ee":
        df = tmva_helper_ee.run_inference(df, col_name="mva_score")
    hists.append(df.Histo1D((f"{ch}_mva_score", "", *(1000, 0, 1)), "mva_score"))

    ## ILC CUT on visible energy
    df = df.Define("rps_no_muons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, zll_leps)")
    df = df.Define("visibleEnergy", "FCCAnalyses::visibleEnergy(rps_no_muons)")
    hists.append(df.Histo1D((f"{ch}_visibleEnergy", "", *bins_p_mu), "visibleEnergy"))
    #df = df.Filter("visibleEnergy > 10")


    # final histograms
    hists.append(df.Histo1D((f"{ch}_leps_p", "", *bins_p_mu), "leps_p"))
    hists.append(df.Histo1D((f"{ch}_zll_p", "", *bins_p_mu), "zll_p"))
    hists.append(df.Histo1D((f"{ch}_zll_m", "", *bins_m_ll), "zll_m"))
    hists.append(df.Histo1D((f"{ch}_zll_recoil", "", *bins_recoil), "zll_recoil_m"))

    hists.append(df.Histo1D((f"{ch}_cosThetaMiss", "", *bins_cosThetaMiss), "cosTheta_miss"))
    hists.append(df.Histo1D((f"{ch}_acoplanarity", "", *bins_aco), "acoplanarity"))
    hists.append(df.Histo1D((f"{ch}_acolinearity", "", *bins_aco), "acolinearity"))




    ########################
    # Final histograms
    ########################
    hists.append(df.Histo2D((f"{ch}_zll_recoil_m", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m", "zll_category"))
    hists.append(df.Histo1D((f"{ch}_zll_recoil_m_final", "", *(bins_m_fine)), "zll_recoil_m"))

    if ch == "mumu" and ecm == 240: mva_sign = 0.83
    if ch == "mumu" and ecm == 365: mva_sign = 0.66
    if ch == "ee" and ecm == 240: mva_sign = 0.88
    if ch == "ee" and ecm == 365: mva_sign = 0.76
    #mva_sign = 0.83 if ch == "mumu" else 0.88 # max significance
    bins_mva_ = [0, mva_sign, 1]
    bins_mrec_ = list(np.arange(100, 150.5, 0.5))
    bins_mva = array.array('d', bins_mva_)
    bins_mrec = array.array('d', bins_mrec_)
    model = ROOT.RDF.TH2DModel(f"{ch}_recoil_m_mva", "", len(bins_mrec_)-1, bins_mrec, len(bins_mva_)-1, bins_mva)
    hists.append(df.Histo2D(model, "zll_recoil_m", "mva_score"))


    # separate recoil plots
    df_low = df.Filter(f"mva_score[0] < {mva_sign}")
    hists.append(df_low.Histo1D((f"{ch}_zll_recoil_m_mva_low", "", *(bins_recoil_fine)), "zll_recoil_m"))

    df_high = df.Filter(f"mva_score[0] > {mva_sign}")
    hists.append(df_high.Histo1D((f"{ch}_zll_recoil_m_mva_high", "", *(bins_recoil_fine)), "zll_recoil_m"))

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

    hists.append(df.Histo2D((f"{ch}_zll_recoil_m_scaleup", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_scaleup", "zll_category"))
    hists.append(df.Histo2D((f"{ch}_zll_recoil_m_scaledw", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_scaledw", "zll_category"))


    # sqrt uncertainty
    df = df.Define("zll_recoil_sqrtsup", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm + 0.002)(zll)")
    df = df.Define("zll_recoil_sqrtsdw", "FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm - 0.002)(zll)")
    df = df.Define("zll_recoil_m_sqrtsup", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_sqrtsup)[0]")
    df = df.Define("zll_recoil_m_sqrtsdw", "FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_sqrtsdw)[0]")

    hists.append(df.Histo2D(("zll_recoil_m_sqrtsup", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_sqrtsup", "zll_category"))
    hists.append(df.Histo2D(("zll_recoil_m_sqrtsdw", "", *(bins_recoil_fine + bins_cat)), "zll_recoil_m_sqrtsdw", "zll_category"))

    return hists


def build_graph(df, dataset):

    hists = []

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

    if "HZZ" in dataset: # remove H(ZZ) invisible decays from HZZ
        df = df.Define("hzz_invisible", "FCCAnalyses::is_hzz_invisible(Particle, Particle1)")
        df = df.Filter("!hzz_invisible")
    if "p8_ee_WW_ecm" in dataset: # remove muons/electrons from inclusive WW
        df = df.Define("ww_leptonic", "FCCAnalyses::is_ww_leptonic(Particle, Particle1)")
        df = df.Filter("!ww_leptonic")

    build_graph_ll(df, hists, dataset, "mumu")
    build_graph_ll(df, hists, dataset, "ee")
    return hists, weightsum