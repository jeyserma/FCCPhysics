
import ROOT
import array
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB

ecm = 240 # 240 365
treemaker = False

# list of all processes
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


processListSignall = {

    f'wzp6_ee_qqH_Hss_ecm{ecm}':           {'fraction':fraction},

}

fractionSig, fractionBkg = 0.5, 0.1
if ecm == 365:
    fractionSig, fractionBkg = 0.5, 1

processListTraining = {
    f'wzp6_ee_qqH_Hbb_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_Hcc_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_Hss_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_Hgg_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_Haa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_HZa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_HWW_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_HZZ_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_qqH_Hmumu_ecm{ecm}':         {'fraction':fractionSig},
    f'wzp6_ee_qqH_Htautau_ecm{ecm}':       {'fraction':fractionSig},
    f'wz3p6_ee_qqH_Hinv_ecm{ecm}':         {'fraction':fractionSig},

    f'wzp6_ee_ssH_Hbb_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_Hcc_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_Hss_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_Hgg_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_Haa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_HZa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_HWW_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_HZZ_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ssH_Hmumu_ecm{ecm}':         {'fraction':fractionSig},
    f'wzp6_ee_ssH_Htautau_ecm{ecm}':       {'fraction':fractionSig},
    f'wz3p6_ee_ssH_Hinv_ecm{ecm}':         {'fraction':fractionSig},

    f'wzp6_ee_ccH_Hbb_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_Hcc_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_Hss_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_Hgg_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_Haa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_HZa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_HWW_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_HZZ_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_ccH_Hmumu_ecm{ecm}':         {'fraction':fractionSig},
    f'wzp6_ee_ccH_Htautau_ecm{ecm}':       {'fraction':fractionSig},
    f'wz3p6_ee_ccH_Hinv_ecm{ecm}':         {'fraction':fractionSig},


    f'wzp6_ee_bbH_Hbb_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_Hcc_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_Hss_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_Hgg_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_Haa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_HZa_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_HWW_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_HZZ_ecm{ecm}':           {'fraction':fractionSig},
    f'wzp6_ee_bbH_Hmumu_ecm{ecm}':         {'fraction':fractionSig},
    f'wzp6_ee_bbH_Htautau_ecm{ecm}':       {'fraction':fractionSig},
    f'wz3p6_ee_bbH_Hinv_ecm{ecm}':         {'fraction':fractionSig},
    
    f'p8_ee_WW_ecm{ecm}':                  {'fraction':fractionBkg},
    f'wz3p6_ee_uu_ecm{ecm}':               {'fraction':fractionBkg},
    f'wz3p6_ee_dd_ecm{ecm}':               {'fraction':fractionBkg},
    f'wz3p6_ee_cc_ecm{ecm}':               {'fraction':fractionBkg},
    f'wz3p6_ee_ss_ecm{ecm}':               {'fraction':fractionBkg},
    f'wz3p6_ee_bb_ecm{ecm}':               {'fraction':fractionBkg},
}

processListTrainingLL = {

    f'wz3p6_ee_uu_ecm{ecm}':               {'fraction':fractionBkg},
    f'wz3p6_ee_dd_ecm{ecm}':               {'fraction':fractionBkg},
    f'wz3p6_ee_cc_ecm{ecm}':               {'fraction':fractionBkg},
}

processList = processListSignal | processListBkg
#processList = processListSignal

if treemaker:
    processList = processListTraining


#processList = {

#    f'p8_ee_WW_ecm{ecm}':                  {'fraction':fraction},
#}


inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "utils.h", "utils_hadronic.h"]


# output directory
if treemaker:
    outputDir   = f"/ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/treemaker/ecm{ecm}"
    outputDir   = f"/ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/treemaker/ecm{ecm}_Sept2026"
else:
    outputDir   = f"output/h_zh_hadronic/histmaker/ecm{ecm}"
    outputDir   = f"output/h_zh_hadronic/histmaker/ecm{ecm}_Sept2026_new_jer5pct"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 64

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 if ecm == 240 else 3e6


# define histograms
bins_recoil = (20000, 0, 200) # 10 MeV bins 
bins_recoil_fine = (20000, 120, 140) # 1 MeV bins 
bins_cosThetaMiss = (10000, 0, 1)

bins_theta = (500, 0, 5)
bins_phi = (500, -5, 5)
bins_aco = (1000, -10, 10)

bins_count = (50, 0, 50)
bins_pdgid = (60, -30, 30)
bins_charge = (10, -5, 5)
bins_iso = (500, 0, 5)
bins_dR = (1000, 0, 10)

bins_cat = (10, 0, 10)
bins_resolution = (10000, 0.95, 1.05)

bins_m_fine = (100, 120, 130) # 100 MeV bins
bins_cos_abs = (100, 0, 1)

ROOT.EnableImplicitMT(nCPUS) # hack to deal correctly with TMVAHelperXGB  # bdt_model_new_0p1_WZqq bdt_model_new_0p1 bdt_model_0p1_inv bdt_model_0p1
## original: output/h_zh_hadronic/training/bdt_model_WW_Zg_thrust_reduced.root
## before: /ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/bdt_model_final.root
tmva_helper = TMVAHelperXGB(f"/ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/bdt_model_ecm{ecm}.root", "bdt_model") # read the XGBoost training
#tmva_helper = TMVAHelperXGB("output/h_zh_hadronic/training/bdt_model_WW_Zg_thrust_reduced_withInv.root", "bdt_model") # read the XGBoost training





def exclusive_clustering(df, njets):
    if njets==0: # inclusive
        df = df.Define("clustered_jets_N0", "JetClustering::clustering_kt(0.6, 0, 5, 1, 0)(pseudo_jets)")
    else:
        df = df.Define(f"clustered_jets_N{njets}", f"JetClustering::clustering_ee_kt(2, {njets}, 0, 10)(pseudo_jets)")
        #df = df.Define(f"clustered_jets_N{njets}", f"JetClustering::clustering_ee_genkt(0.5, 2, {njets}, 0, 10, 0)(pseudo_jets)")
    df = df.Define(f"jets_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets_N{njets})")
    df = df.Define(f"njets_N{njets}", f"jets_N{njets}.size()")
    df = df.Define(f"jetconstituents_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets_N{njets})")
    df = df.Define(f"jets_e_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_e(jets_N{njets})")
    df = df.Define(f"jets_px_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_px(jets_N{njets})")
    df = df.Define(f"jets_py_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_py(jets_N{njets})")
    df = df.Define(f"jets_pz_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_pz(jets_N{njets})")
    df = df.Define(f"jets_m_N{njets}", f"FCCAnalyses::JetClusteringUtils::get_m(jets_N{njets})")
    #df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets})")
    #df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp_smeared(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets}, 1.0, 0.04, 1.0, rdfentry_)")
    #df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp_smeared(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets}, 1.01, 0.04, 1.01, rdfentry_)") # JER 1%
    df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp_smeared(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets}, 1.01, 0.04, 1.05, rdfentry_)") # JER 5%
    #df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp_smeared(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets}, 1.01, 0.04, 1.10, rdfentry_)") # JES 1%
    #df = df.Define(f"jets_rp_N{njets}", f"FCCAnalyses::jets2rp_smeared(jets_px_N{njets}, jets_py_N{njets}, jets_pz_N{njets}, jets_e_N{njets}, jets_m_N{njets}, 1.02, 0.04, 1.10, rdfentry_)") # JES 2%
    df = df.Define(f"jets_rp_cand_N{njets}", f"FCCAnalyses::select_jets(jets_rp_N{njets}, jetconstituents_N{njets}, {njets}, ReconstructedParticles)") # reduces potentially the jet multiplicity
    df = df.Define(f"njets_cand_N{njets}", f"jets_rp_cand_N{njets}.size()")
    return df

def define_variables(df, njets):
    #if njets == 0:
    #df = df.Filter(f"if(njets_cand_N{njets} < 2) std::cout << \"FILTER={njets} \" << njets_cand_N{njets}  <<  std::endl; return njets_cand_N{njets} >= 2") #### tricky need to fix it

    df = df.Define(f"zbuilder_result_N{njets}", f"FCCAnalyses::resonanceBuilder_mass_recoil_hadronic(91.2, 125, 0.0, ecm)(jets_rp_cand_N{njets})")
    df = df.Define(f"zqq_N{njets}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_N{njets}[0]}}") # the Z
    df = df.Define(f"zqq_jets_N{njets}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_N{njets}[1],zbuilder_result_N{njets}[2]}}") # the leptons
    df = df.Define(f"zqq_m_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_mass(zqq_N{njets})[0]")
    df = df.Define(f"zqq_p_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_p(zqq_N{njets})[0]")
    df = df.Define(f"zqq_recoil_N{njets}", f"FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zqq_N{njets})")
    df = df.Define(f"zqq_recoil_m_N{njets}", f"FCCAnalyses::ReconstructedParticle::get_mass(zqq_recoil_N{njets})[0]")


    return df


def veto_selection_leptonic(df, flavor):

    if flavor == "muon":
        df = df.Define(f"leps_{flavor}", "FCCAnalyses::ReconstructedParticle::sel_p(20)(muons_all)")
    else:
        df = df.Define(f"leps_{flavor}", "FCCAnalyses::ReconstructedParticle::sel_p(20)(electrons_all)")

    df = df.Define(f"leps_{flavor}_q", f"FCCAnalyses::ReconstructedParticle::get_charge(leps_{flavor})")
    df = df.Define(f"leps_{flavor}_no", f"FCCAnalyses::ReconstructedParticle::get_n(leps_{flavor})")
    df = df.Define(f"leps_{flavor}_iso", f"FCCAnalyses::coneIsolation(0.01, 0.5)(leps_{flavor}, ReconstructedParticles)")
    df = df.Define(f"leps_{flavor}_sel_iso", f"FCCAnalyses::sel_iso(0.25)(leps_{flavor}, leps_{flavor}_iso)") # 0.25


    # do not veto the H-ll candidates, it suppresses HZZ
    # tighter selection in leptonic channel, can safely remove this here
    df = df.Define(f"zbuilder_result_H{flavor}", f"FCCAnalyses::resonanceBuilder_mass_recoil(125, 91.2, 0.4, ecm, false)(leps_{flavor}, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define(f"zll_H{flavor}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_H{flavor}[0]}}") # the Z
    df = df.Define(f"zll_H{flavor}_m", f"FCCAnalyses::ReconstructedParticle::get_mass(zll_H{flavor})[0]")
    df = df.Define(f"zll_leps_H{flavor}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_H{flavor}[1],zbuilder_result_H{flavor}[2]}}") # the leptons
    df = df.Define(f"zll_leps_dummy_{flavor}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{}}") # the leptons
    df = df.Define(f"leps_to_remove_{flavor}", f"return (zll_H{flavor}_m > (125-3) && zll_H{flavor}_m < (125+3)) ? zll_leps_H{flavor} : zll_leps_dummy_{flavor}")
    df = df.Define(f"leps_good_{flavor}", f"FCCAnalyses::ReconstructedParticle::remove(leps_{flavor}, leps_to_remove_{flavor})")


    df = df.Define(f"zbuilder_result_{flavor}", f"FCCAnalyses::resonanceBuilder_mass_recoil(91.2, 125, 0.4, ecm, false)(leps_{flavor}, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")
    df = df.Define(f"zll_{flavor}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_{flavor}[0]}}") # the Z
    df = df.Define(f"zll_leps_{flavor}", f"ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>{{zbuilder_result_{flavor}[1],zbuilder_result_{flavor}[2]}}") # the leptons
    df = df.Define(f"zll_m_{flavor}", f"FCCAnalyses::ReconstructedParticle::get_mass(zll_{flavor})[0]")
    df = df.Define(f"zll_p_{flavor}", f"FCCAnalyses::ReconstructedParticle::get_p(zll_{flavor})[0]")
    df = df.Define(f"zll_recoil_{flavor}", f"FCCAnalyses::ReconstructedParticle::recoilBuilder(ecm)(zll_{flavor})")
    df = df.Define(f"zll_recoil_m_{flavor}", f"FCCAnalyses::ReconstructedParticle::get_mass(zll_recoil_{flavor})[0]")

    sel_leps = f"leps_{flavor}_no >= 2 && leps_{flavor}_sel_iso.size() > 0 && abs(Sum(leps_{flavor}_q)) < leps_{flavor}_q.size()"
    sel_mll = f"zll_m_{flavor} > 86 && zll_m_{flavor} < 96"
    sel_pll = f"zll_p_{flavor} > 20 && zll_p_{flavor} < 70"
    #sel_recoil = f"zll_recoil_m_{flavor} < 133 && zll_recoil_m_{flavor} > 123" # tighter selection on recoil --> check impact in leptonic
    #sel_recoil = f"zll_recoil_m_{flavor} < 140 && zll_recoil_m_{flavor} > 120" # tighter selection on recoil --> check impact in leptonic
    sel_recoil = f"zll_recoil_m_{flavor} < 150 && zll_recoil_m_{flavor} > 100" # tighter selection on recoil --> check impact in leptonic
    df = df.Filter(f"!({sel_leps} && {sel_mll} && {sel_pll} && {sel_recoil})")
    return df






def build_graph_zqq(df, dataset):

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

    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0")) ## all events




    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Alias("Electron0", "Electron#0.index")
    df = df.Alias("Photon0", "Photon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    df = df.Define("photons_all", "FCCAnalyses::ReconstructedParticle::get(Photon0, ReconstructedParticles)")

    ## orthogonality w.r.t. leptonic channels
    df = veto_selection_leptonic(df, "muon")
    df = veto_selection_leptonic(df, "electron")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1"))




    #### remove isolated photons from clustering
    df = df.Define("photons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(photons_all)")
    if ecm == 240:
        df = df.Define("photons", "FCCAnalyses::sel_range(40, 95, false)(photons_all, photons_all_p)")
    else:
        df = df.Define("photons", "FCCAnalyses::sel_range(20, 170, false)(photons_all, photons_all_p)")
    df = df.Define("rps_no_photons", "FCCAnalyses::ReconstructedParticle::remove(ReconstructedParticles, photons)")


    #### remove isolated muons from clustering
    # ensure orthogonality with leptonic channel
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    if ecm == 240:
        df = df.Define("muons", "FCCAnalyses::sel_range(40, 95, false)(muons_all, muons_all_p)")
    else:
        df = df.Define("muons", "FCCAnalyses::sel_range(20, 170, false)(muons_all, muons_all_p)")

    #### remove isolated electrons from clustering
    # ensure orthogonality with leptonic channel
    df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    if ecm == 240:
        df = df.Define("electrons", "FCCAnalyses::sel_range(40, 95, false)(electrons_all, electrons_all_p)")
    else:
        df = df.Define("electrons", "FCCAnalyses::sel_range(20, 170, false)(electrons_all, electrons_all_p)")

    ### remove isolated electrons from clustering
    # ensure orthogonality with leptonic channel
    df = df.Define("rps_no_photons_muons", "FCCAnalyses::ReconstructedParticle::remove(rps_no_photons, muons)")
    df = df.Define("rps_no_photons_muons_electrons", "FCCAnalyses::ReconstructedParticle::remove(rps_no_photons_muons, electrons)")


    # define PF candidates collection by removing the muons
    if ecm == 240:
        #df = df.Alias("rps_sel", "ReconstructedParticles") ## this was enabled for 240 GeV ### TO BE CHECKED
        df = df.Alias("rps_sel", "rps_no_photons_muons_electrons")
    else:
        df = df.Alias("rps_sel", "rps_no_photons_muons_electrons")
    df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_sel)")
    df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_sel)")
    df = df.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_sel)")
    df = df.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_sel)")
    df = df.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_sel)")
    df = df.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_sel)")
    df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")


    # perform possible clusterings
    
    df = exclusive_clustering(df, 6)
    df = exclusive_clustering(df, 4)
    df = exclusive_clustering(df, 2)
    df = exclusive_clustering(df, 0)

    df = define_variables(df, 6)
    df = define_variables(df, 4)
    df = define_variables(df, 2)
    df = define_variables(df, 0)


    df = df.Define("zqq", "std::vector<Vec_rp> r = {zqq_N0, zqq_N2, zqq_N4, zqq_N6}; return r;")
    df = df.Define("zqq_jets", "std::vector<Vec_rp> r = {jets_rp_cand_N0, jets_rp_cand_N2, jets_rp_cand_N4, jets_rp_cand_N6}; return r;") # wrong, must be jets, not RP
    df = df.Define("zqq_m", "Vec_f r = {zqq_m_N0, zqq_m_N2, zqq_m_N4, zqq_m_N6}; return r;")
    df = df.Define("zqq_p", "Vec_f r = {zqq_p_N0, zqq_p_N2, zqq_p_N4, zqq_p_N6}; return r;")
    df = df.Define("zqq_recoil_m", "Vec_f r = {zqq_recoil_m_N0, zqq_recoil_m_N2, zqq_recoil_m_N4, zqq_recoil_m_N6}; return r;")
    df = df.Define("njets", "Vec_i r = {(int)njets_cand_N0, (int)njets_cand_N2, (int)njets_cand_N4, (int)njets_cand_N6}; return r;")
    df = df.Define("njets_target", "Vec_i r = {0, 2, 4, 6}; return r;")
    df = df.Define("best_clustering_idx", "FCCAnalyses::best_clustering_idx(zqq_m, zqq_p, zqq_recoil_m, njets, njets_target, ecm)")

    hists.append(df.Histo1D(("best_clustering_idx_nosel", "", *(15, -5, 10)), "best_clustering_idx"))

    ## njets for inclusive clustering
    hists.append(df.Histo1D(("njets_inclusive", "", *bins_count), "njets_cand_N0"))
    df_incl = df.Filter("best_clustering_idx == 0")
    hists.append(df_incl.Histo1D(("njets_inclusive_sel", "", *bins_count), "njets_cand_N0"))




    df = df.Filter("best_clustering_idx >= 0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2")) ## after clustering

    df = df.Define("zqq_best", "zqq[best_clustering_idx]")
    df = df.Define("zqq_jets_best", "zqq_jets[best_clustering_idx]")
    df = df.Define("zqq_m_best", "zqq_m[best_clustering_idx]")
    df = df.Define("zqq_recoil_m_best", "zqq_recoil_m[best_clustering_idx]")
    df = df.Define("zqq_p_best", "zqq_p[best_clustering_idx]")
    df = df.Define("zqq_jets_p_best", "FCCAnalyses::ReconstructedParticle::get_p(zqq_jets_best)")
    df = df.Define("zqq_jets_theta_best", "FCCAnalyses::ReconstructedParticle::get_theta(zqq_jets_best)")
    df = df.Define("zqq_jets_costheta_best", "Vec_f ret; for(auto & theta: zqq_jets_theta_best) ret.push_back(std::abs(cos(theta))); return ret;")

    df = df.Define("z_theta", "FCCAnalyses::ReconstructedParticle::get_theta(zqq_best)")
    df = df.Define("z_costheta", "std::abs(cos(z_theta[0]))")
    hists.append(df.Histo1D(("z_costheta_nosel", "", *bins_cos_abs), "z_costheta"))
    

    hists.append(df.Histo1D(("zqq_m_best_nosel", "", *(400, 0, 400)), "zqq_m_best"))
    hists.append(df.Histo1D(("zqq_p_best_nosel", "", *(400, 0, 400)), "zqq_p_best"))
    hists.append(df.Histo1D(("zqq_recoil_m_best_nosel", "", *(400, 0, 400)), "zqq_recoil_m_best"))
    hists.append(df.Histo2D(("zqq_recoil_m_mqq_nosel", "", *((25, 100, 150)+(35, 50, 120))), "zqq_recoil_m_best", "zqq_m_best"))

    ## jet kinematics
    df = df.Define("leading_idx", "(zqq_jets_p_best[0] > zqq_jets_p_best[1]) ? 0 : 1")
    df = df.Define("subleading_idx", "(zqq_jets_p_best[0] > zqq_jets_p_best[1]) ? 1 : 0")
    df = df.Define("leading_jet_p", "zqq_jets_p_best[leading_idx]")
    df = df.Define("subleading_jet_p", "zqq_jets_p_best[subleading_idx]")
    df = df.Define("leading_jet_costheta", "zqq_jets_costheta_best[leading_idx]")
    df = df.Define("subleading_jet_costheta", "zqq_jets_costheta_best[subleading_idx]")
    
    hists.append(df.Histo1D(("leading_jet_p", "", *(400, 0, 400)), "leading_jet_p"))
    hists.append(df.Histo1D(("subleading_jet_p", "", *(400, 0, 400)), "subleading_jet_p"))
    hists.append(df.Histo1D(("leading_jet_costheta", "", *bins_cos_abs), "leading_jet_costheta"))
    hists.append(df.Histo1D(("subleading_jet_costheta", "", *bins_cos_abs), "subleading_jet_costheta"))

    ## attempt to reconstruct WW with 4 jets
    
    df = df.Define("pairs_WW_N4", "FCCAnalyses::pair_WW_N4(jets_rp_cand_N4)")
    df = df.Define("W1", "pairs_WW_N4[0]")
    df = df.Define("W2", "pairs_WW_N4[1]")
    df = df.Define("W1_m", "W1.M()")
    df = df.Define("W2_m", "W2.M()")
    df = df.Define("W1_p", "W1.P()")
    df = df.Define("W2_p", "W2.P()")
    df = df.Define("W1_costheta", "std::abs(W1.Theta())")
    df = df.Define("W2_costheta", "std::abs(W1.Theta())") # should be W2
    df = df.Define("delta_mWW", "std::sqrt((W1_m-78)*(W1_m-78) + (W2_m-78)*(W2_m-78))")

    hists.append(df.Histo1D(("delta_mWW_nosel", "", *(1000, 0, 100)), "delta_mWW"))
    hists.append(df.Histo1D(("W1_m_nosel", "", *(400, 0, 400)), "W1_m"))
    hists.append(df.Histo1D(("W2_m_nosel", "", *(400, 0, 400)), "W2_m"))
    hists.append(df.Histo1D(("W1_p_nosel", "", *(400, 0, 400)), "W1_p"))
    hists.append(df.Histo1D(("W2_p_nosel", "", *(400, 0, 400)), "W2_p"))
    hists.append(df.Histo1D(("W1_costheta_nosel", "", *bins_cos_abs), "W1_costheta"))
    hists.append(df.Histo1D(("W2_costheta_nosel", "", *bins_cos_abs), "W2_costheta"))
    ##
    '''
    df = df.Define("jets_rp_cand_N4_prime", "FCCAnalyses::energyReconstructFourJet_rp(240, jets_rp_cand_N4)")
    df = df.Define("pairs_WW_N4_prime", "FCCAnalyses::pair_WW_N4(jets_rp_cand_N4_prime)")
    df = df.Define("W1_prime", "pairs_WW_N4_prime[0]")
    df = df.Define("W2_prime", "pairs_WW_N4_prime[1]")
    df = df.Define("W1_m_prime", "W1_prime.M()")
    df = df.Define("W2_m_prime", "W2_prime.M()")
    df = df.Define("W1_p_prime", "W1_prime.M()")
    df = df.Define("W2_p_prime", "W2_prime.M()")
    df = df.Define("delta_mWW_prime", "std::sqrt((W1_m_prime-78)*(W1_m_prime-78) + (W2_m_prime-78)*(W2_m_prime-78))")

    hists.append(df.Histo1D(("delta_mWW_nosel_prime", "", *(1000, 0, 100)), "delta_mWW_prime"))
    hists.append(df.Histo1D(("W1_m_nosel_prime", "", *(400, 0, 400)), "W1_m_prime"))
    hists.append(df.Histo1D(("W2_m_nosel_prime", "", *(400, 0, 400)), "W2_m_prime"))
    hists.append(df.Histo1D(("W1_p_nosel_prime", "", *(400, 0, 400)), "W1_p_prime"))
    hists.append(df.Histo1D(("W2_p_nosel_prime", "", *(400, 0, 400)), "W2_p_prime"))
    '''
    ##############################################################
    ## START SELECTION
    ##############################################################



    ## CUT ON m(qq)
    
    hists.append(df.Histo1D(("zqq_m_best_nOne", "", *(400, 0, 400)), "zqq_m_best"))
    ##df = df.Filter("zqq_m_best > 60 && zqq_m_best < 110") # tighter at 80, but aa/mumu problem
    if ecm == 240:
        df = df.Filter("zqq_m_best > 20 && zqq_m_best < 140") ## loose
    else:
        df = df.Filter("zqq_m_best > 20 && zqq_m_best < 200") ## loose
    #df = df.Filter("zqq_m_best > 40 && zqq_m_best < 140") ## loose
    #df = df.Filter("zqq_m_best > 80 && zqq_m_best < 110") ## aggressive 80-60 seems not too much of a difference --> because we fit it!
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))


    ## CUT ON p(qq)
    hists.append(df.Histo1D(("zqq_p_best_nOne", "", *(400, 0, 400)), "zqq_p_best"))
    ##df = df.Filter("zqq_p_best < 60 && zqq_p_best > 30")
    #df = df.Filter("zqq_p_best < 60") # see significance
    if ecm == 240:
        df = df.Filter("zqq_p_best < 90 && zqq_p_best > 20")
    else:
        df = df.Filter("zqq_p_best < 160 && zqq_p_best > 60")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))


    ## cut on cos(qq)
    hists.append(df.Histo1D(("z_costheta_nOne", "", *bins_cos_abs), "z_costheta"))
    #df = df.Filter("z_costheta < 0.95") # tight cut on pqq
    df = df.Filter("z_costheta < 0.85") # after loose cut on pqq
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))

    ## acolinearity
    df = df.Define("acolinearity", "FCCAnalyses::acolinearity(zqq_jets_best)")
    hists.append(df.Histo1D(("acolinearity_nOne", "", *bins_aco), "acolinearity"))
    #df = df.Filter("acolinearity > 0.35") # tight cut on pqq
    df = df.Filter("acolinearity > 0.35") # sigificance 0.45, but keep it
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))


    ## acoplanarity
    df = df.Define("acoplanarity", "FCCAnalyses::acoplanarity(zqq_jets_best)")
    hists.append(df.Histo1D(("acoplanarity_nOne", "", *bins_aco), "acoplanarity"))
    #df = df.Filter("acoplanarity < 5") ## NO CUT
    #hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))




    hists.append(df.Histo2D(("W1_m_W2_m_nOne", "", *((150, 0, 150)+(150, 0, 150))), "W1_m", "W2_m"))
    hists.append(df.Histo1D(("delta_mWW_nOne", "", *(1000, 0, 100)), "delta_mWW"))
    hists.append(df.Histo1D(("W1_m_nOne", "", *(400, 0, 400)), "W1_m"))
    hists.append(df.Histo1D(("W2_m_nOne", "", *(400, 0, 400)), "W2_m"))


    df =df.Filter("delta_mWW > 6")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))


    hists.append(df.Histo1D(("W1_p_nOne", "", *(400, 0, 400)), "W1_p"))
    hists.append(df.Histo1D(("W2_p_nOne", "", *(400, 0, 400)), "W2_p"))


    ####
    ## CUT 7: cos theta(miss)
    ####
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("cosThetaMiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df.Histo1D(("cosThetaMiss_nOne", "", *bins_cosThetaMiss), "cosThetaMiss"))
    
    df = df.Filter("cosThetaMiss < .995")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))

    ## thrust
    df = df.Define("rps_charged_idx", "RP_q != 0") # select only charged particles
    df = df.Define("rps_charged", "ReconstructedParticles[rps_charged_idx]")
    df = df.Define("rps_charged_p", "FCCAnalyses::ReconstructedParticle::get_p(rps_charged)")
    df = df.Define("rps_charged_n", "FCCAnalyses::ReconstructedParticle::get_n(rps_charged)")
    df = df.Define("RP_px_q", "RP_px[rps_charged_idx]")
    df = df.Define("RP_py_q", "RP_py[rps_charged_idx]")
    df = df.Define("RP_pz_q", "RP_pz[rps_charged_idx]")
    df = df.Define("max_p_idx", "int idx=-1; float max=-1; for(int i=0; i<rps_charged_n; i++) if(rps_charged_p[i]>max) {max= rps_charged_p[i]; idx=i;}; return idx")
    
    
    df = df.Define("thrust", "FCCAnalyses::Algorithms::calculate_thrust()(RP_px, RP_py, RP_pz)")
    df = df.Define("thrust_magn", "thrust[0]")
    df = df.Define("thrust_costheta", "abs(thrust[3])")

    
    #df = df.Define("thrust", "FCCAnalyses::minimize_thrust_mc(RP_px_q[max_p_idx], RP_py_q[max_p_idx], RP_pz_q[max_p_idx])(RP_px_q, RP_py_q, RP_pz_q)")
    #df = df.Define("thrust_magn", "thrust[0]")
    #df = df.Define("thrust_costheta", "abs(cos(thrust[5]))")

    hists.append(df.Histo1D(("thrust_magn", "", *(1000, 0, 1)), "thrust_magn"))
    hists.append(df.Histo1D(("thrust_costheta", "", *bins_cos_abs), "thrust_costheta"))

    if ecm == 365:
        df = df.Filter("thrust_magn < 0.85")
        hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))

    df = tmva_helper.run_inference(df, col_name="mva_score")
    hists.append(df.Histo1D(("mva_score", "", *(1000, 0, 1)), "mva_score"))

    #df = df.Filter("mva_score[0] > 0.4")
    ## final histograms
    hists.append(df.Histo1D(("zqq_recoil_m", "", *(400, 0, 400)), "zqq_recoil_m_best"))
    hists.append(df.Histo1D(("zqq_m", "", *(400, 0, 400)), "zqq_m_best"))
    #hists.append(df.Histo2D(("zqq_recoil_m_mqq", "", *((50, 100, 150)+(60, 60, 120))), "zqq_recoil_m_best", "zqq_m_best"))
    hists.append(df.Histo3D(("zqq_recoil_m_mqq_pqq", "", *((50, 100, 150)+(100, 40, 140)+(40, 20, 60))), "zqq_recoil_m_best", "zqq_m_best", "zqq_p_best"))
    hists.append(df.Histo2D(("zqq_recoil_m_mqq", "", *((50, 100, 150)+(100, 40, 140))), "zqq_recoil_m_best", "zqq_m_best"))
    hists.append(df.Histo2D(("zqq_recoil_m_pqq", "", *((50, 100, 150)+(40, 20, 60))), "zqq_recoil_m_best", "zqq_p_best"))

    bins_mva_ = [0, 0.33, 0.66, 1]
    
    if ecm == 240:
        bins_mrec_ = list(range(100, 151, 1))
        bins_mva_ = [0, 0.75, 1] # 0.75 = max significance
        bins_mva_ = [0, 0.70, 1] # 0.75 = max significance
    else:
        bins_mrec_ = list(range(100, 201, 1))
        bins_mva_ = [0, 0.95, 1] # 0.75 = max significance
    bins_m_ = list(range(40, 141, 1)) 
    bins_mva = array.array('d', bins_mva_)
    bins_mrec = array.array('d', bins_mrec_)
    bins_m = array.array('d', bins_m_)
    model = ROOT.RDF.TH3DModel("zqq_recoil_m_mqq_mva", "", len(bins_mrec_)-1, bins_mrec, len(bins_m_)-1, bins_m, len(bins_mva_)-1, bins_mva)
    hists.append(df.Histo3D(model, "zqq_recoil_m_best", "zqq_m_best", "mva_score"))

    hists.append(df.Histo2D(("zqq_recoil_m_mva", "", *((50, 100, 150)+(2, 0, 1))), "zqq_recoil_m_best", "mva_score"))
    #hists.append(df.Histo3D(("zqq_recoil_m_mqq_mva", "", *((50, 100, 150)+(100, 40, 140)+(2, 0, 1))), "zqq_recoil_m_best", "zqq_m_best", "mva_score"))
    return hists, weightsum, df


if treemaker:
    class RDFanalysis():
        def analysers(df):
            hists, weightsum, df = build_graph_zqq(df, "")
            return df

        # define output branches to be saved
        def output():
            branchList = ["zqq_recoil_m_best", "zqq_p_best", "zqq_m_best", "leading_jet_costheta", "subleading_jet_costheta", "leading_jet_p", "subleading_jet_p", "acolinearity", "acoplanarity", "W1_p", "W2_p", "W1_m", "W2_m", "z_costheta", "thrust_magn", "W1_costheta", "W2_costheta"]
            return branchList

else:
    def build_graph(df, dataset):
        hists, weightsum, df = build_graph_zqq(df, dataset)
        return hists, weightsum
