
import ROOT
import array
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB

ecm = 240 # 240 365
treemaker = False

# list of all processes
fraction = 0.1


if ecm == 240:
    processList = {
        #'wzp6_ee_nunuH_HWW_ecm240':         {'fraction':1},
        #'wz3p6_ee_gammagamma_ecm240':       {'fraction':fraction},
        #'wz3p6_ee_tautau_ecm240':           {'fraction':fraction},
        #'wz3p6_ee_mumu_ecm240':             {'fraction':fraction},
        #'wz3p6_ee_ee_Mee_30_150_ecm240':    {'fraction':fraction},
        #'wz3p6_ee_nunu_ecm240':             {'fraction':fraction},


        'wzp6_ee_eeH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HWW_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HWW_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HWW_ecm240':           {'fraction':1},

        'wzp6_ee_nunuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HZZ_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':1},

        'wz3p6_ee_uu_ecm240':               {'fraction':fraction},
        'wz3p6_ee_dd_ecm240':               {'fraction':fraction},
        'wz3p6_ee_cc_ecm240':               {'fraction':fraction},
        'wz3p6_ee_ss_ecm240':               {'fraction':fraction},
        'wz3p6_ee_bb_ecm240':               {'fraction':fraction},

        'p8_ee_ZZ_ecm240':             {'fraction':fraction},
        'p8_ee_WW_ecm240':             {'fraction':fraction},
    }

process = {
        'wzp6_ee_nunuH_HWW_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HWW_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HWW_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HWW_ecm240':           {'fraction':1},

        'wzp6_ee_nunuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HZZ_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':1},
}


processListOld = {
        #'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':1},
        #'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':1},
        #'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':1},
        #'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':1},
        #'wzp6_ee_tautauH_HZZ_ecm240':       {'fraction':1},

        'wzp6_ee_tautauH_HWW_ecm240':       {'fraction':1},
        'wzp6_ee_ccH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HWW_ecm240':           {'fraction':1},
}

processListdd = {
        'wzp6_ee_nunuH_HWW_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HWW_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HWW_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HWW_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HWW_ecm240':           {'fraction':1},

        'wzp6_ee_nunuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_eeH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_tautauH_HZZ_ecm240':       {'fraction':1},
        'wzp6_ee_mumuH_HZZ_ecm240':         {'fraction':1},
        'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':1},
}

inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "../h_zh/utils.h", "utils.h"]


# output directory
if treemaker:
    outputDir   = f"/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/treemaker/ecm{ecm}/"
else:
    outputDir   = f"output/h_zz_ww/histmaker/ecm{ecm}/preSelPlots"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32 ## good for Z/g -- 21 seconds fraction = 0.01 (100 files per DS) also good for WW/ZZ (30 min for 25% of stats)
#nCPUS       = 8 ## WW(40 files) ZZ(20 files)
#nCPUS       = 256 ## WW(40 files) ZZ(20 files)

nCPUs = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 if ecm == 240 else 3e6


# define histograms
bins_m = (400, 0, 400) # 100 MeV bins
bins_p_mu = (2000, 0, 200) # 100 MeV bins
bins_m_ll = (2000, 0, 200) # 100 MeV bins
bins_p_ll = (200, 0, 200) # 1 GeV bins
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

bins_merge = (50000, 0, 50000)
bins_merge_sq = (500, 0, 500)
bins_score = (100, 0, 1)


ROOT.EnableImplicitMT(nCPUS) # hack to deal correctly with TMVAHelperXGB  # bdt_model_new_0p1_WZqq bdt_model_new_0p1 bdt_model_0p1_inv bdt_model_0p1
## original: output/h_zh_hadronic/training/bdt_model_WW_Zg_thrust_reduced.root
## before: /ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/bdt_model_final.root
#tmva_helper = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/zh/hadronic/training/bdt_model_ecm240.root", "bdt_model") # read the XGBoost training
#tmva_helper = TMVAHelperXGB("output/h_zh_hadronic/training/bdt_model_WW_Zg_thrust_reduced_withInv.root", "bdt_model") # read the XGBoost training



################################################################################
## load modules and files for jet clustering and flavor tagging
################################################################################

# files from /eos/experiment/fcc/ee/jet_flavour_tagging/winter2023/wc_pt_7classes_12_04_2023
weaver_preproc = "fccee_flavtagging_edm4hep_wc.json"
weaver_model = "fccee_flavtagging_edm4hep_wc.onnx"

from addons.ONNXRuntime.jetFlavourHelper import JetFlavourHelper
from addons.FastJet.jetClusteringHelper import ExclusiveJetClusteringHelper
from examples.FCCee.weaver.config import collections, njets

################################################################################




def build_graph_hwwzz(df, dataset):

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

    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut0")) ## all events
    if "HWW" in dataset: # remove muons/electrons from inclusive WW ## ONLY HADRONIC
        df = df.Define("ww_hadronic", "FCCAnalyses::is_ww_hadronic(Particle, Particle1)")
        df = df.Filter("ww_hadronic")

    if "HZZ" in dataset: # remove muons/electrons from inclusive ZZ ## ONLY HADRONIC
        df = df.Define("zz_hadronic", "FCCAnalyses::is_zz_hadronic(Particle, Particle1)")
        df = df.Filter("zz_hadronic")


    #df = df.Define("is_hadronic", "FCCAnalyses::is_zz_hadronic(Particle, Particle1)")
    #df = df.Define("is_hadronic", "FCCAnalyses::is_ww_hadronic(Particle, Particle1)")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut1")) ## hadronic events





    # muons
    df = df.Alias("Muon0", "Muon#0.index")
    df = df.Define("muons_all", "FCCAnalyses::ReconstructedParticle::get(Muon0, ReconstructedParticles)")
    df = df.Define("muons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(muons_all)")
    df = df.Define("muons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons_all)")
    df = df.Define("muons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons_all)")
    df = df.Define("muons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons_all)")
    df = df.Define("muons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(muons_all)")

    df = df.Define("muons", "FCCAnalyses::ReconstructedParticle::sel_p(10)(muons_all)")
    df = df.Define("muons_p", "FCCAnalyses::ReconstructedParticle::get_p(muons)")
    df = df.Define("muons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(muons)")
    df = df.Define("muons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(muons)")
    df = df.Define("muons_q", "FCCAnalyses::ReconstructedParticle::get_charge(muons)")
    df = df.Define("muons_no", "FCCAnalyses::ReconstructedParticle::get_n(muons)")

    
    # electrons
    df = df.Alias("Electron0", "Electron#0.index")
    df = df.Define("electrons_all", "FCCAnalyses::ReconstructedParticle::get(Electron0, ReconstructedParticles)")
    df = df.Define("electrons_all_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons_all)")
    df = df.Define("electrons_all_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons_all)")
    df = df.Define("electrons_all_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons_all)")
    df = df.Define("electrons_all_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons_all)")
    df = df.Define("electrons_all_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons_all)")

    df = df.Define("electrons", "FCCAnalyses::ReconstructedParticle::sel_p(10)(electrons_all)")
    df = df.Define("electrons_p", "FCCAnalyses::ReconstructedParticle::get_p(electrons)")
    df = df.Define("electrons_theta", "FCCAnalyses::ReconstructedParticle::get_theta(electrons)")
    df = df.Define("electrons_phi", "FCCAnalyses::ReconstructedParticle::get_phi(electrons)")
    df = df.Define("electrons_q", "FCCAnalyses::ReconstructedParticle::get_charge(electrons)")
    df = df.Define("electrons_no", "FCCAnalyses::ReconstructedParticle::get_n(electrons)")

    ## VETO leptons
    df = df.Filter("muons_no == 0 && electrons_no == 0")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut2"))



    ## MISSING MOMENTUM
    df = df.Define("missingEnergy_rp", "FCCAnalyses::missingEnergy(ecm, ReconstructedParticles)")
    df = df.Define("missingEnergy", "missingEnergy_rp[0].energy")
    df = df.Define("missingEnergy_p", "FCCAnalyses::ReconstructedParticle::get_p(missingEnergy_rp)[0]")
    df = df.Define("cosThetaMiss", "FCCAnalyses::get_cosTheta_miss(missingEnergy_rp)")
    hists.append(df.Histo1D(("missingEnergy_nOne", "", *bins_m), "missingEnergy_p"))
    df = df.Filter("missingEnergy_p < 20")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut3"))


    ## VISIBLE ENERGY
    df = df.Define("visibleEnergy", "FCCAnalyses::visibleEnergy(ReconstructedParticles)")
    hists.append(df.Histo1D(("visibleEnergy_nOne", "", *bins_m), "visibleEnergy"))
    df = df.Filter("visibleEnergy > 115")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut5"))
    

    ## CLUSTERING
    df = df.Alias("rps_sel", "ReconstructedParticles")
    df = df.Define("RP_px", "FCCAnalyses::ReconstructedParticle::get_px(rps_sel)")
    df = df.Define("RP_py", "FCCAnalyses::ReconstructedParticle::get_py(rps_sel)")
    df = df.Define("RP_pz","FCCAnalyses::ReconstructedParticle::get_pz(rps_sel)")
    df = df.Define("RP_e", "FCCAnalyses::ReconstructedParticle::get_e(rps_sel)")
    df = df.Define("RP_m", "FCCAnalyses::ReconstructedParticle::get_mass(rps_sel)")
    df = df.Define("RP_q", "FCCAnalyses::ReconstructedParticle::get_charge(rps_sel)")
    df = df.Define("pseudo_jets", "FCCAnalyses::JetClusteringUtils::set_pseudoJets(RP_px, RP_py, RP_pz, RP_e)")

    df = df.Define("clustered_jets", "JetClustering::clustering_ee_kt(2, 6, 0, 10)(pseudo_jets)")
    df = df.Define("jets", "FCCAnalyses::JetClusteringUtils::get_pseudoJets(clustered_jets)")
    df = df.Define("njets", "jets.size()")
    df = df.Define("jetconstituents", "FCCAnalyses::JetClusteringUtils::get_constituents(clustered_jets)")
    df = df.Define("jets_e", "FCCAnalyses::JetClusteringUtils::get_e(jets)")
    df = df.Define("jets_px", "FCCAnalyses::JetClusteringUtils::get_px(jets)")
    df = df.Define("jets_py", "FCCAnalyses::JetClusteringUtils::get_py(jets)")
    df = df.Define("jets_pz", "FCCAnalyses::JetClusteringUtils::get_pz(jets)")
    df = df.Define("jets_m", "FCCAnalyses::JetClusteringUtils::get_m(jets)")

    hists.append(df.Histo1D(("njets_nOne", "", *bins_count), "njets"))
    df = df.Filter("njets == 6")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut4"))






    ## VISIBLE MASS -- REDUNDANT
    #df = df.Define("visibleMass", "FCCAnalyses::visibleMass(ReconstructedParticles)")
    #hists.append(df.Histo1D(("visibleMass_nOne", "", *bins_m), "visibleMass"))
    #df = df.Filter("visibleMass > 10")
    #hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))




    ## cuts on dmerge
    df = df.Define("dmerge_01", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 0)")
    df = df.Define("dmerge_12", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 1)")
    df = df.Define("dmerge_23", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 2)")
    df = df.Define("dmerge_34", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 3)")
    df = df.Define("dmerge_45", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 4)")
    df = df.Define("dmerge_56", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 5)")
    df = df.Define("dmerge_67", "FCCAnalyses::JetClusteringUtils::get_exclusive_dmerge(clustered_jets, 6)")

    df = df.Define("sq_dmerge_01", "sqrt(dmerge_01)")
    df = df.Define("sq_dmerge_12", "sqrt(dmerge_12)")
    df = df.Define("sq_dmerge_23", "sqrt(dmerge_23)")
    df = df.Define("sq_dmerge_34", "sqrt(dmerge_34)")
    df = df.Define("sq_dmerge_45", "sqrt(dmerge_45)")
    df = df.Define("sq_dmerge_56", "sqrt(dmerge_56)")
    df = df.Define("sq_dmerge_67", "sqrt(dmerge_67)")

    hists.append(df.Histo1D(("dmerge_01", "", *bins_merge), "dmerge_01"))
    hists.append(df.Histo1D(("sq_dmerge_01", "", *bins_merge_sq), "sq_dmerge_01"))
    hists.append(df.Histo1D(("dmerge_12", "", *bins_merge), "dmerge_12"))
    hists.append(df.Histo1D(("sq_dmerge_12", "", *bins_merge_sq), "sq_dmerge_12"))
    hists.append(df.Histo1D(("dmerge_23", "", *bins_merge), "dmerge_23"))
    hists.append(df.Histo1D(("sq_dmerge_23", "", *bins_merge_sq), "sq_dmerge_23"))

    df = df.Filter("sq_dmerge_23 > 50")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut6"))


    hists.append(df.Histo1D(("dmerge_34", "", *bins_merge), "dmerge_34"))
    hists.append(df.Histo1D(("sq_dmerge_34", "", *bins_merge_sq), "sq_dmerge_34"))

    df = df.Filter("sq_dmerge_34 > 25")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut7"))

    hists.append(df.Histo1D(("dmerge_45", "", *bins_merge), "dmerge_45"))
    hists.append(df.Histo1D(("sq_dmerge_45", "", *bins_merge_sq), "sq_dmerge_45"))

    df = df.Filter("sq_dmerge_45 > 15")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut8"))

    hists.append(df.Histo1D(("dmerge_56", "", *bins_merge), "dmerge_56"))
    hists.append(df.Histo1D(("sq_dmerge_56", "", *bins_merge_sq), "sq_dmerge_56"))

    df = df.Filter("sq_dmerge_56 > 10")
    hists.append(df.Histo1D(("cutFlow", "", *bins_count), "cut9"))

    hists.append(df.Histo1D(("dmerge_67", "", *bins_merge), "dmerge_67"))
    hists.append(df.Histo1D(("sq_dmerge_67", "", *bins_merge_sq), "sq_dmerge_67"))



    ## flavor tagging
    df = df.Define("jetconstituents_clusters", "JetConstituentsUtils::build_constituents_cluster(ReconstructedParticles, jetconstituents)")
    jetFlavourHelper = JetFlavourHelper(collections, "jets", "jetconstituents_clusters", "")
    df = jetFlavourHelper.define(df) # define variables
    df = jetFlavourHelper.inference(weaver_preproc, weaver_model, df) # run inference


    return hists, weightsum, df





if treemaker:
    class RDFanalysis():
        def analysers(df):
            hists, weightsum, df = build_graph_hwwzz(df, "")
            return df

        # define output branches to be saved
        def output():
            branchList = ["jets_e", "jets_px", "jets_py", "jets_pz", "jets_m", "visibleEnergy", "visibleMass", "missingEnergy_p", "missingEnergy", "cosThetaMiss", "dmerge_01", "dmerge_12", "dmerge_23", "dmerge_34", "dmerge_45", "dmerge_56", "dmerge_67", "recojet_isB", "recojet_isC", "recojet_isS", "recojet_isG", "recojet_isU", "recojet_isD", "recojet_isTAU"]
            return branchList


else:
    def build_graph(df, dataset):
        hists, weightsum, df = build_graph_hwwzz(df, dataset)
        return hists, weightsum
