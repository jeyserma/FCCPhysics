
import ROOT
import array
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB

ecm = 240 # 240 365
treemaker = False
fraction = 0.01


processList = {
    'wzp6_ee_ccH_HWW_ecm240':           {'fraction':1},
    'wzp6_ee_bbH_HWW_ecm240':           {'fraction':1},
    'wzp6_ee_qqH_HWW_ecm240':           {'fraction':1},
    'wzp6_ee_ssH_HWW_ecm240':           {'fraction':1},

    'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':1},
    'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':1},
    'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':1},
    'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':1},
}

inputDir = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "../h_zh/utils.h", "utils.h"]


if treemaker:
    outputDir   = f"/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/treemaker/ecm{ecm}/jet_pairing"
else:
    outputDir   = f"output/h_zz_ww/histmaker/ecm{ecm}/jet_pairing/"

nCPUS       = 32

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 if ecm == 240 else 3e6

# define histograms

# define histograms
bins_m = (400, 0, 400) # 100 MeV bins
bins_chi2 = (1000, 0, 10000)

def build_graph_hwwzz(df, dataset):

    hists = []

    df = df.Define("ecm", "240" if ecm == 240 else "365")
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")

    if "HWW" in dataset: # remove muons/electrons from inclusive WW ## ONLY HADRONIC
        df = df.Define("ww_hadronic", "FCCAnalyses::is_ww_hadronic(Particle, Particle1)")
        df = df.Filter("ww_hadronic")

    if "HZZ" in dataset: # remove muons/electrons from inclusive ZZ ## ONLY HADRONIC
        df = df.Define("zz_hadronic", "FCCAnalyses::is_zz_hadronic(Particle, Particle1)")
        df = df.Filter("zz_hadronic")


    ## JET CLUSTERING
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

    df = df.Define("init_tlv", "TLorentzVector ret; ret.SetPxPyPzE(0, 0, 0, ecm); return ret;")
    df = df.Define("jets_tlv", "FCCAnalyses::makeLorentzVectors(jets_px, jets_py, jets_pz, jets_e)")






    #df = df.Define("paired_jets_WW", "FCCAnalyses::pairing_WW_ZZ(jets_tlv, 80.385)")
    #df = df.Define("paired_jets_ZZ", "FCCAnalyses::pairing_WW_ZZ(jets_tlv, 91.19)")
    
    df = df.Define("paired_jets_WW", "FCCAnalyses::pairing_test(jets_tlv, 80.385)")
    df = df.Define("paired_jets_ZZ", "FCCAnalyses::pairing_test(jets_tlv, 91.19)")
    #df = df.Define("paired_jets_ZZ", "FCCAnalyses::pairing_ZZ(jets_tlv, 91.19)")
    #df = df.Define("paired_jets_ZZ", "FCCAnalyses::pairing_ZZ_flavor(jets_tlv, recojet_isB, recojet_isC, recojet_isS)")

    df = df.Define("W1_WW", "jets_tlv[paired_jets_WW[0]] + jets_tlv[paired_jets_WW[1]]")
    df = df.Define("W2_WW", "jets_tlv[paired_jets_WW[2]] + jets_tlv[paired_jets_WW[3]]")
    df = df.Define("Z_WW", "jets_tlv[paired_jets_WW[4]] + jets_tlv[paired_jets_WW[5]]")
    df = df.Define("H_WW", "W1_WW+W2_WW")

    df = df.Define("Z1_ZZ", "jets_tlv[paired_jets_ZZ[0]] + jets_tlv[paired_jets_ZZ[1]]")
    df = df.Define("Z2_ZZ", "jets_tlv[paired_jets_ZZ[2]] + jets_tlv[paired_jets_ZZ[3]]")
    df = df.Define("Z_ZZ", "jets_tlv[paired_jets_ZZ[4]] + jets_tlv[paired_jets_ZZ[5]]")
    df = df.Define("H_ZZ", "Z1_ZZ+Z2_ZZ")

    df = df.Define("chi2_WW", "paired_jets_WW[6]")
    df = df.Define("chi2_ZZ", "paired_jets_ZZ[6]")

    hists.append(df.Histo1D(("chi2_WW", "", *bins_chi2), "chi2_WW"))
    hists.append(df.Histo1D(("chi2_ZZ", "", *bins_chi2), "chi2_ZZ"))

    df = df.Define("W1_WW_m", "(float)W1_WW.M()")
    df = df.Define("W2_WW_m", "(float)W2_WW.M()")
    df = df.Define("Z_WW_m", "(float)Z_WW.M()")
    df = df.Define("Z_WW_p", "(float)Z_WW.P()")
    df = df.Define("Z_WW_theta", "(float)Z_WW.Theta()")
    df = df.Define("H_WW_m", "(float)H_WW.M()")
    df = df.Define("H_WW_p", "(float)H_WW.P()")
    df = df.Define("H_WW_theta", "(float)H_WW.Theta()")

    hists.append(df.Histo1D(("W1_WW_m", "", *bins_m), "W1_WW_m"))
    hists.append(df.Histo1D(("W2_WW_m", "", *bins_m), "W2_WW_m"))
    hists.append(df.Histo1D(("Z_WW_m", "", *bins_m), "Z_WW_m"))
    hists.append(df.Histo1D(("Z_WW_p", "", *bins_m), "Z_WW_p"))
    hists.append(df.Histo1D(("H_WW_m", "", *bins_m), "H_WW_m"))
    hists.append(df.Histo1D(("H_WW_p", "", *bins_m), "H_WW_p"))

    df = df.Define("Z1_ZZ_m", "(float)Z1_ZZ.M()")
    df = df.Define("Z2_ZZ_m", "(float)Z2_ZZ.M()")
    df = df.Define("Z_ZZ_m", "(float)Z_ZZ.M()")
    df = df.Define("Z_ZZ_p", "(float)Z_ZZ.P()")
    df = df.Define("Z_ZZ_theta", "(float)Z_ZZ.Theta()")
    df = df.Define("H_ZZ_m", "(float)H_ZZ.M()")
    df = df.Define("H_ZZ_p", "(float)H_ZZ.P()")
    df = df.Define("H_ZZ_theta", "(float)H_ZZ.Theta()")

    hists.append(df.Histo1D(("Z1_ZZ_m", "", *bins_m), "Z1_ZZ_m"))
    hists.append(df.Histo1D(("Z2_ZZ_m", "", *bins_m), "Z2_ZZ_m"))
    hists.append(df.Histo1D(("Z_ZZ_m", "", *bins_m), "Z_ZZ_m"))
    hists.append(df.Histo1D(("Z_ZZ_p", "", *bins_m), "Z_ZZ_p"))
    hists.append(df.Histo1D(("H_ZZ_m", "", *bins_m), "H_ZZ_m"))
    hists.append(df.Histo1D(("H_ZZ_p", "", *bins_m), "H_ZZ_p"))


    return hists, weightsum, df












if treemaker:
    class RDFanalysis():
        def analysers(df):
            hists, weightsum, df = build_graph_hwwzz(df, "")
            return df

        # define output branches to be saved
        def output():

            vars_zz_ww = ["W1_WW_m", "W2_WW_m", "Z_WW_m", "Z_WW_p", "Z_WW_theta", "H_WW_m", "H_WW_p", "H_WW_theta", "Z1_ZZ_m", "Z2_ZZ_m", "Z_ZZ_m", "Z_ZZ_p", "Z_ZZ_theta", "H_ZZ_m", "H_ZZ_p", "H_ZZ_theta", "sum_score_B_Z1", "sum_score_B_Z2", "sum_score_C_Z1", "sum_score_C_Z2", "sum_score_S_Z1", "sum_score_S_Z2", "sum_score_U_Z1", "sum_score_U_Z2", "sum_score_D_Z1", "sum_score_D_Z2", "sum_score_TAU_Z1", "sum_score_TAU_Z2", "sum_score_B_W1", "sum_score_B_W2", "sum_score_C_W1", "sum_score_C_W2", "sum_score_S_W1", "sum_score_S_W2", "sum_score_U_W1", "sum_score_U_W2", "sum_score_D_W1", "sum_score_D_W2", "sum_score_TAU_W1", "sum_score_TAU_W2"]

            vars_kinematics = ["W1_WW_m", "W2_WW_m", "Z_WW_m", "Z_WW_p", "Z_WW_theta", "H_WW_m", "H_WW_p", "H_WW_theta", "Z1_ZZ_m", "Z2_ZZ_m", "Z_ZZ_m", "Z_ZZ_p", "Z_ZZ_theta", "H_ZZ_m", "H_ZZ_p", "H_ZZ_theta"]
            
            vars_jets = ["W1_WW_jet1_p", "W1_WW_jet2_p", "W2_WW_jet1_p", "W2_WW_jet2_p", "W1_WW_jet1_theta", "W1_WW_jet2_theta", "W2_WW_jet1_theta", "W2_WW_jet2_theta", "Z1_ZZ_jet1_p", "Z1_ZZ_jet2_p", "Z2_ZZ_jet1_p", "Z2_ZZ_jet2_p", "Z1_ZZ_jet1_theta", "Z1_ZZ_jet2_theta", "Z2_ZZ_jet1_theta", "Z2_ZZ_jet2_theta"]

            vars_flavor = ["W1_WW_jet1_isB", "W1_WW_jet2_isB", "W1_WW_jet1_isC", "W1_WW_jet2_isC", "W1_WW_jet1_isS", "W1_WW_jet2_isS", "W1_WW_jet1_isU", "W1_WW_jet2_isU", "W1_WW_jet1_isD", "W1_WW_jet2_isD", "W1_WW_jet1_isTAU", "W1_WW_jet2_isTAU", "W2_WW_jet1_isB", "W2_WW_jet2_isB", "W2_WW_jet1_isC", "W2_WW_jet2_isC", "W2_WW_jet1_isS", "W2_WW_jet2_isS", "W2_WW_jet1_isU", "W2_WW_jet2_isU", "W2_WW_jet1_isD", "W2_WW_jet2_isD", "W2_WW_jet1_isTAU", "W2_WW_jet2_isTAU", "Z1_ZZ_jet1_isB", "Z1_ZZ_jet2_isB", "Z1_ZZ_jet1_isC", "Z1_ZZ_jet2_isC", "Z1_ZZ_jet1_isS", "Z1_ZZ_jet2_isS", "Z1_ZZ_jet1_isU", "Z1_ZZ_jet2_isU", "Z1_ZZ_jet1_isD", "Z1_ZZ_jet2_isD", "Z1_ZZ_jet1_isTAU", "Z1_ZZ_jet2_isTAU", "Z2_ZZ_jet1_isB", "Z2_ZZ_jet2_isB", "Z2_ZZ_jet1_isC", "Z2_ZZ_jet2_isC", "Z2_ZZ_jet1_isS", "Z2_ZZ_jet2_isS", "Z2_ZZ_jet1_isU", "Z2_ZZ_jet2_isU", "Z2_ZZ_jet1_isD", "Z2_ZZ_jet2_isD", "Z2_ZZ_jet1_isTAU", "Z2_ZZ_jet2_isTAU"]

            vars_zz_ww = vars_kinematics + vars_jets + vars_flavor




            
            #branchList = ["jets_e", "jets_px", "jets_py", "jets_pz", "jets_m", "visibleEnergy", "visibleMass", "missingEnergy_p", "cosThetaMiss", "paired_jets_WW", "paired_jets_ZZ", "dmerge_01", "dmerge_12", "dmerge_23", "dmerge_34", "dmerge_45", "dmerge_56", "dmerge_67", "sq_dmerge_01", "sq_dmerge_12", "sq_dmerge_23", "sq_dmerge_34", "sq_dmerge_45", "sq_dmerge_56", "sq_dmerge_67", "recojet_isB", "recojet_isC", "recojet_isS", "recojet_isG", "recojet_isU", "recojet_isD", "recojet_isTAU"]
            return vars_zz_ww



else:
    def build_graph(df, dataset):
        hists, weightsum, df = build_graph_hwwzz(df, dataset)
        return hists, weightsum
