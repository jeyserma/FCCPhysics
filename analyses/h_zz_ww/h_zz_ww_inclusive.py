
import ROOT
import array
import numpy as np
ROOT.TH1.SetDefaultSumw2(ROOT.kTRUE)
from addons.TMVAHelper.TMVAHelper import TMVAHelperXGB

ecm = 240 # 240 365
treemaker = False

# list of all processes
fraction = 1


if ecm == 240:
    processList = {
        #'wzp6_ee_nunuH_HWW_ecm240':         {'fraction':1}, # no events
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
        #'wz3p6_ee_gammagamma_ecm240':       {'fraction':fraction}, # no events
        #'wz3p6_ee_tautau_ecm240':           {'fraction':fraction},
        #'wz3p6_ee_mumu_ecm240':             {'fraction':fraction},
        #'wz3p6_ee_ee_Mee_30_150_ecm240':    {'fraction':fraction},
        #'wz3p6_ee_nunu_ecm240':             {'fraction':fraction}, # no events

        'p8_ee_ZZ_ecm240':             {'fraction':fraction},
        'p8_ee_WW_ecm240':             {'fraction':fraction},
    }

    processListSig = {
        #'wzp6_ee_ccH_HWW_ecm240':           {'fraction':1},
        #'wzp6_ee_bbH_HWW_ecm240':           {'fraction':1},
        #'wzp6_ee_qqH_HWW_ecm240':           {'fraction':1},
        #'wzp6_ee_ssH_HWW_ecm240':           {'fraction':1},

        'wzp6_ee_ccH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_bbH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_qqH_HZZ_ecm240':           {'fraction':1},
        'wzp6_ee_ssH_HZZ_ecm240':           {'fraction':1},

    }


#processList = processListSig

inputDir = f"/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/treemaker/ecm{ecm}/"
procDict = "/ceph/submit/data/group/fcc/ee/generation/DelphesEvents/winter2023/IDEA/samplesDict.json"

# additional/custom C++ functions
includePaths = ["../../functions/functions.h", "../../functions/functions_gen.h", "../h_zh/utils.h", "utils.h"]


# output directory
if treemaker:
    outputDir   = f"/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/treemaker/ecm{ecm}/stage2_ww/"
else:
    outputDir   = f"output/h_zz_ww/histmaker/ecm{ecm}/stage2_ww/"

# optional: ncpus, default is 4, -1 uses all cores available
nCPUS       = 32 ## good for Z/g -- 21 seconds fraction = 0.01 (100 files per DS)
nCPUS       = 8 ## WW(40 files) ZZ(20 files)
nCPUS       = 64 ## WW(40 files) ZZ(20 files)

# scale the histograms with the cross-section and integrated luminosity
doScale = True
intLumi = 10.8e6 if ecm == 240 else 3e6

# define histograms

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
bins_score_sum = (200, 0, 2)


ROOT.EnableImplicitMT(nCPUS) # hack to deal correctly with TMVAHelperXGB  # bdt_model_new_0p1_WZqq bdt_model_new_0p1 bdt_model_0p1_inv bdt_model_0p1
tmva_helper = TMVAHelperXGB("/ceph/submit/data/group/fcc/ee/analyses/h_zz_ww/training/ww_incl.root", "ww_incl") # read the XGBoost training







def build_graph_hwwzz(df, dataset):

    hists = []

    df = df.Define("ecm", "240" if ecm == 240 else "365")
    df = df.Define("weight", "1.0")
    weightsum = df.Sum("weight")

    #if "HZZ" in dataset:
    #    df = df.Filter("is_hadronic")
    
    #df = df.Filter("visibleEnergy > 200")

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


    df = df.Define("init_tlv", "TLorentzVector ret; ret.SetPxPyPzE(0, 0, 0, ecm); return ret;")
    df = df.Define("jets_tlv", "FCCAnalyses::makeLorentzVectors(jets_px, jets_py, jets_pz, jets_e)")

    for i in range(0, 6):
        df = df.Define(f"jets_px_{i}", f"jets_px[{i}]")
        df = df.Define(f"jets_py_{i}", f"jets_py[{i}]")
        df = df.Define(f"jets_pz_{i}", f"jets_pz[{i}]")
        
        df = df.Define(f"jets_phi_{i}", f"jets_tlv[{i}].Phi()")
        df = df.Define(f"jets_theta_{i}", f"jets_tlv[{i}].Theta()")

        df = df.Define(f"recojet_isB_{i}", f"recojet_isB[{i}]")
        df = df.Define(f"recojet_isC_{i}", f"recojet_isC[{i}]")
        df = df.Define(f"recojet_isS_{i}", f"recojet_isS[{i}]")
        df = df.Define(f"recojet_isU_{i}", f"recojet_isU[{i}]")
        df = df.Define(f"recojet_isD_{i}", f"recojet_isD[{i}]")
        df = df.Define(f"recojet_isTAU_{i}", f"recojet_isTAU[{i}]")







    # jetFlavourHelper adds new columns for each jet flavor (recojet_isB/C/S/...)
    # each column is a vector of probabilities per jet to be that specific flavor
    hists.append(df.Histo1D(("recojet_isB", "", *bins_score), "recojet_isB"))
    hists.append(df.Histo1D(("recojet_isC", "", *bins_score), "recojet_isC"))
    hists.append(df.Histo1D(("recojet_isS", "", *bins_score), "recojet_isS"))
    hists.append(df.Histo1D(("recojet_isG", "", *bins_score), "recojet_isG"))
    hists.append(df.Histo1D(("recojet_isU", "", *bins_score), "recojet_isU"))
    hists.append(df.Histo1D(("recojet_isD", "", *bins_score), "recojet_isD"))
    hists.append(df.Histo1D(("recojet_isTAU", "", *bins_score), "recojet_isTAU"))




    df = df.Define("paired_jets_WW", "FCCAnalyses::pairing_WW_ZZ(jets_tlv, 80.385)")


    df = df.Define("W1_WW", "jets_tlv[paired_jets_WW[0]] + jets_tlv[paired_jets_WW[1]]")
    df = df.Define("W2_WW", "jets_tlv[paired_jets_WW[2]] + jets_tlv[paired_jets_WW[3]]")
    df = df.Define("Z_WW", "jets_tlv[paired_jets_WW[4]] + jets_tlv[paired_jets_WW[5]]")
    df = df.Define("H_WW", "W1_WW+W2_WW")





    df = df.Define("W1_WW_jet1_p", "jets_tlv[paired_jets_WW[0]].P()")
    df = df.Define("W1_WW_jet2_p", "jets_tlv[paired_jets_WW[1]].P()")
    df = df.Define("W2_WW_jet1_p", "jets_tlv[paired_jets_WW[2]].P()")
    df = df.Define("W2_WW_jet2_p", "jets_tlv[paired_jets_WW[3]].P()")
    df = df.Define("W1_WW_jet1_theta", "jets_tlv[paired_jets_WW[0]].Theta()")
    df = df.Define("W1_WW_jet2_theta", "jets_tlv[paired_jets_WW[1]].Theta()")
    df = df.Define("W2_WW_jet1_theta", "jets_tlv[paired_jets_WW[2]].Theta()")
    df = df.Define("W2_WW_jet2_theta", "jets_tlv[paired_jets_WW[3]].Theta()")

    df = df.Define("W1_WW_jet1_isB", "recojet_isB[paired_jets_WW[0]]")
    df = df.Define("W1_WW_jet2_isB", "recojet_isB[paired_jets_WW[1]]")
    df = df.Define("W1_WW_jet1_isC", "recojet_isC[paired_jets_WW[0]]")
    df = df.Define("W1_WW_jet2_isC", "recojet_isC[paired_jets_WW[1]]")
    df = df.Define("W1_WW_jet1_isS", "recojet_isS[paired_jets_WW[0]]")
    df = df.Define("W1_WW_jet2_isS", "recojet_isS[paired_jets_WW[1]]")
    df = df.Define("W1_WW_jet1_isU", "recojet_isU[paired_jets_WW[0]]")
    df = df.Define("W1_WW_jet2_isU", "recojet_isU[paired_jets_WW[1]]")
    df = df.Define("W1_WW_jet1_isD", "recojet_isD[paired_jets_WW[0]]")
    df = df.Define("W1_WW_jet2_isD", "recojet_isD[paired_jets_WW[1]]")
    df = df.Define("W1_WW_jet1_isTAU", "recojet_isTAU[paired_jets_WW[0]]")
    df = df.Define("W1_WW_jet2_isTAU", "recojet_isTAU[paired_jets_WW[1]]")

    df = df.Define("W2_WW_jet1_isB", "recojet_isB[paired_jets_WW[2]]")
    df = df.Define("W2_WW_jet2_isB", "recojet_isB[paired_jets_WW[3]]")
    df = df.Define("W2_WW_jet1_isC", "recojet_isC[paired_jets_WW[2]]")
    df = df.Define("W2_WW_jet2_isC", "recojet_isC[paired_jets_WW[3]]")
    df = df.Define("W2_WW_jet1_isS", "recojet_isS[paired_jets_WW[2]]")
    df = df.Define("W2_WW_jet2_isS", "recojet_isS[paired_jets_WW[3]]")
    df = df.Define("W2_WW_jet1_isU", "recojet_isU[paired_jets_WW[2]]")
    df = df.Define("W2_WW_jet2_isU", "recojet_isU[paired_jets_WW[3]]")
    df = df.Define("W2_WW_jet1_isD", "recojet_isD[paired_jets_WW[2]]")
    df = df.Define("W2_WW_jet2_isD", "recojet_isD[paired_jets_WW[3]]")
    df = df.Define("W2_WW_jet1_isTAU", "recojet_isTAU[paired_jets_WW[2]]")
    df = df.Define("W2_WW_jet2_isTAU", "recojet_isTAU[paired_jets_WW[3]]")


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


    hists.append(df.Histo2D(("H_WW_m_Z_WW_m", "", *(bins_m+bins_m)), "H_WW_m", "Z_WW_m"))


    # sum of the scores for the Vs

    df = df.Define("sum_score_B_W1", "recojet_isB[paired_jets_WW[0]] + recojet_isB[paired_jets_WW[1]]") # onshell V
    df = df.Define("sum_score_B_W2", "recojet_isB[paired_jets_WW[2]] + recojet_isB[paired_jets_WW[3]]") # offshell V

    df = df.Define("sum_score_C_W1", "recojet_isC[paired_jets_WW[0]] + recojet_isC[paired_jets_WW[1]]") # onshell V
    df = df.Define("sum_score_C_W2", "recojet_isC[paired_jets_WW[2]] + recojet_isC[paired_jets_WW[3]]") # offshell V

    df = df.Define("sum_score_S_W1", "recojet_isS[paired_jets_WW[0]] + recojet_isS[paired_jets_WW[1]]") # onshell V
    df = df.Define("sum_score_S_W2", "recojet_isS[paired_jets_WW[2]] + recojet_isS[paired_jets_WW[3]]") # offshell V

    df = df.Define("sum_score_U_W1", "recojet_isU[paired_jets_WW[0]] + recojet_isU[paired_jets_WW[1]]") # onshell V
    df = df.Define("sum_score_U_W2", "recojet_isU[paired_jets_WW[2]] + recojet_isU[paired_jets_WW[3]]") # offshell V

    df = df.Define("sum_score_D_W1", "recojet_isD[paired_jets_WW[0]] + recojet_isD[paired_jets_WW[1]]") # onshell V
    df = df.Define("sum_score_D_W2", "recojet_isD[paired_jets_WW[2]] + recojet_isD[paired_jets_WW[3]]") # offshell V

    df = df.Define("sum_score_TAU_W1", "recojet_isTAU[paired_jets_WW[0]] + recojet_isTAU[paired_jets_WW[1]]") # onshell V
    df = df.Define("sum_score_TAU_W2", "recojet_isTAU[paired_jets_WW[2]] + recojet_isTAU[paired_jets_WW[3]]") # offshell V

    hists.append(df.Histo1D(("sum_score_B_W1", "", *bins_score_sum), "sum_score_B_W1"))
    hists.append(df.Histo1D(("sum_score_B_W2", "", *bins_score_sum), "sum_score_B_W2"))
    hists.append(df.Histo1D(("sum_score_C_W1", "", *bins_score_sum), "sum_score_C_W1"))
    hists.append(df.Histo1D(("sum_score_C_W2", "", *bins_score_sum), "sum_score_C_W2"))
    hists.append(df.Histo1D(("sum_score_S_W1", "", *bins_score_sum), "sum_score_S_W1"))
    hists.append(df.Histo1D(("sum_score_S_W2", "", *bins_score_sum), "sum_score_S_W2"))
    hists.append(df.Histo1D(("sum_score_U_W1", "", *bins_score_sum), "sum_score_U_W1"))
    hists.append(df.Histo1D(("sum_score_U_W2", "", *bins_score_sum), "sum_score_U_W2"))
    hists.append(df.Histo1D(("sum_score_D_W1", "", *bins_score_sum), "sum_score_D_W1"))
    hists.append(df.Histo1D(("sum_score_D_W2", "", *bins_score_sum), "sum_score_D_W2"))
    hists.append(df.Histo1D(("sum_score_TAU_W1", "", *bins_score_sum), "sum_score_TAU_W1"))
    hists.append(df.Histo1D(("sum_score_TAU_W2", "", *bins_score_sum), "sum_score_TAU_W2"))

    df = tmva_helper.run_inference(df, col_name="mva_score")
    hists.append(df.Histo1D(("mva_score", "", *(1000, 0, 1)), "mva_score"))


    bins_mva_ =  list(np.arange(0, 1, 0.01)) + [1]
    bins_m_ = list(range(0, 401, 1))
    _bins_mva = array.array('d', bins_mva_)
    _bins_m = array.array('d', bins_m_)
    model = ROOT.RDF.TH2DModel("mva_score_1d", "", len(bins_m_)-1, _bins_m, len(bins_mva_)-1, _bins_mva)
    hists.append(df.Histo2D(model, "Z_WW_m", "mva_score"))


    bins_mva_ = [0, 0.22, 1]
    bins_m_ = list(range(0, 401, 1))
    _bins_mva = array.array('d', bins_mva_)
    _bins_m = array.array('d', bins_m_)
    model = ROOT.RDF.TH3DModel("mva_score_2d", "", len(bins_m_)-1, _bins_m, len(bins_m_)-1, _bins_m, len(bins_mva_)-1, _bins_mva)
    hists.append(df.Histo3D(model, "H_WW_m", "Z_WW_m", "mva_score"))



    return hists, weightsum, df










if treemaker:
    class RDFanalysis():
        def analysers(df):
            hists, weightsum, df = build_graph_hwwzz(df, "")
            return df

        # define output branches to be saved
        def output():

            vars_kinematics = ["W1_WW_m", "W2_WW_m", "Z_WW_m", "Z_WW_p", "Z_WW_theta", "H_WW_m", "H_WW_p", "H_WW_theta"]
            
            vars_jets = ["W1_WW_jet1_p", "W1_WW_jet2_p", "W2_WW_jet1_p", "W2_WW_jet2_p", "W1_WW_jet1_theta", "W1_WW_jet2_theta", "W2_WW_jet1_theta", "W2_WW_jet2_theta"]

            vars_flavor = ["W1_WW_jet1_isB", "W1_WW_jet2_isB", "W1_WW_jet1_isC", "W1_WW_jet2_isC", "W1_WW_jet1_isS", "W1_WW_jet2_isS", "W1_WW_jet1_isU", "W1_WW_jet2_isU", "W1_WW_jet1_isD", "W1_WW_jet2_isD", "W1_WW_jet1_isTAU", "W1_WW_jet2_isTAU", "W2_WW_jet1_isB", "W2_WW_jet2_isB", "W2_WW_jet1_isC", "W2_WW_jet2_isC", "W2_WW_jet1_isS", "W2_WW_jet2_isS", "W2_WW_jet1_isU", "W2_WW_jet2_isU", "W2_WW_jet1_isD", "W2_WW_jet2_isD", "W2_WW_jet1_isTAU", "W2_WW_jet2_isTAU"]

            vars_zz_ww = vars_kinematics + vars_jets + vars_flavor




            
            #branchList = ["jets_e", "jets_px", "jets_py", "jets_pz", "jets_m", "visibleEnergy", "visibleMass", "missingEnergy_p", "cosThetaMiss", "paired_jets_WW", "paired_jets_ZZ", "dmerge_01", "dmerge_12", "dmerge_23", "dmerge_34", "dmerge_45", "dmerge_56", "dmerge_67", "sq_dmerge_01", "sq_dmerge_12", "sq_dmerge_23", "sq_dmerge_34", "sq_dmerge_45", "sq_dmerge_56", "sq_dmerge_67", "recojet_isB", "recojet_isC", "recojet_isS", "recojet_isG", "recojet_isU", "recojet_isD", "recojet_isTAU"]
            return vars_zz_ww



else:
    def build_graph(df, dataset):
        hists, weightsum, df = build_graph_hwwzz(df, dataset)
        return hists, weightsum
