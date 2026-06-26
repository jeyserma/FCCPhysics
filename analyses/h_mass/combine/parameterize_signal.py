
import sys,copy,array,os,subprocess
import ROOT
import numpy as np
import argparse

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../../python')
import plotter
import functions

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)



parser = argparse.ArgumentParser()
parser.add_argument("--flavor", type=str, help="Flavor (mumu or ee)", default="mumu")
parser.add_argument("--mode", type=str, help="Detector mode", choices=["IDEA", "IDEA_MC", "IDEA_3T", "CLD", "CLD_FullSim", "IDEA_noBES", "IDEA_2E", "IDEA_BES6pct"], default="IDEA")
parser.add_argument("--cat", type=str, help="Category (0, 1, 2 or 3)", choices=["0", "1", "2", "3"], default="0")
parser.add_argument("--ecm", type=str, help="Center-of-mass", choices=["240", "365"], default="240")
parser.add_argument("--tag", type=str, help="Analysis tag for versioning, optional", default="")
parser.add_argument("--singularity", action='store_true', help="Run in singularity mode")
args = parser.parse_args()



def getHist(hName, procs):
    hist = None
    for proc in procs:
        fInName = f"{inputDir}/{proc}.root"
        if os.path.exists(fInName):
            fIn = ROOT.TFile(fInName)
        else:
            print(f"ERROR: input file {fInName} not found")
            quit()
        fIn.ls()
        h = fIn.Get(hName)
        h.SetDirectory(0)
        if hist == None:
            hist = h
        else:
            hist.Add(h)
        fIn.Close()
    return hist

def doSignal(normYields = True):

    global h_obs

    mHs = [124.9, 124.95, 125.0, 125.05, 125.1]
    mHs = [124.95, 125.0, 125.05]
    #mHs = [124.975, 125.0, 125.025]
    #mHs = [124.975, 125.0, 125.025]
    #mHs = [124.95, 124.975, 125.0, 125.025, 125.05]
    if flavor == "mumu":
        procs = ["p_wzp6_ee_mumuH_mH-lower-100MeV_ecm240", "p_wzp6_ee_mumuH_mH-lower-50MeV_ecm240", "p_wzp6_ee_mumuH_ecm240", "p_wzp6_ee_mumuH_mH-higher-50MeV_ecm240", "p_wzp6_ee_mumuH_mH-higher-100MeV_ecm240"]
        if args.ecm == "240":
            procs = ["wzp6_ee_mumuH_mH-lower-50MeV_ecm240", "wzp6_ee_mumuH_ecm240", "wzp6_ee_mumuH_mH-higher-50MeV_ecm240"]
        else:
            procs = ["wz3p6_ee_mumuH_mH-lower-50MeV_ecm365", "wz3p6_ee_mumuH_ecm365", "wz3p6_ee_mumuH_mH-higher-50MeV_ecm365"]
        #procs = ["wz3p6_ee_mumuH_mH-lower-25MeV_ecm240", "wzp6_ee_mumuH_ecm240", "wz3p6_ee_mumuH_mH-higher-25MeV_ecm240"]
        #procs = ["wzp6_ee_mumuH_mH-lower-50MeV_ecm240", "wz3p6_ee_mumuH_mH-lower-25MeV_ecm240", "wzp6_ee_mumuH_ecm240", "wz3p6_ee_mumuH_mH-higher-25MeV_ecm240", "wzp6_ee_mumuH_mH-higher-50MeV_ecm240"]
    if flavor == "ee":
        procs = ["p_wzp6_ee_eeH_mH-lower-100MeV_ecm240", "p_wzp6_ee_eeH_mH-lower-50MeV_ecm240", "p_wzp6_ee_eeH_ecm240", "p_wzp6_ee_eeH_mH-higher-50MeV_ecm240", "p_wzp6_ee_eeH_mH-higher-100MeV_ecm240"]
        procs = ["wzp6_ee_eeH_mH-lower-50MeV_ecm240", "wzp6_ee_eeH_ecm240", "wzp6_ee_eeH_mH-higher-50MeV_ecm240"]
        #procs = ["wz3p6_ee_eeH_mH-lower-25MeV_ecm240", "wzp6_ee_eeH_ecm240", "wz3p6_ee_eeH_mH-higher-25MeV_ecm240"]
        if args.ecm == "240":
            procs = ["wzp6_ee_eeH_mH-lower-50MeV_ecm240", "wzp6_ee_eeH_ecm240", "wzp6_ee_eeH_mH-higher-50MeV_ecm240"]
        else:
            procs = ["wz3p6_ee_eeH_mH-lower-50MeV_ecm365", "wz3p6_ee_eeH_ecm365", "wz3p6_ee_eeH_mH-higher-50MeV_ecm365"]
    recoilmass = w_tmp.var("zll_recoil_m")
    MH = w_tmp.var("MH")

    param_yield, param_mh, param_mean, param_mean_offset, param_mean_gt, param_mean_gt_offset, param_sigma, param_sigma_gt, param_alpha_1, param_alpha_2, param_n_1, param_n_2, param_cb_1, param_cb_2 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
    param_yield_err, param_mean_err, param_mean_offset_err, param_sigma_err, param_mean_gt_err, param_mean_gt_offset_err, param_sigma_gt_err, param_alpha_1_err, param_alpha_2_err, param_n_1_err, param_n_2_err, param_cb_1_err, param_cb_2_err  = [], [], [], [], [], [], [], [], [], [], [], [], []

    hist_norm = getHist(f"{flavor}_{hName}", [procs[1]])
    hist_norm.Scale(lumiScale)
    hist_norm = hist_norm.ProjectionX("hist_zh_norm", cat_idx_min, cat_idx_max)
    yield_norm = hist_norm.Integral()

    tmp = hist_norm.Clone()
    print(hist_norm.GetNbinsX() / nBins)
    tmp = tmp.Rebin(int(hist_norm.GetNbinsX() / nBins))
    yMax = tmp.GetMaximum()

    # recoil mass plot settings
    cfg = {
 
        'logy'              : False,
        'logx'              : False,

        'xmin'              : 120,
        'xmax'              : 140,
        'ymin'              : 0,
        'ymax'              : yMax,

        'xtitle'            : "m_{rec} (GeV)",
        'ytitle'            : "Events",

        'topRight'          : topRight,
        'topLeft'           : topLeft,

        'ratiofraction'     : 0.3,
        'ytitleR'           : "Pull",
        'yminR'             : -3.5,
        'ymaxR'             : 3.5,
    }



    for i, proc in enumerate(procs):

        if mode == "IDEA_3T":
            proc += "_3T"
        if mode == "CLD":
            proc += "_CLD"
        if mode == "CLD_FullSim":
            proc += "_CLD_FullSim"
            proc = proc.replace("-50MeV", "").replace("mumuH_mH", "mumuH-mH")
        if mode == "IDEA_noBES":
            proc = proc.replace("_ecm240", "_noBES_ecm240")
        if mode == "IDEA_2E" and flavor == "ee":
            proc += "_E2"

        mH = mHs[i]
        mH_ = ("%.3f" % mH).replace(".", "p")
        print("Do mH=%.3f" % mH)

        hist_zh = getHist(f"{flavor}_{hName}", [proc])
        hist_zh.Scale(lumiScale)
        hist_zh = hist_zh.ProjectionX("hist_zh_%s" % mH_, cat_idx_min, cat_idx_max)
        if normYields: hist_zh.Scale(yield_norm/hist_zh.Integral())
        rdh_zh = ROOT.RooDataHist("rdh_zh_%s" % mH_, "rdh_zh", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(hist_zh))
        yield_zh = rdh_zh.sum(False)
        if mH == 125.0 and h_obs == None: h_obs = hist_zh.Clone("h_obs") # take 125.0 GeV to add to observed (need to add background later as well)


        ### fit parameter configuration of 2CBG
        # the gt mean is an offset w.r.t. the CB means (=mean_gt_offset)

        # IDEA
        if args.ecm == "240":
            if cat == 0 and flavor == 'mumu' and (mode == 'IDEA' or mode == 'IDEA_2E' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.818, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)

            if cat == 0 and flavor == 'mumu' and mode == 'IDEA_3T':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.818, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)

            if cat == 0 and flavor == 'mumu' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.818, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)

            if cat == 0 and flavor == 'mumu' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.525, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.26, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.85, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.59, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.55, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.57, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.42, 0, 1)
        
            if cat == 0 and flavor == 'mumu' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.04, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.2, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.37, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.15, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.15, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.8, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.24, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.52, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.39, 0, 1)


            if cat == 0 and flavor == 'ee' and (mode == 'IDEA' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 0 and flavor == 'ee' and mode == 'IDEA_3T':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 0 and flavor == 'ee' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 0 and flavor == 'ee' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36244, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.565, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.16025, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.765, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.716, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 33, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.64, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5714, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.375, 0, 1)

            if cat == 0 and flavor == 'ee' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.04, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.2, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.37, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.15, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.15, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.8, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.24, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.52, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.39, 0, 1)

            if cat == 0 and flavor == 'ee' and mode == 'IDEA_2E':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)



            if cat == 1 and flavor == 'mumu' and (mode == 'IDEA' or mode == 'IDEA_2E' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0930, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4464, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.860, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1733, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.37, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 4.02, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.369, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.450, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.434, 0, 1)

            if cat == 1 and flavor == 'mumu' and mode == 'IDEA_3T':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.5, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)

            if cat == 1 and flavor == 'mumu' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.818, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)
       
            if cat == 1 and flavor == 'mumu' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.525, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.26, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.85, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.59, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.55, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.57, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.42, 0, 1)
                
            if cat == 1 and flavor == 'mumu' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.28, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.525, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.26, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 0.987, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.59, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 12, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.29, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.57, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)

            if cat == 1 and flavor == 'ee' and (mode == 'IDEA' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.1121, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4344, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.688, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1788, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.74, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.90, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.33, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.589, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.570, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.369, 0, 1)

            if cat == 1 and flavor == 'ee' and mode == 'IDEA_3T':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 1 and flavor == 'ee' and mode == 'IDEA_2E':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 1 and flavor == 'ee' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 1 and flavor == 'ee' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36244, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.565, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.16025, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.765, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.716, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 33, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.64, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5714, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.375, 0, 1)

            if cat == 1 and flavor == 'ee' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.04, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.2, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.37, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.15, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.15, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.8, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.24, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.52, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.39, 0, 1)



            if cat == 2 and flavor == 'mumu' and (mode == 'IDEA' or mode == 'IDEA_2E' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0886, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4170, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.670, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.21988, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.96, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.242, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.26, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.66, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5013, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.26, 0, 1)

            if cat == 2 and flavor == 'mumu' and mode == 'IDEA_3T':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.5, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)

            if cat == 2 and flavor == 'mumu' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.818, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)

            if cat == 2 and flavor == 'mumu' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.525, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.26, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.85, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.59, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.55, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.57, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.42, 0, 1)

            if cat == 2 and flavor == 'mumu' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.04, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.2, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.37, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.15, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.15, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.8, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.24, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.52, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.39, 0, 1)

            if cat == 2 and flavor == 'ee' and (mode == 'IDEA' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 2 and flavor == 'ee' and mode == 'IDEA_3T':
                #mean = ROOT.RooRealVar("mean_%s" % mH_, '', 1.25090e+02, mH-1., mH+1.)
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0876, -2, 2)
                #mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.1256)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))
                #sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 4.08196e-01, 0, 1)
                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.3917) # fixed
                #alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -2.00592e-01, -10, 0)
                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1435)
                #alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.05919e+00, 0, 10)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.066)
                #n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.62, -10, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 4.66)
                #n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.25675e-02, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.0028)

                #mean_gt = ROOT.RooRealVar("mean_gt_%s" % mH_, '', 1.25338e+02, recoilMin, recoilMax)
                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.5)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))
                
                #sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 8.30603e-01, 0, 2)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.735) # fixed  
                
                #cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 4.94921e-01 , 0, 1)
                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5448)
                #cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 3.86757e-01 , 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.3780)

            if cat == 2 and flavor == 'ee' and mode == 'IDEA_2E':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 2 and flavor == 'ee' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 2 and flavor == 'ee' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36244, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.565, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.16025, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.765, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.716, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 33, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.64, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5714, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.375, 0, 1)

            if cat == 2 and flavor == 'ee' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.04, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.2, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.37, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.15, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.15, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.8, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.24, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.52, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.39, 0, 1)





            if cat == 3 and flavor == 'mumu' and (mode == 'IDEA' or mode == 'IDEA_2E' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0886, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4170, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.670, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.21988, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.96, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.242, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.26, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.66, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5013, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.26, 0, 1)

            if cat == 3 and flavor == 'mumu' and mode == 'IDEA_3T':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 3 and flavor == 'mumu' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.0919, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.4338, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.818, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.2074, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.13, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.55, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 1.39, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.334, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.4861, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.4132, 0, 1)

            if cat == 3 and flavor == 'mumu' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.525, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.26, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.85, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.59, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.55, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.57, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.42, 0, 1)

            if cat == 3 and flavor == 'mumu' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.28, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.525, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.26, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 0.987, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.59, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 12, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.29, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.57, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)

            if cat == 3 and flavor == 'ee' and (mode == 'IDEA' or mode == 'IDEA_BES6pct'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.153, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.52, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 1.01, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.19, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 4.133, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.48, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.00013, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.71, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.562, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.368, 0, 1)

            if cat == 3 and flavor == 'ee' and mode == 'IDEA_3T':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 3 and flavor == 'ee' and (mode == 'CLD' or mode == 'CLD_FullSim'):
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 3 and flavor == 'ee' and mode == 'IDEA_MC':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.36244, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.565, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.16025, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.765, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 2.716, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 33, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.64, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5714, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.375, 0, 1)

            if cat == 3 and flavor == 'ee' and mode == 'IDEA_2E':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.126, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.46, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.832, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.1721, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.9, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 3.38, -10, 10)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 0.1, -10, 10)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.55, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.556, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.372, 0, 1)

            if cat == 3 and flavor == 'ee' and mode == 'IDEA_noBES':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.05, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 0.28, 0, 1)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 0.525, 0, 2)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.26, -10, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 0.987, 0, 10)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.59, -20, 20)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 12, -20, 20)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.29, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.57, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)

        if args.ecm == "365":
            if cat == 0 and flavor == 'mumu' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.3, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.5, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 1.6, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.6, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 1.24, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.32, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 10, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.47, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.5, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.25, 0, 1)


            if cat == 0 and flavor == 'ee' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.3, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.5, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 1.6, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.6, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 1.24, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.32, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 10, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.47, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.64, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.24, 0, 1)







            if cat == 1 and flavor == 'mumu' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.4, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.4, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 2.0, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.7, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 3.8, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.2, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 15, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.7, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.7, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)



            if cat == 1 and flavor == 'ee' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.3, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.55, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 1.4, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.7, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 1.5, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.3, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 45, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.7, -2, 2)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.7, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)



            if cat == 2 and flavor == 'mumu' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.4, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.6, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 2.8, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.4, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 1.32, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.4, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 30, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.7, -5, 5)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.7, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)





            if cat == 2 and flavor == 'ee' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.4, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.6, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 2.8, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.4, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 1.32, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.4, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 30, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.7, -5, 5)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.7, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)







            if cat == 3 and flavor == 'mumu' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.4, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.6, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 2.8, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.4, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 1.32, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.4, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 30, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.7, -5, 5)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.7, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.28, 0, 1)


            if cat == 3 and flavor == 'ee' and mode == 'IDEA':
                mean_slope = ROOT.RooRealVar("mean_slope_%s" % mH_, '', mH)
                mean_offset = ROOT.RooRealVar("mean_offset_%s" % mH_, '', 0.4, -2, 2)
                mean = ROOT.RooFormulaVar("mean_%s" % mH_, "@0+@1", ROOT.RooArgList(mean_slope, mean_offset))

                sigma = ROOT.RooRealVar("sigma_%s" % mH_, '', 1.6, 0, 5)
                sigma_gt = ROOT.RooRealVar("sigma_gt_%s" % mH_, '', 2.8, 0, 5)

                alpha_1 = ROOT.RooRealVar("alpha_1_%s" % mH_, '', -0.4, -20, 0)
                alpha_2 = ROOT.RooRealVar("alpha_2_%s" % mH_, '', 1.32, 0, 20)
                n_1 = ROOT.RooRealVar("n_1_%s" % mH_, '', 1.4, -50, 50)
                n_2 = ROOT.RooRealVar("n_2_%s" % mH_, '', 30, -50, 50)

                mean_gt_offset = ROOT.RooRealVar("mean_gt_offset_%s" % mH_, '', 0.7, -5, 5)
                mean_gt = ROOT.RooFormulaVar("mean_gt_%s" % mH_, "@0+@1", ROOT.RooArgList(mean, mean_gt_offset))

                cb_1 = ROOT.RooRealVar("cb_1_%s" % mH_, '', 0.52, 0, 1)
                cb_2 = ROOT.RooRealVar("cb_2_%s" % mH_, '', 0.36, 0, 1)


        # construct the 2CBG and perform the fit: pdf = cb_1*cbs_1 + cb_2*cbs_2 + gauss (cb_1 and cb_2 are the fractions, floating)
        cbs_1 = ROOT.RooCBShape("CrystallBall_1_%s" % mH_, "CrystallBall_1", recoilmass, mean, sigma, alpha_1, n_1) # first CrystallBall
        cbs_2 = ROOT.RooCBShape("CrystallBall_2_%s" % mH_, "CrystallBall_2", recoilmass, mean, sigma, alpha_2, n_2) # second CrystallBall
        gauss = ROOT.RooGaussian("gauss_%s" % mH_, "gauss", recoilmass, mean_gt, sigma_gt) # the Gauss

        sig = ROOT.RooAddPdf("sig_%s" % mH_, '', ROOT.RooArgList(cbs_1, cbs_2, gauss), ROOT.RooArgList(cb_1, cb_2)) # half of both CB functions
        sig_norm = ROOT.RooRealVar("sig_%s_norm" % mH_, '', yield_zh, 0, 1e8) # fix normalization
        sig_fit = ROOT.RooAddPdf("zh_model_%s" % mH_, '', ROOT.RooArgList(sig), ROOT.RooArgList(sig_norm))
        sig_fit.fitTo(rdh_zh, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.SumW2Error(sumw2err))

        cb1__ = cb_1.getVal()
        cb2__ = cb_2.getVal()

        # do plotting
        cfg['ymax'] = 1.2*yMax
        plotter.cfg = cfg

        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(rline=0)
        dummyB.GetXaxis().SetTitleOffset(4.*dummyB.GetXaxis().GetTitleOffset()) # hack label
        dummyT.GetYaxis().SetTitleOffset(1.2*dummyT.GetYaxis().GetTitleOffset()) # hack label

        ## TOP PAD ##
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        dummyT.Draw("HIST")

        plt = recoilmass.frame()
        plt.SetTitle("ZH signal")
        rdh_zh.plotOn(plt, ROOT.RooFit.Binning(nBins)) # , ROOT.RooFit.Normalization(yield_zh, ROOT.RooAbsReal.NumEvent)
        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
        chisq = plt.chiSquare()
        sig_fit.paramOn(plt, ROOT.RooFit.Format("NELU", ROOT.RooFit.AutoPrecision(2)), ROOT.RooFit.Layout(0.45, 0.9, 0.9))
        histpull = plt.pullHist()
        plt.Draw("SAME")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatex(0.2, 0.88, label)
        latex.DrawLatex(0.2, 0.82, "#chi^{2} = %.3f" % chisq)
 
        plotter.auxRatio()
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()  

        ## BOTTOM PAD ##
        canvas.cd()
        padB.Draw()
        padB.SetFillStyle(0)
        padB.cd()
        dummyB.Draw("HIST")
        dummyL.Draw("SAME")

        plt = recoilmass.frame()
        plt.addPlotable(histpull, "P")
        plt.Draw("SAME")

        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()
        canvas.SaveAs("%s/fit_mH%s.png" % (outDir, mH_))
        canvas.SaveAs("%s/fit_mH%s.pdf" % (outDir, mH_))

        del dummyB
        del dummyT
        del padT
        del padB
        del canvas

        cfg['ymax'] = 3.*yMax
        #cfg['ymax'] = 10000
        plotter.cfg = cfg
        canvas = plotter.canvas()
        canvas.SetGrid()
        dummy = plotter.dummy()
        dummy.Draw("HIST")
        plt = w_tmp.var("zll_recoil_m").frame()
        colors = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kGreen, ROOT.kCyan] 

        leg = ROOT.TLegend(.50, 0.7, .95, .90)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.SetMargin(0.15)

        cbs_1.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Normalization(cb1__*yield_zh, ROOT.RooAbsReal.NumEvent))
        cbs_2.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Normalization(cb2__*yield_zh, ROOT.RooAbsReal.NumEvent))
        gauss.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kCyan), ROOT.RooFit.Normalization((1.-cb1__-cb2__)*yield_zh, ROOT.RooAbsReal.NumEvent))
        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.Normalization(yield_zh, ROOT.RooAbsReal.NumEvent))

        # define TGraphs for legend
        tmp1 = ROOT.TGraph()
        tmp1.SetPoint(0, 0, 0)
        tmp1.SetLineColor(ROOT.kBlack)
        tmp1.SetLineWidth(3)
        tmp1.Draw("SAME")
        leg.AddEntry(tmp1, "Total PDF", "L")

        tmp2 = ROOT.TGraph()
        tmp2.SetPoint(0, 0, 0)
        tmp2.SetLineColor(ROOT.kRed)
        tmp2.SetLineWidth(3)
        tmp2.Draw("SAME")
        leg.AddEntry(tmp2, "CB1", "L")

        tmp3 = ROOT.TGraph()
        tmp3.SetPoint(0, 0, 0)
        tmp3.SetLineColor(ROOT.kBlue)
        tmp3.SetLineWidth(3)
        tmp3.Draw("SAME")
        leg.AddEntry(tmp3, "CB2", "L")

        tmp4 = ROOT.TGraph()
        tmp4.SetPoint(0, 0, 0)
        tmp4.SetLineColor(ROOT.kCyan)
        tmp4.SetLineWidth(3)
        tmp4.Draw("SAME")
        leg.AddEntry(tmp4, "Gauss", "L")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatex(0.2, 0.92, label)

        plt.Draw("SAME")
        leg.Draw()
        plotter.aux()
        canvas.Modify()
        canvas.Update()
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()
        canvas.Draw()
        canvas.SaveAs("%s/fit_mH%s_decomposition.png" % (outDir, mH_))
        canvas.SaveAs("%s/fit_mH%s_decomposition.pdf" % (outDir, mH_))

        # import
        getattr(w_tmp, 'import')(rdh_zh)
        getattr(w_tmp, 'import')(sig_fit)

        param_mh.append(mH)
        #param_mean.append(mean.getVal())
        param_mean_offset.append(mean_offset.getVal())
        param_sigma.append(sigma.getVal())
        param_mean_gt.append(mean_gt.getVal())
        param_mean_gt_offset.append(mean_gt_offset.getVal())
        param_sigma_gt.append(sigma_gt.getVal())
        param_alpha_1.append(alpha_1.getVal())
        param_alpha_2.append(alpha_2.getVal())
        param_n_1.append(n_1.getVal())
        param_n_2.append(n_2.getVal())
        param_yield.append(sig_norm.getVal())
        param_cb_1.append(cb_1.getVal())
        param_cb_2.append(cb_2.getVal())

        #param_mean_err.append(mean.getError())
        param_mean_offset_err.append(mean_offset.getError())
        param_sigma_err.append(sigma.getError())
        param_mean_gt_err.append(0)
        param_mean_gt_offset_err.append(mean_gt_offset.getError())
        param_sigma_gt_err.append(sigma.getError())
        param_alpha_1_err.append(alpha_1.getError())
        param_alpha_2_err.append(alpha_2.getError())
        param_n_1_err.append(n_1.getError())
        param_n_2_err.append(n_2.getError())
        param_yield_err.append(sig_norm.getError())
        param_cb_1_err.append(cb_1.getError())
        param_cb_2_err.append(cb_2.getError())


    ##################################
    # plot all fitted signals
    ##################################
    cfg['xmin'] = 124
    cfg['xmax'] = 130
    cfg['ymax'] = yMax*2.5
    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")

    plt = w_tmp.var("zll_recoil_m").frame()
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kGreen, ROOT.kCyan]
    for i, mH in enumerate(mHs):
        mH_ = ("%.3f" % mH).replace(".", "p")
        sig_fit = w_tmp.pdf("zh_model_%s" % mH_)
        # need to re-normalize the pdf, as the pdf is normalized to 1
        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[i]), ROOT.RooFit.Normalization(yield_zh, ROOT.RooAbsReal.NumEvent))


    plt.Draw("SAME")

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(13)
    latex.DrawLatex(0.2, 0.92, label)


    plotter.aux()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs("%s/fit_all.png" % (outDir))
    canvas.SaveAs("%s/fit_all.pdf" % (outDir))


    
    # export values
    coeff_mean = np.polyfit(param_mh, param_mean_offset, 1)
    coeff_mean_gt = np.polyfit(param_mh, param_mean_gt, 1)
    fOut = open("%s/coeff.txt" % outDir, "w")

    idx = 1 # take values at central mass 125 GeV (not average using np.average(param_mean_offset))
    fOut.write(str("1.0\n")) # slope
    fOut.write(str(param_mean_offset[idx]) + "\n")
    fOut.write(str("1.0\n"))
    fOut.write(str(param_mean_offset[idx]) + "\n")
    fOut.write(str(param_mean_gt_offset[idx]) + "\n")
    fOut.write(str(param_sigma[idx]) + "\n")
    fOut.write(str(param_sigma_gt[idx]) + "\n")
    fOut.write(str(param_alpha_1[idx]) + "\n")
    fOut.write(str(param_alpha_2[idx]) + "\n")
    fOut.write(str(param_n_1[idx]) + "\n")
    fOut.write(str(param_n_2[idx]) + "\n")
    fOut.write(str(param_cb_1[idx]) + "\n")
    fOut.write(str(param_cb_2[idx]) + "\n")
    fOut.close()


if __name__ == "__main__":

    if not args.singularity:
        cmd = ' '.join(sys.argv)
        cmd_tot = f"singularity exec --bind /work:/work /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif bash -c 'python {cmd} --singularity'"
        os.system(cmd_tot)
        quit()


    ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
    sumw2err = ROOT.kTRUE

    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Fumili2")
    #ROOT.Math.MinimizerOptions.SetMinimizerAlgorithm("Simplex") # Migrad Minimize Simplex Fumili2
    ROOT.Math.MinimizerOptions.PrintDefault("Minuit2")
    ROOT.Math.MinimizerOptions.SetDefaultPrecision(1e-15)
    ROOT.Math.MinimizerOptions.SetDefaultMaxIterations (200)
    #ROOT.Math.MinimizerOptions.PrintDefault()


    mode = args.mode # detector mode
    flavor = args.flavor
    cat = int(args.cat)
    ecm = args.ecm
    tag = args.tag
    flavorLabel = "#mu^{#plus}#mu^{#minus}" if flavor == "mumu" else "e^{#plus}e^{#minus}"
    lumiScale = 1.0 # the fit coeff do slighly depend on the normalization, so the lumi chosen should be close to the real one

    topRight = f"#sqrt{{s}} = {ecm} GeV, 1 ab^{{#minus1}}"
    topLeft = "#bf{FCC-ee} #scale[0.7]{#it{Internal}}"
    label = f"{flavorLabel}, category {cat}"
    #fIn = ROOT.TFile("output_ZH_mass_%s_ecm%s%s_%s.root"%(flavor, args.ecm, "_gen" if mode == "IDEA_MC" else "", args.tag))
    inputDir = f"output/h_mass/histmaker/ecm{ecm}/"
    outDir = f"/work/submit/jaeyserm/public_html/fccee/h_mass/{tag}/combine/{mode}/base_parametric/{flavor}_cat{cat}_ecm{ecm}/"
    if not os.path.exists(outDir): os.makedirs(outDir)
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/h_mass/index.php {outDir}")

    hName = f"zll_recoil_m_cat"

    if cat == 0: cat_idx_min, cat_idx_max = 0, 5
    else: cat_idx_min, cat_idx_max = cat, cat

    nBins = 250 # total number of bins, for plotting
    recoilMin = 120
    recoilMax = 140
    h_obs = None # should hold the data_obs = sum of signal and backgrounds

    recoilmass = ROOT.RooRealVar("zll_recoil_m", "Recoil mass (GeV)", 125, recoilMin, recoilMax)
    MH = ROOT.RooRealVar("MH", "Higgs mass (GeV)", 125, 124.95, 125.05) # name Higgs mass as MH to be compatible with combine

    # define temporary output workspace
    w_tmp = ROOT.RooWorkspace("w_tmp", "workspace")
    w = ROOT.RooWorkspace("w", "workspace") # final workspace for combine

    getattr(w_tmp, 'import')(recoilmass)
    getattr(w_tmp, 'import')(MH)

    yield_norm = -1
    yMax = -1

    doSignal()

    # delete workspaces to avoid segfault
    del w_tmp
    del w
