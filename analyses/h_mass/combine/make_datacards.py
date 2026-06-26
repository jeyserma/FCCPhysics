
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



def getHist(hName, procs):
    hist = None
    for proc in procs:
        fInName = f"{inputDir}/{proc}.root"
        if os.path.exists(fInName):
            fIn = ROOT.TFile(fInName)
        else:
            print(f"ERROR: input file {fInName} not found")
            quit()
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
    global yield_nom
    global yMax
    
    mHs = [124.9, 124.95, 125.0, 125.05, 125.1]
    mHs = [124.95, 125.0, 125.05]
    if flavor == "mumu":
        procs = ["p_wzp6_ee_mumuH_mH-lower-100MeV_ecm240", "p_wzp6_ee_mumuH_mH-lower-50MeV_ecm240", "p_wzp6_ee_mumuH_ecm240", "p_wzp6_ee_mumuH_mH-higher-50MeV_ecm240", "p_wzp6_ee_mumuH_mH-higher-100MeV_ecm240"]
        if args.ecm == "240":
            procs = ["wzp6_ee_mumuH_mH-lower-50MeV_ecm240", "wzp6_ee_mumuH_ecm240", "wzp6_ee_mumuH_mH-higher-50MeV_ecm240"]
        else:
            procs = ["wz3p6_ee_mumuH_mH-lower-50MeV_ecm365", "wz3p6_ee_mumuH_ecm365", "wz3p6_ee_mumuH_mH-higher-50MeV_ecm365"]
    if flavor == "ee":
        procs = ["wzp6_ee_eeH_mH-lower-100MeV_ecm240", "wzp6_ee_eeH_mH-lower-50MeV_ecm240", "wzp6_ee_eeH_ecm240", "wzp6_ee_eeH_mH-higher-50MeV_ecm240", "wzp6_ee_eeH_mH-higher-100MeV_ecm240"]
        procs = ["wzp6_ee_eeH_mH-lower-50MeV_ecm240", "wzp6_ee_eeH_ecm240", "wzp6_ee_eeH_mH-higher-50MeV_ecm240"]
        if args.ecm == "240":
            procs = ["wzp6_ee_eeH_mH-lower-50MeV_ecm240", "wzp6_ee_eeH_ecm240", "wzp6_ee_eeH_mH-higher-50MeV_ecm240"]
        else:
            procs = ["wz3p6_ee_eeH_mH-lower-50MeV_ecm365", "wz3p6_ee_eeH_ecm365", "wz3p6_ee_eeH_mH-higher-50MeV_ecm365"]

        if args.mode == "CLD_FullSim":
            mHs = [125.0, 125.05]
            procs = ["wzp6_ee_eeH_ecm240", "wzp6_ee_eeH_mH-higher-50MeV_ecm240"]

    recoilmass = w_tmp.var("zll_recoil_m")
    MH = w_tmp.var("MH")

    param_yield, param_mh, param_mean, param_mean_gt, param_mean_offset, param_mean_gt_offset, param_sigma, param_sigma_gt, param_alpha_1, param_alpha_2, param_n_1, param_n_2, param_cb_1, param_cb_2 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
    param_yield_err, param_mean_err, param_mean_offset_err, param_sigma_err, param_mean_gt_err, param_mean_gt_offset_err, param_sigma_gt_err, param_alpha_1_err, param_alpha_2_err, param_n_1_err, param_n_2_err, param_cb_1_err, param_cb_2_err  = [], [], [], [], [], [], [], [], [], [], [], [], []

    hist_norm = getHist(f"{flavor}_{hName}", [procs[1]])
    hist_norm.Scale(lumiScale)
    hist_norm = hist_norm.ProjectionX("hist_zh_norm", cat_idx_min, cat_idx_max)
    yield_nom = hist_norm.Integral()

    tmp = hist_norm.Clone()
    tmp = tmp.Rebin(int(hist_norm.GetNbinsX() / nBins))
    yMax = 1.25*tmp.GetMaximum()

    # recoil mass plot settings
    cfg = {
 
        'logy'              : False,
        'logx'              : False,
    
        'xmin'              : 120,
        'xmax'              : 140,
        'ymin'              : 0,
        'ymax'              : yMax,
        
        'xtitle'            : "m_{rec} (GeV)",
        'ytitle'            : "Events / 0.2 GeV",
        
        'topRight'          : topRight,
        'topLeft'           : topLeft,

        'ratiofraction'     : 0.3,
        'ytitleR'           : "Pull",
        'yminR'             : -3.5,
        'ymaxR'             : 3.5,
    }

    ptrs = {}
    fitFunctions = {}
    garbage = [] # need to store the variables for memory issues

    ## Build model
    ## linear functions for mean and mean_gt
    ## constants for all the rest

    # import values
    coeff = np.loadtxt(f"{outDir}/coeff.txt".replace(f"lumi{lumiLabel}", "base_parametric")) # take the coefficients from the baseline (scaled to 1 ab-1)
    param_mean0_ = float(coeff[0])
    param_mean1_ = float(coeff[1])
    param_mean_gt0_ = float(coeff[2])
    param_mean_gt1_ = float(coeff[3])
    param_mean_gt_offset0_ = float(coeff[4])
    param_sigma_ = float(coeff[5])
    param_sigma_gt_ = float(coeff[6])
    param_alpha_1_ = float(coeff[7])
    param_alpha_2_ = float(coeff[8])
    param_n_1_ = float(coeff[9])
    param_n_2_ = float(coeff[10])
    param_cb_1_ = float(coeff[11])
    param_cb_2_ = float(coeff[12])

    mean_argl, sigma_argl = ROOT.RooArgList("mean_argl"), ROOT.RooArgList("sigma_argl")
    sigma_gt_argl, mean_gt_offset_argl = ROOT.RooArgList("sigma_gt_argl"), ROOT.RooArgList("mean_gt_offset_argl")
    alpha1_argl, n1_argl, cb1_argl = ROOT.RooArgList("alpha1_argl"), ROOT.RooArgList("n1_argl"), ROOT.RooArgList("cb1_argl")
    alpha2_argl, n2_argl, cb2_argl = ROOT.RooArgList("alpha2_argl"), ROOT.RooArgList("n2_argl"), ROOT.RooArgList("cb2_argl")

    mean0 = ROOT.RooRealVar("mean0", "", param_mean0_, 0.5, 1.5) # slope
    mean0.setConstant(ROOT.kTRUE)
    mean1 = ROOT.RooRealVar("mean1", "", param_mean1_, -1, 1) # offset
    mean1.setConstant(ROOT.kTRUE)
    mean_argl.add(mean0)
    mean_argl.add(mean1)

    mean_gt_offset = ROOT.RooRealVar("mean_gt_offset", "", param_mean_gt_offset0_, -2, 2)
    mean_gt_offset.setConstant(ROOT.kTRUE)
    mean_gt_offset_argl.add(mean_gt_offset)

    sigma0 = ROOT.RooRealVar("sigma0", "", param_sigma_, 0, 10) # 0.4335
    sigma0.setConstant(ROOT.kTRUE)
    sigma_argl.add(sigma0)

    sigma_gt0 = ROOT.RooRealVar("sigma_gt0", "", param_sigma_gt_, 0, 10)
    sigma_gt0.setConstant(ROOT.kTRUE)
    sigma_gt_argl.add(sigma_gt0)

    alpha10 = ROOT.RooRealVar("alpha10", "", param_alpha_1_, -5, 5)
    alpha10.setConstant(ROOT.kTRUE)
    alpha1_argl.add(alpha10)
    n10 = ROOT.RooRealVar("n10", "", param_n_1_, -50, 50)
    n10.setConstant(ROOT.kTRUE)
    n1_argl.add(n10)
    cb10 = ROOT.RooRealVar("cb10", "", param_cb_1_, 0, 1)
    cb10.setConstant(ROOT.kTRUE)
    cb1_argl.add(cb10)

    alpha20 = ROOT.RooRealVar("alpha20", "", param_alpha_2_, -5, 5)
    alpha20.setConstant(ROOT.kTRUE)
    alpha2_argl.add(alpha20)
    n20 = ROOT.RooRealVar("n20", "", param_n_2_,  -50, 50)
    n20.setConstant(ROOT.kTRUE)
    n2_argl.add(n20)
    cb20 = ROOT.RooRealVar("cb20", "", param_cb_2_, 0, 1)
    cb20.setConstant(ROOT.kTRUE)
    cb2_argl.add(cb20)

    cats = ROOT.RooCategory("category", "") # for each mass bin, define category
    hists = ROOT.std.map("string, RooDataHist*")() # container holding all RooDataHists
    pdf_tot = ROOT.RooSimultaneous("pdf_tot", "", cats) # total pdf, containing all the categories

    garbage = []
    list_alpha1, list_alpha2, list_n1, list_n2, list_cb1, list_cb2, list_sigma, list_sigma_gt, list_mean, list_mean_offset, list_mean_gt, list_mean_gt_offset, list_norm = [], [], [], [], [], [], [], [], [], [], [], [], []
    for i, proc in enumerate(procs):

        if mode == "IDEA_3T":
            proc += "_3T"
        if mode == "CLD":
            proc += "_CLD"
        if mode == "CLD_FullSim":
            proc += "_CLD_FullSim"
            proc = proc.replace("-50MeV", "").replace("mumuH_mH", "mumuH-mH").replace("eeH_mH", "eeH-mH")
        if mode == "IDEA_noBES":
            proc = proc.replace("_ecm240", "_noBES_ecm240")
        if mode == "IDEA_2E" and flavor == "ee":
            proc += "_E2"

        mH = mHs[i]
        mH_ = (f"{mH:.3f}").replace(".", "p")
        print(f"Do mH={mH:.3f}")

        print(f"{proc}/{hName}")
        hist_zh = getHist(f"{flavor}_{hName}", [proc])
        hist_zh.Scale(lumiScale)
        hist_zh = hist_zh.ProjectionX(f"hist_zh_{mH_}", cat_idx_min, cat_idx_max)
        if normYields: hist_zh.Scale(yield_nom/hist_zh.Integral())

        rdh_zh = ROOT.RooDataHist(f"rdh_zh_{mH_}", "", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(hist_zh))
        rdh_zh.SetName(f"rdh_zh_{mH_}")
        yield_zh = rdh_zh.sum(False)

        catIDx = mH_
        hists.insert(ROOT.std.pair("string, RooDataHist*")(catIDx, rdh_zh)) # does not work with recent ROOT versions?
        cats.defineType(catIDx, i)

        mean = ROOT.RooFormulaVar(f"mean_{mH_}", f"x[1] + x[0]*{mH}", mean_argl)
        garbage.append(mean)
        sigma = ROOT.RooFormulaVar(f"sigma_{mH_}", "x[0]", sigma_argl)

        mean_gt_argl = ROOT.RooArgList("mean_gt_argl")
        mean_gt_argl.add(mean)
        mean_gt_argl.add(mean_gt_offset)
        mean_gt = ROOT.RooFormulaVar(f"mean_gt_{mH_}", f"x[0] + x[1]", mean_gt_argl) # mean gt is as an offset w.r.t. the mean
        sigma_gt = ROOT.RooFormulaVar(f"sigma_gt_{mH_}", "x[0]", sigma_gt_argl)

        alpha1 = ROOT.RooFormulaVar(f"alpha1_{mH_}", "x[0]", alpha1_argl)
        n1 = ROOT.RooFormulaVar(f"n1_{mH_}", "x[0]", n1_argl)
        cb1 = ROOT.RooFormulaVar(f"cb1_{mH_}", "x[0]", cb1_argl)

        alpha2 = ROOT.RooFormulaVar(f"alpha2_{mH_}", "x[0]", alpha2_argl)
        n2 = ROOT.RooFormulaVar(f"n2_{mH_}", "x[0]", n2_argl)
        cb2 = ROOT.RooFormulaVar(f"cb2_{mH_}", "x[0]", cb2_argl)


        # construct the 2CBG pdf = cb_1*cbs_1 + cb_2*cbs_2 + gauss (cb_1 and cb_2 are the fractions, floating)
        cbs1 = ROOT.RooCBShape(f"cbs1_{mH_}", "CrystallBall_1", recoilmass, mean, sigma, alpha1, n1) # first CrystallBall
        cbs2 = ROOT.RooCBShape(f"cbs2_{mH_}", "CrystallBall_2", recoilmass, mean, sigma, alpha2, n2) # second CrystallBall
        gauss = ROOT.RooGaussian(f"gauss_{mH_}", "gauss", recoilmass, mean_gt, sigma_gt) # the Gauss

        argl = ROOT.RooArgList(cbs1, cbs2, gauss)
        argl.setName(f"argl_{mH_}")
        norms_argl = ROOT.RooArgList(cb1, cb2)
        norms_argl.setName(f"norms_argl_{mH_}")
        sig = ROOT.RooAddPdf(f"sig_{mH_}", '', argl, norms_argl) # half of both CB functions
        sig_norm = ROOT.RooRealVar(f"sig_norm_{mH_}", '', yield_zh, 0, 1e8) # fix normalization
        #sig_norm.setConstant(ROOT.kTRUE)

        sig_argl = ROOT.RooArgList(sig)
        sig_argl.setName(f"sig_argl_{mH_}")
        sig_norm_argl = ROOT.RooArgList(sig_norm)
        sig_norm_argl.setName(f"sig_norm_argl_{mH_}")
        pdf_sig = ROOT.RooAddPdf(f"zh_model_{mH_}", '', sig_argl, sig_norm_argl)
        pdf_sigs.append(pdf_sig)

        #getattr(w_tmp, 'import')(rdh_zh)
        #getattr(w_tmp, 'import')(pdf_sig)

        # must store the individual vars for later , to extract the values
        # seems not to work with workspace
        list_alpha1.append(alpha1)
        list_alpha2.append(alpha2)
        list_n1.append(n1)
        list_n2.append(n2)
        list_cb1.append(cb1)
        list_cb2.append(cb2)
        list_mean.append(mean)
        list_mean_offset.append(mean1)
        list_sigma.append(sigma)
        list_mean_gt.append(mean_gt)
        list_mean_gt_offset.append(mean_gt_offset)
        list_sigma_gt.append(sigma_gt)
        list_norm.append(sig_norm)
    
        garbage.append(mean_gt_argl)
        garbage.append(cbs1)
        garbage.append(cbs2)
        garbage.append(gauss)
        garbage.append(argl)
        garbage.append(norms_argl)
        garbage.append(sig)
        garbage.append(sig_argl)
        garbage.append(sig_norm_argl)
        garbage.append(pdf_sig)

        pdf_sig.Print()
        pdf_tot.addPdf(pdf_sig, catIDx)

        if mH == 125.0 and h_obs == None: h_obs = hist_zh.Clone("h_obs") # take 125.0 GeV to add to observed (need to add background later as well)

    rdh_tot = ROOT.RooDataHist("rdh_tot", "", ROOT.RooArgList(recoilmass), cats, hists)

    fitRes = pdf_tot.fitTo(rdh_tot, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.Minimizer("Minimizer", "simplex")) # , ROOT.RooFit.SumW2Error(ROOT.kTRUE)
    fitValid = (fitRes.covQual() == 3 and fitRes.status() == 0)
    print(f"****************************")
    print(f"FIT STATUS")
    print(f"Covariance Quality = {fitRes.covQual()}")
    print(f"Fit status = {fitRes.status()}")
    print(f"****************************")

    getattr(w_tmp, 'import')(pdf_tot) # import after fit, to have fit values in the workspace

    cov = fitRes.covarianceMatrix()
    cov.Print()

    # plot
    plotter.cfg = cfg
    cfg['ytitle'] = f"Events / {20000/nBins:.0f} MeV"
    for i, proc in enumerate(procs):

        mH = mHs[i]
        mH_ = (f"{mH:.3f}").replace(".", "p")
        cfg['ymax'] = yMax

        pdf = pdf_sigs[i]
        rdh_zh = hists[mH_]

        # do plotting
        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(rline=0)
        
        ## TOP PAD ##
        canvas.cd()
        padT.Draw()
        padT.SetGrid()
        padT.cd()
        dummyT.Draw("HIST")

        plt = recoilmass.frame()
        plt.SetTitle("ZH signal")
        rdh_zh.plotOn(plt, ROOT.RooFit.Binning(nBins)) # ROOT.RooFit.Normalization(yield_zh, ROOT.RooAbsReal.NumEvent)

        pdf.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
        chisq = plt.chiSquare()
        pdf.paramOn(plt, ROOT.RooFit.Format("NELU", ROOT.RooFit.AutoPrecision(2)), ROOT.RooFit.Layout(0.45, 0.9, 0.9))

        histpull = plt.pullHist()
        plt.Draw("SAME")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatex(0.2, 0.88, label)
        latex.DrawLatex(0.2, 0.82, f"#chi^{{2}} = {chisq:.3f}")

        plotter.auxRatio()
        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()
        
        ## BOTTOM PAD ##
        canvas.cd()
        padB.Draw()
        padB.cd()
        dummyB.GetXaxis().SetTitleOffset(4.0*dummyB.GetXaxis().GetTitleOffset())
        dummyB.Draw("HIST")

        plt = recoilmass.frame()
        plt.addPlotable(histpull, "P")
        plt.Draw("SAME")
        
        line = ROOT.TLine(120, 0, 140, 0)
        line.SetLineColor(ROOT.kBlue+2)
        line.SetLineWidth(2)
        line.Draw("SAME")

        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()
        canvas.SaveAs(f"{outDir}/fit_mH{mH_}.png")
        canvas.SaveAs(f"{outDir}/fit_mH{mH_}.pdf")

        del dummyB
        del dummyT
        del padT
        del padB
        del canvas

        cb1__ = w_tmp.obj(f"cb1_{mH_}").getVal()
        cb2__ = w_tmp.obj(f"cb2_{mH_}").getVal()

        cfg['ymax'] = 2.5*yMax

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

        cbs_1 = w_tmp.obj(f"cbs1_{mH_}")
        cbs_2 = w_tmp.obj(f"cbs2_{mH_}")
        gauss = w_tmp.obj(f"gauss_{mH_}")
        sig_fit = w_tmp.obj(f"zh_model_{mH_}")

        cbs_1.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Normalization(cb1__*yield_nom, ROOT.RooAbsReal.NumEvent))
        cbs_2.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Normalization(cb2__*yield_nom, ROOT.RooAbsReal.NumEvent))
        gauss.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kCyan), ROOT.RooFit.Normalization((1.-cb1__-cb2__)*yield_nom, ROOT.RooAbsReal.NumEvent))
        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

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

        ROOT.gPad.SetTickx()
        ROOT.gPad.SetTicky()
        ROOT.gPad.RedrawAxis()
        canvas.SaveAs(f"{outDir}/fit_mH{mH_}_decomposition.png")
        canvas.SaveAs(f"{outDir}/fit_mH{mH_}_decomposition.pdf")

        # import
        getattr(w_tmp, 'import')(rdh_zh)
        getattr(w_tmp, 'import')(sig_fit)


        param_mh.append(mH)
        param_mean.append(list_mean[i].getVal())
        param_mean_offset.append(list_mean_offset[i].getVal())
        param_sigma.append(list_sigma[i].getVal())
        param_mean_gt.append(list_mean_gt[i].getVal())
        param_mean_gt_offset.append(list_mean_gt_offset[i].getVal())
        param_sigma_gt.append(list_sigma_gt[i].getVal())
        param_alpha_1.append(list_alpha1[i].getVal())
        param_alpha_2.append(list_alpha2[i].getVal())
        param_n_1.append(list_n1[i].getVal())
        param_n_2.append(list_n2[i].getVal())
        param_yield.append(list_norm[i].getVal())
        param_cb_1.append(list_cb1[i].getVal())
        param_cb_2.append(list_cb2[i].getVal())



    ##################################
    # plot all fitted signals
    ##################################
    cfg['xmin'] = 124
    cfg['xmax'] = 128
    cfg['ymax'] = 2.5*np.average(yMax)
    plotter.cfg = cfg

    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")

    plt = w_tmp.var("zll_recoil_m").frame()
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kBlack, ROOT.kGreen, ROOT.kCyan]
    for i, mH in enumerate(mHs):

        mH_ = (f"{mH:.3f}").replace(".", "p")
        sig_fit = pdf_sigs[i]
        # need to re-normalize the pdf, as the pdf is normalized to 1
        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[i]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

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
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_all.png" )
    canvas.SaveAs(f"{outDir}/fit_all.pdf")


    # make splines, to connect the fit parameters a function of the Higgs mass
    # plot them afterwards
    spline_mean = ROOT.RooSpline1D("spline_mean", "spline_mean", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_mean))
    spline_mean_offset = ROOT.RooSpline1D("spline_mean_offset", "spline_mean_offset", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_mean_offset))
    spline_sigma = ROOT.RooSpline1D("spline_sigma", "spline_sigma", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_sigma))
    spline_mean_gt = ROOT.RooSpline1D("spline_mean_gt", "spline_mean_gt", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_mean_gt))
    spline_mean_gt_offset = ROOT.RooSpline1D("spline_mean_gt_offset", "spline_mean_gt_offset", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_mean_gt_offset))
    spline_sigma_gt = ROOT.RooSpline1D("spline_sigma_gt", "spline_sigma_gt", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_sigma_gt))
    spline_yield = ROOT.RooSpline1D("spline_yield", "spline_yield", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_yield))
    spline_alpha_1 = ROOT.RooSpline1D("spline_alpha_1", "spline_alpha_1", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_alpha_1))
    spline_alpha_2 = ROOT.RooSpline1D("spline_alpha_2", "spline_alpha_2", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_alpha_2))
    spline_n_1 = ROOT.RooSpline1D("spline_n_1", "spline_n_1", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_n_1))
    spline_n_2 = ROOT.RooSpline1D("spline_n_2", "spline_n_2", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_n_2))
    spline_cb_1 = ROOT.RooSpline1D("spline_cb_1", "spline_cb_1", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_cb_1))
    spline_cb_2 = ROOT.RooSpline1D("spline_cb_2", "spline_cb_2", MH, len(param_mh), array.array('d', param_mh), array.array('d', param_cb_2))

    form_mean = ROOT.RooFormulaVar("form_mean", "@0*@1 + @2", ROOT.RooArgList(mean0, MH, mean1))
    form_mean.Print()

    ##################################
    # mean
    ##################################
    graph_mean = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_mean), array.array('d', [0]*len(param_mean)), array.array('d', [0]*len(param_mean)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.999*min(param_mean),
        'ymax'              : 1.001*max(param_mean),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#mu (GeV)",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_mean.plotOn(plt)
    graph_mean.SetMarkerStyle(8)
    graph_mean.SetMarkerColor(ROOT.kBlack)
    graph_mean.SetMarkerSize(1.5)
    graph_mean.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_mean.png")
    canvas.SaveAs(f"{outDir}/fit_mean.pdf")
    
    ##################################
    # mean_gt
    ##################################
    graph_mean_gt = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_mean_gt), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.999*min(param_mean_gt),
        'ymax'              : 1.001*max(param_mean_gt),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#mu_{gt} (GeV)",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_mean_gt.plotOn(plt)
    graph_mean_gt.SetMarkerStyle(8)
    graph_mean_gt.SetMarkerColor(ROOT.kBlack)
    graph_mean_gt.SetMarkerSize(1.5)
    graph_mean_gt.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_mean_gt.png")
    canvas.SaveAs(f"{outDir}/fit_mean_gt.pdf")

    ##################################
    # mean_offset
    ##################################
    print(param_mh)
    print(param_mean_offset)
    print(param_mean_offset_err)
    graph_mean_offset = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_mean_offset), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,
    
        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.999*min(param_mean_offset),
        'ymax'              : 1.001*max(param_mean_offset),
        
        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#mu offset (GeV)",
        
        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_mean_offset.plotOn(plt)
    graph_mean_offset.SetMarkerStyle(8)
    graph_mean_offset.SetMarkerColor(ROOT.kBlack)
    graph_mean_offset.SetMarkerSize(1.5)
    graph_mean_offset.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_mean_offset.png")
    canvas.SaveAs(f"{outDir}/fit_mean_offset.pdf")

    ##################################
    # mean_gt_offset
    ##################################
    graph_mean_gt = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_mean_gt_offset), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.999*min(param_mean_gt_offset),
        'ymax'              : 1.001*max(param_mean_gt_offset),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#mu_{gt} offset (GeV)",
        
        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_mean_gt_offset.plotOn(plt)
    graph_mean_gt.SetMarkerStyle(8)
    graph_mean_gt.SetMarkerColor(ROOT.kBlack)
    graph_mean_gt.SetMarkerSize(1.5)
    graph_mean_gt.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_mean_gt_offset.png")
    canvas.SaveAs(f"{outDir}/fit_mean_gt_offset.pdf")

    ##################################
    # signal yield
    ##################################
    graph_yield = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_yield), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,
    
        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.95*min(param_yield),
        'ymax'              : 1.05*max(param_yield),
        
        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "Events",
        
        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_yield.plotOn(plt)
    graph_yield.SetMarkerStyle(8)
    graph_yield.SetMarkerColor(ROOT.kBlack)
    graph_yield.SetMarkerSize(1.5)
    graph_yield.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_yield.png")
    canvas.SaveAs(f"{outDir}/fit_yield.pdf")


    ##################################
    # sigma 
    ##################################
    graph_sigma = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_sigma), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,
    
        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.95*min(param_sigma),
        'ymax'              : 1.05*max(param_sigma),
        
        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#sigma (GeV)",
        
        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_sigma.plotOn(plt)
    graph_sigma.SetMarkerStyle(8)
    graph_sigma.SetMarkerColor(ROOT.kBlack)
    graph_sigma.SetMarkerSize(1.5)
    graph_sigma.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_sigma.png")
    canvas.SaveAs(f"{outDir}/fit_sigma.pdf")

    ##################################
    # sigma_gt
    ##################################
    graph_sigma_gt = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_sigma_gt), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,
    
        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.95*min(param_sigma_gt),
        'ymax'              : 1.05*max(param_sigma_gt),
        
        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#sigma_{gt} (GeV)",
        
        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_sigma_gt.plotOn(plt)
    graph_sigma_gt.SetMarkerStyle(8)
    graph_sigma_gt.SetMarkerColor(ROOT.kBlack)
    graph_sigma_gt.SetMarkerSize(1.5)
    graph_sigma_gt.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_sigma_gt.png")
    canvas.SaveAs(f"{outDir}/fit_sigma_gt.pdf")

    ##################################
    # alpha_1
    ##################################
    graph_alpha_1 = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_alpha_1), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.8*min(param_alpha_1),
        'ymax'              : 1.2*max(param_alpha_1),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#alpha_{1}",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_alpha_1.plotOn(plt)
    graph_alpha_1.SetMarkerStyle(8)
    graph_alpha_1.SetMarkerColor(ROOT.kBlack)
    graph_alpha_1.SetMarkerSize(1.5)
    graph_alpha_1.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_alpha_1.png")
    canvas.SaveAs(f"{outDir}/fit_alpha_1.pdf")

    ##################################
    # alpha_2
    ##################################
    graph_alpha_2 = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_alpha_2), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.8*min(param_alpha_2),
        'ymax'              : 1.2*max(param_alpha_2),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "#alpha_{2}",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_alpha_2.plotOn(plt)
    graph_alpha_2.SetMarkerStyle(8)
    graph_alpha_2.SetMarkerColor(ROOT.kBlack)
    graph_alpha_2.SetMarkerSize(1.5)
    graph_alpha_2.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_alpha_2.png")
    canvas.SaveAs(f"{outDir}/fit_alpha_2.pdf")


    ##################################
    # n_1
    ##################################
    graph_n_1 = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_n_1), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.8*min(param_n_1),
        'ymax'              : 1.2*max(param_n_1),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "n_{1}",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_n_1.plotOn(plt)
    graph_n_1.SetMarkerStyle(8)
    graph_n_1.SetMarkerColor(ROOT.kBlack)
    graph_n_1.SetMarkerSize(1.5)
    graph_n_1.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_n_1.png")
    canvas.SaveAs(f"{outDir}/fit_n_1.pdf")

    ##################################
    # n_2
    ##################################
    graph_n_2 = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_n_2), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.8*min(param_n_2),
        'ymax'              : 1.2*max(param_n_2),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "n_{2}",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_n_2.plotOn(plt)
    graph_n_2.SetMarkerStyle(8)
    graph_n_2.SetMarkerColor(ROOT.kBlack)
    graph_n_2.SetMarkerSize(1.5)
    graph_n_2.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_n_2.png")
    canvas.SaveAs(f"{outDir}/fit_n_2.pdf")

    ##################################
    # cb_1
    ##################################
    graph_cb_1 = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_cb_1), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.8*min(param_cb_1),
        'ymax'              : 1.2*max(param_cb_1),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "cb_{1}",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_cb_1.plotOn(plt)
    graph_cb_1.SetMarkerStyle(8)
    graph_cb_1.SetMarkerColor(ROOT.kBlack)
    graph_cb_1.SetMarkerSize(1.5)
    graph_cb_1.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_cb_1.png")
    canvas.SaveAs(f"{outDir}/fit_cb_1.pdf")

    ##################################
    # cb_2
    ##################################
    graph_cb_2 = ROOT.TGraphErrors(len(param_mh), array.array('d', param_mh), array.array('d', param_cb_2), array.array('d', [0]*len(mHs)), array.array('d', [0]*len(mHs)))

    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 124.9,
        'xmax'              : 125.1,
        'ymin'              : 0.8*min(param_cb_2),
        'ymax'              : 1.2*max(param_cb_2),

        'xtitle'            : "m_{h} (GeV)",
        'ytitle'            : "cb_{2}",

        'topRight'          : topRight,
        'topLeft'           : topLeft,
    }

    plotter.cfg = cfg
    canvas = plotter.canvas(leftMargin=0.2)
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    dummy.GetXaxis().SetNdivisions(305)

    plt = MH.frame()
    spline_cb_2.plotOn(plt)
    graph_cb_2.SetMarkerStyle(8)
    graph_cb_2.SetMarkerColor(ROOT.kBlack)
    graph_cb_2.SetMarkerSize(1.5)
    graph_cb_2.Draw("SAME P")

    latex.DrawLatex(0.25, 0.92, label)
    plt.Draw("SAME")
    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_cb_2.png")
    canvas.SaveAs(f"{outDir}/fit_cb_2.pdf")

    ##################################
    getattr(w_tmp, 'import')(spline_mean)
    getattr(w_tmp, 'import')(spline_sigma)
    getattr(w_tmp, 'import')(spline_yield)
    getattr(w_tmp, 'import')(spline_alpha_1)
    getattr(w_tmp, 'import')(spline_alpha_2)
    getattr(w_tmp, 'import')(spline_n_1)
    getattr(w_tmp, 'import')(spline_n_2)
    getattr(w_tmp, 'import')(spline_cb_1)
    getattr(w_tmp, 'import')(spline_cb_2)
    getattr(w_tmp, 'import')(spline_mean_gt)
    getattr(w_tmp, 'import')(spline_sigma_gt)

    return hist_norm


def doBackgrounds():

    global h_obs

    recoilmass = w_tmp.var("zll_recoil_m")
    hist_bkg = None
    
    
    if flavor == "mumu":
        if args.ecm == "240":
            procs = ["p8_ee_WW_ecm240", "p8_ee_ZZ_ecm240", "wz3p6_ee_mumu_ecm240", "wz3p6_ee_tautau_ecm240", "wzp6_egamma_eZ_Zmumu_ecm240", "wzp6_gammae_eZ_Zmumu_ecm240", "wzp6_gaga_mumu_60_ecm240", "wzp6_gaga_tautau_60_ecm240", "wzp6_ee_nuenueZ_ecm240"]
        else:
            procs = ["p8_ee_WW_ecm365", "p8_ee_ZZ_ecm365", "wzp6_ee_mumu_ecm365", "wz3p6_ee_tautau_ecm365", "wzp6_gammae_eZ_Zmumu_ecm365", "wzp6_gaga_mumu_60_ecm365", "wzp6_gaga_tautau_60_ecm365", "wzp6_ee_nuenueZ_ecm365"]

    if flavor == "ee":
        if args.ecm == "240":
            procs = ["p8_ee_WW_ecm240", "p8_ee_ZZ_ecm240",  "wz3p6_ee_ee_Mee_30_150_ecm240", "wz3p6_ee_tautau_ecm240", "wzp6_egamma_eZ_Zee_ecm240", "wzp6_gammae_eZ_Zee_ecm240", "wzp6_gaga_ee_60_ecm240", "wzp6_gaga_tautau_60_ecm240", "wzp6_ee_nuenueZ_ecm240"]
        else:
            procs = ["p8_ee_WW_ecm365", "p8_ee_ZZ_ecm365", "wz3p6_ee_ee_Mee_30_150_ecm365", "wz3p6_ee_tautau_ecm365", "wzp6_gaga_ee_60_ecm365", "wzp6_gaga_tautau_60_ecm365", "wzp6_ee_nuenueZ_ecm365"]

    for proc in procs:

        hist = getHist(f"{flavor}_{hName}", [proc])
        hist.Scale(lumiScale)
        hist = hist.ProjectionX(f"hist_{proc}", cat_idx_min, cat_idx_max)
        rdh = ROOT.RooDataHist(f"rdh_{proc}", "rdh", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(hist))

        # add to total background
        if hist_bkg == None: hist_bkg = hist
        else: hist_bkg.Add(hist)

        # add to observed 
        if h_obs == None: h_obs = hist.Clone("h_obs")
        else: h_obs.Add(hist)


    hist_bkg.SetName("total_bkg")
    rdh_bkg = ROOT.RooDataHist("rdh_bkg", "rdh_bkg", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(hist_bkg))
    yield_bkg_ = rdh_bkg.sum(False) 

    tmp = hist_bkg.Clone()
    tmp = tmp.Rebin(int(hist_bkg.GetNbinsX() / nBins))
    yMax_bkg = 1.5*tmp.GetMaximum()


    # construct background as 4th order Bernstein polynomial
    b0 = ROOT.RooRealVar("bern0", "bern_coeff", 1, -2, 2)
    b1 = ROOT.RooRealVar("bern1", "bern_coeff", 1, -10, 10)
    b2 = ROOT.RooRealVar("bern2", "bern_coeff", 1, -10, 10)
    b3 = ROOT.RooRealVar("bern3", "bern_coeff", 1, -10, 10)
    bkg = ROOT.RooBernsteinFast(3)("bkg", "bkg", recoilmass, ROOT.RooArgList(b0, b1, b2))

    bkg_norm = ROOT.RooRealVar('bkg_norm_tmp', 'bkg_norm_tmp', yield_bkg_, 0, 1e6)
    bkg_fit = ROOT.RooAddPdf('bkg_fit', '', ROOT.RooArgList(bkg), ROOT.RooArgList(bkg_norm))
    bkg_fit.fitTo(rdh_bkg, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.SumW2Error(sumw2err))


    ########### PLOTTING ###########
    cfg = {

        'logy'              : False,
        'logx'              : False,

        'xmin'              : 120,
        'xmax'              : 140,
        'ymin'              : 0,
        'ymax'              : yMax_bkg,

        'xtitle'            : "m_{rec} (GeV)",
        'ytitle'            : "Events",

        'topRight'          : topRight,
        'topLeft'           : topLeft,

        'ratiofraction'     : 0.3,
        'ytitleR'           : "Pull",
        'yminR'             : -3.5,
        'ymaxR'             : 3.5,
    }

    plotter.cfg = cfg

    canvas, padT, padB = plotter.canvasRatio()
    dummyT, dummyB, dummyL = plotter.dummyRatio(rline=0)
    dummyB.GetXaxis().SetTitleOffset(4.0*dummyB.GetXaxis().GetTitleOffset())

    ## TOP PAD ##
    canvas.cd()
    padT.Draw()
    padT.cd()
    padT.SetGrid()
    dummyT.Draw("HIST")

    plt = recoilmass.frame()
    rdh_bkg.plotOn(plt, ROOT.RooFit.Binning(nBins))

    bkg_fit.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
    bkg_fit.paramOn(plt, ROOT.RooFit.Format("NELU", ROOT.RooFit.AutoPrecision(2)), ROOT.RooFit.Layout(0.5, 0.9, 0.9))

    histpull = plt.pullHist()
    chisq = plt.chiSquare()
    plt.Draw("SAME")

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(13)
    latex.DrawLatex(0.2, 0.88, label)
    latex.DrawLatex(0.2, 0.82, f"#chi^{{2}} = {chisq:.3f}")

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
    canvas.SaveAs(f"{outDir}/fit_bkg.png")
    canvas.SaveAs(f"{outDir}/fit_bkg.pdf")

    # import background parameterization to the workspace
    bkg_norm.setVal(yield_bkg_) # not constant, as
    b0.setConstant(True) # set as constant
    b1.setConstant(True) # set as constant
    b2.setConstant(True) # set as constant
    b3.setConstant(True) # set as constant
    getattr(w_tmp, 'import')(bkg)
    getattr(w_tmp, 'import')(bkg_norm)

    return hist_bkg


def doBES():

    pct = 1
    if mode == "IDEA_BES6pct":
        pct = 6
    
    if ecm == "365":
        pct = 10

    scale_BES = ROOT.RooRealVar("scale_BES", "BES scale parameter", 0, -1, 1)

    ## only consider variation for 125 GeV
    ## assume variations to be indentical for other mass points
    mH = 125.0
    mH_ = (f"{mH:.3f}").replace(".", "p")

    recoilmass = w_tmp.var("zll_recoil_m")
    MH = w_tmp.var("MH")

    param_mean, param_mean_delta, param_sigma, param_alpha_1, param_alpha_2, param_n_1, param_n_2 = [], [], [], [], [], [], []
    param_yield, param_mean_err, param_sigma_err, param_alpha_1_err, param_alpha_2_err, param_n_1_err, param_n_2_err = [], [], [], [], [], [], []


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

    ## for BES, consider only variations in norm, mean and sigma
    ## assume others to be identical to nominal sample
    MH.setVal(125.0) # evaluate all at 125 GeV
    spline_alpha_1 = w_tmp.obj("spline_alpha_1")
    spline_alpha_2 = w_tmp.obj("spline_alpha_2")
    spline_n_1 = w_tmp.obj("spline_n_1")
    spline_n_2 = w_tmp.obj("spline_n_2")
    spline_cb_1 = w_tmp.obj("spline_cb_1")
    spline_cb_2 = w_tmp.obj("spline_cb_2")
    spline_mean_gt = w_tmp.obj("spline_mean_gt")
    spline_mean = w_tmp.obj("spline_mean")
    spline_sigma_gt = w_tmp.obj("spline_sigma_gt")
    spline_sigma = w_tmp.obj("spline_sigma")

    sigma__ = []
    sigma_gt__ = []

    for s in ["Up", "Down"]:

        if s == "Up": proc = f"wzp6_ee_{flavor}H_BES-higher-{pct}pc_ecm{ecm}"
        if s == "Down": proc = f"wzp6_ee_{flavor}H_BES-lower-{pct}pc_ecm{ecm}"
        if ecm == '365' and flavor == "ee":
            proc = "wz3p6_ee_eeH_ecm365" # no BES for electrons at 365

        if mode == "IDEA_3T":
            if flavor == "mumu":
                proc += "_3T"
            else:
                proc = "wzp6_ee_eeH_ecm240_3T" ## MISSING
        if mode == "CLD":
            if flavor == "mumu":
                proc += "_CLD"
            else:
                proc = "wzp6_ee_eeH_ecm240_CLD" ## MISSING
        if mode == "IDEA_2E":
            if flavor == "ee":
                proc = "wzp6_ee_eeH_ecm240_E2"

        if mode == "IDEA_noBES":
            if flavor == "mumu":
                proc = "wzp6_ee_mumuH_noBES_ecm240"
            else:
                proc = "wzp6_ee_eeH_noBES_ecm240"

        hist_zh = getHist(f"{flavor}_{hName}", [proc])
        hist_zh.Scale(lumiScale)
        hist_zh = hist_zh.ProjectionX(f"hist_zh_%s_BES%s" % (mH_, s), cat_idx_min, cat_idx_max)
        print(hist_zh.Integral(), proc)
        hist_zh.SetName(f"hist_zh_%s_BES%s" % (mH_, s))
        hist_zh.Scale(yield_nom/hist_zh.Integral())
        rdh_zh = ROOT.RooDataHist(f"rdh_zh_%s_BES%s" % (mH_, s), "rdh_zh", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(hist_zh))
        yield_zh = rdh_zh.sum(False)


        mean = ROOT.RooRealVar(f"mean_{mH_}_BES{s}", '', spline_mean.getVal())
        sigma = ROOT.RooRealVar(f"sigma_{mH_}_BES{s}", '', spline_sigma.getVal(), 0, 5) # float
        alpha_1 = ROOT.RooRealVar(f"alpha_1_{mH_}_BES{s}", '', spline_alpha_1.getVal())
        alpha_2 = ROOT.RooRealVar(f"alpha_2_{mH_}_BES{s}", '', spline_alpha_2.getVal())
        n_1 = ROOT.RooRealVar(f"n_1_{mH_}_BES{s}", '', spline_n_1.getVal())
        n_2 = ROOT.RooRealVar(f"n_2_{mH_}_BES{s}", '', spline_n_2.getVal())
        mean_gt = ROOT.RooRealVar(f"mean_gt_{mH_}_BES{s}", '', spline_mean_gt.getVal())
        sigma_gt = ROOT.RooRealVar(f"sigma_gt_{mH_}_BES{s}", '', spline_sigma_gt.getVal(), 0, 5) # float
        cb_1 = ROOT.RooRealVar(f"cb_1_{mH_}_BES{s}", '', spline_cb_1.getVal())
        cb_2 = ROOT.RooRealVar(f"cb_2_{mH_}_BES{s}", '', spline_cb_2.getVal())

        cbs_1 = ROOT.RooCBShape(f"CrystallBall_1_{mH_}_BES{s}", "CrystallBall_1", recoilmass, mean, sigma, alpha_1, n_1)
        cbs_2 = ROOT.RooCBShape(f"CrystallBall_2_{mH_}_BES{s}", "CrystallBall_2", recoilmass, mean, sigma, alpha_2, n_2)
        gauss = ROOT.RooGaussian(f"gauss_{mH_}_BES{s}", "gauss", recoilmass, mean_gt, sigma_gt)

        sig = ROOT.RooAddPdf(f"sig_{mH_}_BES{s}", '', ROOT.RooArgList(cbs_1, cbs_2, gauss), ROOT.RooArgList(cb_1, cb_2)) # half of both CB functions
        sig_norm = ROOT.RooRealVar(f"sig_{mH_}_BES{s}_norm", '', yield_zh, 0, 1e6) # fix normalization
        sig_fit = ROOT.RooAddPdf(f"zh_model_{mH_}_BES{s}", '', ROOT.RooArgList(sig), ROOT.RooArgList(sig_norm))
        sig_fit.fitTo(rdh_zh, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.SumW2Error(sumw2err))

        sigma__.append(sigma.getVal())
        sigma_gt__.append(sigma_gt.getVal())

        # do plotting
        cfg['ymax']= yMax
        plotter.cfg = cfg

        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(rline=0)
        dummyB.GetXaxis().SetTitleOffset(4.0*dummyB.GetXaxis().GetTitleOffset())

        ## TOP PAD ##
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        dummyT.Draw("HIST")

        plt = recoilmass.frame()
        plt.SetTitle("ZH signal")
        rdh_zh.plotOn(plt, ROOT.RooFit.Binning(nBins))

        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
        chisq = plt.chiSquare()
        sig_fit.paramOn(plt, ROOT.RooFit.Format("NELU", ROOT.RooFit.AutoPrecision(2)), ROOT.RooFit.Layout(0.25, 0.9, 0.9))

        histpull = plt.pullHist()
        plt.Draw("SAME")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatex(0.2, 0.88, label)
        latex.DrawLatex(0.2, 0.82, f"#chi^{{2}} = {chisq:.3f}")

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
        canvas.SaveAs(f"{outDir}/fit_mH{mH_}_BES{s}.png")
        canvas.SaveAs(f"{outDir}/fit_mH{mH_}_BES{s}.pdf")

        del dummyB
        del dummyT
        del padT
        del padB
        del canvas

        # import
        getattr(w_tmp, 'import')(rdh_zh)
        getattr(w_tmp, 'import')(sig_fit)



    # plot all fitted signals
    cfg['ymax'] = yMax*2.5
    cfg['xmin'] = 124
    cfg['xmax'] = 127
    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()
    dummy.Draw("HIST")
    plt = w_tmp.var("zll_recoil_m").frame()
    colors = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue]

    sig_fit = w_tmp.pdf(f"zh_model_{mH_}_BESUp")
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[0]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))
    sig_fit = w_tmp.pdf(f"zh_model_{mH_}")
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[1]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))
    sig_fit = w_tmp.pdf(f"zh_model_{mH_}_BESDown")
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[2]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

    plt.Draw("SAME")

    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/fit_mH{mH_}_BES.png")
    canvas.SaveAs(f"{outDir}/fit_mH{mH_}_BES.pdf")

    # construct BES uncertainty
    # nominals, w/o the BES uncertainty
    spline_mean = w_tmp.obj("spline_mean")
    spline_sigma = w_tmp.obj("spline_sigma")
    spline_yield = w_tmp.obj("spline_yield")
    MH.setVal(125.0) # evaluate all at 125 GeV
    sigma__nominal = spline_sigma.getVal()
    sigma_gt__nominal = spline_sigma_gt.getVal()

    # sigma param
    #delta = 0.5*(abs(sigma__nominal-sigma__[0]) + abs(sigma__nominal-sigma__[1]))
    delta = 0.5*abs(sigma__[0] - sigma__[1])
    sig_sigma_BES_ = (delta)/sigma__nominal # 1 sigma value  such that (1+bkg_norm_BES)*sigma__nominal = sigma__nominal+delta
    sig_sigma_BES = ROOT.RooRealVar('sig_sigma_BES', 'sig_sigma_BES', sig_sigma_BES_) # constant
    getattr(w_tmp, 'import')(sig_sigma_BES)
    print(sigma__nominal, delta, sig_sigma_BES_)

    # sigma_gt param
    #delta = 0.5*(abs(sigma_gt__nominal-sigma_gt__[0]) + abs(sigma_gt__nominal-sigma_gt__[1]))
    delta = 0.5*abs(sigma_gt__[0] - sigma_gt__[1])
    sig_sigma_gt_BES_ = (delta)/sigma_gt__nominal # 1 sigma value  such that (1+bkg_norm_BES)*sigma__nominal = sigma__nominal+delta
    sig_sigma_gt_BES = ROOT.RooRealVar('sig_sigma_gt_BES', 'sig_sigma_gt_BES', sig_sigma_gt_BES_) # constant
    getattr(w_tmp, 'import')(sig_sigma_gt_BES)
    print(sigma_gt__nominal, delta, sig_sigma_BES_)


def doSQRTS():

    scale_BES = ROOT.RooRealVar("scale_SQRTS", "SQRTS scale parameter", 0, -1, 1)

    ## only consider variation for 125 GeV
    ## assume variations to be indentical for other mass points
    if ecm == "240":
        proc = "wzp6_ee_%sH_ecm240" % (flavor)
    else:
        proc = "wz3p6_ee_%sH_ecm365" % (flavor)

    if mode == "IDEA_3T":
        proc += "_3T"
    if mode == "CLD":
        proc += "_CLD"
    if mode == "IDEA_noBES":
        proc = proc.replace("_ecm240", "_noBES_ecm240")
    if mode == "IDEA_2E" and flavor == "ee":
        proc += "_E2"

    mH = 125.0
    mH_ = ("%.3f" % mH).replace(".", "p")

    recoilmass = w_tmp.var("zll_recoil_m")
    MH = w_tmp.var("MH")

    param_mean, param_mean_delta, param_sigma, param_alpha_1, param_alpha_2, param_n_1, param_n_2 = [], [], [], [], [], [], []
    param_yield, param_mean_err, param_sigma_err, param_alpha_1_err, param_alpha_2_err, param_n_1_err, param_n_2_err = [], [], [], [], [], [], []

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


    ## for SQRTS, consider only variations in norm, CB mean and CB sigma
    ## assume others to be identical to nominal sample
    MH.setVal(125.0) # evaluate all at 125 GeV
    spline_alpha_1 = w_tmp.obj("spline_alpha_1")
    spline_alpha_2 = w_tmp.obj("spline_alpha_2")
    spline_n_1 = w_tmp.obj("spline_n_1")
    spline_n_2 = w_tmp.obj("spline_n_2")
    spline_cb_1 = w_tmp.obj("spline_cb_1")
    spline_cb_2 = w_tmp.obj("spline_cb_2")
    spline_mean_gt = w_tmp.obj("spline_mean_gt")
    spline_mean = w_tmp.obj("spline_mean")
    spline_sigma_gt = w_tmp.obj("spline_sigma_gt")
    spline_sigma = w_tmp.obj("spline_sigma")

    mean__ = []
    mean_gt__ = []

    for s in ["Up", "Down"]:

        if s == "Up": s_ = "up"
        if s == "Down": s_ = "dw"

        hist_zh = getHist(f"{flavor}_{hName}_sqrts{s_}", [proc])
        hist_zh.Scale(lumiScale)
        hist_zh = hist_zh.ProjectionX("hist_zh_%s_SQRTS%s" % (mH_, s), cat_idx_min, cat_idx_max)
        hist_zh.SetName("hist_zh_%s_BES%s" % (mH_, s))
        hist_zh.Scale(yield_nom/hist_zh.Integral())
        rdh_zh = ROOT.RooDataHist("rdh_zh_%s_SQRTS%s" % (mH_, s), "rdh_zh", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(hist_zh))
        yield_zh = rdh_zh.sum(False)

        mean = ROOT.RooRealVar("mean_%s_SQRTS%s" % (mH_, s), '', spline_mean.getVal(), mH-5., mH+5.) # float
        sigma = ROOT.RooRealVar("sigma_%s_SQRTS%s" % (mH_, s), '', spline_sigma.getVal())
        alpha_1 = ROOT.RooRealVar("alpha_1_%s_SQRTS%s" % (mH_, s), '', spline_alpha_1.getVal())
        alpha_2 = ROOT.RooRealVar("alpha_2_%s_SQRTS%s" % (mH_, s), '', spline_alpha_2.getVal())
        n_1 = ROOT.RooRealVar("n_1_%s_SQRTS%s" % (mH_, s), '', spline_n_1.getVal())
        n_2 = ROOT.RooRealVar("n_2_%s_SQRTS%s" % (mH_, s), '', spline_n_2.getVal())
        mean_gt = ROOT.RooRealVar("mean_gt_%s_SQRTS%s" % (mH_, s), '', spline_mean_gt.getVal(), mH-5., mH+5.) # float
        sigma_gt = ROOT.RooRealVar("sigma_gt_%s_SQRTS%s" % (mH_, s), '', spline_sigma_gt.getVal())   

        cbs_1 = ROOT.RooCBShape("CrystallBall_1_%s_SQRTS%s" % (mH_, s), "CrystallBall_1", recoilmass, mean, sigma, alpha_1, n_1)
        cbs_2 = ROOT.RooCBShape("CrystallBall_2_%s_SQRTS%s" % (mH_, s), "CrystallBall_2", recoilmass, mean, sigma, alpha_2, n_2)
        gauss = ROOT.RooGaussian("gauss_%s_SQRTS%s" % (mH_, s), "gauss", recoilmass, mean_gt, sigma_gt)
        cb_1 = ROOT.RooRealVar("cb_1_%s_SQRTS%s" % (mH_, s), '', spline_cb_1.getVal())
        cb_2 = ROOT.RooRealVar("cb_2_%s_SQRTS%s" % (mH_, s), '', spline_cb_2.getVal())

        sig = ROOT.RooAddPdf("sig_%s_SQRTS%s" % (mH_, s), '', ROOT.RooArgList(cbs_1, cbs_2, gauss), ROOT.RooArgList(cb_1, cb_2)) # half of both CB functions
        sig_norm = ROOT.RooRealVar("sig_%s_SQRTS%s_norm" % (mH_, s), '', yield_zh, 0, 1e6) # fix normalization
        sig_fit = ROOT.RooAddPdf("zh_model_%s_SQRTS%s" % (mH_, s), '', ROOT.RooArgList(sig), ROOT.RooArgList(sig_norm))
        sig_fit.fitTo(rdh_zh, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.SumW2Error(sumw2err))

        mean__.append(mean.getVal())
        mean_gt__.append(mean_gt.getVal())

        print("----->", mean.getVal())

        # do plotting
        cfg['ymax']= yMax
        plotter.cfg = cfg

        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(rline=0)
        dummyB.GetXaxis().SetTitleOffset(4.0*dummyB.GetXaxis().GetTitleOffset())

        ## TOP PAD ##
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        dummyT.Draw("HIST")

        plt = recoilmass.frame()
        plt.SetTitle("ZH signal")
        rdh_zh.plotOn(plt, ROOT.RooFit.Binning(nBins))

        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
        chisq = plt.chiSquare()
        sig_fit.paramOn(plt, ROOT.RooFit.Format("NELU", ROOT.RooFit.AutoPrecision(2)), ROOT.RooFit.Layout(0.25, 0.9, 0.9))

        histpull = plt.pullHist()
        plt.Draw("SAME")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatex(0.2, 0.88, label)
        latex.DrawLatex(0.2, 0.82, f"#chi^{{2}} = {chisq:.3f}")

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
        canvas.SaveAs("%s/fit_mH%s_SQRTS%s.png" % (outDir, mH_, s))
        canvas.SaveAs("%s/fit_mH%s_SQRTS%s.pdf" % (outDir, mH_, s))

        del dummyB
        del dummyT
        del padT
        del padB
        del canvas

        # import
        getattr(w_tmp, 'import')(rdh_zh)
        getattr(w_tmp, 'import')(sig_fit)


    # plot all fitted signals
    cfg['ymax'] = yMax*2.5
    cfg['xmin'] = 124
    cfg['xmax'] = 127
    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()

    dummy.Draw("HIST")

    plt = w_tmp.var("zll_recoil_m").frame()
    colors = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue]

    sig_fit = w_tmp.pdf("zh_model_%s_SQRTSUp" % mH_)  # "sig_%s_BES%s_norm" % (mH_, s)
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[0]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

    sig_fit = w_tmp.pdf("zh_model_%s" % mH_)
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[1]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

    sig_fit = w_tmp.pdf("zh_model_%s_SQRTSDown" % mH_)
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[2]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

    plt.Draw("SAME")

    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs("%s/fit_mH%s_SQRTS.png" % (outDir, mH_))
    canvas.SaveAs("%s/fit_mH%s_SQRTS.pdf" % (outDir, mH_))

    # construct SQRTS uncertainty

    # nominals, w/o the BES uncertainty
    MH.setVal(125.0) # evaluate all at 125 GeV
    mean__nominal = spline_mean.getVal()
    mean_gt__nominal = spline_mean_gt.getVal()

    # mean param
    #delta = 0.5*(abs(mean__nominal-mean__[0]) + abs(mean__nominal-mean__[1]))
    delta = 0.5*abs(mean__[0] - mean__[1])
    sig_mean_SQRTS_ = (delta)/mean__nominal # 1 sigma value  such that (1+bkg_norm_BES)*mean__nominal = mean__nominal+delta
    sig_mean_SQRTS = ROOT.RooRealVar('sig_mean_SQRTS', 'sig_mean_SQRTS', sig_mean_SQRTS_) # constant
    getattr(w_tmp, 'import')(sig_mean_SQRTS)
    print(mean__nominal, delta, sig_mean_SQRTS_)


    # mean_gt param
    #delta = 0.5*(abs(mean_gt__nominal-mean_gt__[0]) + abs(mean_gt__nominal-mean_gt__[1]))
    delta = 0.5*abs(mean_gt__[0] - mean_gt__[1])
    sig_mean_gt_SQRTS_ = (delta)/mean_gt__nominal # 1 sigma value  such that (1+bkg_norm_BES)*mean__nominal = mean__nominal+delta
    sig_mean_gt_SQRTS = ROOT.RooRealVar('sig_mean_gt_SQRTS', 'sig_mean_gt_SQRTS', sig_mean_gt_SQRTS_) # constant
    getattr(w_tmp, 'import')(sig_mean_gt_SQRTS)
    print(mean_gt__nominal, delta, sig_mean_gt_SQRTS_)


def doLEPSCALE():

    scale_BES = ROOT.RooRealVar("scale_LEPSCALE", "LEPSCALE scale parameter", 0, -1, 1)

    ## only consider variation for 125 GeV
    ## assume variations to be indentical for other mass points
    if ecm == "240":
        proc = f"wzp6_ee_{flavor}H_ecm240"
    else:
        proc = f"wz3p6_ee_{flavor}H_ecm365"

    if mode == "IDEA_3T":
        proc += "_3T"
    if mode == "CLD":
        proc += "_CLD"
    if mode == "IDEA_noBES":
        proc = proc.replace("_ecm240", "_noBES_ecm240")
    if mode == "IDEA_2E" and flavor == "ee":
        proc += "_E2"


    mH = 125.0
    mH_ = (f"{mH:.3f}").replace(".", "p")

    recoilmass = w_tmp.var("zll_recoil_m")
    MH = w_tmp.var("MH")

    param_mean, param_mean_delta, param_sigma, param_alpha_1, param_alpha_2, param_n_1, param_n_2 = [], [], [], [], [], [], []
    param_yield, param_mean_err, param_sigma_err, param_alpha_1_err, param_alpha_2_err, param_n_1_err, param_n_2_err = [], [], [], [], [], [], []

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


    ## for LEPSCALE, consider only variations in norm, CB mean and CB sigma
    ## assume others to be identical to nominal sample
    MH.setVal(125.0) # evaluate all at 125 GeV
    spline_alpha_1 = w_tmp.obj("spline_alpha_1")
    spline_alpha_2 = w_tmp.obj("spline_alpha_2")
    spline_n_1 = w_tmp.obj("spline_n_1")
    spline_n_2 = w_tmp.obj("spline_n_2")
    spline_cb_1 = w_tmp.obj("spline_cb_1")
    spline_cb_2 = w_tmp.obj("spline_cb_2")
    spline_mean_gt = w_tmp.obj("spline_mean_gt")
    spline_mean = w_tmp.obj("spline_mean")
    spline_sigma_gt = w_tmp.obj("spline_sigma_gt")
    spline_sigma = w_tmp.obj("spline_sigma")

    mean__ = []

    for s in ["Up", "Down"]:

        if s == "Up": s_ = "up"
        if s == "Down": s_ = "dw"

        hist_zh = getHist(f"{flavor}_{hName}_scale{s_}", [proc])
        hist_zh.Scale(lumiScale)
        hist_zh = hist_zh.ProjectionX("hist_zh_%s_LEPSCALE%s" % (mH_, s), cat_idx_min, cat_idx_max)   
        hist_zh.SetName("hist_zh_%s_LEPSCALE%s" % (mH_, s))
        hist_zh.Scale(yield_nom/hist_zh.Integral())
        rdh_zh = ROOT.RooDataHist("rdh_zh_%s_LEPSCALE%s" % (mH_, s), "rdh_zh", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(hist_zh))
        yield_zh = rdh_zh.sum(False)

        mean = ROOT.RooRealVar("mean_%s_LEPSCALE%s" % (mH_, s), '', spline_mean_gt.getVal(), mH-5., mH+5.)
        sigma = ROOT.RooRealVar("sigma_%s_LEPSCALE%s" % (mH_, s), '', spline_sigma.getVal())
        alpha_1 = ROOT.RooRealVar("alpha_1_%s_LEPSCALE%s" % (mH_, s), '', spline_alpha_1.getVal())
        alpha_2 = ROOT.RooRealVar("alpha_2_%s_LEPSCALE%s" % (mH_, s), '', spline_alpha_2.getVal())
        n_1 = ROOT.RooRealVar("n_1_%s_LEPSCALE%s" % (mH_, s), '', spline_n_1.getVal())
        n_2 = ROOT.RooRealVar("n_2_%s_LEPSCALE%s" % (mH_, s), '', spline_n_2.getVal())
        mean_gt = ROOT.RooRealVar("mean_gt_%s_LEPSCALE%s" % (mH_, s), '', spline_mean_gt.getVal())
        sigma_gt = ROOT.RooRealVar("sigma_gt_%s_LEPSCALE%s" % (mH_, s), '', spline_sigma_gt.getVal())

        cbs_1 = ROOT.RooCBShape("CrystallBall_1_%s_LEPSCALE%s" % (mH_, s), "CrystallBall_1", recoilmass, mean, sigma, alpha_1, n_1)
        cbs_2 = ROOT.RooCBShape("CrystallBall_2_%s_LEPSCALE%s" % (mH_, s), "CrystallBall_2", recoilmass, mean, sigma, alpha_2, n_2)
        gauss = ROOT.RooGaussian("gauss_%s_LEPSCALE%s" % (mH_, s), "gauss", recoilmass, mean_gt, sigma_gt)
        cb_1 = ROOT.RooRealVar("cb_1_%s_LEPSCALE%s" % (mH_, s), '', spline_cb_1.getVal())
        cb_2 = ROOT.RooRealVar("cb_2_%s_LEPSCALE%s" % (mH_, s), '', spline_cb_2.getVal())
                
        sig = ROOT.RooAddPdf("sig_%s_LEPSCALE%s" % (mH_, s), '', ROOT.RooArgList(cbs_1, cbs_2, gauss), ROOT.RooArgList(cb_1, cb_2)) # half of both CB functions
        sig_norm = ROOT.RooRealVar("sig_%s_LEPSCALE%s_norm" % (mH_, s), '', yield_zh, 0, 1e6) # fix normalization
        sig_fit = ROOT.RooAddPdf("zh_model_%s_LEPSCALE%s" % (mH_, s), '', ROOT.RooArgList(sig), ROOT.RooArgList(sig_norm))
        sig_fit.fitTo(rdh_zh, ROOT.RooFit.Extended(ROOT.kTRUE), ROOT.RooFit.SumW2Error(sumw2err))

        mean__.append(mean.getVal())

        # do plotting
        cfg['ymax']= yMax
        plotter.cfg = cfg

        canvas, padT, padB = plotter.canvasRatio()
        dummyT, dummyB, dummyL = plotter.dummyRatio(rline=0)
        dummyB.GetXaxis().SetTitleOffset(4.0*dummyB.GetXaxis().GetTitleOffset())

        ## TOP PAD ##
        canvas.cd()
        padT.Draw()
        padT.cd()
        padT.SetGrid()
        dummyT.Draw("HIST")

        plt = recoilmass.frame()
        plt.SetTitle("ZH signal")
        rdh_zh.plotOn(plt, ROOT.RooFit.Binning(nBins))

        sig_fit.plotOn(plt, ROOT.RooFit.LineColor(ROOT.kRed))
        chisq = plt.chiSquare()  
        sig_fit.paramOn(plt, ROOT.RooFit.Format("NELU", ROOT.RooFit.AutoPrecision(2)), ROOT.RooFit.Layout(0.25, 0.9, 0.9))

        histpull = plt.pullHist()
        plt.Draw("SAME")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextColor(1)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatex(0.2, 0.88, label)
        latex.DrawLatex(0.2, 0.82, f"#chi^{{2}} = {chisq:.3f}")

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
        canvas.SaveAs("%s/fit_mH%s_LEPSCALE%s.png" % (outDir, mH_, s))
        canvas.SaveAs("%s/fit_mH%s_LEPSCALE%s.pdf" % (outDir, mH_, s))

        del dummyB
        del dummyT
        del padT
        del padB
        del canvas

        # import
        getattr(w_tmp, 'import')(rdh_zh)
        getattr(w_tmp, 'import')(sig_fit)
        #getattr(w_tmp, 'import')(sig_norm) # already imported with sig_fit


    # plot all fitted signals
    cfg['ymax'] = yMax*2.5
    cfg['xmin'] = 124
    cfg['xmax'] = 127
    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    dummy = plotter.dummy()

    dummy.Draw("HIST")

    plt = w_tmp.var("zll_recoil_m").frame()
    colors = [ROOT.kRed, ROOT.kBlack, ROOT.kBlue]

    sig_fit = w_tmp.pdf("zh_model_%s_LEPSCALEUp" % mH_)  # "sig_%s_BES%s_norm" % (mH_, s)
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[0]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

    sig_fit = w_tmp.pdf("zh_model_%s" % mH_)
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[1]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

    sig_fit = w_tmp.pdf("zh_model_%s_LEPSCALEDown" % mH_)
    sig_fit.plotOn(plt, ROOT.RooFit.LineColor(colors[2]), ROOT.RooFit.Normalization(yield_nom, ROOT.RooAbsReal.NumEvent))

    plt.Draw("SAME")

    plotter.aux()
    canvas.Modify()
    canvas.Update()
    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()
    canvas.Draw()
    canvas.SaveAs("%s/fit_mH%s_LEPSCALE.png" % (outDir, mH_))
    canvas.SaveAs("%s/fit_mH%s_LEPSCALE.pdf" % (outDir, mH_))

    # construct LEPSCALE uncertainty
    # nominals, w/o the BES uncertainty
    spline_mean = w_tmp.obj("spline_mean")
    spline_sigma = w_tmp.obj("spline_sigma")
    MH.setVal(125.0) # evaluate all at 125 GeV
    mean__nominal = spline_mean.getVal()
    sigma__nominal = spline_sigma.getVal()

    # mean param
    #delta = 0.5*(abs(mean__nominal-mean__[0]) + abs(mean__nominal-mean__[1]))
    delta = 0.5*abs(mean__[0] - mean__[1])
    sig_mean_LEPSCALE_ = (delta)/mean__nominal # 1 sigma value  such that (1+bkg_norm_BES)*mean__nominal = mean__nominal+delta
    sig_mean_LEPSCALE = ROOT.RooRealVar('sig_mean_LEPSCALE', 'sig_mean_LEPSCALE', sig_mean_LEPSCALE_) # constant
    getattr(w_tmp, 'import')(sig_mean_LEPSCALE)
    print(mean__nominal, delta, sig_mean_LEPSCALE_)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--flavor", type=str, help="Flavor (mumu or ee)", default="mumu")
    parser.add_argument("--mode", type=str, help="Detector mode", choices=["IDEA", "IDEA_MC", "IDEA_3T", "CLD", "CLD_FullSim", "IDEA_noBES", "IDEA_2E", "IDEA_BES6pct"], default="IDEA")
    parser.add_argument("--cat", type=str, help="Category (0, 1, 2 or 3)", choices=["0", "1", "2", "3"], default="0")
    parser.add_argument("--lumi", type=float, help="Luminosity scale", default=1.0)
    parser.add_argument("--ecm", type=str, help="Center-of-mass", choices=["240", "365"], default="240")
    parser.add_argument("--tag", type=str, help="Analysis tag for versioning, optional", default="")
    parser.add_argument("--singularity", action='store_true', help="Run in singularity mode")
    args = parser.parse_args()



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
    lumiScale = args.lumi

    # re-normalize lumis, by default histograms are scaled to 10.8 and 3.12 ab-1
    if ecm == "240" or ecm == 240:
        lumiScale /= 10.8
    elif ecm == "365" or ecm == 365:
        lumiScale /= 3.12
    else:
        print(f"Center-of-masss {ecm} not supported")
        quit()
    lumiStr = str(int(args.lumi)) if args.lumi.is_integer() else str(args.lumi)
    lumiLabel = lumiStr.replace(".", "p")

    topRight = f"#sqrt{{s}} = {ecm} GeV, {lumiStr} ab^{{#minus1}}"
    topLeft = f"#bf{{FCC-ee}} #scale[0.7]{{#it{{Internal}}}}"
    label = f"{flavorLabel}, category {cat}"
    inputDir = f"output/h_mass/histmaker/ecm{ecm}/"
    outDir = f"/work/submit/jaeyserm/public_html/fccee/h_mass/{tag}/combine/{mode}/lumi{lumiLabel}/{flavor}_cat{cat}_ecm{ecm}/"
    runDir = f"output/h_mass/combine/{tag}/{mode}/lumi{lumiLabel}/{flavor}_cat{cat}_ecm{ecm}/"
    if not os.path.exists(outDir): os.makedirs(outDir)
    if not os.path.exists(runDir): os.makedirs(runDir)
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/h_mass/index.php {outDir}")

    parId = f"{flavor}_cat{cat}_ecm{ecm}"
    hName = f"zll_recoil_m_cat"
    if cat == 0: cat_idx_min, cat_idx_max = 0, 5
    else: cat_idx_min, cat_idx_max = cat, cat

    if not os.path.exists(runDir): os.makedirs(runDir)
    if not os.path.exists(outDir): os.makedirs(outDir)

    nBins = 250 # total number of bins, for plotting
    recoilMin = 120
    recoilMax = 140
    h_obs = None # should hold the data_obs = sum of signal and backgrounds

    recoilmass = ROOT.RooRealVar("zll_recoil_m", "m_{rec} (GeV)", 125, recoilMin, recoilMax)
    MH = ROOT.RooRealVar("MH", "Higgs mass (GeV)", 125, 124.95, 125.05)

    pdf_sigs = []

    # define temporary output workspace
    w_tmp = ROOT.RooWorkspace("w_tmp", "workspace")
    w = ROOT.RooWorkspace("w", "workspace") # final workspace for combine

    getattr(w_tmp, 'import')(recoilmass)
    getattr(w_tmp, 'import')(MH)

    yield_nom = -1
    yMax = -1

    hist_sig = doSignal()
    hist_bkg = doBackgrounds()
    doSyst = True
    if doSyst:

        doBES() # 1 or 6 pct BES variation
        doSQRTS()
        doLEPSCALE()

        # systematic strenghts
        BES = ROOT.RooRealVar('BES_ecm%s'%ecm, 'BES', 0, -5, 5) # BES uncertainty parameter
        #ISR = ROOT.RooRealVar('ISR', 'ISR', 0, -5, 5) # BES uncertainty parameter
        ISR = ROOT.RooRealVar('ISR_ecm%s'%ecm, 'ISR', 0) # BES uncertainty parameter
        SQRTS = ROOT.RooRealVar('SQRTS_ecm%s'%ecm, 'SQRTS', 0, -5, 5) # SQRTS uncertainty parameter
        LEPSCALE = ROOT.RooRealVar('LEPSCALE_%s_ecm%s'%("MU" if flavor=="mumu" else "EL", ecm), 'LEPSCALE', 0, -5, 5) # LEPSCALE uncertainty parameter

        '''
        was before
        sig_mean_ISR = ROOT.RooRealVar('sig_mean_ISR_%s_cat%d'%(flavor,cat), 'sig_mean_ISR', 0)
        sig_sigma_ISR = ROOT.RooRealVar('sig_sigma_ISR_%s_cat%d'%(flavor,cat), 'sig_sigma_ISR', 0)
        sig_norm_ISR = ROOT.RooRealVar('sig_norm_ISR_%s_cat%d'%(flavor,cat), 'sig_norm_ISR', 0)
        sig_n_1_ISR = ROOT.RooRealVar('sig_n_1_ISR_%s_cat%d'%(flavor,cat), 'sig_n_1_ISR', 0)
        sig_n_2_ISR = ROOT.RooRealVar('sig_n_2_ISR_%s_cat%d'%(flavor,cat), 'sig_n_2_ISR', 0)
        sig_mean_gt_ISR = ROOT.RooRealVar('sig_mean_gt_ISR_%s_cat%d'%(flavor,cat), 'sig_mean_gt_ISR', 0)
        
        sig_mean_ISR = ROOT.RooRealVar('sig_mean_ISR', 'sig_mean_ISR', 0)
        sig_sigma_ISR = ROOT.RooRealVar('sig_sigma_ISR', 'sig_sigma_ISR', 0)
        sig_norm_ISR = ROOT.RooRealVar('sig_norm_ISR', 'sig_norm_ISR', 0)
        sig_n_1_ISR = ROOT.RooRealVar('sig_n_1_ISR'), 'sig_n_1_ISR', 0)
        sig_n_2_ISR = ROOT.RooRealVar('sig_n_2_ISR', 'sig_n_2_ISR', 0)
        sig_mean_gt_ISR = ROOT.RooRealVar('sig_mean_gt_ISR', 'sig_mean_gt_ISR', 0)
        '''
        # BES
        #sig_norm_BES = w_tmp.obj("sig_norm_BES")
        #sig_norm_BES.SetName("sig_norm_BES_%s_cat%d"%(flavor,cat))
        #sig_mean_BES = w_tmp.obj("sig_mean_BES")
        #sig_mean_BES.SetName("sig_mean_BES_%s_cat%d"%(flavor,cat))
        sig_sigma_BES = w_tmp.obj("sig_sigma_BES")
        #sig_sigma_BES.SetName("sig_sigma_BES_%s_cat%d"%(flavor,cat))
        sig_sigma_gt_BES = w_tmp.obj("sig_sigma_gt_BES")
        #sig_sigma_gt_BES.SetName("sig_sigma_gt_BES_%s_cat%d"%(flavor,cat))

        # SQRTS
        #sig_norm_SQRTS = w_tmp.obj("sig_norm_SQRTS")
        sig_mean_SQRTS = w_tmp.obj("sig_mean_SQRTS")
        #sig_mean_SQRTS.SetName("sig_mean_SQRTS_%s_cat%d"%(flavor,cat))
        sig_mean_gt_SQRTS = w_tmp.obj("sig_mean_gt_SQRTS")
        #sig_mean_gt_SQRTS.SetName("sig_mean_gt_SQRTS_%s_cat%d"%(flavor,cat))
        
        # LEPSCALE
        #sig_norm_LEPSCALE = w_tmp.obj("sig_norm_LEPSCALE")
        sig_mean_LEPSCALE = w_tmp.obj("sig_mean_LEPSCALE")
        #sig_mean_LEPSCALE.SetName("sig_mean_LEPSCALE_%s_cat%d"%(flavor,cat))
        #sig_sigma_LEPSCALE = w_tmp.obj("sig_sigma_LEPSCALE")
        #sig_sigma_LEPSCALE.SetName("sig_sigma_LEPSCALE_%s_cat%d"%(flavor,cat))


        # build signal model, taking into account all uncertainties
        spline_mean = w_tmp.obj("spline_mean")
        spline_sigma = w_tmp.obj("spline_sigma")
        spline_yield = w_tmp.obj("spline_yield")
        spline_alpha_1 = w_tmp.obj("spline_alpha_1")
        spline_alpha_2 = w_tmp.obj("spline_alpha_2")
        spline_n_1 = w_tmp.obj("spline_n_1")
        spline_n_2 = w_tmp.obj("spline_n_2")
        spline_cb_1 = w_tmp.obj("spline_cb_1")
        spline_cb_2 = w_tmp.obj("spline_cb_2")
        spline_mean_gt = w_tmp.obj("spline_mean_gt")
        spline_sigma_gt = w_tmp.obj("spline_sigma_gt")


        #sig_mean = ROOT.RooFormulaVar("sig_mean", "@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)*(1+@7*@8)", ROOT.RooArgList(spline_mean, BES, sig_mean_BES, ISR, sig_mean_ISR, SQRTS, sig_mean_SQRTS, LEPSCALE, sig_mean_LEPSCALE))
        sig_mean = ROOT.RooFormulaVar("sig_mean", "@0*(1+@1*@2)*(1+@3*@4)", ROOT.RooArgList(spline_mean, LEPSCALE, sig_mean_LEPSCALE, SQRTS, sig_mean_SQRTS))
        sig_sigma = ROOT.RooFormulaVar("sig_sigma", "@0*(1+@1*@2)", ROOT.RooArgList(spline_sigma, BES, sig_sigma_BES))
        sig_alpha_1 = ROOT.RooFormulaVar("sig_alpha_1", "@0", ROOT.RooArgList(spline_alpha_1))
        sig_alpha_2 = ROOT.RooFormulaVar("sig_alpha_2", "@0", ROOT.RooArgList(spline_alpha_2))
        sig_n_1 = ROOT.RooFormulaVar("sig_n_1", "@0", ROOT.RooArgList(spline_n_1))
        sig_n_2 = ROOT.RooFormulaVar("sig_n_2", "@0", ROOT.RooArgList(spline_n_2))
        sig_cb_1 = ROOT.RooFormulaVar("sig_cb_1", "@0", ROOT.RooArgList(spline_cb_1))
        sig_cb_2 = ROOT.RooFormulaVar("sig_cb_2", "@0", ROOT.RooArgList(spline_cb_2))
        sig_mean_gt = ROOT.RooFormulaVar("sig_mean_gt", "@0*(1+@1*@2)", ROOT.RooArgList(spline_mean_gt, SQRTS, sig_mean_SQRTS))
        sig_sigma_gt = ROOT.RooFormulaVar("sig_sigma_gt", "@0*(1+@1*@2)*(1+@3*@4)", ROOT.RooArgList(spline_sigma_gt, BES, sig_sigma_gt_BES, SQRTS, sig_mean_gt_SQRTS))
        sig_norm = ROOT.RooFormulaVar("sig_norm", "@0", ROOT.RooArgList(spline_yield))

    else:
        spline_mean = w_tmp.obj("spline_mean")
        spline_sigma = w_tmp.obj("spline_sigma")
        spline_yield = w_tmp.obj("spline_yield")
        spline_alpha_1 = w_tmp.obj("spline_alpha_1")
        spline_alpha_2 = w_tmp.obj("spline_alpha_2")
        spline_n_1 = w_tmp.obj("spline_n_1")
        spline_n_2 = w_tmp.obj("spline_n_2")
        spline_cb_1 = w_tmp.obj("spline_cb_1")
        spline_cb_2 = w_tmp.obj("spline_cb_2")
        spline_mean_gt = w_tmp.obj("spline_mean_gt")
        spline_sigma_gt = w_tmp.obj("spline_sigma_gt")

        sig_mean = ROOT.RooFormulaVar("sig_mean", "@0", ROOT.RooArgList(spline_mean))
        sig_sigma = ROOT.RooFormulaVar("sig_sigma", "@0", ROOT.RooArgList(spline_sigma))
        sig_alpha_1 = ROOT.RooFormulaVar("sig_alpha_1", "@0", ROOT.RooArgList(spline_alpha_1))
        sig_alpha_2 = ROOT.RooFormulaVar("sig_alpha_2", "@0", ROOT.RooArgList(spline_alpha_2))
        sig_n_1 = ROOT.RooFormulaVar("sig_n_1", "@0", ROOT.RooArgList(spline_n_1))
        sig_n_2 = ROOT.RooFormulaVar("sig_n_2", "@0", ROOT.RooArgList(spline_n_2))
        sig_cb_1 = ROOT.RooFormulaVar("sig_cb_1", "@0", ROOT.RooArgList(spline_cb_1))
        sig_cb_2 = ROOT.RooFormulaVar("sig_cb_2", "@0", ROOT.RooArgList(spline_cb_2))
        sig_mean_gt = ROOT.RooFormulaVar("sig_mean_gt", "@0", ROOT.RooArgList(spline_mean_gt))
        sig_sigma_gt = ROOT.RooFormulaVar("sig_sigma_gt", "@0", ROOT.RooArgList(spline_sigma_gt))
        sig_norm = ROOT.RooFormulaVar("sig_norm", "@0", ROOT.RooArgList(spline_yield))

    # construct final signal pdf
    cbs_1 = ROOT.RooCBShape("CrystallBall_1", "CrystallBall_1", recoilmass, sig_mean, sig_sigma, sig_alpha_1, sig_n_1)
    cbs_2 = ROOT.RooCBShape("CrystallBall_2", "CrystallBall_2", recoilmass, sig_mean, sig_sigma, sig_alpha_2, sig_n_2)
    gauss = ROOT.RooGaussian("gauss", "gauss", recoilmass, sig_mean_gt, sig_sigma_gt)
    sig = ROOT.RooAddPdf("sig", "sig", ROOT.RooArgList(cbs_1, cbs_2, gauss), ROOT.RooArgList(sig_cb_1, sig_cb_2))

    #getattr(w, 'import')(sig_norm)
    getattr(w, 'import')(sig)

    # construct background model 66672.5686008
    bkg_yield = w_tmp.obj("bkg_norm_tmp").getVal()
    #bkg_norm = ROOT.RooRealVar('bkg_%s_ecm%s_norm'%(flavor, ecm), 'bkg_norm', bkg_yield) #, 0, 1e6) # nominal background yield, floating bkg_ee_ecm365_norm
    bkg_norm = ROOT.RooRealVar('bkg_norm', 'bkg_norm', bkg_yield) #, 0, 1e6) # nominal background yield (automatically done by Combine with pdfName_norm, floating)
    bkg_norm.setVal(bkg_yield) # not constant!
    bkg = w_tmp.obj("bkg")
    getattr(w, 'import')(bkg, ROOT.RooFit.RenameAllVariablesExcept(parId, "zll_recoil_m"))
    #getattr(w, 'import')(bkg_norm)
    
    #bkg_norm.Print()
    data_obs = ROOT.RooDataHist("data_obs", "data_obs", ROOT.RooArgList(recoilmass), ROOT.RooFit.Import(h_obs))
    getattr(w, 'import')(data_obs)

    w.writeToFile(f"{runDir}/datacard.root")

    '''
    fOut = ROOT.TFile("%s/datacard.root" % runDir, "UPDATE")
    hist_sig.SetName("hist_sig")
    hist_bkg.SetName("hist_bkg")
    hist_sig.Write()
    hist_bkg.Write()
    fOut.Close()
    '''
    w.Print()
    
    poi = ROOT.RooArgSet(MH)
    w.defineSet("POI", poi)

    del w
    del w_tmp

    if ecm == "240" and flavor == "mumu": bkg_id=1
    elif ecm == "240" and flavor == "ee": bkg_id=2
    elif ecm == "365" and flavor == "mumu": bkg_id=3
    elif ecm == "365" and flavor == "ee": bkg_id=4

    # make datacard
    with open(f"{os.path.dirname(os.path.realpath(__file__))}/datacard_template.txt", 'r') as file:
        dc = file.read()
        dc = dc.replace("$rate_sig", f"{yield_nom}")
        dc = dc.replace("$rate_bkg", f"{bkg_yield}")
        dc = dc.replace("$flavor", f"{flavor}")
        dc = dc.replace("$ecm", f"{ecm}")
        dc = dc.replace("$cat", f"{cat}")
        dc = dc.replace("$bkg_id", f"{bkg_id}")

    with open(f"{runDir}/datacard.txt", 'w') as file:
        file.write(dc)
    
    # remove lepton scale uncertainty for other flavor
    if flavor == "mumu":
        cmd = "sed -i '/LEPSCALE_EL/d' datacard.txt"
    else:
        cmd = "sed -i '/LEPSCALE_MU/d' datacard.txt"
    subprocess.call(cmd, shell=True, cwd=runDir)

    # build the Combine workspace based on the datacard, save it to ws.root
    cmd = "text2workspace.py datacard.txt -o ws.root -v 10  --X-allow-no-background"
    subprocess.call(cmd, shell=True, cwd=runDir)
