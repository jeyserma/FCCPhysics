import sys
import os
import math
import copy
import array

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../python')
import plotter



def getHist(hName, procs, rebin=1):
    hist = None
    for proc in procs:
        if "Hinv" in proc:
            proc = proc.replace("wzp6", "wz3p6")
        fIn = ROOT.TFile(f"{inputDir}/{proc}.root")
        h = fIn.Get(hName)
        h.SetDirectory(0)
        if hist == None:
            hist = h
        else:
            hist.Add(h)
        fIn.Close()
    hist.Rebin(rebin)
    #hist.Scale(3./10.8 * 1.05/1.39)
    return hist

def makeCutFlow(hName="cutFlow", cuts=[], labels=[], sig_scale=1.0):

    totEntries = 1 + len(bkgs)
    #leg = ROOT.TLegend(.5, 1.0-totEntries*0.06, .92, .90)
    leg = ROOT.TLegend(.45, 0.99-(len(bkgs)+2)*0.055, .95, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists_yields = []
    significances = []
    h_sig = getHist(hName, sigs)

    hists_yields.append(copy.deepcopy(h_sig))
    h_sig.Scale(sig_scale)
    h_sig.SetLineColor(sig_color)
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    if sig_scale != 1:
        leg.AddEntry(h_sig, f"{sig_legend} (#times {int(sig_scale)})", "L")
    else:
        leg.AddEntry(h_sig, sig_legend, "L")

    # Get all bkg histograms
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(bkgs):
        h_bkg = getHist(hName, bgks_cfg[bkg])

        if h_bkg_tot == None: h_bkg_tot = h_bkg.Clone("h_bkg_tot")
        else: h_bkg_tot.Add(h_bkg)
        
        h_bkg.SetFillColor(bkgs_colors[i])
        h_bkg.SetLineColor(ROOT.kBlack)
        h_bkg.SetLineWidth(1)
        h_bkg.SetLineStyle(1)

        leg.AddEntry(h_bkg, bkgs_legends[i], "F")
        st.Add(h_bkg)
        hists_yields.append(h_bkg)

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetLineWidth(2)

    for i,cut in enumerate(cuts):
        nsig = h_sig.GetBinContent(i+1) / sig_scale ## undo scaling
        nbkg = 0
        for j,histProc in enumerate(hists_yields):
            nbkg = nbkg + histProc.GetBinContent(i+1)
            print(histProc.GetBinContent(i+1))
        if (nsig+nbkg) == 0:
            print(f"Cut {cut} zero yield sig+bkg")
            s = -1
        else:
            s = nsig / (nsig + nbkg)**0.5
        print(i, cut, s)
        significances.append(s)

    ########### PLOTTING ###########
    cfg = {
        'logy'              : True,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : len(cuts),
        'ymin'              : 1e6,
        'ymax'              : 1e10 ,

        'xtitle'            : "",
        'ytitle'            : "Events",

        'topRight'          : f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}",
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    canvas.SetTicks()
    dummy = plotter.dummy(len(cuts))
    dummy.GetXaxis().SetLabelSize(0.8*dummy.GetXaxis().GetLabelSize())
    dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
    for i,label in enumerate(labels): dummy.GetXaxis().SetBinLabel(i+1, label)
    dummy.GetXaxis().LabelsOption("u")
    dummy.Draw("HIST")

    st.Draw("SAME HIST")
    h_bkg_tot.Draw("SAME HIST")
    h_sig.Draw("SAME HIST")
    leg.Draw("SAME")

    plotter.aux()
    canvas.RedrawAxis()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/{hName}.png")
    canvas.SaveAs(f"{outDir}/{hName}.pdf")

    out_orig = sys.stdout
    with open(f"{outDir}{hName}.txt", 'w') as f:
        sys.stdout = f

        formatted_row = '{:<10} {:<25} {:<25} {:<25} {:<25} {:<25}' # adapt to #bkgs
        print(formatted_row.format(*(["Cut", "Significance", "Signal"]+bkgs)))
        print(formatted_row.format(*(["----------"]+["-----------------------"]*5)))
        for i,cut in enumerate(cuts):
            row = ["Cut %d"%i, "%.3f"%significances[i]]
            for j,histProc in enumerate(hists_yields):
                yield_, err = histProc.GetBinContent(i+1), histProc.GetBinError(i+1)
                row.append("%.4e +/- %.2e" % (yield_, err))

            print(formatted_row.format(*row))
    sys.stdout = out_orig


def makeCutFlowHiggsDecays(hName="cutFlow", cuts=[], labels=[]):

    sigs = [[f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays] for y in higgs_decays]

    totEntries = 1 + len(sigs)
    #leg = ROOT.TLegend(.5, 1.0-totEntries*0.06, .92, .90)
    leg = ROOT.TLegend(.6, 0.99-(len(bkgs)+2)*0.07, .95, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)
    leg.SetNColumns(2)


    hists = []
    for i,sig in enumerate(sigs):
        h_decay = higgs_decays[i]
        h_sig = getHist(hName, sig)
        h_sig.Scale(100./h_sig.GetBinContent(1))

        h_sig.SetLineColor(colors[h_decay])
        h_sig.SetLineWidth(2)
        h_sig.SetLineStyle(1)

        leg.AddEntry(h_sig, higgs_decays_labels[h_decay], "L")
        hists.append(h_sig)


    ########### PLOTTING ###########
    cfg = {
        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : len(cuts),
        'ymin'              : 0,
        'ymax'              : 150 ,

        'xtitle'            : "",
        'ytitle'            : "Selection efficiency (%)",

        'topRight'          : f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}",
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",
        }

    plotter.cfg = cfg

    canvas = plotter.canvas()
    canvas.SetGrid()
    canvas.SetTicks()
    dummy = plotter.dummy(len(cuts))
    dummy.GetXaxis().SetLabelSize(0.8*dummy.GetXaxis().GetLabelSize())
    dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
    for i,label in enumerate(labels): dummy.GetXaxis().SetBinLabel(i+1, label)
    dummy.GetXaxis().LabelsOption("u")
    dummy.Draw("HIST")

    for hist in hists:
        hist.Draw("SAME HIST")
    leg.Draw("SAME")

    plotter.aux()
    canvas.RedrawAxis()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/higgsDecays/{hName}.png")
    canvas.SaveAs(f"{outDir}/higgsDecays/{hName}.pdf")


def makePlotHiggsDecays(hName, outName, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True):

    sigs = [[f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays] for y in higgs_decays]

    leg = ROOT.TLegend(.6, 0.99-(len(bkgs)+2)*0.07, .95, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)
    leg.SetNColumns(2)

    hists = []
    for i,sig in enumerate(sigs):
        h_decay = higgs_decays[i]
        h_sig = getHist(hName, sig)
        h_sig.Scale(1./h_sig.Integral())

        h_sig.SetLineColor(colors[h_decay])
        h_sig.SetLineWidth(2)
        h_sig.SetLineStyle(1)

        leg.AddEntry(h_sig, higgs_decays_labels[h_decay], "L")
        hists.append(h_sig)

    '''
    if yMax < 0:
        if logY:
            yMax = math.ceil(max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])*10000)/10.
        else:
            yMax = 1.4*max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])
    '''
    cfg = {

        'logy'              : logY,
        'logx'              : logX,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : xLabel,
        'ytitle'            : yLabel,
            
        'topRight'          : f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}",
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",

    }

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()
    dummy.Draw("HIST") 

    for i,hist in enumerate(hists):
        hist.Draw("SAME HIST")
        print(i, hist.GetBinContent(5))
    leg.Draw("SAME")
    #quit()
    
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}.png")
    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}.pdf")
    canvas.Close()


def makePlot(hName, outName, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9], sig_scale=1):


    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(legPos[0], legPos[1], legPos[2], legPos[3])
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    h_sig = getHist(hName, sigs, rebin)
    if "TH2" in h_sig.ClassName(): h_sig = h_sig.ProjectionX("h_sig")
    h_sig.SetLineColor(sig_color)
    h_sig.SetLineWidth(3)
    h_sig.SetLineStyle(1)
    h_sig.Scale(sig_scale)
    leg.AddEntry(h_sig, sig_legend, "L")

    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(bkgs):

        hist = getHist(hName, bgks_cfg[bkg], rebin)
        if "TH2" in hist.ClassName(): hist = hist.ProjectionX()
        hist.SetName(bkg)
        hist.SetFillColor(bkgs_colors[i])
        hist.SetLineColor(ROOT.kBlack)
        hist.SetLineWidth(1)
        hist.SetLineStyle(1)

        leg.AddEntry(hist, bkgs_legends[i], "F")
        st.Add(hist)
        if h_bkg_tot == None:
            h_bkg_tot = copy.deepcopy(hist)
            h_bkg_tot.SetName("h_bkg_tot")
        else: h_bkg_tot.Add(hist)

    if yMax < 0:
        if logY:
            yMax = math.ceil(max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])*10000)/10.
        else:
            yMax = 1.4*max([h_bkg_tot.GetMaximum(), h_sig.GetMaximum()])

    cfg = {

        'logy'              : logY,
        'logx'              : logX,
        
        'xmin'              : xMin,
        'xmax'              : xMax,
        'ymin'              : yMin,
        'ymax'              : yMax,
            
        'xtitle'            : xLabel,
        'ytitle'            : yLabel,
            
        'topRight'          : f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}",
        'topLeft'           : "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}",

    }

    plotter.cfg = cfg
    canvas = plotter.canvas()
    dummy = plotter.dummy()
    dummy.Draw("HIST") 
    st.Draw("HIST SAME")

    '''
    hTot_err = hTot.Clone("hTot_err")
    hTot_err.SetFillColor(ROOT.kBlack)
    hTot_err.SetMarkerColor(ROOT.kBlack)
    hTot_err.SetFillStyle(3004)
    leg.AddEntry(hTot_err, "Stat. Unc.", "F")
    '''

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetFillColor(0)
    h_bkg_tot.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg_tot.Draw("HIST SAME")
    h_sig.Draw("HIST SAME")
    leg.Draw("SAME")
    
    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"{outDir}/{outName}.png")
    canvas.SaveAs(f"{outDir}/{outName}.pdf")
    canvas.Close()


def significance(hName, xMin=-10000, xMax=10000, reverse=False):

    h_sig = getHist(hName, sigs)
    sig_tot = h_sig.Integral()

    bkgs_procs = []
    for bkg in bkgs:
        bkgs_procs.extend(bgks_cfg[bkg])

    h_bkg = getHist(hName, bkgs_procs)
    x, y, l = [], [], []

    for i in range(1, h_sig.GetNbinsX()+1):
        if reverse:
            iStart = 1
            iEnd = i
        else:
            iStart = i
            iEnd = h_sig.GetNbinsX()+1
        center = h_sig.GetBinCenter(i)
        if center > xMax or center < xMin:
            continue
        sig = h_sig.Integral(iStart, iEnd)
        bkg = h_bkg.Integral(iStart, iEnd)
        significance = sig / (sig + bkg)**0.5
        sig_loss = sig / sig_tot
        print(f"{i} {center:.5f} {significance:.5f} {sig_loss:.5f}")
        x.append(center)
        y.append(significance)
        l.append(sig_loss)


    graph_sign = ROOT.TGraph(len(x), array.array('d', x), array.array('d', y))
    graph_l = ROOT.TGraph(len(x), array.array('d', x), array.array('d', l))
    max_y = max(y)
    max_index = y.index(max_y)
    max_x = x[max_index]
    max_l = l[max_index]

    canvas = ROOT.TCanvas("", "", 800, 800)
    graph_sign.SetMarkerStyle(20)
    graph_sign.SetMarkerColor(ROOT.kBlue)
    graph_sign.GetXaxis().SetRangeUser(xMin, xMax)
    graph_sign.Draw("AP")
    canvas.Update()

    # Add a marker for the maximum point
    max_marker = ROOT.TMarker(max_x, max_y, 20)
    max_marker.SetMarkerColor(ROOT.kRed)
    max_marker.SetMarkerSize(1.5)
    max_marker.Draw()



    rightmax = 1.0
    print(ROOT.gPad.GetUymin())
    print(ROOT.gPad.GetUymax())
    scale =  ROOT.gPad.GetUymax()/rightmax
    rightmin = ROOT.gPad.GetUymin()/ROOT.gPad.GetUymax()
    graph_l.Scale(scale)
    graph_l.SetLineColor(ROOT.kRed)
    graph_l.SetLineWidth(2)
    graph_l.Draw("SAME L")


    axis_r = ROOT.TGaxis(ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),rightmin,rightmax,510,"+L")
    axis_r.SetLineColor(ROOT.kRed)
    axis_r.SetLabelColor(ROOT.kRed)
    axis_r.Draw()


    # Add a text box to indicate the maximum value
    text = ROOT.TLatex()
    text.SetTextSize(0.03)
    text.SetTextColor(ROOT.kBlack)
    text.DrawLatexNDC(0.1, 0.92, f"Max: x = {max_x}, y = {max_y:.5f}, signal loss = {max_l:.5f}")


    suffix = "_reverse" if reverse else ""
    canvas.SaveAs(f"{outDir}/significance/{hName}{suffix}.png")
    canvas.SaveAs(f"{outDir}/significance/{hName}{suffix}.pdf")
    canvas.Close()


if __name__ == "__main__":

    ecm, lumi = 240, "10.8"
    #ecm, lumi = 365, "3"

    inputDir = f"output/h_zh_hadronic/histmaker/"
    outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zh_hadronic/plots_ecm{ecm}/"

    z_decays = ["qq", "bb", "cc", "ss", "ee", "mumu", "tautau", "nunu"]
    z_decays = ["qq", "bb", "cc", "ss"]
    #z_decays = ["ss"]
    higgs_decays = ["bb", "cc", "ss", "gg", "mumu", "tautau", "ZZ", "WW", "Za", "aa", "inv"]
    higgs_decays_labels = {"bb": "H#rightarrowb#bar{b}", "cc": "H#rightarrowc#bar{c}", "ss": "H#rightarrows#bar{s}", "gg": "H#rightarrowgg", "mumu": "H#rightarrow#mu#bar{#mu}", "tautau": "H#rightarrow#tau#bar{#tau}", "ZZ": "H#rightarrowZZ*", "WW": "H#rightarrowWW*", "Za": "H#rightarrowZ#gamma", "aa": "H#rightarrow#gamma#gamma", "inv": "H#rightarrowInv"}
    colors = {"bb": ROOT.kBlack, "cc": ROOT.kBlue , "ss": ROOT.kRed, "gg": ROOT.kGreen+1, "mumu": ROOT.kOrange, "tautau": ROOT.kCyan, "ZZ": ROOT.kGray, "WW": ROOT.kGray+2, "Za": ROOT.kGreen+2, "aa": ROOT.kRed+2, "inv": ROOT.kBlue+2}

    #higgs_decays = ["mumu", "tautau", "ZZ", "WW", "Za", "aa", "inv"]
    #higgs_decays = ["tautau", "ZZ", "WW", "Za"]
    #higgs_decays = ["mumu", "aa", "inv"]

    bkgs = ["WW", "ZZ", "Zgamma"]
    bkgs_legends = ["WW", "ZZ", "Z/#gamma^{*} #rightarrow f#bar{f}+#gamma(#gamma)"]
    # colors from https://github.com/mpetroff/accessible-color-cycles
    bkgs_colors = [ROOT.TColor.GetColor("#f89c20"), ROOT.TColor.GetColor("#5790fc"), ROOT.TColor.GetColor("#964a8b")]
    sig_color = ROOT.TColor.GetColor("#e42536") #

    #z_decays = ["qq", "bb", "cc", "ss", "ee", "mumu", "tautau", "nunu"]
    #higgs_decays = ["bb", "cc", "gg", "ss", "mumu", "ZZ", "WW", "Za", "aa"] ## add tau

    if ecm == 240:
        sig_legend = "ZH(#gamma#gamma)"

        sigs = [f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays for y in higgs_decays]


        bgks_cfg = { 
            "WW"      : ['p8_ee_WW_ecm240'],
            "ZZ"      : ['p8_ee_ZZ_ecm240'],
            "Zgamma"    : ['p8_ee_Zqq_ecm240']
        }

        cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9"]
        labels = ["All events", "Veto leptonic", "Clustering", "m(qq)", "p(qq)", "cos(qq)", "acol", "acop", "WW", "cos(#thetamiss)"]

        makeCutFlow("cutFlow", cuts, labels, 100.)
        makeCutFlowHiggsDecays("cutFlow", cuts, labels)


    # significance
    if True and ecm == 240:
        significance("zqq_m_best_nOne", 50, 130)
        significance("zqq_m_best_nOne", 50, 130, reverse=True)

        significance("zqq_p_best_nOne", 0, 100)
        significance("zqq_p_best_nOne", 0, 100, reverse=True)

        significance("z_costheta_nOne", 0, 1)
        significance("z_costheta_nOne", 0, 1, reverse=True)
        
        significance("acolinearity_nOne", 0, 3.1)
        significance("acolinearity_nOne", 0, 3.1, reverse=True)
        
        significance("acoplanarity_nOne", 0, 3.1)
        significance("acoplanarity_nOne", 0, 3.1, reverse=True)
        
        significance("delta_mWW_nOne", 0, 20)
        significance("cosThetaMiss_nOne", 0.95, 1, reverse=True)

        significance("mva_score", 0, 0.9)
        significance("mva_score", 0, 1, reverse=True)

    makePlotHiggsDecays("best_clustering_idx_nosel", "best_clustering_idx_nosel", xMin=-2, xMax=5, yMin=1e-5, yMax=1e3, xLabel="best_clustering_idx_nosel", yLabel="Events", logY=True)
    makePlotHiggsDecays("zqq_m_best_nosel", "zqq_m_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_m_best_nosel", yLabel="Events", logY=True)
    makePlotHiggsDecays("zqq_p_best_nosel", "zqq_p_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_p_best_nosel", yLabel="Events", logY=True)
    makePlotHiggsDecays("zqq_recoil_m_best_nosel", "zqq_recoil_m_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_recoil_m_best_nosel", yLabel="Events", logY=True)
    makePlotHiggsDecays("z_costheta_nosel", "z_costheta_nosel", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="cos(qq) (GeV)", yLabel="Events", logY=True)
    makePlotHiggsDecays("cosThetaMiss_nOne", "cosThetaMiss_nOne", xMin=0.95, xMax=1, yMin=1e-5, yMax=10, xLabel="cosThetaMiss_nOne", yLabel="Events", logY=True)

    makePlotHiggsDecays("njets_inclusive", "njets_inclusive", xMin=0, xMax=15, yMin=1e-5, yMax=1e3, xLabel="njets_inclusive", yLabel="Events", logY=True)
    makePlotHiggsDecays("njets_inclusive_sel", "njets_inclusive_sel", xMin=0, xMax=15, yMin=1e-5, yMax=1e3, xLabel="njets_inclusive_sel", yLabel="Events", logY=True)
    makePlotHiggsDecays("delta_mWW_nOne", "delta_mWW_nOne", xMin=0, xMax=50, yMin=1e-5, yMax=1e3, xLabel="delta_mWW_nOne", yLabel="Events", logY=True)
    #makePlotHiggsDecays("mva_score", "mva_score", xMin=0, xMax=1, yMin=1e-5, yMax=1e3, xLabel="mva_score", yLabel="Events", logY=True)

    makePlotHiggsDecays("leading_jet_p", "leading_jet_p", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="leading_jet_p", yLabel="Events", logY=True)
    makePlotHiggsDecays("subleading_jet_p", "subleading_jet_p", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="subleading_jet_p", yLabel="Events", logY=True)
    makePlotHiggsDecays("leading_jet_costheta", "leading_jet_costheta", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="leading_jet_costheta", yLabel="Events", logY=True)
    makePlotHiggsDecays("subleading_jet_costheta", "subleading_jet_costheta", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="subleading_jet_costheta", yLabel="Events", logY=True)

    makePlotHiggsDecays("acoplanarity_nOne", "acoplanarity_nOne", xMin=0, xMax=3.14, yMin=1e-5, yMax=10, xLabel="acoplanarity_best", yLabel="Events", logY=True)
    makePlotHiggsDecays("acolinearity_nOne", "acolinearity_nOne", xMin=0, xMax=3.14, yMin=1e-5, yMax=10, xLabel="acolinearity_best", yLabel="Events", logY=True)



    makePlot("zqq_m_best_nosel", "zqq_m_best_nosel", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_p_best_nosel", "zqq_p_best_nosel", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p(qq) (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_recoil_m_best_nosel", "zqq_recoil_m_best_nosel", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="Recoil z(qq) (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("best_clustering_idx_nosel", "best_clustering_idx_nosel", xMin=-2, xMax=5, yMin=1e-1, yMax=-1, xLabel="best_clustering_idx_nosel", yLabel="Events", logY=True, rebin=1)
    makePlot("z_costheta_nosel", "z_costheta_nosel", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos(qq) (GeV)", yLabel="Events", logY=True, rebin=1)


    makePlot("delta_mWW_nosel", "delta_mWW_nosel", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="delta_mWW_nosel", yLabel="Events", logY=True, rebin=1)
    makePlot("delta_mWW_nOne", "delta_mWW_nOne", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="delta_mWW_nOne", yLabel="Events", logY=True, rebin=1)
    makePlot("W1_m_nOne", "W1_m_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W1 (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("W2_m_nOne", "W2_m_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W2 (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("W1_p_nOne", "W1_p_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p(qq) W1 (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("W2_p_nOne", "W2_p_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p(qq) W2 (GeV)", yLabel="Events", logY=True, rebin=1)

    #makePlot("w1_prime", "w1_prime", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W1 (GeV)", yLabel="Events", logY=True, rebin=1)
    #makePlot("w2_prime", "w2_prime", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W2 (GeV)", yLabel="Events", logY=True, rebin=1)


    makePlot("leading_jet_p", "leading_jet_p", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="leading_jet_p", yLabel="Events", logY=True, rebin=1)
    makePlot("subleading_jet_p", "subleading_jet_p", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="subleading_jet_p", yLabel="Events", logY=True, rebin=1)
    makePlot("leading_jet_costheta", "leading_jet_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="leading_jet_costheta", yLabel="Events", logY=True, rebin=1)
    makePlot("subleading_jet_costheta", "subleading_jet_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="subleading_jet_costheta", yLabel="Events", logY=True, rebin=1)

    makePlot("njets_inclusive", "njets_inclusive", xMin=0, xMax=15, yMin=1e-1, yMax=-1, xLabel="njets_inclusive", yLabel="Events", logY=True)
    makePlot("njets_inclusive_sel", "njets_inclusive_sel", xMin=0, xMax=15, yMin=1e-1, yMax=-1, xLabel="njets_inclusive_sel", yLabel="Events", logY=True)

    makePlot("acoplanarity_nOne", "acoplanarity_nOne", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="acoplanarity_best", yLabel="Events", logY=True, rebin=1)
    makePlot("acolinearity_nOne", "acolinearity_nOne", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="acolinearity_best", yLabel="Events", logY=True, rebin=1)


    makePlot("zqq_m_best_nOne", "zqq_m_best_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_p_best_nOne", "zqq_p_best_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p(qq) (GeV)", yLabel="Events", logY=True, rebin=1)

    makePlot("z_costheta_nOne", "z_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos(qq) (GeV)", yLabel="Events", logY=True, rebin=1)


    makePlot("thrust_magn", "thrust_magn", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="thrust_magn", yLabel="Events", logY=True, rebin=1)
    makePlot("thrust_costheta", "thrust_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="MVA score", yLabel="Events", logY=True, rebin=1)

    # final recoil plot
    makePlot("zqq_recoil_m", "zqq_recoil_m", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="Final Recoil z(qq) (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_m", "zqq_m", xMin=00, xMax=200, yMin=1e-1, yMax=-1, xLabel="Final m(qq) (GeV)", yLabel="Events", logY=True, rebin=1)
    
    
    makePlot("mva_score", "mva_score", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="MVA score", yLabel="Events", logY=True, rebin=1)