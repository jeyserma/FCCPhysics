import sys
import os
import math
import copy
import array
import argparse

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

sys.path.insert(0, f'{os.path.dirname(os.path.realpath(__file__))}/../../python')
import plotter


parser = argparse.ArgumentParser()
parser.add_argument("--cat", type=str, help="Category (qq, mumu, ee)", default="qq")
parser.add_argument("--ecm", type=int, help="Center-of-mass energy", default=240)
args = parser.parse_args()


def getHist(hName, procs, rebin=1):
    hist = None
    for proc in procs:
        fInName = f"{inputDir}/{proc}.root"
        if os.path.exists(fInName):
            fIn = ROOT.TFile(fInName)
        elif os.path.exists(fInName.replace("wzp6", "wz3p6")):
            fIn = ROOT.TFile(fInName.replace("wzp6", "wz3p6"))
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
    hist.Rebin(rebin)
    #hist.Scale(3./10.8 * 1.05/1.39)
    return hist

def makeCutFlow(hName="cutFlow", cuts=[], labels=[], sig_scale=1.0, yMin=1e6, yMax=1e10):

    leg = ROOT.TLegend(.65, 0.99-(len(procs))*0.06, .99, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists_yields = []
    significances = []
    h_sig = getHist(hName, procs_cfg[procs[0]])

    hists_yields.append(copy.deepcopy(h_sig))
    h_sig.Scale(sig_scale)
    h_sig.SetLineColor(procs_colors[procs[0]])
    h_sig.SetLineWidth(4)
    h_sig.SetLineStyle(1)
    if sig_scale != 1:
        leg.AddEntry(h_sig, f"{procs_labels[procs[0]]} (#times {int(sig_scale)})", "L")
    else:
        leg.AddEntry(h_sig, procs_labels[procs[0]], "L")

    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(procs[1:]):
        h_bkg = getHist(hName, procs_cfg[bkg])

        if h_bkg_tot == None: h_bkg_tot = h_bkg.Clone("h_bkg_tot")
        else: h_bkg_tot.Add(h_bkg)
        
        h_bkg.SetFillColor(procs_colors[bkg])
        h_bkg.SetLineColor(ROOT.kBlack)
        h_bkg.SetLineWidth(1)
        h_bkg.SetLineStyle(1)

        leg.AddEntry(h_bkg, procs_labels[bkg], "F")
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
        'ymin'              : yMin,
        'ymax'              : yMax ,

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
    dummy.GetXaxis().SetLabelSize(0.75*dummy.GetXaxis().GetLabelSize())
    dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
    for i,label in enumerate(labels): dummy.GetXaxis().SetBinLabel(i+1, label)
    dummy.GetXaxis().LabelsOption("u")
    dummy.Draw("HIST")

    st.Draw("SAME HIST")
    h_bkg_tot.Draw("SAME HIST")
    h_sig.Draw("SAME ][ HIST")
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

        formatted_row = '{:<10} {:<25} ' + ' '.join(['{:<25}']*len(procs))
        print(formatted_row.format(*(["Cut", "Significance"]+procs)))
        print(formatted_row.format(*(["----------"]+["-----------------------"]*(len(procs)+1))))
        for i,cut in enumerate(cuts):
            row = ["Cut %d"%i, "%.3f"%significances[i]]
            for j,histProc in enumerate(hists_yields):
                yield_, err = histProc.GetBinContent(i+1), histProc.GetBinError(i+1)
                row.append("%.4e +/- %.2e" % (yield_, err))

            print(formatted_row.format(*row))
    sys.stdout = out_orig


def acceptance(hName="", lastBin=5, z_decays=[]):
    
    for z_decay in z_decays:
        sigs = [f'wzp6_ee_{x}H_H{y}_ecm{ecm}' for x in z_decay for y in h_decays]
        h = getHist(hName, sigs)
        print(h.Integral())
        eff = h.GetBinContent(lastBin) / h.GetBinContent(1)
        print(z_decay, f"{100.*eff:.2f}")


def makeCutFlowHiggsDecays(hName="cutFlow", outName="", cuts=[], cut_labels=[], yMin=0, yMax=150, xMax=-1, z_decays=[], h_decays=[], h_decays_labels=[], h_decays_colors=[]):
    if outName == "":
        outName = hName

    #sigs = [[f'wzp6_ee_{x}H_H{y}_ecm{ecm}' for x in z_decays] for y in h_decays]
    sigs = [[f'wzp6_ee_{x}H_H{yy}_ecm{ecm}' for x in z_decays for yy in y.split(",")] for y in h_decays]
    #leg = ROOT.TLegend(.2, .925-(len(sigs)/4+1)*0.07, .95, .925)
    #leg.SetBorderSize(0)
    #leg.SetFillStyle(0)
    #leg.SetTextSize(0.03)
    #leg.SetMargin(0.25)
    #leg.SetNColumns(4)

    leg = ROOT.TLegend(.65, .925-(len(sigs)/4+1)*0.1, .95, .925)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.25)


    hists = []
    hist_tot = None
    eff_final, eff_final_err = [], []
    for i,sig in enumerate(sigs):
        h_decay = h_decays[i]
        print(hName, sig)
        h_sig = getHist(hName, sig)
        h_sig.Scale(100./h_sig.GetBinContent(1))

        h_sig.SetLineColor(h_decays_colors[h_decay])
        h_sig.SetLineWidth(3)
        h_sig.SetLineStyle(1)

        leg.AddEntry(h_sig, h_decays_labels[h_decay], "L")
        hists.append(h_sig)

        eff_final.append(h_sig.GetBinContent(len(cuts)))
        eff_final_err.append(h_sig.GetBinError(len(cuts)))
        
        if hist_tot == None:
            hist_tot = h_sig.Clone("h_tot")
        else:
            hist_tot.Add(h_sig)

    hist_tot.Scale(1./len(sigs))
    eff_avg = sum(eff_final) / float(len(eff_final))
    eff_avg, eff_avg_err = hist_tot.GetBinContent(len(cuts)), hist_tot.GetBinError(len(cuts))
    eff_min, eff_max = eff_avg-min(eff_final), max(eff_final)-eff_avg

    ########### PLOTTING ###########
    cfg = {
        'logy'              : False,
        'logx'              : False,
        
        'xmin'              : 0,
        'xmax'              : len(cuts),
        'ymin'              : yMin,
        'ymax'              : yMax ,

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
    for i,label in enumerate(cut_labels): dummy.GetXaxis().SetBinLabel(i+1, label)
    dummy.GetXaxis().LabelsOption("u")
    if xMax > 0:
        dummy.GetXaxis().SetRangeUser(0, xMax)
    dummy.Draw("HIST")

    txt = ROOT.TLatex()
    txt.SetTextSize(0.04)
    txt.SetTextColor(1)
    txt.SetTextFont(42)
    txt.SetNDC()
    #txt.DrawLatex(0.2, 0.2, f"Avg. eff.: {eff_avg:.2f}#pm{eff_avg_err:.2f} %")
    #txt.DrawLatex(0.2, 0.15, f"Min/max: {eff_min:.2f}/{eff_max:.2f}")
    txt.DrawLatex(0.2, 0.2, f"Avg. eff.: {eff_avg:.2f}^{{#plus{eff_max:.2f}}}_{{#minus{eff_min:.2f}}} %")
    txt.Draw("SAME")

    for hist in hists:
        hist.Draw("SAME ][ HIST")
    leg.Draw("SAME")

    plotter.aux()
    canvas.RedrawAxis()
    canvas.Modify()
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}.png")
    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}.pdf")

    out_orig = sys.stdout
    
    with open(f"{outDir}/higgsDecays/{outName}.txt", 'w') as f:
        sys.stdout = f

        formatted_row = '{:<10} ' + ' '.join(['{:<18}']*len(hists))
        print(formatted_row.format(*(["Cut"]+h_decays)))
        print(formatted_row.format(*(["----------"]+["----------------"]*len(hists))))
        for i,cut in enumerate(cuts):
            row = ["Cut %d"%i]
            for j,histProc in enumerate(hists):
                yield_, err = hists[j].GetBinContent(i+1), hists[j].GetBinError(i+1)
                row.append("%.2f +/- %.2f" % (yield_, err))
            print(formatted_row.format(*row))
        print("\n")
        print(f"Average: {eff_avg:.3f} +/- {eff_avg_err:.3f}")
        print(f"Min/max: {eff_min:.3f}/{eff_max:.3f}")
    sys.stdout = out_orig
    del canvas

    ## make final efficiency plot eff_final, eff_final_err
    sel_eff_avg = eff_avg 
    sel_eff_tot_err = 0
    if args.cat == "mumu" and args.ecm == 240:
        xMin, xMax = 68, 78
        #xMin, xMax = 68-16, 78 # cos theta miss
        if len(z_decays) > 1:
            xMin, xMax = 0, 5
    if args.cat == "ee" and args.ecm == 240:
        xMin, xMax = 62, 68
        #xMin, xMax = 62-16, 68 # cos theta miss
        if len(z_decays) > 1:
            xMin, xMax = 0, 5
    if args.cat == "qq":
        if args.ecm == 240:
            xMin, xMax = 55, 85
        else:
            xMin, xMax = 45, 75
    h_pulls = ROOT.TH2F("pulls", "pulls", (xMax-xMin)*10, xMin, xMax, len(sigs)+1, 0, len(sigs)+1)
    g_pulls = ROOT.TGraphErrors(len(sigs)+1)

    ip = 0 # counter for TGraph
    g_pulls.SetPoint(ip, eff_avg, float(ip) + 0.5)
    g_pulls.SetPointError(ip, eff_avg_err, 0.)
    h_pulls.GetYaxis().SetBinLabel(ip+1, "Average")
    ip += 1

    for i,sig in enumerate(sigs):
        h_decay = h_decays[i]
        g_pulls.SetPoint(ip, eff_final[i], float(ip) + 0.5)
        g_pulls.SetPointError(ip, eff_final_err[i], 0.)
        h_pulls.GetYaxis().SetBinLabel(ip+1, h_decays_labels[h_decay])
        ip += 1

    canvas = ROOT.TCanvas("c", "c", 800, 800)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.1)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetFillStyle(4000) # transparency?
    canvas.SetGrid(1, 0)
    canvas.SetTickx(1)

    #xTitle = "Selection efficiency Z(#mu^{#plus}#mu^{#minus})H  (%)"
    #if flavor == "ee": xTitle = "Selection efficiency Z(e^{#plus}e^{#minus})H  (%)"
    xTitle = "Selection efficiency (%)"

    h_pulls.GetXaxis().SetTitleSize(0.04)
    h_pulls.GetXaxis().SetLabelSize(0.035)
    h_pulls.GetXaxis().SetTitle(xTitle)
    h_pulls.GetXaxis().SetTitleOffset(1)
    h_pulls.GetYaxis().SetLabelSize(0.055)
    h_pulls.GetYaxis().SetTickLength(0)
    h_pulls.GetYaxis().LabelsOption('v')
    h_pulls.SetNdivisions(506, 'XYZ')
    h_pulls.Draw("HIST 0")

    if args.cat == "mumu" and args.ecm == 240:
        sel_avg, sel_err = eff_avg, 0.68
    elif args.cat == "ee" and args.ecm == 240:
        sel_avg, sel_err = eff_avg, 0.81
    elif args.cat == "qq" and args.ecm == 240:
        sel_avg, sel_err = eff_avg, 0.38
    else:
        pass

    maxx = len(sigs)+1
    line = ROOT.TLine(eff_avg, 0, eff_avg, maxx)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    line.Draw("SAME")




    shade = ROOT.TGraph()
    shade.SetPoint(0, sel_avg-sel_err, 0)
    shade.SetPoint(1, sel_avg+sel_err, 0)
    shade.SetPoint(2, sel_avg+sel_err, maxx)
    shade.SetPoint(3, sel_avg-sel_err, maxx)
    shade.SetPoint(4, sel_avg-sel_err, 0)
    #shade.SetFillStyle(3013)
    shade.SetFillColor(16)
    shade.SetFillColorAlpha(16, 0.35);
    shade.Draw("SAME F")

    g_pulls.SetMarkerSize(1.2)
    g_pulls.SetMarkerStyle(20)
    g_pulls.SetLineWidth(2)
    g_pulls.Draw('P0 SAME')

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(30) # 0 special vertical aligment with subscripts
    latex.DrawLatex(0.95, 0.925, f"#sqrt{{s}} = {ecm} GeV, {lumi} ab^{{#minus1}}")

    latex.SetTextAlign(13)
    latex.SetTextFont(42)
    latex.SetTextSize(0.045)
    latex.DrawLatex(0.15, 0.96, "#bf{FCC-ee} #scale[0.7]{#it{Simulation}}")

    #latex.SetTextAlign(13)
    #latex.SetTextFont(42)
    #latex.SetTextSize(0.045)
    #latex.SetTextColor(ROOT.kGray+1)
    #latex.DrawLatex(0.2, 0.17, "avg. #pm 0.1 %")

    txt = ROOT.TLatex()
    txt.SetTextSize(0.04)
    txt.SetTextColor(1)
    txt.SetTextFont(42)
    txt.SetNDC()
    #txt.DrawLatex(0.2, 0.2, f"Avg. eff.: {eff_avg:.2f}#pm{eff_avg_err:.2f} %")
    #txt.DrawLatex(0.2, 0.15, f"Min/max: {eff_min:.2f}/{eff_max:.2f}")
    txt.DrawLatex(0.2, 0.2, f"Avg. eff.: {eff_avg:.2f}^{{#plus{eff_max:.2f}}}_{{#minus{eff_min:.2f}}} %")
    txt.Draw("SAME")

    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}_finalSelection.png")
    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}_finalSelection.pdf")


def makePlotHiggsDecays(hName, outName="", xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="", yLabel="Events", logX=False, logY=True, rebin=1, xLabels=[], h_decays=[], h_decays_labels=[], h_decays_colors=[], legPos="right"):
    if outName == "":
        outName = hName

    #sigs = [[f'wzp6_ee_{x}H_H{yy}_ecm{ecm}' for x in z_decays for yy in y.split(",")] for y in h_decays]
    #leg = ROOT.TLegend(.2, .925-(len(sigs)/4+1)*0.07, .95, .925)
    #leg.SetBorderSize(0)
    #leg.SetFillStyle(0)
    #leg.SetTextSize(0.03)
    #leg.SetMargin(0.25)
    #leg.SetNColumns(4)

    if outName in ['njets_inclusive_c', 'best_clustering_idx_nosel_c']:
        sigs = [[yy for x in z_decays for yy in y.split(",")] for y in h_decays]
    else:
        sigs = [[f'wzp6_ee_{x}H_H{yy}_ecm{ecm}' for x in z_decays for yy in y.split(",")] for y in h_decays]
    if legPos == "left":
        leg = ROOT.TLegend(.2, .925-(len(sigs)/4+1)*0.1, .5, .925)
    else:
        leg = ROOT.TLegend(.65, .925-(len(sigs)/4+1)*0.1, .95, .925)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.25)

    hists = []
    for i,sig in enumerate(sigs):
        h_decay = h_decays[i]
        h_sig = getHist(hName, sig)
        h_sig.Rebin(rebin)
        if h_sig.Integral() > 0:
            h_sig.Scale(1./h_sig.Integral())

        h_sig.SetLineColor(h_decays_colors[h_decay])
        h_sig.SetLineWidth(3)
        h_sig.SetLineStyle(1)

        leg.AddEntry(h_sig, h_decays_labels[h_decay], "L")
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
    dummy = plotter.dummy(1 if len(xLabels) == 0 else len(xLabels))
    if len(xLabels) > 0:
        dummy.GetXaxis().SetLabelSize(0.8*dummy.GetXaxis().GetLabelSize())
        dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
        for i,label in enumerate(xLabels): dummy.GetXaxis().SetBinLabel(i+1, label)
        if not "njets_" in hName:
            dummy.GetXaxis().LabelsOption("u")
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


def makePlotSignalRatios(hName, outName="", xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", rebin=1):

    if outName == "":
        outName = hName

    sigs = [[f'wzp6_ee_{x}H_H{y}_ecm{ecm}' for x in z_decays] for y in h_decays]
    leg = ROOT.TLegend(.2, .925-(len(sigs)/4+1)*0.07, .95, .925)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.25)
    leg.SetNColumns(4)

    hists = []
    h_tot = None
    for i,sig in enumerate(sigs):
        h_decay = h_decays[i]
        h_sig = getHist(hName, sig)
        h_sig.Rebin(rebin)
        #h_sig.Scale(1./h_sig.Integral())

        h_sig.SetLineColor(h_decays_colors[h_decay])
        h_sig.SetLineWidth(2)
        h_sig.SetLineStyle(1)

        leg.AddEntry(h_sig, h_decays_labels[h_decay], "L")
        hists.append(h_sig)
        
        if h_tot == None:
            h_tot = h_sig.Clone("h_tot")
        else:
            h_tot.Add(h_sig)

    h_tot.Scale(1./h_tot.Integral())

    for h in hists:
        h.Scale(1./h.Integral())
        h.Divide(h_tot)

    cfg = {

        'logy'              : False,
        'logx'              : False,

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
    leg.Draw("SAME")

    canvas.SetGrid()
    canvas.Modify()
    canvas.Update()

    plotter.aux()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()

    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}.png")
    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}.pdf")
    canvas.Close()


def makePlot(hName, outName="", xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9], sig_scale=1, xLabels=[]):

    if logY:
        #ROOT.TGaxis.SetMaxDigits(2)
        ROOT.TGaxis.SetExponentOffset(0, 0, "y")
    else:
        #ROOT.TGaxis.SetMaxDigits(2)
        ROOT.TGaxis.SetExponentOffset(-0.07, 0.01, "y")

    if outName == "":
        outName = hName
    
    alt = hName
    if hName == "zqq_mva_low" or hName == "recoil_mva_low" or hName == "zqq_mva_high" or hName == "recoil_mva_high":
        hName = "zqq_recoil_m_mqq_mva"

    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(.65, 0.99-(len(procs))*0.06, .99, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    h_sig = getHist(hName, procs_cfg[procs[0]], rebin=rebin)
    if alt == "zqq_mva_low":
        h_sig = h_sig.ProjectionY("sig", 1, h_sig.GetXaxis().GetNbins(), 1, 1)
    if alt == "recoil_mva_low":
        h_sig = h_sig.ProjectionX("sig", 1, h_sig.GetYaxis().GetNbins(), 1, 1)
    if alt == "zqq_mva_high":
        h_sig = h_sig.ProjectionY("sig", 1, h_sig.GetXaxis().GetNbins(), 2, 2)
    if alt == "recoil_mva_high":
        h_sig = h_sig.ProjectionX("sig", 1, h_sig.GetYaxis().GetNbins(), 2, 2)
    h_sig.SetLineColor(procs_colors[procs[0]])
    h_sig.SetLineWidth(3)
    h_sig.SetLineStyle(1)
    h_sig.Scale(sig_scale)
    if sig_scale != 1:
        leg.AddEntry(h_sig, f"{procs_labels[procs[0]]} (#times {int(sig_scale)})", "L")
    else:
        leg.AddEntry(h_sig, procs_labels[procs[0]], "L")

    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,bkg in enumerate(procs[1:]):
        h_bkg = getHist(hName, procs_cfg[bkg], rebin=rebin)
        if alt == "zqq_mva_low":
            h_bkg = h_bkg.ProjectionY(bkg, 1, h_bkg.GetXaxis().GetNbins(), 1, 1)
        if alt == "recoil_mva_low":
            h_bkg = h_bkg.ProjectionX(bkg, 1, h_bkg.GetYaxis().GetNbins(), 1, 1)
        if alt == "zqq_mva_high":
            h_bkg = h_bkg.ProjectionY(bkg, 1, h_bkg.GetXaxis().GetNbins(), 2, 2)
        if alt == "recoil_mva_high":
            h_bkg = h_bkg.ProjectionX(bkg, 1, h_bkg.GetYaxis().GetNbins(), 2, 2)


        if h_bkg_tot == None: h_bkg_tot = h_bkg.Clone("h_bkg_tot")
        else: h_bkg_tot.Add(h_bkg)
        
        h_bkg.SetFillColor(procs_colors[bkg])
        h_bkg.SetLineColor(ROOT.kBlack)
        h_bkg.SetLineWidth(1)
        h_bkg.SetLineStyle(1)

        leg.AddEntry(h_bkg, procs_labels[bkg], "F")
        st.Add(h_bkg)


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
    dummy = plotter.dummy(1 if len(xLabels) == 0 else len(xLabels))
    if len(xLabels) > 0:
        dummy.GetXaxis().SetLabelSize(0.8*dummy.GetXaxis().GetLabelSize())
        dummy.GetXaxis().SetLabelOffset(1.3*dummy.GetXaxis().GetLabelOffset())
        for i,label in enumerate(xLabels): dummy.GetXaxis().SetBinLabel(i+1, label)
        dummy.GetXaxis().LabelsOption("u")
    dummy.Draw("HIST")
    st.Draw("HIST SAME")

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetFillColor(0)
    h_bkg_tot.SetLineWidth(2)
    #hTot_err.Draw("E2 SAME")
    h_bkg_tot.Draw("HIST SAME")
    h_sig.Draw("HIST SAME")
    #if hName == "best_clustering_idx_nosel":
    #    h_sig.Draw("HIST [ SAME")
    #else:
    #    h_sig.Draw("HIST ][ SAME")
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

    h_sig = getHist(hName, procs_cfg[procs[0]])
    sig_tot = h_sig.Integral()

    bkgs_procs = []
    for i,bkg in enumerate(procs[1:]):
        bkgs_procs.extend(procs_cfg[bkg])

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
        if (sig+bkg) <= 0 or sig_tot <= 0:
            significance = 0
            sig_loss = 0
        else:
            significance = sig / (sig + bkg)**0.5
            sig_loss = sig / sig_tot
        #significance = sig / (sig + bkg)**0.5
        #sig_loss = sig / sig_tot
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

    ecm = args.ecm
    cat = args.cat
    lumi = "10.8" if ecm == 240 else "3.12"


    z_decays_all = ["qq", "bb", "cc", "ss", "ee", "mumu", "tautau", "nunu"]
    z_decays = ["qq", "bb", "cc", "ss", "ee", "mumu", "tautau", "nunu"]
    h_decays = ["bb", "cc", "gg", "ss", "mumu", "tautau", "ZZ", "WW", "Za", "aa", "inv"]
    h_decays_no_mumu = ["bb", "cc", "gg", "ss", "tautau", "ZZ", "WW", "Za", "aa", "inv"]
    h_decays_rare_bsm = ["uu", "dd", "ee", "taue", "taumu", "bs", "bd", "cu", "sd"]
    h_decays_rare_bsm = ["uu", "dd", "ee", "taue", "taumu"]
    

    colors_hex = [
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#009E73",  # bluish green
    "#F0E442",  # yellow
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    ]
    lc = [ROOT.TColor.GetColor(c) for c in colors_hex]

    h_decays_labels = {"bb": "H#rightarrowb#bar{b}", "cc": "H#rightarrowc#bar{c}", "ss": "H#rightarrows#bar{s}", "gg": "H#rightarrowgg", "mumu": "H#rightarrow#mu^{#plus}#mu^{#minus}", "tautau": "H#rightarrow#tau^{#plus}#tau^{#minus}", "ZZ": "H#rightarrowZZ*", "WW": "H#rightarrowWW*", "Za": "H#rightarrowZ#gamma", "aa": "H#rightarrow#gamma#gamma", "inv": "H#rightarrowInv"}
    h_decays_colors = {"bb": ROOT.kBlack, "cc": ROOT.kBlue , "ss": ROOT.kRed, "gg": ROOT.kGreen+1, "mumu": ROOT.kOrange, "tautau": ROOT.kCyan, "ZZ": ROOT.kGray, "WW": ROOT.kGray+2, "Za": ROOT.kGreen+2, "aa": ROOT.kRed+2, "inv": ROOT.kBlue+2}


    # separate aggregated plots for Higgs decays and backgrounds
    h_decays_agg = ["bb,cc,ss", "gg", "WW,ZZ", "tautau", "Za", "aa,mumu,inv"]
    h_decays_agg_a = ["bb,cc,ss", "gg", "WW,ZZ"]
    h_decays_agg_b = ["tautau", "Za", "aa,mumu,inv"]
    h_decays_agg_c = ["p8_ee_WW_ecm240,p8_ee_WW_mumu_ecm240,p8_ee_WW_ee_ecm240", "p8_ee_ZZ_ecm240", "wz3p6_ee_tautau_ecm240,wz3p6_ee_mumu_ecm240,wz3p6_ee_ee_Mee_30_150_ecm240,wz3p6_ee_uu_ecm240,wz3p6_ee_dd_ecm240,wz3p6_ee_cc_ecm240,wz3p6_ee_ss_ecm240,wz3p6_ee_bb_ecm240,wz3p6_ee_nunu_ecm240"]
    h_decays_labels_agg = {"bb,cc,ss": "H#rightarrowbb/cc/ss", "gg": "H#rightarrowgg", "WW,ZZ": "H#rightarrowWW*/ZZ*", "tautau": "H#rightarrow#tau^{#plus}#tau^{#minus}", "Za": "H#rightarrowZ#gamma", "aa,mumu,inv": "H#rightarrow#gamma#gamma/#mu^{#plus}#mu^{#minus}/Inv"}
    h_decays_colors_agg = {"bb,cc,ss": lc[0], "gg": lc[1], "WW,ZZ": lc[2], "tautau": lc[3], "Za": lc[4], "aa,mumu,inv": lc[5]}
    h_decays_colors_agg = {"bb,cc,ss": ROOT.TColor.GetColor("#a9c593"), "gg": ROOT.TColor.GetColor("#a15f86"), "WW,ZZ":  ROOT.TColor.GetColor("#3d4c74"), "tautau": ROOT.TColor.GetColor("#a9c593"), "Za": ROOT.TColor.GetColor("#a15f86"), "aa,mumu,inv":  ROOT.TColor.GetColor("#3d4c74")}
    h_decays_colors_agg = {"bb,cc,ss": ROOT.TColor.GetColor("#E69F00"), "gg": ROOT.TColor.GetColor("#56B4E9"), "WW,ZZ":  ROOT.TColor.GetColor("#009E73"), "tautau": ROOT.TColor.GetColor("#F0E442"), "Za": ROOT.TColor.GetColor("#0072B2"), "aa,mumu,inv":  ROOT.TColor.GetColor("#D55E00")}
    
    h_decays_labels_agg["p8_ee_WW_ecm240,p8_ee_WW_mumu_ecm240,p8_ee_WW_ee_ecm240"] = "WW"
    h_decays_labels_agg["p8_ee_ZZ_ecm240"] = "ZZ"
    h_decays_labels_agg["wz3p6_ee_tautau_ecm240,wz3p6_ee_mumu_ecm240,wz3p6_ee_ee_Mee_30_150_ecm240,wz3p6_ee_uu_ecm240,wz3p6_ee_dd_ecm240,wz3p6_ee_cc_ecm240,wz3p6_ee_ss_ecm240,wz3p6_ee_bb_ecm240,wz3p6_ee_nunu_ecm240"] = "Z/#gamma^{*}"
    h_decays_colors_agg["p8_ee_WW_ecm240,p8_ee_WW_mumu_ecm240,p8_ee_WW_ee_ecm240"] = ROOT.TColor.GetColor("#a9c593")
    h_decays_colors_agg["p8_ee_ZZ_ecm240"] = ROOT.TColor.GetColor("#a15f86")
    h_decays_colors_agg["wz3p6_ee_tautau_ecm240,wz3p6_ee_mumu_ecm240,wz3p6_ee_ee_Mee_30_150_ecm240,wz3p6_ee_uu_ecm240,wz3p6_ee_dd_ecm240,wz3p6_ee_cc_ecm240,wz3p6_ee_ss_ecm240,wz3p6_ee_bb_ecm240,wz3p6_ee_nunu_ecm240"] = ROOT.TColor.GetColor("#3d4c74")


    '''
            "ZH"        : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZqqH"      : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZmumuH"    : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZnunuH"    : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZeeH"      : ROOT.TColor.GetColor("#e42536"),#e42536
        "WW"        : ROOT.TColor.GetColor("#a9c593"),#f89c20
        "ZZ"        : ROOT.TColor.GetColor("#a15f86"),#5790fc
        "Zgamma"    : ROOT.TColor.GetColor("#3d4c74"),#964a8b
        "Zqqgamma"  : ROOT.TColor.GetColor("#3d4c74"),#964a8b
        "Rare"      : ROOT.TColor.GetColor("#9c9ca1") #9c9ca1
    '''

    #h_decays_labels_agg = {"bb,cc,ss": "bb,cc,ss", "gg": "gg", "WW,ZZ": "WW,ZZ", "tautau": "tau", "Za": "Za", "aa,mumu,inv": "aa,mumu,inv"}
    #h_decays_colors_agg = {"bb,cc,ss": lc[0], "gg": lc[1], "WW,ZZ": lc[2], "tautau": lc[3], "Za": lc[4], "aa,mumu,inv": lc[5]}

    h_decays_rare_bsm_labels = {"uu": "H#rightarrowu#bar{u}", "dd": "H#rightarrowd#bar{d}", "ee": "H#rightarrowee", "taue": "H#rightarrow#taue", "taumu": "H#rightarrow#tau#mu", "bs": "H#rightarrowbs", "bd": "H#rightarrowbd", "cu": "H#rightarrowcu", "sd": "H#rightarrowsd", "inv": "H#rightarrowInv"}
    h_decays_rare_bsm_colors = {"uu": ROOT.kBlack, "dd": ROOT.kBlue , "ee": ROOT.kRed, "taue": ROOT.kGreen+1, "taumu": ROOT.kOrange, "bs": ROOT.kGray, "bd": ROOT.kGray+2, "cu": ROOT.kGreen+2, "sd": ROOT.kRed+2}


    procs_cfg_240 = {
        "ZH"        : [f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays for y in h_decays],
        "ZqqH"      : [f'wzp6_ee_{x}H_H{y}_ecm240' for x in ["qq", "bb", "cc", "ss"] for y in h_decays],
        "ZmumuH"    : [f'wzp6_ee_mumuH_H{y}_ecm240' for y in h_decays],
        "ZnunuH"    : [f'wzp6_ee_nunuH_H{y}_ecm240' for y in h_decays],
        "ZeeH"      : [f'wzp6_ee_{x}H_H{y}_ecm240' for x in ["ee"] for y in h_decays],
        "WW"        : ['p8_ee_WW_ecm240', 'p8_ee_WW_mumu_ecm240', 'p8_ee_WW_ee_ecm240'], # p8_ee_WW_mumu_ecm240 p8_ee_WW_ecm240
        "ZZ"        : ['p8_ee_ZZ_ecm240'],
        "Zgamma"    : ['wz3p6_ee_tautau_ecm240', 'wz3p6_ee_mumu_ecm240', 'wz3p6_ee_ee_Mee_30_150_ecm240', 'wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_nunu_ecm240'],
        #"Zgamma"    : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240', 'wz3p6_ee_nunu_ecm240'],
        #"Zgamma"    : ['p8_ee_Zqq_ecm240'],
        "Zqqgamma"  : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240'],
        "Rare"      : ['wzp6_egamma_eZ_Zmumu_ecm240', 'wzp6_gammae_eZ_Zmumu_ecm240', 'wzp6_gaga_mumu_60_ecm240', 'wzp6_egamma_eZ_Zee_ecm240', 'wzp6_gammae_eZ_Zee_ecm240', 'wzp6_gaga_ee_60_ecm240', 'wzp6_gaga_tautau_60_ecm240', 'wzp6_ee_nuenueZ_ecm240'],
    }

    procs_cfg_365 = {
        "ZH"        : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in z_decays for y in h_decays],
        #"ZqqH"      : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in ["qq", "bb", "cc", "ss"] for y in h_decays],
        "ZqqH"      : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in ["qq", "bb", "cc", "ss"] for y in h_decays],
        "ZmumuH"    : [f'wzp6_ee_mumuH_H{y}_ecm365' for y in h_decays],
        "ZeeH"      : [f'wzp6_ee_{x}H_H{y}_ecm365' for x in ["ee"] for y in h_decays],
        "WW"        : ['p8_ee_WW_ecm365', 'p8_ee_WW_mumu_ecm365', 'p8_ee_WW_ee_ecm365'], # p8_ee_WW_mumu_ecm365 p8_ee_WW_ecm365
        "ZZ"        : ['p8_ee_ZZ_ecm365'],
        "Zgamma"    : ['wz3p6_ee_tautau_ecm365', 'wz3p6_ee_mumu_ecm365', 'wz3p6_ee_ee_Mee_30_150_ecm365', 'wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_nunu_ecm365'],
        #"Zgamma"    : ['wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365', 'wz3p6_ee_nunu_ecm365'],
        #"Zgamma"    : ['p8_ee_Zqq_ecm365'],
        "Zqqgamma"  : ['wz3p6_ee_uu_ecm365', 'wz3p6_ee_dd_ecm365', 'wz3p6_ee_cc_ecm365', 'wz3p6_ee_ss_ecm365', 'wz3p6_ee_bb_ecm365'],
        "Rare"      : ['wzp6_egamma_eZ_Zmumu_ecm365', 'wzp6_gammae_eZ_Zmumu_ecm365', 'wzp6_gaga_mumu_60_ecm365', 'wzp6_egamma_eZ_Zee_ecm365', 'wzp6_gammae_eZ_Zee_ecm365', 'wzp6_gaga_ee_60_ecm365', 'wzp6_gaga_tautau_60_ecm365', 'wzp6_ee_nuenueZ_ecm365'],
    }
    
    if ecm == 240:
        procs_cfg = procs_cfg_240
    if ecm == 365:
        procs_cfg = procs_cfg_365

    # cubehelix colors
    #3d4c74
    #73528c
    #a15f86
    #b67b76
    #b2a275
    #a9c593
    #b1dcc2


    procs_labels = {
        "ZH"        : "ZH",
        "ZqqH"      : "Z(qq)H",
        "ZmumuH"    : "Z(#mu^{+}#mu^{#minus})H",
        "ZnunuH"    : "Z(#nu#nu)H",
        "ZeeH"      : "Z(e^{+}e^{#minus})H",
        "WW"        : "WW",
        "ZZ"        : "ZZ",
        "Zgamma"    : "Z/#gamma^{*}", #  #rightarrow f#bar{f}+#gamma(#gamma)
        "Zqqgamma"  : "Z/#gamma^{*}", #  #rightarrow f#bar{f}+#gamma(#gamma)
        "Rare"      : "Rare"
    }

    # colors from https://github.com/mpetroff/accessible-color-cycles
    procs_colors = {
        "ZH"        : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZqqH"      : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZmumuH"    : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZnunuH"    : ROOT.TColor.GetColor("#e42536"),#e42536
        "ZeeH"      : ROOT.TColor.GetColor("#e42536"),#e42536
        "WW"        : ROOT.TColor.GetColor("#a9c593"),#f89c20
        "ZZ"        : ROOT.TColor.GetColor("#a15f86"),#5790fc
        "Zgamma"    : ROOT.TColor.GetColor("#3d4c74"),#964a8b
        "Zqqgamma"  : ROOT.TColor.GetColor("#3d4c74"),#964a8b
        "Rare"      : ROOT.TColor.GetColor("#9c9ca1") #9c9ca1
    }


    if cat == "qq":

        inputDir = f"output/h_zh_hadronic/histmaker/ecm{ecm}"
        outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zh_hadronic/plots_ecm{ecm}_paper/"

        z_decays = ["qq", "bb", "cc", "ss"]
        procs = ["ZqqH", "WW", "ZZ", "Zgamma", "Rare"] # first must be signal

        if ecm == 240:
            cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8"]
            cut_labels = ["All events", "Veto leptonic", "Clustering", "20 < m_{qq} < 140", "20 < p_{qq} < 90", "cos(qq) < 0.85", "acol(qq) > 0.35", "WW pair mass", "|cos#theta_{miss}| < 0.995"]
        else:
            cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9"]
            cut_labels = ["All events", "Veto leptonic", "Clustering", "20 < m_{qq} < 200", "60 < p_{qq} < 160", "cos(qq) < 0.85", "acol(qq) > 0.35", "WW pair mass", "|cos#theta_{miss}| < 0.995", "|T| < 0.85"]

        makeCutFlowHiggsDecays("cutFlow", outName="cutFlow_agg", cuts=cuts, cut_labels=cut_labels, yMin=60, yMax=120, z_decays=z_decays, h_decays=h_decays_agg, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)
        makeCutFlow("cutFlow", cuts, cut_labels, 100., yMin=1e7 if ecm==240 else 1e9, yMax=1e10 if ecm==240 else 1e9)

        makeCutFlowHiggsDecays("cutFlow", cuts=cuts, cut_labels=cut_labels, yMin=40, yMax=150, z_decays=z_decays, h_decays=h_decays, h_decays_labels=h_decays_labels, h_decays_colors=h_decays_colors)

        #makeCutFlowHiggsDecays("cutFlow", outName="cutFlow_vvH", cuts=cuts, cut_labels=cut_labels, yMin=0, yMax=150, z_decays=["nunu"], h_decays=h_decays, h_decays_labels=h_decays_labels, h_decays_colors=h_decays_colors)

        # significance 
        if False:
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

            significance("thrust_magn", 0.8, 1, reverse=True)
    

            significance("mva_score", 0, 1)
            significance("mva_score", 0, 1, reverse=True)


        # signal plots for paper
        makePlot("zqq_mva_low", xMin=60, xMax=120, yMin=1e2, yMax=-1, xLabel="Jet pair invariant mass (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("zqq_mva_high", xMin=60, xMax=120, yMin=1e2, yMax=-1, xLabel="Jet pair invariant mass (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("zqq_mva_low", outName="zqq_mva_low_noLog", xMin=60, xMax=120, yMin=0, yMax=4000e3, xLabel="Jet pair invariant mass (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=10)
        makePlot("zqq_mva_high", outName="zqq_mva_high_noLog", xMin=60, xMax=120, yMin=0, yMax=800e3, xLabel="Jet pair invariant mass (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=10)

        makePlot("recoil_mva_low", xMin=100, xMax=150, yMin=1e2, yMax=-1, xLabel="Recoil mass (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("recoil_mva_high", xMin=100, xMax=150, yMin=1e2, yMax=-1, xLabel="Recoil mass (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("recoil_mva_low", outName="recoil_mva_low_noLog", xMin=100, xMax=150, yMin=0, yMax=5000e3, xLabel="Recoil mass (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=10)
        makePlot("recoil_mva_high", outName="recoil_mva_high_noLog", xMin=100, xMax=150, yMin=0, yMax=800e3, xLabel="Recoil mass (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=10)

        makePlotHiggsDecays("best_clustering_idx_nosel", outName="best_clustering_idx_nosel_a",  xMin=-1, xMax=4, yMin=0, yMax=1.2, xLabels=["No pairs", "Inclusive", "Exclusive N=2", "Exclusive N=4", "Exclusive N=6"], yLabel="Events", logY=False, h_decays=h_decays_agg_a, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)
        makePlotHiggsDecays("best_clustering_idx_nosel", outName="best_clustering_idx_nosel_b",  xMin=-1, xMax=4, yMin=0, yMax=1.2, xLabels=["No pairs", "Inclusive", "Exclusive N=2", "Exclusive N=4", "Exclusive N=6"], yLabel="Events", logY=False, h_decays=h_decays_agg_b, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)
        makePlotHiggsDecays("best_clustering_idx_nosel", outName="best_clustering_idx_nosel_c",  xMin=-1, xMax=4, yMin=0, yMax=1.2, xLabels=["No pairs", "Inclusive", "Exclusive N=2", "Exclusive N=4", "Exclusive N=6"], yLabel="Events", logY=False, h_decays=h_decays_agg_c, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)

        makePlotHiggsDecays("njets_inclusive", outName="njets_inclusive_a", xMin=0, xMax=11, yMin=0, yMax=0.7, yLabel="Events (normalized)", xLabel="Number of jets (inclusive)", xLabels=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], logY=False, h_decays=h_decays_agg_a, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)
        makePlotHiggsDecays("njets_inclusive", outName="njets_inclusive_b", xMin=0, xMax=11, yMin=0, yMax=0.7, yLabel="Events (normalized)", xLabel="Number of jets (inclusive)", xLabels=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], logY=False, h_decays=h_decays_agg_b, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)
        makePlotHiggsDecays("njets_inclusive", outName="njets_inclusive_c", xMin=0, xMax=11, yMin=0, yMax=0.7, yLabel="Events (normalized)", xLabel="Number of jets (inclusive)", xLabels=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"], logY=False, h_decays=h_decays_agg_c, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)


        makePlotHiggsDecays("zqq_m_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_m_best_nosel", yLabel="Events", logY=True)
        makePlotHiggsDecays("zqq_p_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_p_best_nosel", yLabel="Events", logY=True)
        makePlotHiggsDecays("zqq_recoil_m_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_recoil_m_best_nosel", yLabel="Events", logY=True)
        makePlotHiggsDecays("z_costheta_nosel", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="cos(qq) (GeV)", yLabel="Events", logY=True)
        makePlotHiggsDecays("cosThetaMiss_nOne", xMin=0.95, xMax=1, yMin=1e-5, yMax=10, xLabel="cosThetaMiss_nOne", yLabel="Events", logY=True)

        #makePlotHiggsDecays("njets_inclusive", xMin=0, xMax=15, yMin=1e-5, yMax=1e3, xLabel="njets_inclusive", yLabel="Events", logY=True)
        makePlotHiggsDecays("njets_inclusive_sel", xMin=0, xMax=15, yMin=1e-5, yMax=1e3, xLabel="njets_inclusive_sel", yLabel="Events", logY=True)
        makePlotHiggsDecays("delta_mWW_nOne", xMin=0, xMax=50, yMin=1e-5, yMax=1e3, xLabel="delta_mWW_nOne", yLabel="Events", logY=True)
        makePlotHiggsDecays("mva_score", xMin=0, xMax=1, yMin=1e-5, yMax=1e1, xLabel="MVA score", yLabel="Events", logY=True, rebin=2)
        makePlotHiggsDecays("mva_score", outName="mva_score_agg", xMin=0, xMax=1, yMin=1e-5, yMax=1e1, xLabel="MVA score", yLabel="Events", logY=True, rebin=2, h_decays=h_decays_agg, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)
        makePlotHiggsDecays("mva_score", outName="mva_score_agg_noLog", xMin=0, xMax=1, yMin=0, yMax=0.02, xLabel="MVA score", yLabel="Events", logY=False, rebin=2, h_decays=h_decays_agg, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg, legPos="left")

        makePlotHiggsDecays("leading_jet_p", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="leading_jet_p", yLabel="Events", logY=True)
        makePlotHiggsDecays("subleading_jet_p", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="subleading_jet_p", yLabel="Events", logY=True)
        makePlotHiggsDecays("leading_jet_costheta", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="leading_jet_costheta", yLabel="Events", logY=True)
        makePlotHiggsDecays("subleading_jet_costheta", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="subleading_jet_costheta", yLabel="Events", logY=True)

        makePlotHiggsDecays("acoplanarity_nOne", xMin=0, xMax=3.14, yMin=1e-5, yMax=10, xLabel="acoplanarity_best", yLabel="Events", logY=True)
        makePlotHiggsDecays("acolinearity_nOne", xMin=0, xMax=3.14, yMin=1e-5, yMax=10, xLabel="acolinearity_best", yLabel="Events", logY=True)

        makePlotHiggsDecays("thrust_magn", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="thrust_magn", yLabel="Events", logY=True)
        
        makePlotHiggsDecays("W1_m_nOne", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="W1_m_nOne", yLabel="Events", logY=True)
        makePlotHiggsDecays("W2_m_nOne", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="W2_m_nOne", yLabel="Events", logY=True)
        makePlotHiggsDecays("W1_p_nOne", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="W1_p_nOne", yLabel="Events", logY=True)
        makePlotHiggsDecays("W2_p_nOne", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="W2_p_nOne", yLabel="Events", logY=True)
        

        makePlot("zqq_m_best_nosel", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("zqq_p_best_nosel", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("zqq_recoil_m_best_nosel", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="Recoil qq (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("best_clustering_idx_nosel", xMin=-1, xMax=4, yMin=1e3, yMax=-1, xLabel="", xLabels=["No pairs", "Inclusive", "Exclusive N=2", "Exclusive N=4", "Exclusive N=6"], yLabel="Events", logY=True, rebin=1)
        makePlot("z_costheta_nosel", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos(#theta_{qq})", yLabel="Events", logY=True, rebin=1)


        makePlot("delta_mWW_nosel", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="delta_mWW_nosel", yLabel="Events", logY=True, rebin=1)
        makePlot("delta_mWW_nOne", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="delta_mWW_nOne", yLabel="Events", logY=True, rebin=1)
        makePlot("W1_m_nOne", xMin=20, xMax=140, yMin=0, yMax=6000e3, xLabel="W mass pair 1 (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=100)
        makePlot("W2_m_nOne", xMin=20, xMax=140, yMin=0, yMax=-1, xLabel="m_{qq} W2 (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=100)
        makePlot("W1_p_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p_{qq} W1 (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("W2_p_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p_{qq} W2 (GeV)", yLabel="Events", logY=True, rebin=1)

        #makePlot("w1_prime", "w1_prime", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W1 (GeV)", yLabel="Events", logY=True, rebin=1)
        #makePlot("w2_prime", "w2_prime", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W2 (GeV)", yLabel="Events", logY=True, rebin=1)


        makePlot("leading_jet_p", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="Leading jet momentum (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("subleading_jet_p", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="Subleading jet momentum (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("leading_jet_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="Leading jet cos(#theta)", yLabel="Events", logY=True, rebin=1)
        makePlot("subleading_jet_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="Subleading jet cos(#theta)", yLabel="Events", logY=True, rebin=1)

        makePlot("njets_inclusive", xMin=0, xMax=14, yMin=1e0, yMax=1e10, xLabel="Inclusive jet multiplicity", yLabel="Events", logY=True)
        makePlot("njets_inclusive_sel", xMin=0, xMax=14, yMin=1e0, yMax=1e10, xLabel="njets_inclusive_sel", yLabel="Events", logY=True)

        makePlot("acoplanarity_nOne", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="Acoplanarity (rad)", yLabel="Events", logY=True, rebin=1)
        makePlot("acolinearity_nOne", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="Acolinearity (rad)", yLabel="Events", logY=True, rebin=1)


        makePlot("zqq_m_best_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("zqq_p_best_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)

        makePlot("z_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos(#theta_{qq})", yLabel="Events", logY=True, rebin=1)


        makePlot("thrust_magn", xMin=0.6, xMax=0.9, yMin=0, yMax=4000e3, xLabel="Thrust magnitude", yLabel="Events", logY=False, rebin=2, sig_scale=100)
        makePlot("thrust_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="MVA score", yLabel="Events", logY=True, rebin=1)

        # final recoil plot
        makePlot("zqq_recoil_m", xMin=100, xMax=150, yMin=1e2, yMax=-1, xLabel="Recoil qq (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("zqq_m", xMin=40, xMax=140, yMin=1e0, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot("zqq_recoil_m", outName="zqq_recoil_m_noLog", xMin=100, xMax=150, yMin=0, yMax=-1, xLabel="Recoil qq (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=50)
        makePlot("zqq_m", outName="zqq_m_noLog", xMin=40, xMax=140, yMin=0, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=50)



        

        makePlot("mva_score", xMin=0, xMax=1, yMin=1e2, yMax=1e8, xLabel="MVA score", yLabel="Events", logY=True, rebin=1)
        makePlot("mva_score", outName="mva_score_noLog", xMin=0, xMax=1, yMin=0, yMax=5e5, xLabel="MVA score", yLabel="Events", logY=False, rebin=1, sig_scale=50)



    if cat == "ee" or cat == "mumu":

        inputDir = f"output/h_zh_leptonic/histmaker/ecm{ecm}/"
        outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zh_leptonic/plots_{cat}_ecm{ecm}_paper"


        z_decays = [cat] # only mumu or ee
        #z_decays = ["qq", "bb", "cc", "ss", "ee", "tautau", "nunu"] # only mumu or ee
        procs = [f"Z{cat}H", "WW", "ZZ", "Zgamma", "Rare"] # first must be signal
        #procs = ["ZnunuH", "WW", "ZZ", "Zgamma", "Rare"] # first must be signal
        #procs = ["ZqqH", "WW", "ZZ", "Zgamma", "Rare"] # first must be signal

        cuts = ["cut0", "cut1", "cut2", "cut3", "cut4"] # , "cut6"
        cut_labels = ["All events", "#geq 1 #mu^{#pm} + ISO", "#geq 2 #mu^{#pm} + OS", "86 < m_{#mu^{+}#mu^{#minus}} < 96", "20 < p_{#mu^{+}#mu^{#minus}} < 70", "100 < m_{rec} < 150", "|cos#theta_{miss}| < 0.98"]

        makeCutFlow(f"{cat}_cutFlow", cuts, cut_labels, sig_scale=100., yMin=1e4, yMax=1e10)
        #makeCutFlowHiggsDecays(f"{cat}_cutFlow", outName="cutFlow", cuts=cuts, cut_labels=cut_labels, yMin=40, yMax=150, z_decays=z_decays, h_decays=h_decays, h_decays_labels=h_decays_labels, h_decays_colors=h_decays_colors)
        #makeCutFlowHiggsDecays(f"{cat}_cutFlow", outName="cutFlow_all", cuts=cuts, cut_labels=cut_labels, yMin=0, yMax=150, z_decays=z_decays_all, h_decays=h_decays, h_decays_labels=h_decays_labels, h_decays_colors=h_decays_colors)
        #makeCutFlowHiggsDecays(f"{cat}_cutFlow", outName="cutFlow_others", cuts=cuts, cut_labels=cut_labels, yMin=0, yMax=150, z_decays=[n for n in z_decays_all if n != cat], h_decays=h_decays, h_decays_labels=h_decays_labels, h_decays_colors=h_decays_colors)

        #makeCutFlowHiggsDecays(f"{cat}_cutFlow", outName="cutFlow_rare_BSM", cuts=cuts, cut_labels=cut_labels, yMin=40, yMax=150, h_decays=h_decays_rare_bsm, h_decays_labels=h_decays_rare_bsm_labels, h_decays_colors=h_decays_rare_bsm_colors)
        #makeCutFlowHiggsDecays(f"{cat}_cutFlow", outName=f"{cat}_cutFlow_agg", cuts=cuts, cut_labels=cut_labels, yMin=60, yMax=120, z_decays=z_decays, h_decays=h_decays_agg, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg)

        #acceptance(hName=f"{cat}_cutFlow", z_decays=[["ee"], ["mumu"], ["tautau"], ["qq", "bb", "cc", "ss"], ["nunu"]])
        #quit()
        # significance
        if False:
            significance(f"{cat}_cosThetaMiss_nOne", 0.95, 1, reverse=True)
            significance(f"{cat}_mva_score", 0, 0.99)
            significance(f"{cat}_mva_score", 0, 0.99, reverse=True)
            significance(f"{cat}_zll_p_nOne", 0, 100)
            significance(f"{cat}_zll_p_nOne", 0, 100, reverse=True)
            significance(f"{cat}_zll_m_nOne", 50, 150)
            significance(f"{cat}_zll_m_nOne", 50, 150, reverse=True)

        makePlotHiggsDecays(f"{cat}_zll_m_nOne", outName="zll_m_nOne", xMin=0, xMax=150, yMin=1e-5, yMax=1e3, xLabel="m(ll) (GeV)", yLabel="Events", logY=True)
        makePlotHiggsDecays(f"{cat}_zll_p_nOne", outName="zll_p_nOne", xMin=0, xMax=150, yMin=1e-5, yMax=1e3, xLabel="p(ll) (GeV)", yLabel="Events", logY=True)
        makePlotHiggsDecays(f"{cat}_zll_recoil_m_final", outName="zll_recoil_m_final", xMin=120, xMax=130, yMin=1e-5, yMax=100, xLabel="Recoil (GeV)", yLabel="Events", logY=True)
        makePlotHiggsDecays(f"{cat}_cosThetaMiss_nOne", outName="cosThetaMiss_nOne", xMin=0.95, xMax=1, yMin=1e-5, yMax=1, xLabel="|cos#theta_{miss}|", yLabel="Events", logY=True, rebin=2)
        makePlotHiggsDecays(f"{cat}_mva_score", outName="mva_score", xMin=0, xMax=1, yMin=1e-4, yMax=1e1, xLabel="MVA score", yLabel="Events", logY=True, rebin=10)
        makePlotHiggsDecays(f"{cat}_mva_score", outName="mva_score_noLog", xMin=0, xMax=1, yMin=0, yMax=0.12, xLabel="MVA score", yLabel="Events", logY=False, rebin=10, h_decays=h_decays_agg, h_decays_labels=h_decays_labels_agg, h_decays_colors=h_decays_colors_agg, legPos="left")
        makePlotHiggsDecays(f"{cat}_acoplanarity", outName="acoplanarity", xMin=0, xMax=5, yMin=1e-5, yMax=1e3, xLabel="m(ll) (GeV)", yLabel="Events", logY=True)
        makePlotHiggsDecays(f"{cat}_acolinearity", outName="acolinearity", xMin=0, xMax=5, yMin=1e-5, yMax=1e3, xLabel="p(ll) (GeV)", yLabel="Events", logY=True)
        #makePlotSignalRatios(f"{cat}_mva_score", outName="mva_score_signalRatio", xMin=0.75, xMax=1, yMin=0.5, yMax=1.5, xLabel="MVA score", yLabel="Signal/nominal", rebin=5)

        makePlot(f"{cat}_zll_m_nOne", "zll_m_nOne", xMin=0, xMax=150, yMin=1e-3, yMax=-1, xLabel="m(ll) (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot(f"{cat}_zll_p_nOne", "zll_p_nOne", xMin=0, xMax=150, yMin=1e-3, yMax=-1, xLabel="p(ll) (GeV)", yLabel="Events", logY=True, rebin=1)
        makePlot(f"{cat}_zll_recoil_nOne", "zll_recoil_nOne", xMin=100, xMax=150, yMin=0, yMax=-1, xLabel="Recoil (GeV)", yLabel="Events", logY=False, rebin=1)
        makePlot(f"{cat}_zll_recoil_m_final", "zll_recoil_m_final", xMin=120, xMax=130, yMin=0, yMax=-1, xLabel="Recoil (GeV)", yLabel="Events", logY=False, rebin=1)
        makePlot(f"{cat}_mva_score", "mva_score", xMin=0, xMax=1, yMin=1e1 if ecm==240 else 1e-1, yMax=1e5, xLabel="MVA score", yLabel="Events", logY=True, rebin=4)
        makePlot(f"{cat}_mva_score", "mva_score_noLog", xMin=0, xMax=1, yMin=0, yMax=5e3, xLabel="MVA score", yLabel="Events", logY=False, rebin=4)
        makePlot(f"{cat}_zll_recoil_m_mva_low", "zll_recoil_m_mva_low", xMin=100, xMax=150, yMin=0, yMax=25e3, xLabel="Recoil mass (GeV)", yLabel="Events", logY=False, rebin=4)
        makePlot(f"{cat}_zll_recoil_m_mva_high", "zll_recoil_m_mva_high", xMin=120, xMax=150, yMin=0, yMax=3000, xLabel="Recoil mass (GeV)", yLabel="Events", logY=False, rebin=1)
        
        makePlot(f"{cat}_zll_recoil", "zll_recoil", xMin=120, xMax=140, yMin=0, yMax=700, xLabel="Recoil mass (GeV)", yLabel="Events", logY=False, rebin=10, sig_scale=5)
        # 16e3  
        
        ## MVA plots  
        makePlot(f"{cat}_acolinearity", "zll_acolinearity", xMin=0, xMax=1.2, yMin=0, yMax=25e3, xLabel="Acolinearity (rad)", yLabel="Events", logY=False, rebin=1, sig_scale=10)
        makePlot(f"{cat}_zll_p", "zll_p", xMin=20, xMax=70, yMin=0, yMax=16e3, xLabel="Lepton pair momentum (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=10)
        
        
        '''
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
        '''