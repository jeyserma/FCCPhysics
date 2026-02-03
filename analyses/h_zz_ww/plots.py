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

def makeCutFlow(hName="cutFlow", cuts=[], labels=[], sig_scales=[], yMin=1e6, yMax=1e10):

    leg = ROOT.TLegend(.5, 0.99-(len(procs))*0.06, .99, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists_yields = []
    significances = []





    sigs = []
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,proc in enumerate(procs):
        isSig = False
        if "SIGNAL:" in proc:
            isSig = True
            proc = proc.replace("SIGNAL:", "")
        h = getHist(hName, procs_cfg[proc])
        hists_yields.append(h)
        if isSig:
            h.SetLineColor(procs_colors[proc])
            h.SetLineWidth(3)
            h.SetLineStyle(1)
            h.Scale(sig_scales[i])
            if sig_scales[i] != 1:
                leg.AddEntry(h, f"{procs_labels[proc]} (#times {int(sig_scales[i])})", "L")
            else:
                leg.AddEntry(h, procs_labels[proc], "L")

            sigs.append(h)
        else:
            if h_bkg_tot == None: h_bkg_tot = h.Clone("h_bkg_tot")
            else: h_bkg_tot.Add(h)
            h.SetFillColor(procs_colors[proc])
            h.SetLineColor(ROOT.kBlack)
            h.SetLineWidth(1)
            h.SetLineStyle(1)
            leg.AddEntry(h, procs_labels[proc], "F")
            st.Add(h)
            

    h_bkg_tot.SetLineColor(ROOT.kBlack)
    h_bkg_tot.SetLineWidth(2)

    for i,cut in enumerate(cuts):
        nsig = sigs[0].GetBinContent(i+1) / sig_scales[0] ## undo scaling
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
    for h_sig in sigs:
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


def makeCutFlowHiggsDecays(hName="cutFlow", outName="", cuts=[], cut_labels=[], yMin=0, yMax=150, z_decays=[], h_decays=[], h_decays_labels=[], h_decays_colors=[]):

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
    hist_tot = None
    eff_final, eff_final_err = [], []
    for i,sig in enumerate(sigs):
        h_decay = h_decays[i]
        h_sig = getHist(hName, sig)
        h_sig.Scale(100./h_sig.GetBinContent(1))

        h_sig.SetLineColor(h_decays_colors[h_decay])
        h_sig.SetLineWidth(2)
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
    dummy.Draw("HIST")

    txt = ROOT.TLatex()
    txt.SetTextSize(0.04)
    txt.SetTextColor(1)
    txt.SetTextFont(42)
    txt.SetNDC()
    txt.DrawLatex(0.2, 0.2, f"Avg. eff.: {eff_avg:.2f}#pm{eff_avg_err:.2f} %")
    txt.DrawLatex(0.2, 0.15, f"Min/max: {eff_min:.2f}/{eff_max:.2f}")
    txt.Draw("SAME")

    for hist in hists:
        hist.Draw("SAME HIST")
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
    if args.cat == "mumu":
        xMin, xMax = 68, 74
        if len(z_decays) > 1:
            xMin, xMax = 0, 5
        xMin, xMax = 0, 100
    if args.cat == "ee":
        xMin, xMax = 58, 66
        if len(z_decays) > 1:
            xMin, xMax = 0, 5
    if args.cat == "qq":
        xMin, xMax = 55, 85
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

    maxx = len(sigs)+1
    line = ROOT.TLine(eff_avg, 0, eff_avg, maxx)
    line.SetLineColor(ROOT.kGray)
    line.SetLineWidth(2)
    #line.Draw("SAME")


    shade = ROOT.TGraph()
    shade.SetPoint(0, eff_avg-eff_avg_err, 0)
    shade.SetPoint(1, eff_avg+eff_avg_err, 0)
    shade.SetPoint(2, eff_avg+eff_avg_err, maxx)
    shade.SetPoint(3, eff_avg-eff_avg_err, maxx)
    shade.SetPoint(4, eff_avg-eff_avg_err, 0)
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
    txt.DrawLatex(0.2, 0.2, f"Avg. eff.: {eff_avg:.2f}#pm{eff_avg_err:.2f} %")
    txt.DrawLatex(0.2, 0.15, f"Min/max: {eff_min:.2f}/{eff_max:.2f}")
    txt.Draw("SAME")

    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}_finalSelection.png")
    canvas.SaveAs(f"{outDir}/higgsDecays/{outName}_finalSelection.pdf")


def makePlotHiggsDecays(hName, outName="", xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="", yLabel="Events", logX=False, logY=True, rebin=1, xLabels=[]):

    if outName == "":
        outName = hName

    sigs = [[f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays] for y in h_decays]
    leg = ROOT.TLegend(.2, .925-(len(sigs)/4+1)*0.07, .95, .925)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.25)
    leg.SetNColumns(4)

    hists = []
    for i,sig in enumerate(sigs):
        h_decay = h_decays[i]
        h_sig = getHist(hName, sig)
        h_sig.Rebin(rebin)
        if h_sig.Integral() > 0:
            h_sig.Scale(1./h_sig.Integral())

        h_sig.SetLineColor(h_decays_colors[h_decay])
        h_sig.SetLineWidth(2)
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

    sigs = [[f'wzp6_ee_{x}H_H{y}_ecm240' for x in z_decays] for y in h_decays]
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

    if outName == "":
        outName = hName

    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(.55, 0.99-(len(procs))*0.06, .99, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    '''
    h_sig = getHist(hName, procs_cfg[procs[0]])
    h_sig.SetLineColor(procs_colors[procs[0]])
    h_sig.SetLineWidth(3)
    h_sig.SetLineStyle(1)
    h_sig.Scale(sig_scale)
    if sig_scale != 1:
        leg.AddEntry(h_sig, f"{procs_labels[procs[0]]} (#times {int(sig_scale)})", "L")
    else:
        leg.AddEntry(h_sig, procs_labels[procs[0]], "L")
    '''
    sigs = []
    st = ROOT.THStack()
    st.SetName("stack")
    h_bkg_tot = None
    for i,proc in enumerate(procs):
        isSig = False
        if "SIGNAL:" in proc:
            isSig = True
            proc = proc.replace("SIGNAL:", "")
        h = getHist(hName, procs_cfg[proc])

        if isSig:
            h.SetLineColor(procs_colors[proc])
            h.SetLineWidth(3)
            h.SetLineStyle(1)
            leg.AddEntry(h, procs_labels[proc], "L")
            sigs.append(h)
        else:
            if h_bkg_tot == None: h_bkg_tot = h.Clone("h_bkg_tot")
            else: h_bkg_tot.Add(h)
            h.SetFillColor(procs_colors[proc])
            h.SetLineColor(ROOT.kBlack)
            h.SetLineWidth(1)
            h.SetLineStyle(1)
            leg.AddEntry(h, procs_labels[proc], "F")
            st.Add(h)



    if yMax < 0:
        if logY:
            yMax = math.ceil(max([h_bkg_tot.GetMaximum(), max([h.GetMaximum() for h in sigs])])*10000)/10.
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
    for h_sig in sigs:
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



def compareHists(procHists, outName, xMin=0, xMax=100, yMin=1, yMax=1e5, xLabel="xlabel", yLabel="Events", logX=False, logY=True, rebin=1, legPos=[0.3, 0.75, 0.9, 0.9], norm=True):

    st = ROOT.THStack()
    st.SetName("stack")

    leg = ROOT.TLegend(.55, 0.99-(len(procHists))*0.1, .99, .90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    hists = []
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue]
    for i,procHist in enumerate(procHists):
        proc, hName = procHist.split(':')[0], procHist.split(':')[1]
        h = getHist(hName, procs_cfg[proc])
        if norm:
            h.Scale(1./h.Integral())
        h.SetLineColor(colors[i])
        h.SetLineWidth(3)
        h.SetLineStyle(1)
        leg.AddEntry(h, f"{proc} {hName}", "L")
        hists.append(h)



    if yMax < 0:
        if logY:
            yMax = math.ceil(max([h.GetMaximum() for h in hists])*10000)/10.
        else:
            yMax = 1.3*max([h.GetMaximum() for h in hists])

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
    dummy = plotter.dummy(1)
    dummy.Draw("HIST")

    for h in hists:
        h.Draw("HIST SAME")
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


def significance(hName, xMin=-10000, xMax=10000, reverse=False, norm=False):

    h_sig = getHist(hName, procs_cfg[procs[0]])
    
  

    bkgs_procs = []
    for i,bkg in enumerate(procs[1:]):
        bkgs_procs.extend(procs_cfg[bkg])

    h_bkg = getHist(hName, bkgs_procs)
    if norm:
        h_sig.Scale(1./h_sig.Integral())
        h_bkg.Scale(1./h_bkg.Integral())
    
    sig_tot = h_sig.Integral()
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
    lumi = "10.8" if ecm == 240 else "3"



    procs_cfg_240 = {
        "ZHWW"        : ['wzp6_ee_ccH_HWW_ecm240', 'wzp6_ee_bbH_HWW_ecm240', 'wzp6_ee_qqH_HWW_ecm240', 'wzp6_ee_ssH_HWW_ecm240', 'wzp6_ee_nunuH_HWW_ecm240', 'wzp6_ee_eeH_HWW_ecm240', 'wzp6_ee_tautauH_HWW_ecm240', 'wzp6_ee_mumuH_HWW_ecm240'],
        "ZqqHWW"        : ['wzp6_ee_ccH_HWW_ecm240', 'wzp6_ee_bbH_HWW_ecm240', 'wzp6_ee_qqH_HWW_ecm240', 'wzp6_ee_ssH_HWW_ecm240'],
        "ZqqHZZ"        : ['wzp6_ee_ccH_HZZ_ecm240', 'wzp6_ee_bbH_HZZ_ecm240', 'wzp6_ee_qqH_HZZ_ecm240', 'wzp6_ee_ssH_HZZ_ecm240'],
        "ZHZZ"        : ['wzp6_ee_nunuH_HZZ_ecm240', 'wzp6_ee_eeH_HZZ_ecm240', 'wzp6_ee_tautauH_HZZ_ecm240', 'wzp6_ee_mumuH_HZZ_ecm240', 'wzp6_ee_ccH_HZZ_ecm240', 'wzp6_ee_bbH_HZZ_ecm240', 'wzp6_ee_qqH_HZZ_ecm240', 'wzp6_ee_ssH_HZZ_ecm240'],
        "WW"        : ['p8_ee_WW_ecm240'],
        "ZZ"        : ['p8_ee_ZZ_ecm240'],
        "Zgamma"    : ['wz3p6_ee_uu_ecm240', 'wz3p6_ee_dd_ecm240', 'wz3p6_ee_cc_ecm240', 'wz3p6_ee_ss_ecm240', 'wz3p6_ee_bb_ecm240'],
    }

    
    if ecm == 240:
        procs_cfg = procs_cfg_240
    if ecm == 365:
        procs_cfg = procs_cfg_365


    procs_labels = {
        "ZHWW"      : "ZH(WW)",
        "ZqqHWW"      : "Z(qq)H(W(qq)W(qq))",
        "ZqqHZZ"      : "Z(qq)H(Z(qq)Z(qq))",
        "ZHZZ"      : "ZHZZ",
        "WW"        : "WW",
        "ZZ"        : "ZZ",
        "Zgamma"    : "Z/#gamma^{*} #rightarrow f#bar{f}+#gamma(#gamma)",
    }

    # colors from https://github.com/mpetroff/accessible-color-cycles
    procs_colors = {
        "ZHWW"      : ROOT.TColor.GetColor("#e42536"),
        "ZqqHWW"      : ROOT.TColor.GetColor("#e42536"),
        "ZqqHZZ"      : ROOT.TColor.GetColor("#9c9ca1"),
        "WW"        : ROOT.TColor.GetColor("#f89c20"),
        "ZZ"        : ROOT.TColor.GetColor("#5790fc"),
        "Zgamma"    : ROOT.TColor.GetColor("#964a8b"),
        "ZHZZ"      : ROOT.TColor.GetColor("#9c9ca1") 
    }



    inputDir = f"output/h_zz_ww/histmaker/ecm{ecm}/stage2/"
    inputDir = f"output/h_zz_ww/histmaker/ecm{ecm}/preSelPlots"
    outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zz_ww/plots_ecm{ecm}/"

    procs = ["SIGNAL:ZqqHWW", "SIGNAL:ZqqHZZ", "WW", "ZZ", "Zgamma"] # first must be signal
    procs = [ "SIGNAL:ZqqHZZ","SIGNAL:ZqqHWW", "WW", "ZZ", "Zgamma"] # first must be signal
    #procs = ["ZqqHZZ", "ZHWW", "WW", "ZZ", "Zgamma"] # first must be signal


    #procs = ["ZqqHWW", "ZqqHZZ"]
    #significance("mva_score", 0, 1, norm=True)
    #significance("mva_score", 0, 1, reverse=True, norm=True)
    
    #procs = ["ZqqHWW", "WW", "ZZ", "Zgamma"]
    #significance("mva_score", 0, 1)
    
    #procs = ["ZqqHWW", "WW", "ZZ", "Zgamma"]
    #significance("mva_score_ww", 0, 1)
    #procs = ["ZqqHZZ", "WW", "ZZ", "Zgamma"]
    #significance("mva_score_zz", 0, 1)
    #quit()

    if True:
        inputDir = f"output/h_zz_ww/histmaker/ecm{ecm}/jet_pairing/"
        outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zz_ww/plots_ecm{ecm}/jet_pairing"
        compareHists(['ZqqHWW:W1_WW_m', 'ZqqHZZ:Z1_ZZ_m'], "V1_m", xMax=125, yMin=0, yMax=-1, xLabel="Invariant mass Z/W 1 (GeV)", yLabel="Events", logY=False, rebin=1)
        compareHists(['ZqqHWW:W2_WW_m', 'ZqqHZZ:Z2_ZZ_m'], "V2_m", xMax=125, yMin=0, yMax=-1, xLabel="Invariant mass Z/W 2 (GeV)", yLabel="Events", logY=False, rebin=1)
        compareHists(['ZqqHWW:Z_WW_m', 'ZqqHZZ:Z_ZZ_m'], "Z_m", xMax=125, yMin=0, yMax=-1, xLabel="Invariant mass associated Z (GeV)", yLabel="Events", logY=False, rebin=1)
        compareHists(['ZqqHWW:H_WW_m', 'ZqqHZZ:H_ZZ_m'], "H_m", xMin=100, xMax=150, yMin=0, yMax=-1, xLabel="Invariant mass Higgs (GeV)", yLabel="Events", logY=False, rebin=1)
        compareHists(['ZqqHWW:chi2_WW', 'ZqqHZZ:chi2_ZZ'], "chi2", xMin=0, xMax=1000, yMin=0, yMax=-1, xLabel="chi2 min", yLabel="Events", logY=False, rebin=1)
        quit()
    
    ### PRESELECTION PLOTS
    inputDir = f"output/h_zz_ww/histmaker/ecm{ecm}/preSelPlots"
    
    makePlot("missingEnergy_nOne", xMin=0, xMax=120, yMin=1e-1, yMax=-1, xLabel="Missing momentum (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("visibleEnergy_nOne", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Visible energy (GeV)", yLabel="Events", logY=True, rebin=1)
    ##makePlot("visibleMass_nOne", xMin=120, xMax=250, yMin=1e-1, yMax=-1, xLabel="Visible mass (GeV)", yLabel="Events", logY=True, rebin=1)

    makePlot("sq_dmerge_01", xMin=0, xMax=500, yMin=1e-1, yMax=-1, xLabel="#sqrt{d_{01}}", yLabel="Events", logY=True, rebin=1)
    makePlot("sq_dmerge_12", xMin=0, xMax=300, yMin=1e-1, yMax=-1, xLabel="#sqrt{d_{12}}", yLabel="Events", logY=True, rebin=1)
    makePlot("sq_dmerge_23", xMin=0, xMax=150, yMin=1e-1, yMax=-1, xLabel="#sqrt{d_{23}}", yLabel="Events", logY=True, rebin=1)
    makePlot("sq_dmerge_34", xMin=0, xMax=100, yMin=1e-1, yMax=-1, xLabel="#sqrt{d_{34}}", yLabel="Events", logY=True, rebin=1)
    makePlot("sq_dmerge_45", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="#sqrt{d_{45}}", yLabel="Events", logY=True, rebin=1)
    makePlot("sq_dmerge_56", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="#sqrt{d_{56}}", yLabel="Events", logY=True, rebin=1)
    makePlot("sq_dmerge_67", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="#sqrt{d_{67}}", yLabel="Events", logY=True, rebin=1)


    cuts = ["cut0", "cut1", "cut2", "cut3", "cut4", "cut5", "cut6", "cut7", "cut8", "cut9"]
    cut_labels = ["All events", "Signal hadronic", "Veto leptons", "p_{miss} < 20", "E_{vis} > 115", "6 jets", "#sqrt{d_{23}} > 50", "#sqrt{d_{34}} > 25", "#sqrt{d_{45}} > 15", "#sqrt{d_{56}} > 10"]



    makeCutFlow("cutFlow", cuts, cut_labels, sig_scales=[100, 1000])

    '''
    



    # significance
    if True:
        significance("visibleEnergy_nOne", 100, 250)
        significance("visibleEnergy_nOne", 100, 250, reverse=True)
        significance("missingEnergy_nOne", 0, 100, reverse=True)
        
        significance("sq_dmerge_23", 0, 100)
        significance("sq_dmerge_34", 0, 100)
    
        


    makePlot("njets_nOne", xMin=0, xMax=15, yMin=1e-1, yMax=-1, xLabel="njets_nOne", yLabel="Events", logY=True, rebin=1)

    #makePlot("W1_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="W1_m", yLabel="Events", logY=True, rebin=1)
    #makePlot("W2_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="W2_m", yLabel="Events", logY=True, rebin=1)
    #makePlot("Z_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Z_m", yLabel="Events", logY=True, rebin=1)
    #makePlot("H_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="H_m", yLabel="Events", logY=True, rebin=1)
    '''

    makePlot("mva_score", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="mva_score", yLabel="MVA score", logY=True, rebin=1)


    # pairing
    outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zz_ww/plots_ecm{ecm}/pairing_WW"
    makePlot("W1_WW_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="W1_WW_m", yLabel="Events", logY=True, rebin=1)
    makePlot("W2_WW_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="W2_WW_m", yLabel="Events", logY=True, rebin=1)
    makePlot("Z_WW_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Z_WW_m", yLabel="Events", logY=True, rebin=1)
    makePlot("Z_WW_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Z_WW_p", yLabel="Events", logY=True, rebin=1)
    makePlot("H_WW_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="H_WW_m", yLabel="Events", logY=True, rebin=1)
    makePlot("H_WW_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="H_WW_p", yLabel="Events", logY=True, rebin=1)

    outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zz_ww/plots_ecm{ecm}/pairing_ZZ"
    makePlot("Z1_ZZ_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Z1_ZZ_m", yLabel="Events", logY=True, rebin=1)
    makePlot("Z2_ZZ_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Z2_ZZ_m", yLabel="Events", logY=True, rebin=1)
    makePlot("Z_ZZ_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Z_ZZ_m", yLabel="Events", logY=True, rebin=1)
    makePlot("Z_ZZ_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="Z_ZZ_p", yLabel="Events", logY=True, rebin=1)
    makePlot("H_ZZ_m", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="H_ZZ_m", yLabel="Events", logY=True, rebin=1)
    makePlot("H_ZZ_p", xMin=0, xMax=250, yMin=1e-1, yMax=-1, xLabel="H_ZZ_p", yLabel="Events", logY=True, rebin=1)


    '''
    outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zz_ww/plots_ecm{ecm}/jet_tagging"
    compareHists(['ZqqHWW:sum_score_B_V1', 'ZqqHZZ:sum_score_B_V1'], "sum_score_B_V1", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_B_V1", yLabel="Events", logY=False, rebin=1)
    compareHists(['ZqqHWW:sum_score_B_V2', 'ZqqHZZ:sum_score_B_V2'], "sum_score_B_V2", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_B_V2", yLabel="Events", logY=False, rebin=1)
    
    compareHists(['ZqqHWW:sum_score_C_V1', 'ZqqHZZ:sum_score_C_V1'], "sum_score_C_V1", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_C_V1", yLabel="Events", logY=False, rebin=1)
    compareHists(['ZqqHWW:sum_score_C_V2', 'ZqqHZZ:sum_score_C_V2'], "sum_score_C_V2", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_C_V2", yLabel="Events", logY=False, rebin=1)
    
    compareHists(['ZqqHWW:sum_score_S_V1', 'ZqqHZZ:sum_score_S_V1'], "sum_score_S_V1", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_S_V1", yLabel="Events", logY=False, rebin=1)
    compareHists(['ZqqHWW:sum_score_S_V2', 'ZqqHZZ:sum_score_S_V2'], "sum_score_S_V2", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_S_V2", yLabel="Events", logY=False, rebin=1)
    
    compareHists(['ZqqHWW:sum_score_U_V1', 'ZqqHZZ:sum_score_U_V1'], "sum_score_U_V1", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_U_V1", yLabel="Events", logY=False, rebin=1)
    compareHists(['ZqqHWW:sum_score_U_V2', 'ZqqHZZ:sum_score_U_V2'], "sum_score_U_V2", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_U_V2", yLabel="Events", logY=False, rebin=1)
    
    compareHists(['ZqqHWW:sum_score_D_V1', 'ZqqHZZ:sum_score_D_V1'], "sum_score_D_V1", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_D_V1", yLabel="Events", logY=False, rebin=1)
    compareHists(['ZqqHWW:sum_score_D_V2', 'ZqqHZZ:sum_score_D_V2'], "sum_score_D_V2", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_D_V2", yLabel="Events", logY=False, rebin=1)

    compareHists(['ZqqHWW:sum_score_TAU_V1', 'ZqqHZZ:sum_score_TAU_V1'], "sum_score_TAU_V1", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_TAU_V1", yLabel="Events", logY=False, rebin=1)
    compareHists(['ZqqHWW:sum_score_TAU_V2', 'ZqqHZZ:sum_score_TAU_V2'], "sum_score_TAU_V2", xMin=0, xMax=2, yMin=0, yMax=-1, xLabel="sum_score_TAU_V2", yLabel="Events", logY=False, rebin=1)
    '''

    outDir = f"/home/submit/jaeyserm/public_html/fccee/h_zz_ww/plots_ecm{ecm}/dmerge"
    makePlot("dmerge_01", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_01", yLabel="Events", logY=True, rebin=1)
    makePlot("dmerge_12", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_12", yLabel="Events", logY=True, rebin=1)
    makePlot("dmerge_23", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_23", yLabel="Events", logY=True, rebin=1)
    makePlot("dmerge_34", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_34", yLabel="Events", logY=True, rebin=1)
    makePlot("dmerge_45", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_45", yLabel="Events", logY=True, rebin=1)
    makePlot("dmerge_56", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_56", yLabel="Events", logY=True, rebin=1)
    makePlot("dmerge_67", xMin=0, xMax=50000, yMin=1e-1, yMax=-1, xLabel="dmerge_67", yLabel="Events", logY=True, rebin=1)






    quit()


    makePlotHiggsDecays("best_clustering_idx_nosel", xMin=-1, xMax=4, yMin=0, yMax=1.5, xLabels=["No pairs", "Inclusive", "Exclusive N=2", "Exclusive N=4", "Exclusive N=6"], yLabel="Events", logY=False)

    makePlotHiggsDecays("zqq_m_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_m_best_nosel", yLabel="Events", logY=True)
    makePlotHiggsDecays("zqq_p_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_p_best_nosel", yLabel="Events", logY=True)
    makePlotHiggsDecays("zqq_recoil_m_best_nosel", xMin=0, xMax=200, yMin=1e-5, yMax=10, xLabel="zqq_recoil_m_best_nosel", yLabel="Events", logY=True)
    makePlotHiggsDecays("z_costheta_nosel", xMin=0, xMax=1, yMin=1e-5, yMax=10, xLabel="cos(qq) (GeV)", yLabel="Events", logY=True)
    makePlotHiggsDecays("cosThetaMiss_nOne", xMin=0.95, xMax=1, yMin=1e-5, yMax=10, xLabel="cosThetaMiss_nOne", yLabel="Events", logY=True)

    makePlotHiggsDecays("njets_inclusive", xMin=0, xMax=15, yMin=1e-5, yMax=1e3, xLabel="njets_inclusive", yLabel="Events", logY=True)
    makePlotHiggsDecays("njets_inclusive_sel", xMin=0, xMax=15, yMin=1e-5, yMax=1e3, xLabel="njets_inclusive_sel", yLabel="Events", logY=True)
    makePlotHiggsDecays("delta_mWW_nOne", xMin=0, xMax=50, yMin=1e-5, yMax=1e3, xLabel="delta_mWW_nOne", yLabel="Events", logY=True)
    makePlotHiggsDecays("mva_score", xMin=0, xMax=1, yMin=1e-5, yMax=1e1, xLabel="MVA score", yLabel="Events", logY=True, rebin=2)

    
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
    makePlot("best_clustering_idx_nosel", xMin=-1, xMax=4, yMin=1e2, yMax=-1, xLabel="", xLabels=["No pairs", "Inclusive", "Exclusive N=2", "Exclusive N=4", "Exclusive N=6"], yLabel="Events", logY=True, rebin=1)
    makePlot("z_costheta_nosel", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos(#theta_{qq})", yLabel="Events", logY=True, rebin=1)


    makePlot("delta_mWW_nosel", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="delta_mWW_nosel", yLabel="Events", logY=True, rebin=1)
    makePlot("delta_mWW_nOne", xMin=0, xMax=50, yMin=1e-1, yMax=-1, xLabel="delta_mWW_nOne", yLabel="Events", logY=True, rebin=1)
    makePlot("W1_m_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m_{qq} W1 (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("W2_m_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m_{qq} W2 (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("W1_p_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p_{qq} W1 (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("W2_p_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p_{qq} W2 (GeV)", yLabel="Events", logY=True, rebin=1)

    #makePlot("w1_prime", "w1_prime", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W1 (GeV)", yLabel="Events", logY=True, rebin=1)
    #makePlot("w2_prime", "w2_prime", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m(qq) W2 (GeV)", yLabel="Events", logY=True, rebin=1)


    makePlot("leading_jet_p", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="Leading jet momentum (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("subleading_jet_p", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="Subleading jet momentum (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("leading_jet_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="Leading jet cos(#theta)", yLabel="Events", logY=True, rebin=1)
    makePlot("subleading_jet_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="Subleading jet cos(#theta)", yLabel="Events", logY=True, rebin=1)

    makePlot("njets_inclusive", xMin=0, xMax=15, yMin=1e-1, yMax=-1, xLabel="Inclusive jet multiplicity", yLabel="Events", logY=True)
    makePlot("njets_inclusive_sel", xMin=0, xMax=15, yMin=1e-1, yMax=-1, xLabel="njets_inclusive_sel", yLabel="Events", logY=True)

    makePlot("acoplanarity_nOne", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="Acoplanarity (rad)", yLabel="Events", logY=True, rebin=1)
    makePlot("acolinearity_nOne", xMin=0, xMax=4, yMin=1e-1, yMax=-1, xLabel="Acolinearity (rad)", yLabel="Events", logY=True, rebin=1)


    makePlot("zqq_m_best_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_p_best_nOne", xMin=0, xMax=200, yMin=1e-1, yMax=-1, xLabel="p_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)

    makePlot("z_costheta_nOne", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="cos(#theta_{qq})", yLabel="Events", logY=True, rebin=1)


    makePlot("thrust_magn", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="Thrust magnitude", yLabel="Events", logY=True, rebin=1)
    makePlot("thrust_costheta", xMin=0, xMax=1, yMin=1e-1, yMax=-1, xLabel="MVA score", yLabel="Events", logY=True, rebin=1)

    # final recoil plot
    makePlot("zqq_recoil_m", xMin=100, xMax=150, yMin=1e2, yMax=-1, xLabel="Recoil qq (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_m", xMin=40, xMax=140, yMin=1e0, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=True, rebin=1)
    makePlot("zqq_recoil_m", outName="zqq_recoil_m_noLog", xMin=100, xMax=150, yMin=0, yMax=-1, xLabel="Recoil qq (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=50)
    makePlot("zqq_m", outName="zqq_m_noLog", xMin=40, xMax=140, yMin=0, yMax=-1, xLabel="m_{qq} (GeV)", yLabel="Events", logY=False, rebin=1, sig_scale=50)


    makePlot("mva_score", xMin=0, xMax=1, yMin=1e1, yMax=-1, xLabel="MVA score", yLabel="Events", logY=True, rebin=1)
    makePlot("mva_score", outName="mva_score_noLog", xMin=0, xMax=1, yMin=0, yMax=5e5, xLabel="MVA score", yLabel="Events", logY=False, rebin=1, sig_scale=50)

