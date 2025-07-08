
import sys,argparse,copy
import os
import ROOT
import math
from collections import OrderedDict

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)



def compare(inputFiles, histNames, legends, outName="out", yMin=0, yMax=-1):
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kOrange, ROOT.kBlue, ROOT.kGreen, ROOT.kGray+1, ROOT.kMagenta]

    hists = []
    for i,fin in enumerate(inputFiles):
        fIn = ROOT.TFile(fin)
        h = copy.deepcopy(fIn.Get(histNames[i]))
        fIn.Close()
        hists.append(h)

    c = ROOT.TCanvas("", "", 800, 800)
    c.SetTopMargin(0.055)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.11)

    hists[0].GetXaxis().SetTitleFont(43)
    hists[0].GetXaxis().SetTitleSize(32)
    hists[0].GetXaxis().SetLabelFont(43)
    hists[0].GetXaxis().SetLabelSize(28)

    hists[0].GetYaxis().SetTitleFont(43)
    hists[0].GetYaxis().SetTitleSize(32)
    hists[0].GetYaxis().SetLabelFont(43)
    hists[0].GetYaxis().SetLabelSize(28)


    leg = ROOT.TLegend(0.2, 0.9-0.04*(len(inputFiles)+1), 0.55, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)

    for i,h in enumerate(hists):
        h.SetLineColor(colors[i])
        h.SetLineWidth(2)
        h.SetLineStyle(1)
        h.SetFillStyle(0)
        leg.AddEntry(h, legends[i], "L")

        if i == 0:
            h.GetYaxis().SetRangeUser(yMin, yMax if yMax > 0 else h.GetMaximum()*1.05)
            h.GetXaxis().SetTitle("cos(#theta)")
            h.GetYaxis().SetTitle("Material budget x/X_{0} (%)")
            h.Draw("HIST")
        else:
            h.Draw("SAME HIST")

    leg.Draw()
    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()
    c.SaveAs(f"{outDir}/{outName}.png")
    c.SaveAs(f"{outDir}/{outName}.pdf")


if __name__ == '__main__':

    outDir = "/home/submit/jaeyserm/public_html/fccee/material_budget/"

    inputFiles = []
    histNames = []
    legends = []


    inputFiles.append("/home/submit/jaeyserm/public_html/fccee/material_budget//x0_grouped.root")
    inputFiles.append("/home/submit/jaeyserm/public_html/fccee/material_budget//delphes_x0.root")
    legends.append("FullSim")
    legends.append("FastSim (Delphes)")

    # plot beampipe
    histNames = ["Beampipe", "Beampipe"]
    compare(inputFiles, histNames, legends, outName="beampipe", yMin=0, yMax=5)

    histNames = ["Vertex", "Vertex"]
    compare(inputFiles, histNames, legends, outName="vertex", yMin=0, yMax=10)

