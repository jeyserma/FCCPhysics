from __future__ import print_function
import argparse

import sys, os
import json

sys.path.append(os.path.expandvars("$FCCSW") + "/Examples/scripts")
from plotstyle import FCCStyle

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

def main():

    parser = argparse.ArgumentParser(description="Material Plotter")
    parser.add_argument("--fname", "-f", dest="fname", type=str, help="name of file to read")
    parser.add_argument("--outDir", "-o", dest="outDir", type=str, default="/home/submit/jaeyserm/public_html/fccee/material_budget/")
    parser.add_argument("--param", "-p", dest="param", type=str, help="Parameter to plot (x0, lambda, depth)", choices=["x0", "lambda", "depth"], default="x0")
    parser.add_argument("--ignoreMats", "-i", dest="ignoreMats", nargs="+", default=[], help="List of materials that should be ignored")
    parser.add_argument("-y", "--ymax", type=int, default=-1, help="Maximum y-axis")
    args = parser.parse_args()

    angleMin = 0
    angleMax = 0.99
    angleDef = "cosTheta"
    angleBinning = 0.01

    material_grouping = {}
    material_grouping['Beampipe'] = ["Gold", "AlBeMet162", "LiquidNDecane", "beam"]
    material_grouping['Vertex'] = ["Aluminum", "CarbonFleece", "CarbonFiberVtx", "CarbonFiberVtx66D", "RohacellVtx66D", "RohacellVtx", "KaptonVtx", "GlueEcobond45", "Water", "Silicon", "PCB", "PEEK"]
    colors_grouping = {'Beampipe': ROOT.kRed , 'Vertex': ROOT.kOrange}

    ############################################################################
    angleBinsN = (int)((angleMax - angleMin) / angleBinning) # number of bins

    f = ROOT.TFile.Open(args.fname, "read")
    tree = f.Get("materials")
    histDict = {}

    material_grouping_hists = {}
    material_hists = {}
    for group in material_grouping.keys():
        material_grouping_hists[group] = ROOT.TH1F(group, "", angleBinsN, angleMin, angleMax)

    def material2group(material, material_grouping):
        for mg in material_grouping.keys():
            if material in material_grouping[mg]:
                return mg
        print(f"Group for material {material} not found")
        return None


    for angleBinning, entry in enumerate(tree):
        nMat = entry.nMaterials
        for i in range(nMat):
            material = entry.material.at(i)
            group = material2group(material, material_grouping)
            print(angleBinning, material, group)

            if material in args.ignoreMats:
                continue

            if material not in material_hists.keys():
                material_hists[material] = ROOT.TH1F(material, "", angleBinsN, angleMin, angleMax)

            hs = material_hists[material]
            val = hs.GetBinContent(angleBinning+1) # original value
            if args.param == "x0":
                val += entry.nX0.at(i) * 100.0
            elif args.param == "lambda":
                val += entry.nLambda.at(i)
            elif args.param == "depth":
                val += entry.matDepth.at(i)
            hs.SetBinContent(angleBinning+1, val)

            if group != None:
                material_grouping_hists[group].SetBinContent(angleBinning+1, material_grouping_hists[group].GetBinContent(angleBinning+1)+val)


    axis_titles = ["Material budget x/X_{0} [%] ", "Number of #lambda", "Material depth [cm]"]

    if args.param == "x0":
        ytitle = "Material budget x/X_{0} (%)"
    elif args.param == "lambda":
        ytitle = "Number of #lambda"
    elif args.param == "depth":
        ytitle = "Material depth (cm)"

    if angleDef == "eta":
        xtitle = "#eta"
    elif angleDef == "theta":
        xtitle = "#theta"
    elif angleDef == "thetaRad":
        xtitle = "#theta [rad]"
    elif angleDef == "cosTheta":
        xtitle = "cos(#theta)"

    dummy = ROOT.TH1D("dummy", "", angleBinsN, angleMin, angleMax)
    dummy.GetXaxis().SetTitle(xtitle)
    dummy.GetYaxis().SetTitle(ytitle)

    dummy.GetXaxis().SetTitleFont(43)
    dummy.GetXaxis().SetTitleSize(32)
    dummy.GetXaxis().SetLabelFont(43)
    dummy.GetXaxis().SetLabelSize(28)

    dummy.GetYaxis().SetTitleFont(43)
    dummy.GetYaxis().SetTitleSize(32)
    dummy.GetYaxis().SetLabelFont(43)
    dummy.GetYaxis().SetLabelSize(28)

    # plot all materials
    legend = ROOT.TLegend(0.2, 0.6, 0.5, 0.94)
    legend.SetLineColor(0)
    ths = ROOT.THStack()
    for i, material in enumerate(material_hists.keys()):
        linecolor = 1
        fillcolor = FCCStyle.fillcolors[i if i < 7 else 0]
        match material:
            case "CarbonFiber":
                fillcolor = FCCStyle.fillcolors[0]
            case "CarbonFoam":
                fillcolor = FCCStyle.fillcolors[0]
            case "CarbonFleece":
                fillcolor = ROOT.kBlack
            case "Rohacell":
                fillcolor = FCCStyle.fillcolors[4]
            case "Silicon":
                fillcolor = FCCStyle.fillcolors[2]
            case "Aluminum":
                fillcolor = FCCStyle.fillcolors[1]
            case "Kapton":
                fillcolor = FCCStyle.fillcolors[3]
            case "GlueEcobond45":
                fillcolor = FCCStyle.fillcolors[6]
            case "Water":
                fillcolor = FCCStyle.fillcolors[5]
            case "PCB":
                fillcolor = ROOT.kGreen

        material_hists[material].SetLineColor(linecolor)
        material_hists[material].SetFillColor(fillcolor)
        material_hists[material].SetLineWidth(1)
        material_hists[material].SetFillStyle(1001)
        ths.Add(material_hists[material])

    for i, material in enumerate(reversed(material_hists.keys())):
        legend.AddEntry(material_hists[material], material, "f")

    c = ROOT.TCanvas("", "", 800, 800)
    c.SetTopMargin(0.055)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.11)
    dummy.GetYaxis().SetRangeUser(0, args.ymax if args.ymax > 0 else ths.GetStack().Last().GetMaximum()*1.05)
    dummy.Draw("HIST")
    ths.Draw("SAME HIST")

    legend.SetTextSize(0.04)
    legend.Draw()

    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()
    c.SaveAs(f"{args.outDir}/{args.param}.pdf")
    c.SaveAs(f"{args.outDir}/{args.param}.png")
    fOut = ROOT.TFile(f"{args.outDir}/{args.param}.root", "RECREATE")
    for h in material_hists.values():
        h.Write()
    fOut.Close()
    del ths, legend, c



    # plot grouped materials
    if angleDef == "theta":
        legend = ROOT.TLegend(0.5, 0.9-0.04*(len(material_grouping)+1), 0.9, 0.9)
    else:
        legend = ROOT.TLegend(0.2, 0.9-0.04*(len(material_grouping)+1), 0.55, 0.9)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    legend.SetMargin(0.2)
    legend.SetHeader("")

    ths = ROOT.THStack()
    for group in material_grouping.keys():
        material_grouping_hists[group].SetLineColor(ROOT.kBlack)
        material_grouping_hists[group].SetFillColor(colors_grouping[group])
        material_grouping_hists[group].SetLineWidth(1)
        material_grouping_hists[group].SetFillStyle(1001)
        ths.Add(material_grouping_hists[group])
        legend.AddEntry(material_grouping_hists[group], group, "F")

    c = ROOT.TCanvas("", "", 800, 800)
    c.SetTopMargin(0.055)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.11)
    dummy.GetYaxis().SetRangeUser(0, args.ymax if args.ymax > 0 else ths.GetStack().Last().GetMaximum()*1.05)
    dummy.Draw("HIST")
    ths.Draw("SAME HIST")
    legend.Draw()

    ROOT.gPad.SetTicks()
    ROOT.gPad.RedrawAxis()
    c.SaveAs(f"{args.outDir}/{args.param}_grouped.pdf")
    c.SaveAs(f"{args.outDir}/{args.param}_grouped.png")
    fOut = ROOT.TFile(f"{args.outDir}/{args.param}_grouped.root", "RECREATE")
    for h in material_grouping_hists.values():
        h.Write()
    fOut.Close()


if __name__ == "__main__":
    main()
