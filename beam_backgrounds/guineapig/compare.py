import os
import glob
import copy

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)


baseDir = "/home/submit/jaeyserm/public_html/fccee/pairs_check"

inputs = {}

inputs["FCCee_Z_4IP_04may23_FCCee_Z128_noTracking"] = {
    "legend": "Default window",
    "color": ROOT.kBlack
}

inputs["FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_doubleZWindow"] = {
    "legend": "Double window in Z",
    "color": ROOT.kRed
}
inputs["FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfYWindow"] = {
    "legend": "Half window in Y",
    "color": ROOT.kGreen+1
}

inputs["FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_reducedWindow"] = {
    "legend": "Reduced window",
    "color": ROOT.kRed
}
inputs["FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfWindow"] = {
    "legend": "Half window in X, Y and Z",
    "color": ROOT.kBlue
}

inputs["FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfZWindow"] = {
    "legend": "Half window in Z",
    "color": ROOT.kBlue
}


zoom_cfg = [0.6, 0.5, 4]

outName = "output"
zoom = False
to_plot = ["FCCee_Z_4IP_04may23_FCCee_Z128_noTracking", "FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfZWindow", "FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_doubleZWindow", "FCCee_Z_4IP_04may23_FCCee_Z128_noTracking_halfYWindow"]


for i, vec in enumerate(["x", "y", "z"]):
    c1 = ROOT.TCanvas("c1", "Histograms", 1000, 1000)
    c1.SetTopMargin(0.055)
    c1.SetRightMargin(0.05)
    c1.SetLeftMargin(0.12)
    c1.SetBottomMargin(0.11)
    
    leg = ROOT.TLegend(0.2, 0.8, 0.8, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetMargin(0.2)
    leg.SetNColumns(2)

    hists = []
    for j, plot in enumerate(to_plot):
        print(plot)
        fIn = ROOT.TFile(f"{baseDir}/{plot}/output.root")
        h = copy.deepcopy(fIn.Get(f"hist_{vec}"))
        fIn.Close()

        leg.AddEntry(h, inputs[plot]['legend'], "L")
        h.SetLineColor(inputs[plot]['color'])
        h.SetLineWidth(2)
        h.Scale(1./h.Integral())
        hists.append(h)
    
    yMax = max([h.GetMaximum() for h in hists])
    for j,h in enumerate(hists):
        if j==0:
            h.Draw("HIST")
            h.GetYaxis().SetRangeUser(0, 1.3*yMax)
            if zoom:
                h.GetXaxis().SetRangeUser(-zoom_cfg[i], zoom_cfg[i])
        else:
            h.Draw("HIST SAME")

    leg.Draw()
    suffix = '_zoom' if zoom else ''
    c1.SaveAs(f"{baseDir}/{outName}_{vec}{suffix}.png")
    del c1



