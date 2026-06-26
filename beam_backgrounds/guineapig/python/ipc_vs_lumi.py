#!/usr/bin/env python3

import json
import argparse
import ROOT
from array import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

def read_json(path, base_dir):
    with open(f"{path}/summary.json") as f:
        return json.load(f)


def get_nested(d, key):
    for k, val in d.items():
        if key in val:
            return float(val[key])
    print(f"Key {key} not found in json!")
    return None



def make_graph_errors(xvals, yvals, yerrs, name="g"):
    # sort simultaneously by x
    pairs = sorted(zip(xvals, yvals, yerrs), key=lambda p: p[0])

    x_sorted  = array('d', [p[0] for p in pairs])
    y_sorted  = array('d', [p[1] for p in pairs])
    ey_sorted = array('d', [p[2] for p in pairs])

    # no horizontal errors
    ex_sorted = array('d', [0.0] * len(x_sorted))

    g = ROOT.TGraphErrors(
        len(x_sorted),
        x_sorted,
        y_sorted,
        ex_sorted,
        ey_sorted
    )

    g.SetName(name)

    return g

def fit_graph_linear(graph, name="fit"):
    x = graph.GetX()
    n = graph.GetN()

    xmin = min(x[i] for i in range(n))
    xmax = max(x[i] for i in range(n))

    fit = ROOT.TF1(name, "[0] + [1]*x", xmin, xmax)

    graph.Fit(fit, "QSN")


    # make an independent copy for drawing
    fit_draw = fit.Clone(name + "_draw")
    fit_draw.SetLineColor(graph.GetMarkerColor())
    fit_draw.SetLineStyle(2)
    fit_draw.SetLineWidth(2)

    slope = fit_draw.GetParameter(1)
    slope_err = fit_draw.GetParError(1)

    return fit_draw, slope, slope_err


def main():

    inputs = ["FCCee_Z_GHC_V25p1", "FCCee_Z_LCC_V105_v2_50ns", "FCCee_Z_LCC_V105_v2_25ns", "FCCee_Z_LCC_V105", "FCCee_Z_GHC_V25p3_4",  "FCCee_Z_4IP_GHC_V24p4", "FCCee_Z_GHC_V23", "FCCee_Z_CDR"]

    #FCCee_ZH_GHC_V25p1
    #FCCee_TOP_GHC_V25p1
    #FCCee_WW_GHC_V25p1

    base_dir = "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc/"
    out_dir = "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc/plots/"

    colors = [
        ROOT.kBlack,
        ROOT.kRed + 1,
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 1,
    ]





    lumi_tot, data_tot, err_tot = [], [], []
    lumi_ll, data_ll, err_ll = [], [], []
    lumi_bh, data_bh, err_bh = [], [], []
    lumi_bw, data_bw, err_bw = [], [], []
    for i,tag in enumerate(inputs):

        data = read_json(f"{base_dir}/{tag}/Z256_2T_grids8/", base_dir)
        lumi_ee = get_nested(data, "lumi_ee") / 1e28
        pairs_n_average = get_nested(data, "pairs_n_average")
        pairs_n_average_BW = get_nested(data, "pairs_n_average_BW")
        pairs_n_average_BH = get_nested(data, "pairs_n_average_BH")
        pairs_n_average_LL = get_nested(data, "pairs_n_average_LL")

        pairs_n_95pctquant = get_nested(data, "pairs_n_95pctquant")
        pairs_n_95pctquant_BW = get_nested(data, "pairs_n_95pctquant_BW")
        pairs_n_95pctquant_BH = get_nested(data, "pairs_n_95pctquant_BH")
        pairs_n_95pctquant_LL = get_nested(data, "pairs_n_95pctquant_LL")

        lumi_tot.append(lumi_ee)
        data_tot.append(pairs_n_average)
        err_tot.append((pairs_n_95pctquant-pairs_n_average)/1.645)

        lumi_ll.append(lumi_ee)
        data_ll.append(pairs_n_average_LL)
        err_ll.append((pairs_n_95pctquant_LL-pairs_n_average_LL)/1.645)

        lumi_bh.append(lumi_ee)
        data_bh.append(pairs_n_average_BH)
        err_bh.append((pairs_n_95pctquant_BH-pairs_n_average_BH)/1.645)

        lumi_bw.append(lumi_ee)
        data_bw.append(pairs_n_average_BW)
        err_bw.append((pairs_n_95pctquant_BW-pairs_n_average_BW)/1.645)


    g_tot = make_graph_errors(lumi_tot, data_tot, err_tot, "g_tot")
    g_ll  = make_graph_errors(lumi_ll,  data_ll,  err_ll,  "g_ll")
    g_bh  = make_graph_errors(lumi_bh,  data_bh,  err_bh,  "g_bh")
    g_bw  = make_graph_errors(lumi_bw,  data_bw,  err_bw,  "g_bw")

    g_tot.SetLineColor(colors[0])
    g_tot.SetMarkerColor(colors[0])
    g_tot.SetMarkerStyle(20)
    g_tot.SetLineWidth(1)
    

    g_ll.SetLineColor(colors[1])
    g_ll.SetMarkerColor(colors[1])
    g_ll.SetMarkerStyle(20)
    g_ll.SetLineWidth(1)
    

    g_bh.SetLineColor(colors[2])
    g_bh.SetMarkerColor(colors[2])
    g_bh.SetMarkerStyle(20)
    g_bh.SetLineWidth(1)
    

    g_bw.SetLineColor(colors[3])
    g_bw.SetMarkerColor(colors[3])
    g_bw.SetMarkerStyle(20)
    g_bw.SetLineWidth(1)



    canvas = ROOT.TCanvas("c", "c", 800, 800)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.05)
    canvas.SetGrid()

    frame = ROOT.TH1F("frame", "", 100, 3, 9)
    frame.GetXaxis().SetTitle("Integrated luminosity (10^{28}cm^{-2})")
    frame.GetYaxis().SetTitle("Number of IPC pairs")
    frame.SetLineColor(0)
    frame.SetMarkerColor(0)
    frame.SetFillStyle(0)
    frame.Draw()

    frame.SetMinimum(0)
    frame.SetMaximum(max(data_tot)*1.6)  

    
    g_tot.Draw("PE SAME")
    g_ll.Draw("PE SAME")
    g_bh.Draw("PE SAME")
    g_bw.Draw("PE SAME")

    fit_tot, slope_tot, slopeerr_tot = fit_graph_linear(g_tot, "fit_tot")
    fit_ll,  slope_ll,  slopeerr_ll  = fit_graph_linear(g_ll,  "fit_ll")
    fit_bh,  slope_bh,  slopeerr_bh  = fit_graph_linear(g_bh,  "fit_bh")
    fit_bw,  slope_bw,  slopeerr_bw  = fit_graph_linear(g_bw,  "fit_bw")

    print("TOT:", slope_tot, "+/-", slopeerr_tot)
    print("LL :", slope_ll,  "+/-", slopeerr_ll)
    print("BH :", slope_bh,  "+/-", slopeerr_bh)
    print("BW :", slope_bw,  "+/-", slopeerr_bw)

    legend = ROOT.TLegend(0.18, 0.95-0.04*(len(inputs)), 0.45, 0.95)
    legend.SetHeader("")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)

    legend.AddEntry(g_tot, f"All pairs (#sigma = {slope_tot*0.1:.2f}#pm{slopeerr_tot*0.1:.2f} mbarn)", "lp")
    legend.AddEntry(g_ll, f"Landau-Lifshitz (LL) (#sigma = {slope_ll*0.1:.2f}#pm{slopeerr_ll*0.1:.2f} mbarn)", "lp")
    legend.AddEntry(g_bh, f"Bethe-Heitler (BH) (#sigma = {slope_bh*0.1:.2f}#pm{slopeerr_bh*0.1:.2f} mbarn) ", "lp")
    legend.AddEntry(g_bw, f"Breit-Wheeler (BW) (#sigma = {slope_bw*0.1:.2f}#pm{slopeerr_bw*0.1:.2f} mbarn)", "lp")

    legend.Draw()

    fit_tot.Draw("SAME")
    fit_ll.Draw("SAME")
    fit_bh.Draw("SAME")
    fit_bw.Draw("SAME")

    # Optional text
    extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}"
    extra_text_right = "91.2 GeV"
    if extra_text_left:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatexNDC(0.16, 0.95, extra_text_left)
    if extra_text_right:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(33)
        latex.DrawLatexNDC(0.95, 0.95, extra_text_right)

    canvas.SetGrid()
    canvas.Modified()
    canvas.Update()

    canvas.SaveAs(f"{out_dir}/ipc_vs_lumi.png")
    canvas.SaveAs(f"{out_dir}/ipc_vs_lumi.pdf")

if __name__ == "__main__":
    main()