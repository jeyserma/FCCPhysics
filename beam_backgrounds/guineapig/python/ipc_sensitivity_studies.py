#!/usr/bin/env python3

import json
import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

def read_json(path, base_dir):
    with open(f"{base_dir}/{path}/Z256_2T_grids8/summary.json") as f:
        return json.load(f)


def get_nested(d, key):
    for k, val in d.items():
        if key in val:
            return float(val[key])
    print(f"Key {key} not found in json!")
    return None


def make_graph(points):
    """
    points = [(x, y), ...]
    """
    points = sorted(points)

    g = ROOT.TGraph(len(points))
    for i, (x, y) in enumerate(points):
        g.SetPoint(i, x, y)

    return g


def main():
    
    observable, name, leg_pos = "lumi_ee", "Metric: luminosity", "right"
    observable, name, leg_pos = "pairs_n_average", "Metric: number of IPC pairs", "right"
    vatiations_group_name = "xing" # sigma beta xing np

    base_dir = "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc_sensitivity_studies/"

    nominal = read_json("FCCee_Z_GHC_V25p1/", base_dir)
    nominal_value = get_nested(nominal, observable)

    # Define variations here
    if vatiations_group_name == "sigma":
        max_var_x = 30
        variations = {
            "#sigma_{x}": {
                -25: "FCCee_Z_GHC_V25p1_sigmaX_25pctDown",
                -10: "FCCee_Z_GHC_V25p1_sigmaX_10pctDown",
                10: "FCCee_Z_GHC_V25p1_sigmaX_10pctUp",
                25: "FCCee_Z_GHC_V25p1_sigmaX_25pctUp",
            },
            "#sigma_{y}": {
                -25: "FCCee_Z_GHC_V25p1_sigmaY_25pctDown",
                -10: "FCCee_Z_GHC_V25p1_sigmaY_10pctDown",
                10: "FCCee_Z_GHC_V25p1_sigmaY_10pctUp",
                25: "FCCee_Z_GHC_V25p1_sigmaY_25pctUp",
            },
            "#sigma_{z}": {
                -25: "FCCee_Z_GHC_V25p1_sigmaZ_25pctDown",
                -10: "FCCee_Z_GHC_V25p1_sigmaZ_10pctDown",
                10: "FCCee_Z_GHC_V25p1_sigmaZ_10pctUp",
                25: "FCCee_Z_GHC_V25p1_sigmaZ_25pctUp",
            },
        }

    if vatiations_group_name == "beta":
        max_var_x = 30
        variations = {
            "#beta_{x}": {
                -25: "FCCee_Z_GHC_V25p1_betaX_25pctDown",
                -10: "FCCee_Z_GHC_V25p1_betaX_10pctDown",
                10: "FCCee_Z_GHC_V25p1_betaX_10pctUp",
                25: "FCCee_Z_GHC_V25p1_betaX_25pctUp",
            },
            "#beta_{y}": {
                -25: "FCCee_Z_GHC_V25p1_betaY_25pctDown",
                -10: "FCCee_Z_GHC_V25p1_betaY_10pctDown",
                10: "FCCee_Z_GHC_V25p1_betaY_10pctUp",
                25: "FCCee_Z_GHC_V25p1_betaY_25pctUp",
            },
        }

    if vatiations_group_name == "xing":
        max_var_x = 15
        variations = {
            "Crossing angle": {
                -10: "FCCee_Z_GHC_V25p1_xing_10pctUp", # inverted
                -5: "FCCee_Z_GHC_V25p1_xing_5pctDown",
                5: "FCCee_Z_GHC_V25p1_xing_5pctUp",
                10: "FCCee_Z_GHC_V25p1_xing_10pctDown", # inverted
            }
        }
    if vatiations_group_name == "np":
        max_var_x = 30
        variations = {
            "Bunch intensity": {
                -25: "FCCee_Z_GHC_V25p1_np_25pctDown",
                -10: "FCCee_Z_GHC_V25p1_np_10pctDown",
                10: "FCCee_Z_GHC_V25p1_np_10pctUp",
                25: "FCCee_Z_GHC_V25p1_np_25pctUp",
            }
        }

    canvas = ROOT.TCanvas("c", "c", 800, 800)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.05)
    canvas.SetGrid()

    frame = ROOT.TH1F("frame", "", 100, -max_var_x, max_var_x)
    frame.GetXaxis().SetTitle("Beam-parameter variation (%)")
    frame.GetYaxis().SetTitle("Relative deviation w.r.t. nominal (%)")
    frame.SetLineColor(0)
    frame.SetMarkerColor(0)
    frame.SetFillStyle(0)
    frame.Draw()

    
    if leg_pos == "left":
        legend = ROOT.TLegend(0.18, 0.88-0.04*(len(variations)+1), 0.45, 0.88)
    if leg_pos == "right":
        legend = ROOT.TLegend(0.50, 0.88-0.04*(len(variations)+1), 0.90, 0.88)
    legend.SetHeader(name)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)

    colors = [
        ROOT.kBlue + 1,
        ROOT.kRed + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
    ]

    graphs = []
    all_rel_devs = []
    for i, (label, files) in enumerate(variations.items()):
        points = []

        for delta, path in files.items():
            data = read_json(path, base_dir)
            value = get_nested(data, observable)

            print(value, nominal_value)
            rel_dev = (value / nominal_value - 1.0)*100.
            points.append((delta, rel_dev))
            all_rel_devs.append(rel_dev)

        g = make_graph(points)
        g.SetLineColor(colors[i % len(colors)])
        g.SetMarkerColor(colors[i % len(colors)])
        g.SetMarkerStyle(20 + i)
        g.SetLineWidth(2)

        g.Draw("LP SAME")
        legend.AddEntry(g, label, "lp")
        graphs.append(g)

    line = ROOT.TLine(-max_var_x, 0, max_var_x, 0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("SAME")

    legend.Draw()

    max_dev = max(abs(x) for x in all_rel_devs) * 1.2
    max_dev = round(max_dev, 3)
    frame.SetMinimum(-max_dev)
    frame.SetMaximum(+max_dev)

    # Optional text
    extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}"
    extra_text_right = "GHC 25.1 (FSR), 91.2 GeV"
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

    canvas.SaveAs(f"{base_dir}/{vatiations_group_name}_{observable}.png")
    canvas.SaveAs(f"{base_dir}/{vatiations_group_name}_{observable}.pdf")

if __name__ == "__main__":
    main()