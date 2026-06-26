#!/usr/bin/env python3

import os
import ROOT
import plotter
import config as gpconfig
import argparse
import array

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "y")


def make_plot(theta_pt, theta_pt_acc, outname, label):

    fraction = 100.* theta_pt_acc.Integral() / theta_pt.Integral()
    print(f"Fraction of events (%): {fraction:.2f}")

    h = theta_pt.Clone(f"{theta_pt.GetName()}_2d_clone")
    h.SetDirectory(0)

    h.GetXaxis().SetRangeUser(-3, 0.2)
    #h.GetYaxis().SetRangeUser(y_range[0], y_range[1])

    #h.GetXaxis().SetTitle(x_title)
    #h.GetYaxis().SetTitle(y_title)
    #h.GetZaxis().SetTitle(z_title)

    # Axis sizes
    h.GetXaxis().SetTitleSize(0.04)
    h.GetYaxis().SetTitleSize(0.04)
    h.GetZaxis().SetTitleSize(0.04)

    h.GetXaxis().SetLabelSize(0.04)
    h.GetYaxis().SetLabelSize(0.04)
    h.GetZaxis().SetLabelSize(0.04)

    # Offsets
    h.GetXaxis().SetTitleOffset(1.15)
    h.GetYaxis().SetTitleOffset(1.35)
    h.GetZaxis().SetTitleOffset(1.20)

    h.GetXaxis().SetLabelOffset(0.01)
    h.GetYaxis().SetLabelOffset(0.01)
    h.GetZaxis().SetLabelOffset(0.01)

    h.GetXaxis().SetNdivisions(510)
    h.GetYaxis().SetNdivisions(510)
    h.GetZaxis().SetNdivisions(510)

    c = ROOT.TCanvas(f"c_{h.GetName()}", "", 900, 800)

    # Slightly wider right margin for color palette
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.18)

    h.Draw("COLZ")

    theta_pt_acc.SetFillColorAlpha(ROOT.kRed, 0.75)
    theta_pt_acc.SetMarkerColorAlpha(ROOT.kRed, 0.75)
    theta_pt_acc.Draw("BOX SAME")

    c.Modified()
    c.Update()


    # Optional text
    extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}"
    extra_text_right = label
    if extra_text_left:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatexNDC(0.14, 0.95, extra_text_left)
    if extra_text_right:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(33)
        latex.DrawLatexNDC(0.82, 0.95, extra_text_right)

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.DrawLatex(0.18, 0.85, f"Acceptance fraction: {fraction:.2f} % ")
    latex.DrawLatex(0.18, 0.80, f"2T, R=13 mm, z=100 mm")

    c.Modified()
    c.Update()
    c.SaveAs(f"{outname}.png")
    c.SaveAs(f"{outname}.pdf")

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--accelerator", type=str, help="Accelerator config", default="FCCee_Z_GHC_V25p1")
    parser.add_argument("--parameter_set", type=str, help="Parameter set", default="Z256_2T_grids8")
    parser.add_argument("--campaign", type=str, help="Campaign (as in config.py)", choices=["ipc", "ipc_studies", "ipc_sensitivity_studies"], default="ipc")
    args = parser.parse_args()

    #FCCee_Z_GHC_V25p1
    #FCCee_Z_LCC_V105_v2_50ns
    #FCCee_Z_LCC_V105_v2_25ns
    #FCCee_Z_LCC_V105
    #FCCee_Z_GHC_V25p3_4
    #FCCee_TOP_GHC_V25p1
    #FCCee_WW_GHC_V25p1
    #FCCee_ZH_GHC_V25p1
    #FCCee_Z_4IP_GHC_V24p4
    #FCCee_Z_CDR
    #FCCee_Z_GHC_V23


    
    accelerator = args.accelerator
    parameter_set = args.parameter_set
    campaign = args.campaign
    cfg = getattr(gpconfig, campaign)
    label = cfg[accelerator][parameter_set]['label']
    
    output_dir = f"/home/submit/jaeyserm/public_html/fccee/guineapig/validation/{campaign}/{accelerator}/{parameter_set}/"
    os.system(f"mkdir -p {output_dir}")
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/guineapig/validation/index.php {output_dir}")
    input_file = f"{output_dir}/output.root"

    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {input_file}")

    theta_pt = f.Get("pairs/pairs_theta_pt")
    theta_pt_acc = f.Get("pairs/pairs_theta_pt_acc")
    make_plot(theta_pt, theta_pt_acc, f"{output_dir}/pairs1_theta_pt_acceptance_overlay", label)

    theta_pt = f.Get("pairs0/pairs_theta_pt")
    theta_pt_acc = f.Get("pairs0/pairs_theta_pt_acc")
    make_plot(theta_pt, theta_pt_acc, f"{output_dir}/pairs0_theta_pt_acceptance_overlay", label)

if __name__ == "__main__":
    main()