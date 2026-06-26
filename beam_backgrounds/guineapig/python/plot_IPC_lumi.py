#!/usr/bin/env python3

import os
import ROOT
import plotter
import config as gpconfig
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "y")


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--accelerator", type=str, help="Accelerator config", default="FCCee_Z_GHC_V23")
    parser.add_argument("--parameter_set", type=str, help="Parameter set", default="Z256_2T_grids8")
    parser.add_argument("--maxfiles", type=int, help="number of files", default=-1)
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
    maxfiles = args.maxfiles
    campaign = args.campaign
    cfg = getattr(gpconfig, campaign)
    label = cfg[accelerator][parameter_set]['label']
    vertical_offset = 1.3
    
    output_dir = f"/home/submit/jaeyserm/public_html/fccee/guineapig/validation/{campaign}/{accelerator}/{parameter_set}/"
    os.system(f"mkdir -p {output_dir}")
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/guineapig/validation/index.php {output_dir}")
    input_file = f"{output_dir}/output.root"

    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {input_file}")

    h_lumi_raw = f.Get("h_lumi")
    
    h_lumi = h_lumi_raw.Clone("h_lumi_noraw")
    for i in range(1, h_lumi.GetNbinsX() + 1):
        h_lumi.SetBinContent(i, h_lumi_raw.GetBinContent(i)/1e24) # convert to powers of 24
    h_lumi_int = h_lumi.GetCumulative()
    h_lumi_int.Scale(1e-4) # convert to 1e28

    c = ROOT.TCanvas("c", "two y axes", 800, 650)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.10)
    c.SetGrid()

    h_lumi.GetYaxis().SetLabelSize(h_lumi.GetYaxis().GetLabelSize()*1.1)
    h_lumi.GetYaxis().SetTitleSize(h_lumi.GetYaxis().GetTitleSize()*1.1)

    h_lumi.GetXaxis().SetLabelSize(h_lumi.GetXaxis().GetLabelSize()*1.1)
    h_lumi.GetXaxis().SetTitleSize(h_lumi.GetXaxis().GetTitleSize()*1.1)

    # Draw left-axis histogram
    h_lumi.SetLineColor(ROOT.kBlue + 1)
    h_lumi.GetYaxis().SetTitle("Luminosity (10^{24}cm^{-2}/step)")
    h_lumi.GetXaxis().SetTitle("Beam crossing step")
    h_lumi.SetLineWidth(2)
    h_lumi.SetMaximum(h_lumi.GetMaximum()*vertical_offset)
    h_lumi.Draw("HIST")


    # Scale h2 to the visible left-axis range
    right_min, right_max = 0.0, h_lumi_int.GetMaximum()*vertical_offset
    left_min  = h_lumi.GetMinimum()
    left_max  = h_lumi.GetMaximum()

    scale = (left_max - left_min) / (right_max - right_min)
    h_lumi_int_scaled = h_lumi_int.Clone("h_lumi_int_scaled")
    for i in range(1, h_lumi_int_scaled.GetNbinsX() + 1):
        y = h_lumi_int.GetBinContent(i)
        h_lumi_int_scaled.SetBinContent(i, left_min + (y - right_min) * scale)

    #h_lumi_int_scaled = h_lumi_int
    h_lumi_int_scaled.SetLineColor(ROOT.kRed + 1)
    h_lumi_int_scaled.SetLineWidth(2)
    h_lumi_int_scaled.Draw("hist same")

    # Force pad update before drawing extra axis
    c.Update()

    # Right vertical axis
    x_right = ROOT.gPad.GetUxmax()
    axis_right = ROOT.TGaxis(
        x_right, left_min,
        x_right, left_max,
        right_min, right_max,
        510, "+L"
    )

    axis_right.SetTitle("Integrated luminosity (10^{28}cm^{-2})")
    axis_right.SetLineWidth(1)
    axis_right.SetLabelFont(42)
    axis_right.SetTitleFont(42)
    axis_right.SetLabelSize(h_lumi.GetYaxis().GetLabelSize())
    axis_right.SetTitleSize(h_lumi.GetYaxis().GetTitleSize())
    axis_right.SetTitleOffset(axis_right.GetTitleOffset()*1.1)
    axis_right.Draw()

    # Optional legend
    leg = ROOT.TLegend(0.15, 0.775, 0.75, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetMargin(0.15)
    leg.SetTextSize(0.045)
    leg.AddEntry(h_lumi, "Luminosity per step", "l")
    leg.AddEntry(h_lumi_int_scaled, f"Integrated luminosity (tot. {h_lumi_int.GetMaximum():.1f}#times10^{{28}}cm^{{-2}})", "l")
    leg.Draw()

    # Optional text
    extra_text_left = "#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}"
    extra_text_right = label
    if extra_text_left:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.045)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatexNDC(0.12, 0.96, extra_text_left)
    if extra_text_right:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.045)
        latex.SetTextFont(42)
        latex.SetTextAlign(33)
        latex.DrawLatexNDC(0.90, 0.96, extra_text_right)



    c.SaveAs(f"{output_dir}/h_lumi_int.png")
    c.SaveAs(f"{output_dir}/h_lumi_int.pdf")

if __name__ == "__main__":
    main()