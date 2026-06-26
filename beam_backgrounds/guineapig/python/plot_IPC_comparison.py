#!/usr/bin/env python3

import os
import ROOT
import plotter
import config as gpconfig
import argparse

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TGaxis.SetExponentOffset(-0.08, 0.01, "y")


def get_hist(file_handle, hist_name, unique_tag=""):
    h = file_handle.Get(hist_name)
    if not h:
        raise RuntimeError(
            f"Histogram '{hist_name}' not found in file '{file_handle.GetName()}'."
        )

    out_name = f"{h.GetName()}_{unique_tag}" if unique_tag else h.GetName()
    h = h.Clone(out_name)
    h.SetDirectory(0)
    return h


def main():

    output_dir = "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/comparison/LCC/"
    input_files = [
        "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc/FCCee_Z_LCC_V105/Z256_2T_grids8/output.root",
        "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc/FCCee_Z_LCC_V105_v2_25ns/Z256_2T_grids8/output.root",
        "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc/FCCee_Z_LCC_V105_v2_50ns/Z256_2T_grids8/output.root",
    ]
    labels = ["LCC v105", "LCC 25ns", "LCC 50ns"]



    os.system(f"mkdir -p {output_dir}")
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/guineapig/validation/index.php {output_dir}")

    # Histograms to overlay from the same file
    hist_configs_1d = [
        {
            "hist_name": "h_lumi",
            "outname": "h_lumi",
            "draw_option": "HIST",
        },
        {
            "hist_name": "pairs0/pairs_type",
            "scale_factor": 100,
            "outname": "pairs_type",
            "x_title": ["Breit-Wheeler (BW)", "Bethe-Heitler (BH)", "Landau-Lifshitz (LL)"],
            "y_title": "Fraction (%)",
            "x_range": (0, 3),
            "y_range": None,
            "logy": False,
            "normalize": True,
            "rebin": 1,
            "draw_option": "hist",
        },
        {
            "hist_name": "pairs0/pairs_n",
            "outname": "pairs_n",
            "x_title": "Pair multiplicity",
            "y_title": "Entries",
            "x_range": (0, 2000),
            "rebin": 5
        },
        {
            "hist_name": "pairs0/pairs_theta",
            "outname": "pairs0_theta",
            "normalize": False,
        },
        {
            "hist_name": "pairs/pairs_theta",
            "outname": "pairs1_theta",
            "normalize": False,
        },
        {
            "hist_name": "pairs0/pairs_phi",
            "outname": "pairs0_phi",
            "normalize": False,
        },
        {
            "hist_name": "pairs/pairs_phi",
            "outname": "pairs1_phi",
            "normalize": False,
        },
        {
            "hist_name": "pairs0/pairs_E_log10",
            "outname": "pairs0_E",
            "normalize": False,
        },
        {
            "hist_name": "pairs/pairs_E_log10",
            "outname": "pairs1_E",
            "normalize": False,
        },
        {
            "hist_name": "pairs0/pairs_p_log10",
            "outname": "pairs0_p",
            "normalize": False,
        },
        {
            "hist_name": "pairs/pairs_p_log10",
            "outname": "pairs1_p",
            "normalize": False,
        },
        {
            "hist_name": "pairs0/pairs_pt_log10",
            "outname": "pairs0_pt",
            "normalize": False,
        },
        {
            "hist_name": "pairs/pairs_pt_log10",
            "outname": "pairs1_pt",
            "normalize": False,
        },
        {
            "hist_name": "pairs0/pairs_x",
            "outname": "pairs0_x",
            "normalize": True,
        },
        {
            "hist_name": "pairs/pairs_x",
            "outname": "pairs1_x",
            "normalize": True,
        },

        
        {
            "hist_name": "pairs0/pairs_y",
            "outname": "pairs0_y",
            "normalize": True,
        },
        {
            "hist_name": "pairs/pairs_y",
            "outname": "pairs1_y",
            "normalize": True,
        },

        
        {
            "hist_name": "pairs0/pairs_z",
            "outname": "pairs0_z",
            "normalize": True,
        },
        {
            "hist_name": "pairs/pairs_z",
            "outname": "pairs1_z",
            "normalize": True,
        },
    ]

    # Open all files
    root_files = []
    for fname in input_files:
        f = ROOT.TFile.Open(fname)
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open ROOT file: {fname}")
        root_files.append(f)

    colors = [
        ROOT.kBlack,
        ROOT.kRed + 1,
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 1,
    ]

    line_styles = [1, 1, 1, 1, 2, 2, 3]
    marker_styles = [0] * 10

    for cfg in hist_configs_1d:
        hist_name = cfg["hist_name"]
        print(f"Plotting {cfg['outname']}")

        hists = []
        for i, f in enumerate(root_files):
            h = get_hist(f, hist_name, unique_tag=f"file{i}")
            hists.append(h)

        plotter.plot_hists_1d(
            hists=hists,
            labels=labels,
            outname=os.path.join(output_dir, cfg["outname"]),
            x_title=cfg.get("x_title", None),
            y_title=cfg.get("y_title", None),
            x_range=cfg.get("x_range", None),
            y_range=cfg.get("y_range", None),
            logy=cfg.get("logy", False),
            logx=cfg.get("logx", False),
            colors=cfg.get("colors", colors),
            line_width=3,
            line_styles=line_styles,
            marker_styles=marker_styles,
            marker_size=0.0,
            canvas_size=(800, 800),
            normalize=cfg.get("normalize", False),
            rebin=cfg.get("rebin", False),
            scale_factor=[cfg.get("scale_factor", 1.0)]*len(hists),
            draw_option=cfg.get("draw_option", "hist"),
            extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
            extra_text_right="", 
        )



if __name__ == "__main__":
    main()