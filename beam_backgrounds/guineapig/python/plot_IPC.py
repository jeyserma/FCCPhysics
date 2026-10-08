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

    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(i) == 0:
            h.SetBinContent(i, -1)
            h.SetBinError(i, 0.0)


    return h


def main():

    input_dir = "/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_studies/FCCee_Z_LCC_V105/Z128_2T_grids1/"
    input_dir = "/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_secondary_grid_studies/FCCee_Z_GHC_V25p1/CFG_TEST_2G_CLOSURE"
    #input_dir = "/ceph/submit/data/group/fcc/ee/beam_backgrounds/guineapig/ipc_primary_grid_studies/FCCee_Z_GHC_V25p1/CFG_GRIDT_64_64_64_1/"

    '''
CFG_GRIDD_16_350_16
CFG_GRIDD_32_700_32
CFG_GRIDD_64_1400_64
CFG_GRIDD_128_2800_128

CFG_GRIDD1_16_175_16
CFG_GRIDD1_32_350_32
CFG_GRIDD1_64_700_64
CFG_GRIDD1_128_1400_128


    '''
    ###################################
    input_dir = input_dir.rstrip('\\').rstrip('/').replace("//", "/")
    parameter_set = input_dir.split("/")[-1]
    accelerator = input_dir.split("/")[-2]
    campaign =  input_dir.split("/")[-3]
    label = gpconfig.acc_cfg[accelerator]['label']

    output_dir = f"/home/submit/jaeyserm/public_html/fccee/guineapig/validation/{campaign}/{accelerator}/{parameter_set}/"
    os.system(f"mkdir -p {output_dir}")
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/guineapig/validation/index.php {output_dir}")
    input_file = f"{output_dir}/output.root"

    # Histograms to overlay from the same file
    hist_configs_1d = [
        {
            "hist_names": ["h_lumi"],
            "labels": [],
            "outname": "h_lumi",
            "x_range": None,
            "y_range": None,
            "logy": False,
            "normalize": False,
            "rebin": 1,
            "draw_option": "HIST",
            "colors": [ROOT.kRed + 1]
        },
        {
            "hist_names": ["pairs0/pairs_t"],
            "labels": [],
            "outname": "pairs0_t",
            "x_range": None,
            "y_range": None,
            "logy": False,
            "normalize": False,
            "rebin": 100,
            "draw_option": "HIST",
            "colors": [ROOT.kRed + 1],
            "x_range": (-150, 150),
        },
        {
            "hist_names": ["pairs/pairs_t"],
            "labels": [],
            "outname": "pairs1_t",
            "x_range": None,
            "y_range": None,
            "logy": False,
            "normalize": False,
            "rebin": 100,
            "draw_option": "HIST",
            "colors": [ROOT.kRed + 1],
            "x_range": (-150, 150),
        },
        {
            "hist_names": ["pairs0/pairs_type"],
            "labels": [],
            "scale_factor": [100],
            "outname": "pairs_type",
            "x_title": ["Breit-Wheeler (BW)", "Bethe-Heitler (BH)", "Landau-Lifshitz (LL)"],
            "y_title": "Fraction (%)",
            "x_range": (0, 3),
            "y_range": None,
            "logy": False,
            "normalize": True,
            "rebin": 1,
            "draw_option": "hist text",
            "colors": [ROOT.kRed + 1]
        },
        {
            "hist_names": ["pairs0/pairs_n", "pairs0/pairs_n_ll", "pairs0/pairs_n_bh", "pairs0/pairs_n_bw"],
            "hist_names": ["pairs0/pairs_n"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "scale_factor": [1, 1, 0.5, 0.15],
            "outname": "pairs_n",
            "x_title": "Pair multiplicity",
            "y_title": "Entries",
            "x_range": (0, 2500),
            "y_range": None,
            "logy": False,
            "normalize": False,
            "rebin": 5
        },
        {
            "hist_names": ["pairs0/pairs_theta", "pairs0/pairs_theta_ll", "pairs0/pairs_theta_bh", "pairs0/pairs_theta_bw"],
            "hist_names": ["pairs0/pairs_theta"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_theta",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_theta", "pairs/pairs_theta_ll", "pairs/pairs_theta_bh", "pairs/pairs_theta_bw"],
            "hist_names": ["pairs/pairs_theta"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_theta",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_theta", "pairs/pairs_theta"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_theta",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_phi", "pairs0/pairs_phi_ll", "pairs0/pairs_phi_bh", "pairs0/pairs_phi_bw"],
            "hist_names": ["pairs0/pairs_phi"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_phi",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_phi", "pairs/pairs_phi_ll", "pairs/pairs_phi_bh", "pairs/pairs_phi_bw"],
            "hist_names": ["pairs/pairs_phi"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_phi",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_phi", "pairs/pairs_phi"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_phi",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_E_log10", "pairs0/pairs_E_log10_ll", "pairs0/pairs_E_log10_bh", "pairs0/pairs_E_log10_bw"],
            "hist_names": ["pairs0/pairs_E_log10"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_E",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_E_log10", "pairs/pairs_E_log10_ll", "pairs/pairs_E_log10_bh", "pairs/pairs_E_log10_bw"],
            "hist_names": ["pairs/pairs_E_log10"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_E",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_E_log10", "pairs/pairs_E_log10"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_E",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_p_log10", "pairs0/pairs_p_log10_ll", "pairs0/pairs_p_log10_bh", "pairs0/pairs_p_log10_bw"],
            "hist_names": ["pairs0/pairs_p_log10"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_p",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_p_log10", "pairs/pairs_p_log10_ll", "pairs/pairs_p_log10_bh", "pairs/pairs_p_log10_bw"],
            "hist_names": ["pairs/pairs_p_log10"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_p",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_p_log10", "pairs/pairs_p_log10"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_p",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_pt_log10", "pairs0/pairs_pt_log10_ll", "pairs0/pairs_pt_log10_bh", "pairs0/pairs_pt_log10_bw"],
            "hist_names": ["pairs0/pairs_pt_log10"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_pt",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_pt_log10", "pairs/pairs_pt_log10_ll", "pairs/pairs_pt_log10_bh", "pairs/pairs_pt_log10_bw"],
            "hist_names": ["pairs/pairs_pt_log10"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_pt",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_pt_log10", "pairs/pairs_pt_log10"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_pt",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_x", "pairs0/pairs_x_ll", "pairs0/pairs_x_bh", "pairs0/pairs_x_bw"],
            "hist_names": ["pairs0/pairs_x"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_x",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_x", "pairs/pairs_x_ll", "pairs/pairs_x_bh", "pairs/pairs_x_bw"],
            "hist_names": ["pairs/pairs_x"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_x",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_x", "pairs/pairs_x"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_x",
            "normalize": True,
        },

        
        {
            "hist_names": ["pairs0/pairs_y", "pairs0/pairs_y_ll", "pairs0/pairs_y_bh", "pairs0/pairs_y_bw"],
            "hist_names": ["pairs0/pairs_y"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_y",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_y", "pairs/pairs_y_ll", "pairs/pairs_y_bh", "pairs/pairs_y_bw"],
            "hist_names": ["pairs/pairs_y"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_y",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_y", "pairs/pairs_y"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_y",
            "normalize": True,
        },

        
        {
            "hist_names": ["pairs0/pairs_z", "pairs0/pairs_z_ll", "pairs0/pairs_z_bh", "pairs0/pairs_z_bw"],
            "hist_names": ["pairs0/pairs_z"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs0_z",
            "normalize": True,
        },
        {
            "hist_names": ["pairs/pairs_z", "pairs/pairs_z_ll", "pairs/pairs_z_bh", "pairs/pairs_z_bw"],
            "hist_names": ["pairs/pairs_z"],
            "labels": ["All pairs", "Landau-Lifshitz (LL)", "Bethe-Heitler (BH)", "Breit-Wheeler (BW)"],
            "outname": "pairs1_z",
            "normalize": True,
        },
        {
            "hist_names": ["pairs0/pairs_z", "pairs/pairs_z"],
            "labels": ["Pairs0", "Pairs1"],
            "outname": "pairs_z",
            "normalize": True,
        },
    ]


    hist_configs_2d = [
        {
            "hist_name": "pairs0/pairs_theta_pt",
            "outname": "pairs0_theta_pt",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs0/pairs_theta_pt_bw",
            "outname": "pairs0_theta_pt_bw",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs0/pairs_theta_pt_bh",
            "outname": "pairs0_theta_pt_bh",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs0/pairs_theta_pt_ll",
            "outname": "pairs0_theta_pt_ll",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs0/pairs_theta_pt_acc",
            "outname": "pairs0_theta_pt_acc",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs/pairs_theta_pt",
            "outname": "pairs1_theta_pt",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs/pairs_theta_pt_bw",
            "outname": "pairs1_theta_pt_bw",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs/pairs_theta_pt_bh",
            "outname": "pairs1_theta_pt_bh",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs/pairs_theta_pt_ll",
            "outname": "pairs1_theta_pt_ll",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs/pairs_theta_pt_acc",
            "outname": "pairs1_theta_pt_acc",
            "draw_option": "COLZ",
            "x_range": (-3, 0.2),
        },
        {
            "hist_name": "pairs0/pairs_xy",
            "outname": "pairs0_xy",
            "draw_option": "COLZ",
        },
        {
            "hist_name": "pairs/pairs_xy",
            "outname": "pairs1_xy",
            "draw_option": "COLZ",
        },
        
        {
            "hist_name": "pairs0/pairs_xz",
            "outname": "pairs0_xz",
            "draw_option": "COLZ",
        },
        {
            "hist_name": "pairs/pairs_xz",
            "outname": "pairs1_xz",
            "draw_option": "COLZ",
        },
        
        {
            "hist_name": "pairs0/pairs_yz",
            "outname": "pairs0_yz",
            "draw_option": "COLZ",
        },
        {
            "hist_name": "pairs/pairs_yz",
            "outname": "pairs1_yz",
            "draw_option": "COLZ",
        },
    ]
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {input_file}")

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
        hist_names = cfg["hist_names"]
        labels = cfg["labels"]

        print(f"Plotting {cfg['outname']}")

        hists = []
        for i, hist_name in enumerate(hist_names):
            h = get_hist(f, hist_name, unique_tag=f"h{i}")
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
            scale_factor=cfg.get("scale_factor", False),
            draw_option=cfg.get("draw_option", "hist"),
            extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
            extra_text_right=label, 
        )

    for cfg in hist_configs_2d:
        hist_name = cfg["hist_name"]

        print(f"Plotting {cfg['outname']}")

        hist = get_hist(f, hist_name)

        plotter.plot_hist_2d(
            hist=hist,
            outname=os.path.join(output_dir, cfg["outname"]),
            x_title=cfg.get("x_title", None),
            y_title=cfg.get("y_title", None),
            z_title=cfg.get("z_title", None),
            x_range=cfg.get("x_range", None),
            y_range=cfg.get("y_range", None),
            z_range=cfg.get("z_range", None),
            canvas_size=(900, 800),
            draw_option=cfg.get("draw_option", "colz"),
            extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
            extra_text_right=label, 
        )

    f.Close()


if __name__ == "__main__":
    main()