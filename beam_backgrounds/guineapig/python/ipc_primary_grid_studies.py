#!/usr/bin/env python3

import json
import os
from array import array

import ROOT


def read_json_value(json_path, key=None):
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if key is None:
        value = data
    else:
        value = data

        for part in key.split("."):
            if isinstance(value, dict):
                if part not in value:
                    raise KeyError(
                        f"Key '{part}' not found while reading '{key}' "
                        f"from {json_path}"
                    )
                value = value[part]

            elif isinstance(value, list):
                try:
                    index = int(part)
                except ValueError as error:
                    raise KeyError(
                        f"Expected a list index, but got '{part}' "
                        f"while reading '{key}'"
                    ) from error

                value = value[index]

            else:
                raise KeyError(
                    f"Cannot access '{part}' in object {value!r} "
                    f"while reading '{key}'"
                )

    if not isinstance(value, (int, float)):
        raise TypeError(
            f"Value at '{key}' in {json_path} is not numerical: {value!r}"
        )

    return float(value)


def make_graphs(
    graph_config,
    output="graphs.pdf",
    json_key=None,
    x_title="x",
    y_title="Value",
    title="",
    draw_points=True,
    log_x=False,
    log_y=False,
    x_range=None,
    y_range=None,
    extra_text_left=None,
    extra_text_right=None,
    legend_label=None,
):

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.05)
    canvas.SetBottomMargin(0.12)
    canvas.SetTopMargin(0.08)

    canvas.SetLogx(log_x)
    canvas.SetLogy(log_y)

    nlegs = len(graph_config.items())
    if legend_label:
        nlegs+=1
    legend = ROOT.TLegend(0.55, 0.50-0.05*nlegs, 0.93, 0.50)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.035)
    if legend_label:
        legend.SetHeader(legend_label)

    colors = [
        ROOT.kBlue + 1,
        ROOT.kRed + 1,
        ROOT.kGreen + 2,
        #ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 2,
        ROOT.kViolet + 1,
        ROOT.kGray + 2,
    ]

    markers = [
        20,
        21,
        22,
        23,
        33,
        34,
        29,
        24,
    ]

    multigraph = ROOT.TMultiGraph()
    graphs = {}

    norm = False
    for index, (graph_name, config) in enumerate(graph_config.items()):
        x_values = config["x"]
        json_files = config["files"]
        curve_json_key = config.get("json_key", json_key)

        if len(x_values) != len(json_files):
            raise ValueError(
                f"Graph '{graph_name}' has {len(x_values)} x values "
                f"but {len(json_files)} JSON files."
            )

        y_values = []

        for x_value, json_file in zip(x_values, json_files):
            value_norm = 1.0
            if 'norm' in config:
                norm = True
                try:
                    value_norm = read_json_value(config["norm"], curve_json_key)
                except (FileNotFoundError, KeyError, TypeError, json.JSONDecodeError) as error:
                    raise RuntimeError(
                        f"Could not read point x={x_value} for graph "
                        f"'{graph_name}' from '{json_file}'"
                    ) from error
            try:
                value = read_json_value(json_file, curve_json_key)
            except (FileNotFoundError, KeyError, TypeError, json.JSONDecodeError) as error:
                raise RuntimeError(
                    f"Could not read point x={x_value} for graph "
                    f"'{graph_name}' from '{json_file}'"
                ) from error

            y_values.append(value/value_norm)
            print(
                f"{graph_name:20s} "
                f"x = {x_value:10g}, y = {value:12g}, "
                f"file = {json_file}"
            )

        graph = ROOT.TGraph(
            len(x_values),
            array("d", map(float, x_values)),
            array("d", y_values),
        )

        graph.SetName(
            "graph_" + "".join(
                character if character.isalnum() else "_"
                for character in graph_name
            )
        )
        graph.SetTitle(graph_name)

        color = colors[index % len(colors)]
        marker = markers[index % len(markers)]

        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetLineWidth(3)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(1.2)

        multigraph.Add(graph, "LP" if draw_points else "L")
        legend.AddEntry(graph, graph_name, "LP" if draw_points else "L")
        graphs[graph_name] = graph

    multigraph.SetTitle(f"{title};{x_title};{y_title}")
    multigraph.Draw("A")

    
    if x_range is not None:
        multigraph.GetXaxis().SetLimits(x_range[0], x_range[1])

    if y_range is not None:
        multigraph.SetMinimum(y_range[0])
        multigraph.SetMaximum(y_range[1])

    multigraph.GetXaxis().SetTitleSize(0.045)
    multigraph.GetYaxis().SetTitleSize(0.045)
    multigraph.GetXaxis().SetLabelSize(0.040)
    multigraph.GetYaxis().SetLabelSize(0.040)
    multigraph.GetXaxis().SetTitleOffset(1.15)
    multigraph.GetYaxis().SetTitleOffset(1.35)

    legend.Draw()

    if norm:
        y_reference = 1.0

        x_min = multigraph.GetXaxis().GetXmin()
        x_max = multigraph.GetXaxis().GetXmax()

        reference_line = ROOT.TLine(x_min, y_reference, x_max, y_reference)
        reference_line.SetLineStyle(3)   # dotted
        reference_line.SetLineWidth(2)
        reference_line.SetLineColor(ROOT.kBlack)
        reference_line.Draw("SAME")


    # Optional text
    if extra_text_left:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(13)
        latex.DrawLatexNDC(0.13, 0.95, extra_text_left)
    if extra_text_right:
        latex = ROOT.TLatex()
        latex.SetTextSize(0.035)
        latex.SetTextFont(42)
        latex.SetTextAlign(33)
        latex.DrawLatexNDC(0.95, 0.95, extra_text_right)

    canvas.RedrawAxis()
    canvas.SaveAs(output)
    canvas.SaveAs(output.replace('.png', '.pdf'))

    # Also create a ROOT file so that the graphs can be edited later.
    #root_output = os.path.splitext(output)[0] + ".root"
    #root_file = ROOT.TFile(root_output, "RECREATE")

    #canvas.Write("canvas")
    #multigraph.Write("multigraph")

    #for graph in graphs.values():
    #    graph.Write()

    #root_file.Close()

    print(f"\nSaved plot to: {output}")
    #print(f"Saved ROOT objects to: {root_output}")

    return canvas, multigraph, graphs


if __name__ == "__main__":

    base_input_dir = "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc_primary_grid_studies/FCCee_Z_GHC_V25p1/"
    output_dir = "/home/submit/jaeyserm/public_html/fccee/guineapig/validation/ipc_primary_grid_studies/FCCee_Z_GHC_V25p1/plots/"

    #####################################
    ## GRID SIZE
    #####################################


    graph_config = {
        "All varied": {
            "x": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            "files": [
                f"{base_input_dir}/CFG_MXYZ_0p5_0p5_0p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_1_1_1/summary.json",
                f"{base_input_dir}/CFG_MXYZ_1p5_1p5_1p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_2_2_2/summary.json",
                f"{base_input_dir}/CFG_MXYZ_2p5_2p5_2p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_3_3_3/summary.json",
                f"{base_input_dir}/CFG_MXYZ_3p5_3p5_3p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
        },
        "Variation x": {
            "x": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            "files": [
                f"{base_input_dir}/CFG_MXYZ_0p5_4_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_1_4_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_1p5_4_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_2_4_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_2p5_4_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_3_4_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_3p5_4_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
        },
        
        "Variation y": {
            "x": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            "files": [
                f"{base_input_dir}/CFG_MXYZ_4_0p5_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_1_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_1p5_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_2_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_2p5_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_3_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_3p5_4/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
        },
        "Variation z": {
            "x": [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0],
            "files": [
                f"{base_input_dir}/CFG_MXYZ_4_4_0p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_1/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_1p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_2/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_2p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_3/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_3p5/summary.json",
                f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_MXYZ_4_4_4/summary.json",
        },
    }

    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_lumi_ee.png",
        json_key="luminosity.lumi_ee",
        x_title="m_{x}, m_{y}, m_{z}",
        y_title="L_{ee}/L_{ee}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV", 
        legend_label="Grid size scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_lumi_eg.png",
        json_key="luminosity.lumi_eg",
         x_title="m_{x}, m_{y}, m_{z}",
        y_title="L_{eg}/L_{eg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_lumi_ge.png",
        json_key="luminosity.lumi_ge",
         x_title="m_{x}, m_{y}, m_{z}",
        y_title="L_{ge}/L_{ge}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_lumi_gg.png",
        json_key="luminosity.lumi_gg",
        x_title="m_{x}, m_{y}, m_{z}",
        y_title="L_{gg}/L_{gg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_upsmax.png",
        json_key="luminosity.upsmax",
        x_title="m_{x}, m_{y}, m_{z}",
        y_title="Y/Ymax",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_pairs_n_average.png",
        json_key="ipcs.pairs_n_average",
        x_title="m_{x}, m_{y}, m_{z}",
        y_title="n_{pairs}/n_{pairs}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_pairs_n_average_LL.png",
        json_key="ipcs.pairs_n_average_LL",
        x_title="m_{x}, m_{y}, m_{z}",
        y_title="n_{pairs,LL}/n_{pairs, LL}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_pairs_n_average_BH.png",
        json_key="ipcs.pairs_n_average_BH",
        x_title="m_{x}, m_{y}, m_{z}",
        y_title="n_{pairs,BH}/n_{pairs, BH}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/size/grid_size_pairs_n_average_BW.png",
        json_key="ipcs.pairs_n_average_BW",
        x_title="m_{x}, m_{y}, m_{z}",
        y_title="n_{pairs,BW}/n_{pairs, BW}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 5.0),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid size scan",
        # log_y=True,
    )



    #####################################
    ## GRID DENSITY
    #####################################




    graph_config = {
        "All varied": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_32_32_32_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_64_64_64_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_128_128_128_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_256_256_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
        },
        "Variation x": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_32_512_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_64_512_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_128_512_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_512_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
        },
        
        "Variation y": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_512_32_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_64_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_128_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_256_512_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
        },
        "Variation z": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_512_512_32_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_64_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_128_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_256_1/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_1/summary.json",
        },
    }

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_lumi_ee.png",
        json_key="luminosity.lumi_ee",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{ee}/L_{ee}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV", 
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_lumi_eg.png",
        json_key="luminosity.lumi_eg",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{eg}/L_{eg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_lumi_ge.png",
        json_key="luminosity.lumi_ge",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{ge}/L_{ge}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_lumi_gg.png",
        json_key="luminosity.lumi_gg",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{gg}/L_{gg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_upsmax.png",
        json_key="luminosity.upsmax",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="Y/Ymax",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_pairs_n_average.png",
        json_key="ipcs.pairs_n_average",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs}/n_{pairs}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_pairs_n_average_LL.png",
        json_key="ipcs.pairs_n_average_LL",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs,LL}/n_{pairs, LL}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_pairs_n_average_BH.png",
        json_key="ipcs.pairs_n_average_BH",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs,BH}/n_{pairs, BH}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt1_pairs_n_average_BW.png",
        json_key="ipcs.pairs_n_average_BW",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs,BW}/n_{pairs, BW}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid density scan, n_{t}=1",
        # log_y=True,
    )




    graph_config = {
        "All varied": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_32_32_32_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_64_64_64_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_128_128_128_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_256_256_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
        },
        "Variation x": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_32_512_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_64_512_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_128_512_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_512_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
        },
        
        "Variation y": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_512_32_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_64_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_128_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_256_512_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
        },
        "Variation z": {
            "x": [32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_512_512_32_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_64_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_128_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_256_4/summary.json",
                f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_512_512_512_4/summary.json",
        },
    }

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_lumi_ee.png",
        json_key="luminosity.lumi_ee",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{ee}/L_{ee}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_lumi_eg.png",
        json_key="luminosity.lumi_eg",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{eg}/L_{eg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_lumi_ge.png",
        json_key="luminosity.lumi_ge",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{ge}/L_{ge}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_lumi_gg.png",
        json_key="luminosity.lumi_gg",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="L_{gg}/L_{gg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_upsmax.png",
        json_key="luminosity.upsmax",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="Y/Ymax",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_pairs_n_average.png",
        json_key="ipcs.pairs_n_average",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs}/n_{pairs}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_pairs_n_average_LL.png",
        json_key="ipcs.pairs_n_average_LL",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs,LL}/n_{pairs, LL}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_pairs_n_average_BH.png",
        json_key="ipcs.pairs_n_average_BH",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs,BH}/n_{pairs, BH}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/density/grid_density_nt4_pairs_n_average_BW.png",
        json_key="ipcs.pairs_n_average_BW",
        x_title="n_{x}, n_{y}, n_{z}",
        y_title="n_{pairs,BW}/n_{pairs, BW}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 550),
        y_range=(0.4, 1.40),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",   
        legend_label="Grid density scan, n_{t}=4",
        # log_y=True,
    )


    #####################################
    ## TIME
    #####################################





    graph_config = {
        "32": {
            "x": [1, 2, 3, 4, 5],
            "files": [
                f"{base_input_dir}/CFG_NT_32_32_32_1/summary.json",
                f"{base_input_dir}/CFG_NT_32_32_32_2/summary.json",
                f"{base_input_dir}/CFG_NT_32_32_32_3/summary.json",
                f"{base_input_dir}/CFG_NT_32_32_32_4/summary.json",
                f"{base_input_dir}/CFG_NT_32_32_32_5/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NT_32_32_32_5/summary.json",
        },
        "64": {
            "x": [1, 2, 3, 4, 5],
            "files": [
                f"{base_input_dir}/CFG_NT_64_64_64_1/summary.json",
                f"{base_input_dir}/CFG_NT_64_64_64_2/summary.json",
                f"{base_input_dir}/CFG_NT_64_64_64_3/summary.json",
                f"{base_input_dir}/CFG_NT_64_64_64_4/summary.json",
                f"{base_input_dir}/CFG_NT_64_64_64_5/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NT_64_64_64_5/summary.json",
        },
        "128": {
            "x": [1, 2, 3, 4, 5],
            "files": [
                f"{base_input_dir}/CFG_NT_128_128_128_1/summary.json",
                f"{base_input_dir}/CFG_NT_128_128_128_2/summary.json",
                f"{base_input_dir}/CFG_NT_128_128_128_3/summary.json",
                f"{base_input_dir}/CFG_NT_128_128_128_4/summary.json",
                f"{base_input_dir}/CFG_NT_128_128_128_5/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NT_128_128_128_5/summary.json",
        },
        "256": {
            "x": [1, 2, 3, 4, 5],
            "files": [
                f"{base_input_dir}/CFG_NT_256_256_256_1/summary.json",
                f"{base_input_dir}/CFG_NT_256_256_256_2/summary.json",
                f"{base_input_dir}/CFG_NT_256_256_256_3/summary.json",
                f"{base_input_dir}/CFG_NT_256_256_256_4/summary.json",
                f"{base_input_dir}/CFG_NT_256_256_256_5/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NT_256_256_256_5/summary.json",
        },
        "512": {
            "x": [1, 2, 3, 4, 5],
            "files": [
                f"{base_input_dir}/CFG_NT_512_512_512_1/summary.json",
                f"{base_input_dir}/CFG_NT_512_512_512_2/summary.json",
                f"{base_input_dir}/CFG_NT_512_512_512_3/summary.json",
                f"{base_input_dir}/CFG_NT_512_512_512_4/summary.json",
                f"{base_input_dir}/CFG_NT_512_512_512_5/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NT_512_512_512_5/summary.json",
        },
    }

    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_lumi_ee.png",
        json_key="luminosity.lumi_ee",
        x_title="n_{t}",
        y_title="L_{ee}/L_{ee}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.8, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV", 
        legend_label="Grid time scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_lumi_eg.png",
        json_key="luminosity.lumi_eg",
        x_title="n_{t}",
        y_title="L_{eg}/L_{eg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.8, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_lumi_ge.png",
        json_key="luminosity.lumi_ge",
        x_title="n_{t}",
        y_title="L_{ge}/L_{ge}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.8, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_lumi_gg.png",
        json_key="luminosity.lumi_gg",
        x_title="n_{t}",
        y_title="L_{gg}/L_{gg}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.8, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_upsmax.png",
        json_key="luminosity.upsmax",
        x_title="n_{t}",
        y_title="Y/Ymax",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_pairs_n_average.png",
        json_key="ipcs.pairs_n_average",
        x_title="n_{t}",
        y_title="n_{pairs}/n_{pairs}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )

    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_pairs_n_average_LL.png",
        json_key="ipcs.pairs_n_average_LL",
        x_title="n_{t}",
        y_title="n_{pairs,LL}/n_{pairs, LL}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_pairs_n_average_BH.png",
        json_key="ipcs.pairs_n_average_BH",
        x_title="n_{t}",
        y_title="n_{pairs,BH}/n_{pairs, BH}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )
    make_graphs(
        graph_config,
        output=f"{output_dir}/time/time_pairs_n_average_BW.png",
        json_key="ipcs.pairs_n_average_BW",
        x_title="n_{t}",
        y_title="n_{pairs,BW}/n_{pairs, BW}^{max}",
        title="",
        draw_points=True,
        x_range=(0.0, 6),
        y_range=(0.2, 1.20),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",  
        legend_label="Grid time scan",
        # log_y=True,
    )




    ############################################


    quit()

    graph_config = {
        "ny=2": {
            "x": [16, 32, 64, 128, 256],
            "files": [
                f"{base_input_dir}/CFG_GRIDD_16_350_16/summary.json",
                f"{base_input_dir}/CFG_GRIDD_32_700_32/summary.json",
                f"{base_input_dir}/CFG_GRIDD_64_1400_64/summary.json",
                f"{base_input_dir}/CFG_GRIDD_128_2800_128/summary.json",
                f"{base_input_dir}/CFG_GRIDD_128_2800_128/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_GRIDD_128_2800_128/summary.json",
        },
        "ny=1": {
            "x": [16, 32, 64, 128, 256],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_16_175_16/summary.json",
                f"{base_input_dir}/CFG_NXYZT_32_350_32/summary.json",
                f"{base_input_dir}/CFG_NXYZT_64_700_64/summary.json",
                f"{base_input_dir}/CFG_NXYZT_128_1400_128/summary.json",
                f"{base_input_dir}/CFG_GRIDD_256_5600_256/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_GRIDD_256_5600_256/summary.json",
        },


    }

    make_graphs(
        graph_config,
        output=f"{output_dir}/grid_density_lumi_ee.png",
        json_key="luminosity.lumi_ee",
        x_title="Grid density",
        y_title="Luminosity",
        title="",
        draw_points=True,
        x_range=(0, 300),
        y_range=(0, 1.1),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",
        log_x=False,
    )


    graph_config = {
        "NO_LEGEND": {
            "x": [1000, 10000, 100000, 1000000],
            "files": [
                f"{base_input_dir}/CFG3_1000/summary.json",
                f"{base_input_dir}/CFG3_10000/summary.json",
                f"{base_input_dir}/CFG3_100000/summary.json",
                f"{base_input_dir}/CFG3_1000000/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG3_1000000/summary.json",
        },
    }

    make_graphs(
        graph_config,
        output=f"{output_dir}/nmacro_lumi_ee.png",
        json_key="luminosity.lumi_ee",
        x_title="Number of macro particles",
        y_title="Luminosity",
        title="",
        draw_points=True,
        x_range=(100, 10000000),
        y_range=(0.8, 1.05),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",
        log_x=True,
    )


    graph_config = {
        "NO_LEGEND": {
            "x": [16, 32, 64, 128, 256, 512],
            "files": [
                f"{base_input_dir}/CFG_NXYZT_256_16_256/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_32_256/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_64_256/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_128_256/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_256_256/summary.json",
                f"{base_input_dir}/CFG_NXYZT_256_512_256/summary.json",
            ],
            "norm": f"{base_input_dir}/CFG_NXYZT_256_512_256/summary.json",
        },
    }

    make_graphs(
        graph_config,
        output=f"{output_dir}/scan_y.png",
        json_key="luminosity.lumi_gg",
        x_title="Number of macro particles",
        y_title="Luminosity",
        title="",
        draw_points=True,
        x_range=(0, 600),
        y_range=(0.0, 1.05),
        extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
        extra_text_right="GHC V25.1 (FSR) 91.2 GeV",
    )


        







