import os
import ROOT
import plotter
import config as gpconfig
import argparse
import numpy as np
import math


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

outDir = "/home/submit/jaeyserm/public_html/fccee/guineapig/visualization/fields/"

def extract(GRID_ID=1):
    INFILE   = "/work/submit/jaeyserm/fccee/FCCAnalyses/FCCPhysics/beam_backgrounds/guineapig/guinea-pig-15122025_dev_time/fields.txt"
    #GRID_ID = 1             # gridsPtr_ index, == "GRID8" in the GP banner
    FIELD   = "total"       # "total", "beam1" or "beam2"

    # ---- geometry from the FIELDGRID header of this grid
    geo = None
    rows = []
    with open(INFILE) as f:
        for line in f:
            if line.startswith("FIELDGRID "):
                p = line.split()
                if int(p[1]) == GRID_ID:
                    geo = {p[k]: float(p[k+1]) for k in range(2, len(p)-1, 2)}
            elif line.startswith("FIELD "):
                p = line.split()
                if int(p[1]) == GRID_ID:
                    rows.append(p[4:12])          # i j x y E1x E1y E2x E2y

    if geo is None:
        raise RuntimeError(f"no FIELDGRID header for grid {GRID_ID}")
    print("geometry:", geo)

    d  = np.array(rows, dtype=float)
    print(f"read {len(d)} entries for grid {GRID_ID}")

    i, j = d[:, 0].astype(int), d[:, 1].astype(int)
    x, y = d[:, 2], d[:, 3]

    if FIELD == "beam1":
        ex, ey = d[:, 4], d[:, 5]
    elif FIELD == "beam2":
        ex, ey = d[:, 6], d[:, 7]
    else:
        ex, ey = d[:, 4] + d[:, 6], d[:, 5] + d[:, 7]
    mag = np.hypot(ex, ey) * 2.0                  # field a particle feels (E+B)

    # ---- envelope: max over all slice pairs, per cell
    i0, j0 = int(i.min()), int(j.min())
    ni, nj = int(i.max()) - i0 + 1, int(j.max()) - j0 + 1
    env = np.zeros((ni, nj))
    np.maximum.at(env, (i - i0, j - j0), mag)

    def pick_unit(vmax):
        """Return (scale, label) so that vmax/scale lands in a readable range."""
        v = abs(vmax)
        if v < 1e3:  return 1.0, "nm"      # < 1 um
        if v < 1e6:  return 1e3, "#mum"    # < 1 mm
        return 1e6, "mm"

    dx, dy = geo["delta_x"], geo["delta_y"]
    xlo, xhi = float(x.min() - 0.5*dx), float(x.max() + 0.5*dx)
    ylo, yhi = float(y.min() - 0.5*dy), float(y.max() + 0.5*dy)

    sx, ux = pick_unit(max(abs(xlo), abs(xhi)))
    sy, uy = pick_unit(max(abs(ylo), abs(yhi)))

    # field: GV/nm -> GV/m  (1 GV/nm = 1e9 GV/m)
    env_plot = env * 1e9

    h = ROOT.TH2D("field_env",
                f";x ({ux});y ({uy});max |E| over slices (GV/m)",
                ni, xlo/sx, xhi/sx,
                nj, ylo/sy, yhi/sy)
    for a in range(ni):
        for b in range(nj):
            h.SetBinContent(a + 1, b + 1, env_plot[a, b])

    c = ROOT.TCanvas("c", "", 900, 800)
    c.SetRightMargin(0.16)
    c.SetLogz()
    h.Draw("COLZ")
    c.SaveAs(f"{outDir}/field_env_grid{GRID_ID}.png")

    print(f"cells {ni} x {nj}   x: {xlo:.4g}..{xhi:.4g} nm   y: {ylo:.4g}..{yhi:.4g} nm")
    print(f"expected half-extent from header: cut_x={geo['cut_x']:.4g}  cut_y={geo['cut_y']:.4g} nm")
    print(f"max |E| = {env.max():.3e}")
    print(f"min |E| = {env.min():.3e}")
    print(f"Ratio max/min = {env.max()/env.min():.3e}")

    fout = ROOT.TFile(f"{outDir}/field_env_grid{GRID_ID}.root", "RECREATE")
    h.Write()
    fout.Close()

    print(f"cells {ni} x {nj}, max |E| = {env.max():.3e}")


    # --- envelope statistics: global max vs max on the outer ring of dumped cells
    b = np.zeros_like(env, dtype=bool)
    b[0, :] = True; b[-1, :] = True          # x-low / x-high edges
    b[:, 0] = True; b[:, -1] = True          # y-low / y-high edges

    xs = np.zeros(ni); xs[i - i0] = x
    ys = np.zeros(nj); ys[j - j0] = y

    max_global = env.max()
    a_max, b_max = np.unravel_index(np.argmax(env), env.shape)
    x_max, y_max = xs[a_max], ys[b_max]

    edges = {
        "x-low ": env[0,  :].max(),
        "x-high": env[-1, :].max(),
        "y-low ": env[:,  0].max(),
        "y-high": env[:, -1].max(),
    }
    max_bnd = env[b].max()

    print(f"global max      : {max_global*1e9:.4g} GV/m  at x={x_max/sx:.4g} {ux}, y={y_max/sy:.4g} {uy}")
    for k, v in edges.items():
        print(f"  max on {k}  : {v*1e9:.4g} GV/m   (max/edge = {max_global/v:.3e})")
    print(f"max on boundary : {max_bnd*1e9:.4g} GV/m")
    print(f"RATIO max/boundary = {max_global/max_bnd:.4e}")


def plot(GRID_ID):

    label = "GHC V25.1 (FSR) 91.2 GeV"
    label = f"Grid {GRID_ID+1}"
    fIn = ROOT.TFile(f"{outDir}/field_env_grid{GRID_ID}.root")
    h = fIn.Get("field_env")

    max_global = h.GetMaximum()
    h.Scale(1./max_global)

    plotter.plot_hist_2d(
            hist=h,
            outname=f"{outDir}/field_env_grid{GRID_ID}",
            x_title=h.GetXaxis().GetTitle(),
            y_title=h.GetYaxis().GetTitle(),
            z_title=h.GetZaxis().GetTitle(),
            x_range=None,
            y_range=None,
            z_range=(0, 1),
            canvas_size=(800, 800),
            draw_option="COLZ",
            extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
            extra_text_right=label, 
            logz=False,
    )

    nx = h.GetNbinsX()
    ny = h.GetNbinsY()


    max_global = h.GetMaximum()

    edges = {
        "x-low":  max(h.GetBinContent(1,  iy) for iy in range(1, ny+1)),
        "x-high": max(h.GetBinContent(nx, iy) for iy in range(1, ny+1)),
        "y-low":  max(h.GetBinContent(ix, 1)  for ix in range(1, nx+1)),
        "y-high": max(h.GetBinContent(ix, ny) for ix in range(1, nx+1)),
    }

    max_bnd = max(edges.values())

    print(f"Global max: {max_global:.4g}")
    for k, v in edges.items():
        print(f"  {k}: {v:.4g} (max/edge = {max_global/v:.3e})")

    print(f"Boundary max: {max_bnd:.4g}")
    print(f"Max/boundary: {max_global/max_bnd:.4e}")

    fIn.Close()

    return max_global, max_bnd


def main():

    #extract(1)

    grids = list(range(0, 8))
    max_vals = []
    max_global_ref = -1
    for ig in grids:
        max_global, max_bnd = plot(ig)
        max_vals.append(max_global/max_bnd)

    x_arr = np.array(list(range(1, 9)), dtype=np.float64) # grids starting from 1
    y_arr = np.array(max_vals, dtype=np.float64)
    graph = ROOT.TGraph(len(grids), x_arr, y_arr)


    # Assuming graph is your TGraph
    xmin = min(graph.GetPointX(i) for i in range(graph.GetN()))
    xmax = max(graph.GetPointX(i) for i in range(graph.GetN()))

    ymin = min(graph.GetPointY(i) for i in range(graph.GetN()))
    ymax = max(graph.GetPointY(i) for i in range(graph.GetN()))

    # Sigmoid function
    sigmoid = ROOT.TF1(
        "sigmoid",
        "[0] + ([1]-[0])/(1+exp(-(x-[2])/[3]))",
        xmin, xmax
    )

    sigmoid.SetParNames("y_min", "y_max", "x0", "width")
    sigmoid.SetLineColor(ROOT.kBlack)

    # Initial parameter estimates
    sigmoid.SetParameters(
        ymin,
        ymax,
        (xmin + xmax) / 2,
        (xmax - xmin) / 10
    )

    # Fit the graph
    fit_result = graph.Fit(sigmoid, "RS")

    label = "GHC V25.1 (FSR) 91.2 GeV"
    plotter.plot_graphs_1d(
            graphs=[sigmoid, graph],
            labels=["Sigmoid fit", "Max E / Max E boundary", ],
            outname=f"{outDir}/grid_max_field_convergence",
            x_title="Grid number",
            y_title="Max E / Max E boundary",
            x_range=None,
            y_range=None,
            canvas_size=(800, 800),
            draw_option=["L", "P"],
            extra_text_left="#bf{FCC-ee}#scale[0.7]{#it{ GuineaPig Simulation}}", 
            extra_text_right=label, 
            legend_pos=(0.30, 0.75, 0.75, 0.88)
    )
    

if __name__ == "__main__":
    main()