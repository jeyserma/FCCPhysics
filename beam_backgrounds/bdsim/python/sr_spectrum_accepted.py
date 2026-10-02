#!/usr/bin/env python3
"""
sr_spectrum_accepted.py
=======================
Energy spectrum of the ACCEPTED SR photons, over a whole production directory of
EDM4hep files.

    python python/sr_spectrum_accepted.py \\
        /ceph/submit/data/group/fcc/ee/beam_backgrounds/bdsim/LCC_v2short_halo_cvmfs_v1.7.7

The companion `sr_spectrum.py` needs the raw BDSIM output to draw the
before-cuts curve, so it only works on a local or `--dryrun` job. Production
directories hold just the skimmed `*_edm4hep.root`, which is what this reads.

Normalisation comes from the files themselves: the per-event `nPrimaries` are
summed across everything read, and

    weight = chargeFraction * bunchIntensity / sum(nPrimaries)

so the y axis is photons per bunch crossing (one beam) per 0.1% bandwidth,
exactly as in sr_spectrum.py. Reading a subset of the directory therefore still
gives the right normalisation -- the weight scales with what was actually read.

Needs key4hep (podio + ROOT).
"""

import argparse
import array
import glob
import math
import os
import sys

import ROOT

from sample_name import resolve_sample

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError

E_MIN, E_MAX, NBINS = 1e-9, 1e-2, 60      # GeV
COL = "#2a78d6"                            # blue, matching sr_spectrum.py
COL_MUTED = "#6b6b6b"
COL_FAINT = "#9a9a9a"


def sci(v):
    if v <= 0:
        return "0"
    exp = int(math.floor(math.log10(v)))
    return f"{v / 10 ** exp:.2f}#times10^{{{exp}}}"


def log_bins():
    lo, hi = math.log10(E_MIN), math.log10(E_MAX)
    return array.array("d", [10 ** (lo + (hi - lo) * i / NBINS) for i in range(NBINS + 1)])


def read_directory(files, max_files=None):
    """
    Fill the energy histogram from the EDM4hep files and sum their primaries.

    The energies are read with a TChain rather than the podio python API: for a
    2000-file production that is the difference between about a minute and most
    of an hour, because the per-particle python loop dominates everything else.
    The primary counts still come through podio, but that is one small read per
    file rather than 270000.
    """
    from podio import root_io

    if max_files:
        files = files[:max_files]

    n_primaries, n_files, runconfig = 0, 0, None
    charge_fraction, bunch_intensity = None, None
    for path in files:
        try:
            reader = root_io.Reader(path)
        except Exception:                        # noqa: BLE001 -- unreadable file
            print(f"  skipping unreadable {os.path.basename(path)}", file=sys.stderr)
            continue
        for event in reader.get("events"):
            try:
                n_primaries += int(event.get_parameter("nPrimaries"))
            except Exception:                    # noqa: BLE001
                pass
        if runconfig is None and "metadata" in reader.categories:
            meta = reader.get("metadata")[0]
            keys = list(meta.parameters)
            if "runConfig" in keys:
                runconfig = meta.get_parameter("runConfig")
            if "chargeFraction" in keys:
                charge_fraction = float(meta.get_parameter("chargeFraction"))
            if "bunchIntensity" in keys:
                bunch_intensity = float(meta.get_parameter("bunchIntensity"))
        n_files += 1

    edges = log_bins()
    # NB: no SetDirectory(0) before the Draw. TTree::Draw resolves ">>h_acc" by
    # name in gDirectory, so a detached histogram is invisible to it and ROOT
    # silently fills a fresh one instead -- leaving this one empty.
    h = ROOT.TH1D("h_acc", "", NBINS, edges)
    chain = ROOT.TChain("events")
    for path in files:
        chain.Add(path)
    expr = ("sqrt(MCParticles.momentum.x*MCParticles.momentum.x"
            "+MCParticles.momentum.y*MCParticles.momentum.y"
            "+MCParticles.momentum.z*MCParticles.momentum.z)>>h_acc")
    n = chain.Draw(expr, "", "goff")
    if n <= 0:
        sys.exit("TChain::Draw returned no entries -- is MCParticles present "
                 "in these files?")
    h = ROOT.gDirectory.Get("h_acc")
    h.SetDirectory(0)

    return h, n_primaries, n_files, runconfig, charge_fraction, bunch_intensity


def main():
    ap = argparse.ArgumentParser(
        description="Accepted SR photon spectrum over a production directory.")
    ap.add_argument("input", help="directory of *_edm4hep.root files (or a glob)")
    ap.add_argument("-o", "--output",
                    default="/home/submit/jaeyserm/public_html/fccee/bdsim//",
                    help="output base directory")
    ap.add_argument("--max-files", type=int, default=100,
                    help="read only the first N files (the normalisation still "
                         "comes out right -- the weight scales with what was read)")
    args = ap.parse_args()


    files = sorted(glob.glob(os.path.join(args.input, "*_edm4hep.root")))
    # tag = full sample name, e.g. LCC_v2short_halo_standalone_G4SR, so samples
    # from different BDSIM builds get separate plot directories
    sample = resolve_sample(args.input)
    tag = sample[3] if sample else os.path.basename(os.path.normpath(args.input))

    if not files:
        sys.exit(f"no *_edm4hep.root found in {args.input}")
    print(f"found {len(files)} files in {args.input}")

    h, n_primaries, n_files, runconfig, cf, bi = read_directory(files, args.max_files)
    if n_primaries <= 0:
        sys.exit("no nPrimaries found in these files; cannot normalise")
    if cf is None or bi is None:
        sys.exit("metadata has no chargeFraction / bunchIntensity; cannot normalise")

    rc = runconfig or (sample[1] if sample else "")

    n_photons = int(h.GetEntries())
    weight = cf * bi / n_primaries

    # per 0.1% bandwidth: dN/dlnE * 1e-3
    ymax = 0.0
    for b in range(1, NBINS + 1):
        dlnE = math.log(h.GetBinLowEdge(b + 1)) - math.log(h.GetBinLowEdge(b))
        h.SetBinContent(b, h.GetBinContent(b) * weight / dlnE * 1e-3)
        h.SetBinError(b, h.GetBinError(b) * weight / dlnE * 1e-3)
        ymax = max(ymax, h.GetBinContent(b))

    col = ROOT.TColor.GetColor(COL)
    h.SetLineColor(col)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(col, 0.18)
    filled = [b for b in range(1, NBINS + 1) if h.GetBinContent(b) > 0]
    if filled:
        h.GetXaxis().SetRange(filled[0], filled[-1])

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetGridColor(ROOT.TColor.GetColor("#dcdcdc"))
    ROOT.gStyle.SetGridStyle(1)

    c = ROOT.TCanvas("c", "", 1100, 760)
    c.SetLogx(); c.SetLogy(); c.SetGridx(); c.SetGridy()
    c.SetLeftMargin(0.115); c.SetRightMargin(0.035)
    c.SetTopMargin(0.150); c.SetBottomMargin(0.125)

    frame = c.DrawFrame(E_MIN, ymax * 1e-5, E_MAX, ymax * 60)
    frame.GetXaxis().SetTitle("photon energy (GeV)")
    frame.GetYaxis().SetTitle("photons per bunch crossing per 0.1% BW")
    frame.GetXaxis().CenterTitle(True); frame.GetYaxis().CenterTitle(True)
    frame.GetXaxis().SetTitleOffset(1.25); frame.GetYaxis().SetTitleOffset(1.45)
    frame.GetXaxis().SetLabelSize(0.033); frame.GetYaxis().SetLabelSize(0.033)
    frame.GetXaxis().SetTitleSize(0.038); frame.GetYaxis().SetTitleSize(0.038)

    h.Draw("hist same")
    h.Draw("hist same ][")

    lat = ROOT.TLatex(); lat.SetNDC(True); lat.SetTextFont(42)
    lat.SetTextAlign(31); lat.SetTextSize(0.030); lat.SetTextColor(col)
    lat.DrawLatex(0.950, 0.795, f"{sci(n_photons * weight)} photons / bunch crossing")
    lat.SetTextSize(0.027); lat.SetTextColor(ROOT.TColor.GetColor(COL_MUTED))
    lat.DrawLatex(0.950, 0.750, f"{n_photons:,} simulated, weight {weight:.4g}")
    lat.DrawLatex(0.950, 0.707, f"{n_files} files, {n_primaries:,} primaries")

    lat.SetTextAlign(11); lat.SetTextColor(ROOT.kBlack)
    lat.SetTextFont(62); lat.SetTextSize(0.038)
    lat.SetTextAlign(21)          # centred on the frame, not left-aligned
    lat.DrawLatex(0.115 + (1 - 0.115 - 0.035) / 2, 0.955,
                  f"Accepted SR photons #minus {tag} {rc}")

    top = ROOT.TGaxis(E_MIN, ymax * 60, E_MAX, ymax * 60,
                      E_MIN * 1e6, E_MAX * 1e6, 510, "-G")
    top.SetLabelFont(42); top.SetLabelSize(0.030)
    top.Draw()
    lat.SetTextAlign(21); lat.SetTextFont(42); lat.SetTextSize(0.033)
    lat.DrawLatex(0.115 + (1 - 0.115 - 0.035) / 2, 0.905, "photon energy (keV)")

    c.RedrawAxis()

    out = f"{args.output}/{tag}"
    os.makedirs(out, exist_ok=True)
    os.system(f"cp {args.output}/index.php {out}/")
    c.SaveAs(f"{out}/sr_spectrum_acceptance.png")
    c.SaveAs(f"{out}/sr_spectrum_acceptance.pdf")

    print(f"runconfig  : {runconfig}")
    print(f"files      : {n_files}")
    print(f"primaries  : {n_primaries:,}")
    print(f"photons    : {n_photons:,}")
    print(f"weight     : {weight:.6g}")
    print(f"yield      : {n_photons / n_primaries:.5f} photons/primary")
    print(f"total      : {n_photons * weight:.4g} photons per bunch crossing")
    print(f"wrote      : {out}")


if __name__ == "__main__":
    main()
