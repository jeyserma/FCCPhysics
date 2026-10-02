#!/usr/bin/env python3
"""
sr_spectrum.py
==============
SR photon energy spectrum at the BDSIM sampler, before and after the acceptance
skim that convert_edm4hep.py applies -- both curves on one set of axes.

Made for a local run or a `submit_new.py --dryrun`, which leaves the raw BDSIM
output in  <storagedir>/<sample name>/tmp/output_12345.root  (the older
/tmp/bdsim/<lattice>/<runconfig>/ layout from submit_old.py works too).  That raw file is never
transferred back from the grid, so this is the only place the *pre-skim* photons
exist: the plot cannot be made from a production sample.

Usage
-----
    # point it at the dryrun directory and let it find the file
    python sr_spectrum.py /ceph/submit/data/group/fcc/ee/beam_backgrounds/bdsim/LCC_v2short_halo_standalone_G4SR/tmp

    # or give the file directly
    python sr_spectrum.py /ceph/submit/data/group/fcc/ee/beam_backgrounds/bdsim/LCC_v2short_halo_standalone_G4SR/tmp/output_12345.root -o halo.png

    # the run config sets the charge fraction; inferred from the path if it can be
    python sr_spectrum.py <file> --runconfig core

Environment
-----------
Needs the key4hep stack, NOT the BDSIM stack:

    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08

Like convert_edm4hep.py, it deliberately avoids pybdsim: the sampler is an
unsplit BDSOutputROOTEventSampler<float>, but ROOT resolves it through the
StreamerInfo in the file, so no BDSIM dictionaries are needed. TTree::Draw is
used for the same reason.

Normalisation
-------------
Each simulated primary stands for  chargeFraction * bunchIntensity / nPrimaries
real positrons, so the y axis is photons per bunch crossing (one beam) per 0.1%
relative bandwidth,  dN/dlnE * 1e-3.  Integrating the curve over dlnE/1e-3
returns the quoted total, which is the closure check in the printout.
"""

import argparse
import array
import glob
import math
import os
import sys

import ROOT

from sample_name import resolve_sample, tried_message

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError      # silence the expected "no dictionary" notices

# ---------------------------------------------------------------- constants --
SAMPLER = "QD0AL"                 # sampler element, 2.4 m upstream of the IP


E_MIN, E_MAX, NBINS = 3e-9, 1e-2, 60

BUNCH_INTENSITY = 2.02e11         # fallback if sample_config cannot be imported
CHARGE_FRACTION = {"halo": 0.01, "core": 0.99}

COL_ALL = "#eb6834"               # orange -- every photon at the sampler
COL_CUT = "#2a78d6"               # blue   -- survives the acceptance cuts
COL_MUTED = "#6b6b6b"


def load_sample_config(config_dir):
    """
    Pick up BUNCH_INTENSITY / CHARGE_FRACTION from the sample's own
    sample_config.py so this script cannot drift from the production values.
    Falls back to the constants above if it is not importable.
    """
    global BUNCH_INTENSITY, CHARGE_FRACTION
    if not config_dir or not os.path.isdir(config_dir):
        return False
    sys.path.insert(0, config_dir)
    try:
        import sample_config
        BUNCH_INTENSITY = sample_config.BUNCH_INTENSITY
        CHARGE_FRACTION = dict(sample_config.CHARGE_FRACTION)
        return True
    except Exception:                       # noqa: BLE001 -- optional convenience
        return False
    finally:
        sys.path.pop(0)


def find_input(path):
    """
    Accept a file, or a directory holding one raw BDSIM output.

    Returns (raw, edm4hep) -- the EDM4hep sibling is None if it is not there.
    """
    if os.path.isfile(path):
        raw = path
    else:
        if not os.path.isdir(path):
            sys.exit(f"no such file or directory: {path}")
        hits = [f for f in sorted(glob.glob(os.path.join(path, "output_*.root")))
                if not f.endswith("_edm4hep.root")]
        if not hits:
            sys.exit(f"no raw BDSIM output_*.root in {path} "
                     "(for a directory of skimmed files use sr_spectrum_accepted.py)")
        if len(hits) > 1:
            print(f"note: {len(hits)} raw files found, using {os.path.basename(hits[0])}")
        raw = hits[0]

    conv = raw.replace(".root", "_edm4hep.root")
    return raw, (conv if os.path.exists(conv) else None)


def read_accepted_edm4hep(path):
    """
    Photon energies (GeV) and the primary count from the converted file.

    Taking the accepted sample from the EDM4hep output rather than re-applying
    the skim here means the blue curve is exactly what the production contains.
    Reimplementing the cuts would be a second copy of convert_edm4hep.py's logic,
    free to drift from it silently.
    """
    from podio import root_io

    reader = root_io.Reader(path)
    energies = []
    n_primaries = 0
    e_cut_min = None
    for event in reader.get("events"):
        for p in event.get("MCParticles"):
            m = p.getMomentum()
            energies.append((m.x * m.x + m.y * m.y + m.z * m.z) ** 0.5)
        try:
            n_primaries += int(event.get_parameter("nPrimaries"))
        except Exception:                      # noqa: BLE001 -- older files
            pass
    if "metadata" in reader.categories:
        meta = reader.get("metadata")[0]
        keys = list(meta.parameters)
        if n_primaries == 0 and "nPrimariesTotal" in keys:
            n_primaries = int(meta.get_parameter("nPrimariesTotal"))
        # The threshold is only used to annotate the plot, but take it from the
        # file rather than hard-coding it: then the line always marks the cut
        # that was actually applied.
        if "cutEnergyMin_GeV" in keys:
            e_cut_min = float(meta.get_parameter("cutEnergyMin_GeV"))
    return energies, n_primaries, e_cut_min


def draw_columns(tree, columns, n_estimate):
    """Read up to 4 branches with TTree::Draw into python lists of floats."""
    tree.SetEstimate(n_estimate)
    n = tree.Draw(":".join(columns), "", "goff")
    if n < 0:
        raise RuntimeError(f"TTree::Draw failed for {columns}")
    out = []
    for i in range(len(columns)):
        buf = tree.GetVal(i)
        out.append([buf[j] for j in range(n)])
    return out, n


def read_all_photons(rootfile):
    """
    Every photon at the sampler, before any cut: (energies in GeV, n_primaries).

    Only energy and partID are read. The accepted subset is NOT recomputed here
    -- it comes from the converted file, so there is exactly one implementation
    of the skim (convert_edm4hep.py) and no second copy to drift from it.
    """
    f = ROOT.TFile.Open(rootfile)
    if not f or f.IsZombie():
        sys.exit(f"cannot open {rootfile}")
    tree = f.Get("Event")
    if not tree:
        sys.exit(f"no Event tree in {rootfile} -- is this a raw BDSIM output?")

    n_primaries = int(tree.GetEntries())
    tree.SetEstimate(n_primaries + 1)
    tree.Draw(f"{SAMPLER}.n", "", "goff")
    nbuf = tree.GetVal(0)
    n_hits = int(sum(nbuf[i] for i in range(n_primaries)))
    if n_hits == 0:
        return [], n_primaries

    (energy, partID), _ = draw_columns(
        tree, [f"{SAMPLER}.energy", f"{SAMPLER}.partID"], n_hits)
    e_all = [energy[i] for i in range(n_hits) if partID[i] == 22]
    f.Close()
    return e_all, n_primaries


def log_bins():
    lo, hi = math.log10(E_MIN), math.log10(E_MAX)
    return array.array("d", [10 ** (lo + (hi - lo) * i / NBINS) for i in range(NBINS + 1)])


def make_hist(name, energies, weight, edges):
    """Fill, then convert each bin to dN/dlnE * 1e-3 = photons per 0.1% BW."""
    h = ROOT.TH1D(name, "", NBINS, edges)
    h.SetDirectory(0)
    for e in energies:
        h.Fill(e, weight)
    for b in range(1, NBINS + 1):
        dlnE = math.log(h.GetBinLowEdge(b + 1)) - math.log(h.GetBinLowEdge(b))
        h.SetBinContent(b, h.GetBinContent(b) / dlnE * 1e-3)
        h.SetBinError(b, h.GetBinError(b) / dlnE * 1e-3)
    return h


def style_hist(h, hex_colour):
    col = ROOT.TColor.GetColor(hex_colour)
    h.SetLineColor(col)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(col, 0.18)
    return col


def sci(v):
    if v <= 0:
        return "0"
    exp = int(math.floor(math.log10(v)))
    return f"{v / 10 ** exp:.1f}#times10^{{{exp}}}"


def median(values):
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def main():
    global SAMPLER
    ap = argparse.ArgumentParser(
        description="SR photon spectrum before/after the acceptance skim, "
                    "from a local or --dryrun BDSIM output.")
    ap.add_argument("input", help="raw BDSIM output_*.root, or the directory holding it")
    ap.add_argument("-o", "--output", default="/home/submit/jaeyserm/public_html/fccee/bdsim//",
                    help="output directory")
    ap.add_argument("--config-dir", default=None,
                    help="the sample's config dir, to import sample_config.py "
                         "(default: inferred from the path)")
    ap.add_argument("--sampler", default=SAMPLER, help=f"sampler name (default {SAMPLER})")
    args = ap.parse_args()

    SAMPLER = args.sampler

    infile, convfile = find_input(args.input)
    sample = resolve_sample(infile, args.config_dir)
    if sample is None:
        sys.exit(tried_message(infile) + "\n  or pass --config-dir")
    lattice, runconfig, cfgdir, tag = sample
    got_cfg = load_sample_config(cfgdir)

    if runconfig not in CHARGE_FRACTION:
        sys.exit(f"runconfig '{runconfig}' from {infile} has no charge fraction; "
                 f"known: {sorted(CHARGE_FRACTION)}.")

    print(f"input      : {infile}")
    print(f"sample     : {tag}  (lattice {lattice})")
    print(f"runconfig  : {runconfig}  (chargeFraction {CHARGE_FRACTION[runconfig]}, "
          f"{'from sample_config.py' if got_cfg else 'built-in fallback'})")

    if convfile is None:
        sys.exit(f"no _edm4hep.root alongside {os.path.basename(infile)}.\n"
                 "The accepted photons are read from the converted file, not "
                 "recomputed here -- run convert.sh (or convert_edm4hep.py) first.")

    e_all, n_primaries = read_all_photons(infile)
    if not e_all:
        sys.exit("no photons at the sampler")

    e_cut, n_prim_conv, e_cut_min = read_accepted_edm4hep(convfile)
    if e_cut_min is None:
        e_cut_min = 2e-6
        print("  note: metadata has no cutEnergyMin_GeV; annotating at 2 keV")
    print(f"accepted   : {os.path.basename(convfile)}  ({len(e_cut)} photons)")
    if n_prim_conv and n_prim_conv != n_primaries:
        print(f"  note: EDM4hep says {n_prim_conv} primaries, raw tree has "
              f"{n_primaries}; using the EDM4hep value for the weight")
        n_primaries = n_prim_conv

    cutflow = [("all photons at sampler", len(e_all)),
               ("accepted (from EDM4hep)", len(e_cut))]

    weight = CHARGE_FRACTION[runconfig] * BUNCH_INTENSITY / n_primaries
    tot_all, tot_cut = len(e_all) * weight, len(e_cut) * weight
    acc = 100.0 * len(e_cut) / len(e_all)

    edges = log_bins()
    h_all = make_hist("h_all", e_all, weight, edges)
    h_cut = make_hist("h_cut", e_cut, weight, edges)
    col_all = style_hist(h_all, COL_ALL)
    col_cut = style_hist(h_cut, COL_CUT)

    # ------------------------------------------------------------- drawing --
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickX(0)
    ROOT.gStyle.SetPadTickY(1)

    c = ROOT.TCanvas("c", "", 1100, 760)
    c.SetLogx()
    c.SetLogy()
    c.SetLeftMargin(0.115)
    c.SetRightMargin(0.035)
    c.SetTopMargin(0.150)
    c.SetBottomMargin(0.125)
    c.SetFillColor(ROOT.kWhite)
    c.SetGridx()
    c.SetGridy()
    ROOT.gStyle.SetGridColor(ROOT.TColor.GetColor("#dcdcdc"))
    ROOT.gStyle.SetGridStyle(1)

    ymax = max(h_all.GetMaximum(), h_cut.GetMaximum())
    # Outside the filled span ROOT would draw the empty bins as a flat line along
    # the frame floor of a log axis, so restrict each histogram to its own range.
    for h in (h_all, h_cut):
        filled = [b for b in range(1, NBINS + 1) if h.GetBinContent(b) > 0]
        if filled:
            h.GetXaxis().SetRange(filled[0], filled[-1])
        for b in filled[0:1] and range(filled[0], filled[-1] + 1):
            if h.GetBinContent(b) <= 0:
                h.SetBinContent(b, ymax * 1e-12)   # interior gap, clipped away

    frame = c.DrawFrame(E_MIN, ymax * 1e-4, E_MAX, ymax * 60)
    frame.GetXaxis().SetTitle("photon energy (GeV)")
    frame.GetYaxis().SetTitle("photons per bunch crossing per 0.1% BW")
    frame.GetXaxis().CenterTitle(True)
    frame.GetYaxis().CenterTitle(True)
    frame.GetXaxis().SetTitleOffset(1.25)
    frame.GetYaxis().SetTitleOffset(1.45)
    frame.GetXaxis().SetLabelSize(0.033)
    frame.GetYaxis().SetLabelSize(0.033)
    frame.GetXaxis().SetTitleSize(0.038)
    frame.GetYaxis().SetTitleSize(0.038)

    h_all.Draw("hist same")
    h_cut.Draw("hist same")
    h_all.Draw("hist same ][")          # redraw outlines over the fills
    h_cut.Draw("hist same ][")

    ylo, yhi = ymax * 1e-4, ymax * 60

    def vline(xval, hex_colour, style):
        ln = ROOT.TLine(xval, ylo, xval, yhi)
        ln.SetLineColor(ROOT.TColor.GetColor(hex_colour))
        ln.SetLineStyle(style)
        ln.SetLineWidth(2)
        ln.Draw()
        return ln

    # No threshold markers: both edges are already where the curves visibly
    # start, so the lines and their labels only added clutter.

    lat = ROOT.TLatex()
    lat.SetNDC(False)
    lat.SetTextFont(42)

    # direct labels, so identity is not carried by colour alone
    lat.SetTextColor(col_all)
    lat.SetTextSize(0.034)
    lat.SetTextAlign(11)
    b = h_all.FindBin(1.4e-7)
    lat.DrawLatex(1.6e-7, h_all.GetBinContent(b) * 2.4, "all photons at the sampler")

    lat.SetTextColor(col_cut)
    lat.SetTextAlign(31)
    b = h_cut.FindBin(2.0e-5)
    lat.DrawLatex(1.4e-4, h_cut.GetBinContent(b) * 0.34, "after acceptance cuts")

    # numbers, upper right where both curves have fallen away
    lat.SetNDC(True)
    lat.SetTextAlign(31)
    lat.SetTextSize(0.030)
    lat.SetTextColor(col_all)
    lat.DrawLatex(0.950, 0.795,
                  f"{sci(tot_all)} photons / bunch crossing    median "
                  f"{median(e_all) * 1e6:.1f} keV")
    lat.SetTextColor(col_cut)
    lat.DrawLatex(0.950, 0.745,
                  f"{sci(tot_cut)} photons / bunch crossing    median "
                  f"{median(e_cut) * 1e6:.1f} keV")
    lat.SetTextColor(ROOT.TColor.GetColor(COL_MUTED))
    lat.SetTextSize(0.027)
    lat.DrawLatex(0.950, 0.700, f"acceptance {acc:.1f} %")

    # title and one provenance line
    lat.SetTextAlign(11)
    lat.SetTextColor(ROOT.kBlack)
    lat.SetTextFont(62)
    lat.SetTextSize(0.038)
    lat.SetTextAlign(21)          # centred on the frame, not left-aligned
    lat.DrawLatex(0.115 + (1 - 0.115 - 0.035) / 2, 0.955,
                  f"SR photon spectrum at the sampler #minus {tag}")
    lat.SetTextFont(42)

    # keV scale on top, the units the audience thinks in
    top = ROOT.TGaxis(E_MIN, yhi, E_MAX, yhi, E_MIN * 1e6, E_MAX * 1e6, 510, "-G")
    top.SetLabelFont(42)
    top.SetLabelSize(0.030)
    top.Draw()

    lat.SetNDC(True)
    lat.SetTextAlign(21)
    lat.SetTextFont(42)
    lat.SetTextSize(0.033)
    lat.SetTextColor(ROOT.kBlack)
    lat.DrawLatex(0.115 + (1 - 0.115 - 0.035) / 2, 0.905, "photon energy (keV)")

    c.RedrawAxis()

    out = f"{args.output}/{tag}"
    os.makedirs(out, exist_ok=True)
    os.system(f"cp {args.output}/index.php {out}/")
    c.SaveAs(f"{out}/sr_spectrum.png")
    c.SaveAs(f"{out}/sr_spectrum.pdf")

    # ------------------------------------------------------------ printout --
    print(f"primaries  : {n_primaries}")
    print(f"weight     : {weight:.6g} photons/BX per simulated photon")
    print("cut flow:")
    for label, count in cutflow:
        print(f"  {label:26s} {count:9d}")
    print(f"acceptance : {acc:.2f} %")
    print(f"totals /BX : all {tot_all:.4g}   after cuts {tot_cut:.4g}")
    print(f"median keV : all {median(e_all) * 1e6:.2f}   after cuts {median(e_cut) * 1e6:.2f}")

    # closure: integrating the plotted curve must return the quoted total
    integral = sum(h_cut.GetBinContent(b) *
                   (math.log(h_cut.GetBinLowEdge(b + 1)) - math.log(h_cut.GetBinLowEdge(b)))
                   / 1e-3 for b in range(1, NBINS + 1))
    print(f"closure    : curve integral {integral:.6g} vs N*w {tot_cut:.6g}")


if __name__ == "__main__":
    main()
