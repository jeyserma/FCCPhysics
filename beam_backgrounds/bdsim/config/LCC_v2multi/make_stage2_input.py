#!/usr/bin/env python3
"""
Bridge between the two stages: read the stage-1 BDSIM output, take everything
that crossed the split plane, and write it out as the stage-2 beam file (one
file per species -- see below).

    source env.sh                       # BDSIM stack (or key4hep, either works)
    python make_stage2_input.py output_12345.root

Writes GMAD/stage2_input.dat in the format

    pdgid  x[m]  xp[rad]  y[m]  yp[rad]  E[GeV]  t[ns]  weight

Design notes
------------
*Every species* crossing the plane is carried across, not just the positrons.
Photons and shower products already produced upstream of the split still have to
traverse the final focus -- some are stopped by the QC2 mask -- so dropping them
would silently remove the soft part of the spectrum.

*One copy only.* The replay factor comes from running stage 2 several times with
different --seed, not from replicating here. Both are computationally identical,
but seeds cost no extra disk and let the passes go to separate batch jobs. Each
stage-2 output then carries nPrimariesTotal = NGENERATE and the K outputs sum to
K * NGENERATE by themselves, which sr_weight() and read_sample() already handle.

*The catch.* Replaying one particle N times does NOT give N independent samples:
the copies share their incoming phase space. Deterministic trajectories -- an
upstream photon that simply drifts -- give literally identical copies and zero
gain. Run with --report to see how much of the dump is in that category.
"""

import argparse
import os
import sys

import numpy as np
import ROOT

import sample_config

M_E_GEV = 0.000510998950   # e+/e- mass, for the E -> Ek conversion

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError    # expected "no dictionary" notices


def find_sampler(rootfile):
    """
    Name of the stage-1 sampler branch.

    Selected by CLASS, not by name: BDSIM writes Primary, Eloss, Histos,
    ApertureImpacts and friends alongside the sampler, and several of those are
    also sampler-shaped. Only the plane samplers are
    BDSOutputROOTEventSampler<...>, and 'Primary' is one too, so it is excluded
    explicitly. Branch names carry a trailing dot, which is stripped.
    """
    f = ROOT.TFile.Open(rootfile)
    tree = f.Get("Event")
    if not tree:
        sys.exit(f"no Event tree in {rootfile}")
    names = []
    for br in tree.GetListOfBranches():
        cls = br.GetClassName() or ""
        name = br.GetName().rstrip(".")
        if "Sampler" not in cls:
            continue
        if name in ("Primary", "PrimaryGlobal"):
            continue
        names.append(name)
    f.Close()
    if not names:
        sys.exit(f"no sampler branch found in {rootfile}")
    if len(names) > 1:
        print(f"note: several sampler branches {names}, using {names[0]}")
    return names[0]


def draw(tree, columns, n_estimate):
    tree.SetEstimate(n_estimate)
    n = tree.Draw(":".join(columns), "", "goff")
    if n < 0:
        raise RuntimeError(f"TTree::Draw failed for {columns}")
    return [np.array([tree.GetVal(i)[j] for j in range(n)]) for i in range(len(columns))], n


def read_stage1(rootfile, sampler):
    f = ROOT.TFile.Open(rootfile)
    tree = f.Get("Event")
    n_primaries = int(tree.GetEntries())

    tree.SetEstimate(n_primaries + 1)
    tree.Draw(f"{sampler}.n", "", "goff")
    nb = tree.GetVal(0)
    n_hits = int(sum(nb[i] for i in range(n_primaries)))
    if n_hits == 0:
        sys.exit("stage-1 sampler is empty")

    (x, y, xp, yp), _ = draw(tree, [f"{sampler}.{v}" for v in ("x", "y", "xp", "yp")], n_hits)
    (zp, energy, T), _ = draw(tree, [f"{sampler}.{v}" for v in ("zp", "energy", "T")], n_hits)
    (partID, weight), _ = draw(tree, [f"{sampler}.{v}" for v in ("partID", "weight")], n_hits)
    f.Close()

    # xp/yp in the sampler are direction cosines; the userfile wants angles.
    # For the tiny angles here (microradians) they agree to far better than
    # needed, but divide by zp anyway so it stays right for scattered tracks.
    with np.errstate(divide="ignore", invalid="ignore"):
        xang = np.where(zp != 0, xp / zp, 0.0)
        yang = np.where(zp != 0, yp / zp, 0.0)

    # Backward-going particles cannot be injected into a forward line.
    forward = zp > 0
    dropped = int((~forward).sum())

    return dict(pdg=partID[forward].astype(int), x=x[forward], xp=xang[forward],
                y=y[forward], yp=yang[forward], E=energy[forward],
                t=T[forward], w=weight[forward]), n_primaries, n_hits, dropped


def main():
    ap = argparse.ArgumentParser(
        description="Replicate a stage-1 phase-space dump into the stage-2 beam file.")
    ap.add_argument("input", help="stage-1 BDSIM output_<seed>.root")
    ap.add_argument("-o", "--output", default=os.path.join("GMAD", "stage2_input.dat"),
                    help="only the directory is used; the files are named per species")
    ap.add_argument("--sampler", default=None, help="sampler branch name (auto by default)")
    ap.add_argument("--report", action="store_true",
                    help="report the species mix and how much of the dump is "
                         "deterministic (gains nothing from replay)")
    args = ap.parse_args()

    sampler = args.sampler or find_sampler(args.input)
    d, n_primaries, n_hits, dropped = read_stage1(args.input, sampler)
    n_part = len(d["pdg"])

    outdir = os.path.dirname(args.output) or "."
    os.makedirs(outdir, exist_ok=True)

    # Split by species and write one file each. BDSIM 1.7.7 applies the
    # REFERENCE beam particle's mass to every line regardless of pdgid, so a
    # mixed file injected with particle="e+" silently loses every photon below
    # 511 keV. One file per species, each run with a matching reference
    # particle, is the only safe way -- see README.md section 4.
    groups = {"e+": d["pdg"] == -11, "gamma": d["pdg"] == 22}
    other = ~(groups["e+"] | groups["gamma"])
    if other.any():
        import collections
        rest = collections.Counter(d["pdg"][other].tolist())
        print(f"WARNING: {int(other.sum())} particles of other species are NOT "
              f"carried across: {dict(rest)}")

    written = {}
    for species, mask in groups.items():
        idx = np.nonzero(mask)[0]
        path = os.path.join(outdir, os.path.basename(
            {"e+": "stage2_input_ep.dat", "gamma": "stage2_input_gamma.dat"}[species]))
        # Ek, not total E: mass is zero for the photon, m_e for the positron.
        mass = 0.0 if species == "gamma" else M_E_GEV
        with open(path, "w") as f:
            for i in idx:
                f.write("%d\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.9g\t%.6g\n" %
                        (d["pdg"][i], d["x"][i], d["xp"][i], d["y"][i], d["yp"][i],
                         max(d["E"][i] - mass, 0.0), d["t"][i], d["w"][i]))
        written[species] = (len(idx), path)

    print(f"stage-1 file        : {args.input}")
    print(f"sampler             : {sampler}")
    print(f"stage-1 primaries   : {n_primaries}")
    print(f"crossed the plane   : {n_hits}" +
          (f"  ({dropped} backward-going dropped)" if dropped else ""))
    for species, (cnt, path) in written.items():
        print(f"  {species:<6}            : {cnt:>8}  ->  {path}")
    print()
    # The converter defaults to sample_config.NGENERATE, because in production
    # stage 1 always runs with it. If stage 1 was run with anything else the
    # default is silently wrong, so say so loudly here.
    expected = sample_config.NGENERATE.get("halo")
    if expected is not None and n_primaries != expected:
        print()
        print(f"  *** stage 1 used {n_primaries} primaries, but sample_config.NGENERATE"
              f'["halo"] is {expected}.')
        print(f"  *** convert_edm4hep.py would normalise to {expected} and be WRONG by"
              f" x{expected / n_primaries:.4g}.")
        print(f"  *** pass --n-primaries {n_primaries} to the converter.")
        print()

    print(f"Each stage-2 run over this dump has nPrimariesTotal = {n_primaries}.")
    print("Run stage 2 K times with different --seed to get a replay factor K;")
    print("the K outputs sum to K * that on their own, so nothing needs to know K.")

    if args.report:
        print()
        uniq, cnt = np.unique(d["pdg"], return_counts=True)
        order = np.argsort(-cnt)
        print("species crossing the split plane:")
        for k in order:
            name = {22: "gamma", -11: "e+", 11: "e-", 2212: "proton",
                    2112: "neutron"}.get(int(uniq[k]), str(int(uniq[k])))
            print("  %-8s %8d  (%5.2f %%)" % (name, cnt[k], 100 * cnt[k] / n_part))
        gam = (d["pdg"] == 22).sum()
        print()
        print("Replay gains information only where the downstream physics is")
        print("stochastic. The %d photons (%.1f %% of the dump) mostly drift the" % (gam, 100 * gam / n_part))
        print("last 6.8 m deterministically, so their copies are near-identical:")
        print("their effective sample size does NOT grow with N, even though the")
        print("raw count does. Treat the soft part of the spectrum accordingly.")


if __name__ == "__main__":
    main()
