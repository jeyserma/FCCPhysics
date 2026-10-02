"""
Convert the BDSIM output ROOT file directly into EDM4hep.

Replaces convert.py (which wrote HEPEvt text). Same sampler, same physics cuts,
same time reference -- only the output format differs. See README.md section 9.

IMPORTANT -- this script runs in the key4hep stack, NOT the BDSIM stack:

    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
    python convert_edm4hep.py --input output_<seed>.root

It deliberately does not use pybdsim. The BDSIM sampler is stored as an unsplit
BDSOutputROOTEventSampler<float>, but ROOT can read it through the StreamerInfo
embedded in the file, so no BDSIM dictionaries (and no BDSIM stack) are needed.
TTree::Draw resolves the members; RDataFrame and uproot cannot (uproot falls back
to a pure-Python AsObjects path that is ~400x slower).

Units written, following the EDM4hep conventions:
    momentum  GeV       vertex  mm        time  ns        mass  GeV
Note that HEPEvt used mm/c for the time; here it is nanoseconds.
"""

import argparse
import os

import sample_config
import numpy as np

import ROOT
from podio import root_io, Frame
import edm4hep

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError   # silence the expected "no dictionary" notices

C_M_PER_NS = 0.299792458          # speed of light, m/ns
SAMPLER = "QD0AL"                 # the sampler element, 2.4 m upstream of the IP
Z_SAMPLER_M = -2.4                # its position with respect to the IP, in metres
DRIFT_TO_S6_M = 2.4 + 6.0         # extrapolation length for the rp cut -> s = +6.0 m


def draw_columns(tree, columns, n_estimate):
    """
    Read `columns` (at most 4) from `tree` with TTree::Draw and return them as
    numpy arrays. Draw is used because the sampler branch is an unsplit object
    that only the StreamerInfo can resolve.
    """
    tree.SetEstimate(n_estimate)
    n = tree.Draw(":".join(columns), "", "goff")
    if n < 0:
        raise RuntimeError(f"TTree::Draw failed for {columns}")
    out = []
    for i in range(len(columns)):
        buf = tree.GetVal(i)
        buf.reshape((n,))
        out.append(np.array(buf, dtype=np.float64, copy=True))
    return out


def sampler_s_position(rootfile):
    """
    Path length from the start of the beamline to the sampler plane, in metres.

    Read from Model.samplerSPosition rather than by looking up the element name
    in Model.componentName: the component names are a vector<string> inside the
    emulated Model class, which ROOT cannot hand back without the BDSIM
    dictionary. The last sampler is taken, matching pybdsim's SamplerData(d, -1)
    used by convert.py. Cross-checked against pybdsim's SEnd for this lattice:
    286.865599 m vs 286.865601 m.
    """
    f = ROOT.TFile.Open(rootfile)
    model = f.Get("Model")
    if not model or model.GetEntries() == 0:
        raise RuntimeError(f"No Model tree in {rootfile}")
    model.SetEstimate(1000000)
    n = model.Draw("Model.samplerSPosition", "", "goff")
    if n < 1:
        raise RuntimeError(f"No sampler found in the Model tree of {rootfile}")
    buf = model.GetVal(0)
    buf.reshape((n,))
    return float(np.array(buf, dtype=np.float64, copy=True)[-1])


def read_photons(rootfile):
    """Read the sampler, apply the skim, return a dict of numpy arrays."""
    f = ROOT.TFile.Open(rootfile)
    tree = f.Get("Event")
    if not tree:
        raise RuntimeError(f"No Event tree in {rootfile}")

    # hits per primary, needed to size the Draw buffers and to keep track of
    # which primary positron each photon came from (see the splitting in main())
    n_primaries = int(tree.GetEntries())
    tree.SetEstimate(n_primaries + 1)
    tree.Draw(f"{SAMPLER}.n", "", "goff")
    nbuf = tree.GetVal(0)
    nbuf.reshape((n_primaries,))
    n_per_primary = np.array(nbuf, dtype=np.int64, copy=True)
    n_hits = int(n_per_primary.sum())
    if n_hits == 0:
        empty = {k: np.array([]) for k in ("px", "py", "pz", "E", "x", "y", "t", "primary")}
        return empty, 0, n_primaries

    c = [f"{SAMPLER}.{v}" for v in ("x", "y", "xp", "yp", "zp", "energy", "p", "T", "partID")]
    x, y, xp, yp = draw_columns(tree, c[0:4], n_hits)
    zp, energy, p, T = draw_columns(tree, c[4:8], n_hits)
    (partID,) = draw_columns(tree, c[8:9], n_hits)

    # ---- skim, identical to convert.py ----
    # photons only; forward-going only; E > 2 keV;
    # r < 18 mm at the sampler; r > 9 mm extrapolated to s = +6.0 m
    with np.errstate(divide="ignore", invalid="ignore"):
        slope_x = np.where(zp != 0, xp / zp, 0.0)
        slope_y = np.where(zp != 0, yp / zp, 0.0)
    r = np.hypot(x, y)
    rp = np.hypot(x + slope_x * DRIFT_TO_S6_M, y + slope_y * DRIFT_TO_S6_M)
    keep = (partID == 22) & (zp > 0) & (energy > 2e-6) & (r < 18e-3) & (rp > 9e-3)

    # index of the primary positron that produced each surviving photon; the
    # sampler is filled in primary order, so this array is sorted
    primary = np.repeat(np.arange(n_primaries, dtype=np.int64), n_per_primary)[keep]

    x, y, xp, yp, zp = x[keep], y[keep], xp[keep], yp[keep], zp[keep]
    energy, p, T = energy[keep], p[keep], T[keep]

    # time of the reference particle at the IP, in ns: it still has
    # (endS + 2.4) m to travel from the start of the beamline
    t0_ns = (sampler_s_position(rootfile) - Z_SAMPLER_M) / C_M_PER_NS

    return {
        "px": xp * p,                    # GeV  (xp/yp/zp are momentum direction cosines)
        "py": yp * p,
        "pz": zp * p,
        "E": energy,                     # GeV
        "x": 1e3 * x,                    # mm
        "y": 1e3 * y,                    # mm
        "t": T - t0_ns,                  # ns, relative to the bunch crossing at the IP
        "primary": primary,
    }, n_hits, n_primaries


WEIGHT_FORMULA = ("weight = chargeFraction * bunchIntensity / "
                  "sum(nPrimaries over all events used); "
                  "counts per bunch crossing, one beam")


def build_metadata(args, rootfile, n_primaries, n_events, n_photons, s_sampler):
    """
    Everything needed to normalise this file, so that the weight can be
    recomputed from the data rather than from an external note.

    The weight is deliberately expressed per *event* rather than per file:
    summing `nPrimaries` over whatever events were actually read gives the right
    denominator for any subset, any ngenerate and any --events-per-file.

    chargeFraction and bunchIntensity come from sample_config.py. They are not
    simulation inputs -- BDSIM never sees them -- they are applied afterwards to
    turn counts into rates.
    """
    return {
        # --- normalisation ---
        "runConfig": args.runconfig,
        # which BDSIM build made the sample (key of BDSIM_SOURCES in submit_new.py)
        "bdsimSource": os.environ.get("BDSIM_SOURCE", "unknown"),
        "chargeFraction": float(sample_config.CHARGE_FRACTION[args.runconfig]),
        "bunchIntensity": float(sample_config.BUNCH_INTENSITY),
        "weightFormula": WEIGHT_FORMULA,
        # --- measured from the file ---
        "nPrimariesTotal": int(n_primaries),
        "nEvents": int(n_events),
        "nPhotons": int(n_photons),
        "samplerName": SAMPLER,
        "samplerSPosition_m": float(s_sampler),
        "zSampler_mm": float(1e3 * Z_SAMPLER_M),
        # --- the skim that was applied ---
        "cutEnergyMin_GeV": 2e-6,
        "cutRadiusMax_m": 18e-3,
        "cutRadiusMinExtrapolated_m": 9e-3,
        "cutExtrapolationLength_m": float(DRIFT_TO_S6_M),
        "cutForwardGoingOnly": 1,
        # --- provenance ---
        "inputFile": os.path.basename(rootfile),
        "mcParticleCollection": args.collection,
        "units": "momentum GeV, vertex mm, time ns, mass GeV",
    }


def write_edm4hep(photons, outfile, collection, category, n_primaries, n_events=1,
                  metadata=None):
    """
    Write the photons as `n_events` Frames.

    n_events=1 (the default) reproduces the HEPEvt granularity: one file, one
    event, every photon in it.

    For n_events>1 the *primary positrons* are divided into that many contiguous
    blocks and each block's photons become one event. Splitting on primaries
    rather than on photons keeps all the photons radiated by a given positron in
    the same event, and makes each output event an independent sample of the
    halo. The event boundary is arbitrary either way -- one file is a fixed
    number of primaries, not a bunch crossing -- and the per-photon weight is
    unaffected, so the normalisation of section 7 does not change.
    """
    px, py, pz = photons["px"], photons["py"], photons["pz"]
    x, y, t = photons["x"], photons["y"], photons["t"]
    primary = photons["primary"]
    z = 1e3 * Z_SAMPLER_M

    n_events = max(1, min(int(n_events), max(1, n_primaries)))
    # primary index at which each event starts; `primary` is sorted, so a binary
    # search turns those into photon-array slice points
    edges = np.linspace(0, n_primaries, n_events + 1).astype(np.int64)
    cuts = np.searchsorted(primary, edges)

    writer = root_io.Writer(outfile)

    # one frame of parameters describing the whole file
    if metadata:
        meta = Frame()
        for key, value in metadata.items():
            meta.put_parameter(key, value)
        writer.write_frame(meta, "metadata")

    for ev in range(n_events):
        lo, hi = int(cuts[ev]), int(cuts[ev + 1])
        particles = edm4hep.MCParticleCollection()
        for i in range(lo, hi):
            mcp = particles.create()
            mcp.setPDG(22)
            mcp.setGeneratorStatus(1)        # stable, to be tracked
            mcp.setCharge(0.0)
            mcp.setMass(0.0)
            mcp.setMomentum(edm4hep.Vector3d(px[i], py[i], pz[i]))
            mcp.setVertex(edm4hep.Vector3d(x[i], y[i], z))
            mcp.setTime(t[i])
        frame = Frame()
        frame.put(particles, collection)
        # Number of primary positrons this event accounts for. Carried on every
        # event so that summing it over whatever events are read gives the exact
        # denominator of the weight, independently of file size or splitting.
        frame.put_parameter("nPrimaries", int(edges[ev + 1] - edges[ev]))
        writer.write_frame(frame, category)
    # podio 1.7's Writer has no explicit finish()/close(); the underlying
    # ROOTWriter is finalised in its destructor, so drop the reference here
    # rather than relying on interpreter shutdown.
    del writer
    return len(px), n_events


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=str, required=True, help="BDSIM output ROOT file")
    parser.add_argument("--output", type=str, default=None,
                        help="EDM4hep output file (default: <input> with _edm4hep.root)")
    parser.add_argument("--collection", type=str, default="MCParticles",
                        help="Name of the MCParticle collection. The default 'MCParticles' is "
                             "what ddsim's EDM4hep reader looks for, so no ddsim flag is needed. "
                             "If it is not found ddsim does NOT fail -- it silently falls back to "
                             "the particle gun and writes a plausible-looking file of gun events.")
    parser.add_argument("--category", type=str, default="events", help="podio Frame category")
    parser.add_argument("--events-per-file", type=int, default=sample_config.EVENTS_PER_FILE,
                        help="Split the photons over this many events, by dividing the primary "
                             "positrons into contiguous blocks (default 1, i.e. one event per "
                             "file, matching the old HEPEvt granularity)")
    parser.add_argument("--runconfig", type=str, required=True,
                        choices=sorted(sample_config.CHARGE_FRACTION),
                        help="Which sample this is; selects the charge fraction from "
                             "sample_config.py and is recorded in the metadata. Adding a new "
                             "sample is just a new key in sample_config.py.")
    args = parser.parse_args()

    outfile = args.output or args.input.replace(".root", "_edm4hep.root")

    photons, n_raw, n_primaries = read_photons(args.input)
    n_events_req = max(1, min(int(args.events_per_file), max(1, n_primaries)))
    metadata = build_metadata(args, args.input, n_primaries, n_events_req,
                              len(photons["px"]), sampler_s_position(args.input))

    n, n_events = write_edm4hep(photons, outfile, args.collection, args.category,
                                n_primaries, args.events_per_file, metadata)
    print(f"{args.input}: {n_raw} sampler hits -> {n} photons from {n_primaries} primaries "
          f"written to {outfile} in {n_events} event(s) "
          f"(collection '{args.collection}', category '{args.category}')")


if __name__ == "__main__":
    main()
