#!/usr/bin/env python3
"""
Read the 'metadata' frame of EDM4hep files and collect it for a whole sample
into a small ROOT file of TParameter / TNamed objects, so the analysis chain can
compute the event weight without knowing how the sample was produced.

    weight = chargeFraction * bunchIntensity / nPrimariesTotal

giving photon counts per bunch crossing for one beam.

Usable as a module or from the command line.

    from extract_metadata import read_metadata
    meta = read_metadata("halo_00123.root")     # {} if the file has no metadata

    python3 extract_metadata.py -o halo_v2_meta.root /path/to/halo_v2/*.root

read_metadata() reads the 'metadata' frame and nothing else: no event frames, no
other Frame categories, no top-level ROOT objects. It does not assume which keys
are present -- whatever convert_edm4hep.py put there comes back -- so adding a
field to build_metadata() needs no change here.

The sample-level part (summing, consistency checks, the weight) is separate:
read_metadata() describes one file, read_sample() describes a sample.

Needs key4hep (podio + ROOT).
"""

import argparse
import contextlib
import glob
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)

# Frame categories that may hold the converter's metadata, most likely first.
# convert_edm4hep.py writes it as 'metadata', but ddsim renames it to 'meta' on
# the way through -- ddsim's own 'metadata' frame holds the CellID encodings --
# so the SR metadata is under a different name before and after detector
# simulation. Categories outside this list are searched too, so a further
# rename does not break the lookup.
METADATA_CATEGORIES = ("metadata", "meta")

# The frame that carries this key is the converter's, whatever it is called.
MARKER_KEY = "nPrimariesTotal"

EVENT_CATEGORY = "events"

# Keys that are per-file counters and must be added up across a sample.
SUMMED = {"nPrimariesTotal", "nEvents", "nPhotons"}

# Keys that legitimately differ from file to file and carry no sample-level
# meaning, so they are neither checked for consistency nor written out.
PER_FILE = {"inputFile"}

# Everything else describes the sample (cuts, sampler position, charge fraction,
# ...) and must agree across files; a mismatch means the file list mixes
# samples, which would silently produce a wrong weight.


@contextlib.contextmanager
def _quiet_root():
    """
    Silence ROOT's own error printing for the duration of a probe.

    Asking a file without a metadata tree for one is an expected miss, but ROOT
    reports it on stderr before raising, which would make a clean "no metadata
    here" look like a failure in the analysis logs.
    """
    previous = ROOT.gErrorIgnoreLevel
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    try:
        yield
    finally:
        ROOT.gErrorIgnoreLevel = previous


def read_metadata(path, strict=False, category=None):
    """
    The parameters of the converter's metadata frame of `path`, as a flat dict.

    Works on both convert_edm4hep.py output and the ddsim output made from it:
    the frame is located by content rather than by name, since ddsim renames the
    category (see METADATA_CATEGORIES). Pass `category` to force one name and
    skip the search.

    Returns {} when the file has no such frame -- that is the normal answer for
    a file produced before metadata was added, not an error. Set strict=True to
    raise instead, which is what you want when a missing weight would go
    unnoticed downstream.
    """
    if not os.path.exists(path):
        if strict:
            raise FileNotFoundError(path)
        return {}

    try:
        with _quiet_root():
            from podio import root_io

            reader = root_io.Reader(path)
            available = [c for c in reader.categories if c != EVENT_CATEGORY]

            if category is not None:
                candidates = [category]
            else:
                # Preferred names first, then anything else the file happens to
                # have, so an unexpected rename still resolves.
                candidates = ([c for c in METADATA_CATEGORIES if c in available] +
                              [c for c in available if c not in METADATA_CATEGORIES])

            for name in candidates:
                frames = reader.get(name)
                if len(frames) == 0:
                    continue
                # Bind the Frame to a name and build the dict while it is alive:
                # reading parameters off an already collected Frame is a
                # segfault, not an exception.
                frame = frames[0]
                keys = list(frame.parameters)
                if category is None and MARKER_KEY not in keys:
                    continue        # ddsim's own 'metadata' frame, not ours
                return {key: frame.get_parameter(key) for key in keys}

            raise KeyError(f"{path}: no metadata frame with '{MARKER_KEY}' "
                           f"(categories: {available})")

    except Exception:               # noqa: BLE001 -- not podio, unreadable, or no metadata
        if strict:
            raise
        return {}


def sr_weight(path, nevents):
    """
    The SR weight for `nevents` events of the sample that `path` belongs to.

    Multiply any photon count by this to get counts per bunch crossing, one beam.

    `path` is any one file of the sample -- only its constants are used, so it
    does not matter which. `nevents` is how many events the analysis actually
    read, across all files. The metadata counts primaries, not events, so the
    per-event primary count from the file converts between them; every job in a
    sample uses the same ngenerate, which is what makes one file enough.

        w = sr_weight("halo_00123.root", chain.GetEntries())
    """
    if nevents <= 0:
        raise ValueError(f"nevents must be positive, got {nevents}")

    meta = read_metadata(path, strict=True)
    primaries_per_event = float(meta["nPrimariesTotal"]) / float(meta["nEvents"])

    return (float(meta["chargeFraction"]) * float(meta["bunchIntensity"]) /
            (nevents * primaries_per_event))


def primaries_of(meta):
    """The number of primaries this file's metadata implies, or None."""
    for key in ("nPrimariesTotal", "nPrimaries"):
        if key in meta:
            return int(meta[key])
    return None


# --------------------------------------------------------------------------
# combining files into a sample
# --------------------------------------------------------------------------

def read_sample(paths, allow_missing=False):
    """
    Fold the per-file metadata of `paths` into one sample-level dict, including
    the weight. Counters in SUMMED are added; everything else must agree across
    files, or a ValueError names the two files that disagree.
    """
    const, const_from, sums = {}, {}, {}
    n_primaries, n_files = 0, 0
    skipped = []

    for path in paths:
        try:
            meta = read_metadata(path, strict=True)
        except Exception as exc:     # noqa: BLE001 -- report whatever ROOT/podio raised
            if not allow_missing:
                raise
            skipped.append((path, str(exc)))
            continue

        for key, value in meta.items():
            if key in PER_FILE:
                continue
            if key in SUMMED:
                sums[key] = sums.get(key, 0) + int(value)
            elif key not in const:
                const[key] = value
                const_from[key] = path
            elif const[key] != value:
                raise ValueError(
                    f"'{key}' differs between files -- the list mixes samples:\n"
                    f"  {const_from[key]}: {const[key]!r}\n"
                    f"  {path}: {value!r}"
                )

        n = primaries_of(meta)
        if n is None:
            raise KeyError(f"{path}: metadata has no primary count, cannot weight")
        n_primaries += n
        n_files += 1

    if n_files == 0:
        raise RuntimeError("no readable files with metadata")
    if n_primaries <= 0:
        raise ValueError(f"nPrimaries summed to {n_primaries}; cannot form a weight")

    out = dict(const)
    out.update({k: int(v) for k, v in sums.items()})
    out["nPrimaries"] = int(n_primaries)
    out["nFiles"] = int(n_files)
    out["nFilesSkipped"] = int(len(skipped))

    if "chargeFraction" in const and "bunchIntensity" in const:
        out["weight"] = (float(const["chargeFraction"]) *
                         float(const["bunchIntensity"]) / n_primaries)

    return out, skipped


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------

def _as_tobject(name, value):
    """The ROOT object matching a Python value's type."""
    if isinstance(value, str):
        # TParameter has no string specialisation; TNamed carries it in the
        # title, read back with GetTitle().
        return ROOT.TNamed(name, value)
    if isinstance(value, bool):                    # before int: bools are ints
        return ROOT.TParameter(bool)(name, value)
    if isinstance(value, int):
        # Long64_t, not int: nPrimaries for a full sample exceeds 2^31.
        # cppyy converts the Python int; ROOT.Long64_t itself does not exist.
        return ROOT.TParameter("Long64_t")(name, value)
    if isinstance(value, float):
        # "double" spelled out: ROOT.TParameter(float) is TParameter<float>,
        # single precision, which mangles bunchIntensity = 2.02e11.
        return ROOT.TParameter("double")(name, value)
    return ROOT.TNamed(name, str(value))


def write_metadata(meta, outfile):
    """Write a flat metadata dict as top-level TParameter / TNamed objects."""
    fout = ROOT.TFile.Open(outfile, "RECREATE")
    if not fout or fout.IsZombie():
        raise IOError(f"cannot open {outfile} for writing")
    for name, value in sorted(meta.items()):
        fout.WriteTObject(_as_tobject(name, value), name)
    fout.Close()


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------

def expand(items):
    """Expand directories and globs, keeping a stable order."""
    files = []
    for item in items:
        if os.path.isdir(item):
            files.extend(sorted(glob.glob(os.path.join(item, "*.root"))))
        elif any(ch in item for ch in "*?["):
            files.extend(sorted(glob.glob(item)))
        else:
            files.append(item)
    return files


def main():
    parser = argparse.ArgumentParser(
        description="Extract the metadata frame of EDM4hep files into a flat "
                    "ROOT file of TParameters.")
    parser.add_argument("inputs", nargs="+",
                        help="EDM4hep files, or directories/globs of them")
    parser.add_argument("-o", "--output",
                        help="output ROOT file. Omit to only print what was found.")
    parser.add_argument("--allow-missing", action="store_true",
                        help="skip files without metadata instead of failing. Use with "
                             "care: a silently skipped file makes the weight too large "
                             "for the events that remain.")
    parser.add_argument("--dump", action="store_true",
                        help="print the metadata of each input file separately and stop")
    args = parser.parse_args()

    files = expand(args.inputs)
    if not files:
        sys.exit("no input files matched")

    if args.dump:
        for path in files:
            meta = read_metadata(path)
            print(f"--- {path}")
            if not meta:
                print("    (no metadata)")
            for key, value in sorted(meta.items()):
                print(f"    {key:28s} {value!r}")
        return

    try:
        meta, skipped = read_sample(files, args.allow_missing)
    except (ValueError, KeyError, RuntimeError) as exc:
        sys.exit(str(exc))

    print(f"sample        : {meta.get('runConfig', '?')}")
    print(f"files         : {meta['nFiles']}" +
          (f"  ({len(skipped)} skipped)" if skipped else ""))
    print(f"nPrimaries    : {meta['nPrimaries']:,}")
    print(f"nPhotons      : {meta.get('nPhotons', 0):,}")
    if "weight" in meta:
        print(f"chargeFraction: {meta['chargeFraction']}")
        print(f"bunchIntensity: {meta['bunchIntensity']:.4g}")
        print(f"weight        : {meta['weight']:.6g}   "
              f"(photons per bunch crossing, one beam)")
    else:
        print("weight        : not computed -- "
              "chargeFraction / bunchIntensity missing from the metadata")

    if args.output:
        write_metadata(meta, args.output)
        print(f"wrote         : {args.output}")

    for path, exc in skipped:
        print(f"  SKIPPED {path}: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
