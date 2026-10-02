"""
sample_name.py
==============
Work out (lattice, runconfig, config directory, sample name) from a run path.

Samples are named by submit_new.py as

    <lattice>_<runconfig>_<bdsim source>[_<suffix>]      e.g. LCC_v2short_halo_standalone_G4SR

both on storage and for the --dryrun directory <storagedir>/<sample name>/tmp/. The older
layouts are still accepted:

    /ceph/.../bdsim/LCC_v2short_halo/        flat, no source
    /tmp/bdsim/LCC_v2short/halo/             nested, from submit_old.py --dryrun

Lattice names contain underscores themselves (LCC_v2short), so the name cannot
be split on "_". Instead the lattice is matched against the existing
config/<lattice>/ directories and the runconfig against their run_<runconfig>.py,
both as the leading part of the name. Whatever follows (source, suffix) stays in
the sample name and is never mistaken for the runconfig.
"""

import os

CONFIG_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            "..", "config"))


def _runconfigs(cfg):
    return sorted((f[len("run_"):-len(".py")] for f in os.listdir(cfg)
                   if f.startswith("run_") and f.endswith(".py")), key=len, reverse=True)


def _match_runconfig(cfg, rest):
    for rc in _runconfigs(cfg):
        if rest == rc or rest.startswith(rc + "_"):
            return rc
    return None


def resolve_sample(path, configdir=None):
    """
    (lattice, runconfig, configdir, name) for a run directory or a file inside it,
    or None if the path does not match any config/<lattice>/.

    `name` is the flat sample name, also for the nested layout, so it can be used
    directly as the output directory for plots.
    """
    path = os.path.normpath(os.path.abspath(path))
    if os.path.isfile(path):
        path = os.path.dirname(path)
    parts = path.split(os.sep)
    if parts[-1] == "tmp":                  # dryrun directory inside the sample
        parts = parts[:-1]

    if configdir:
        cfgs = [os.path.normpath(os.path.abspath(configdir))]
    else:
        cfgs = [os.path.join(CONFIG_ROOT, d) for d in os.listdir(CONFIG_ROOT)
                if os.path.isdir(os.path.join(CONFIG_ROOT, d))]
    cfgs.sort(key=lambda c: -len(os.path.basename(c)))     # LCC_v2short before LCC_v2

    for cfg in cfgs:
        lattice = os.path.basename(cfg)
        # flat: <lattice>_<runconfig>[_...]
        if parts[-1].startswith(lattice + "_"):
            rc = _match_runconfig(cfg, parts[-1][len(lattice) + 1:])
            if rc:
                return lattice, rc, cfg, parts[-1]
        # nested: <lattice>/<runconfig>[_...]
        if len(parts) >= 2 and parts[-2] == lattice:
            rc = _match_runconfig(cfg, parts[-1])
            if rc:
                return lattice, rc, cfg, f"{lattice}_{parts[-1]}"
    return None


def tried_message(path):
    lattices = sorted(d for d in os.listdir(CONFIG_ROOT)
                      if os.path.isdir(os.path.join(CONFIG_ROOT, d)))
    return (f"cannot work out the lattice / runconfig from {path}\n"
            f"  lattices known: {', '.join(lattices)}\n"
            "  expected .../<lattice>_<runconfig>[_<source>][_<suffix>]/ "
            "or .../<lattice>/<runconfig>/")
