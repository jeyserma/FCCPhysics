# BDSIM

Synchrotron-radiation backgrounds at the FCC-ee IP. BDSIM tracks positrons
through the last ~289 m of the beamline and records the SR photons crossing a
sampler 2.4 m before the IP; those photons are the input to ddsim.

Lattice: LCC_v1 (corresponds to LCC V106.2).

```bash
# report on a sample: files, size, weight, what is still missing (submits nothing)
python submit.py --analysis --lattice LCC_v2 --runconfig halo

# dry run (one job, locally, in /tmp/bdsim/<lattice>/<runconfig>/)
python submit.py --dryrun --lattice LCC_v2 --runconfig halo

# production
python submit.py --submit --cms_pool --lattice LCC_v2 --runconfig halo --njobs 2600
```

Output goes to `<storagedir>/<lattice>_<runconfig>/`.

`--analysis` reads `sample_config.py`, counts and sizes what is already on disk,
and reports the weight that set currently gives plus what completing it costs:

```
LCC_v2 / halo
  bunch intensity        2.02e+11
  charge fraction        0.01
  particles in one bunch 2.02e+09
  ngenerate per file     500000

  /ceph/.../LCC_v2_halo/
  files                  4133
  total size             21.18 GB
  mean file size         5.50 MB
  primaries simulated    2.066e+09
  current weight         0.977498

  files for weight 1     4040
  still needed           0
  size at weight 1       20.7 GB
```

It only stats the files, so it needs no key4hep and stays fast on a directory of
any size — which also means the primary count assumes every file used the
`ngenerate` it reports.

`--ngenerate N` overrides `sample_config.py` — useful with `--dryrun` for a
short test run, and with `--analysis` to see what a different job size implies.

> `--dryrun` exits 0 even when the job script fails — `dryrun()` does not check
> `subprocess.run`'s return value. Read the log, not the exit code.

## Configs

| | |
|---|---|
| `config/LCC_v1` | the original. Writes HEPEvt text. Kept to reproduce the 2026 campaign. |
| `config/LCC_v2` | **use this.** Fixed halo sampler, writes EDM4hep. The beam differs from v1 by the sampling method *only* — the energy-spread and dispersion bugs are documented but deliberately left in place (sections 2 and 3 of its README). |

## Two stacks per job

Generation and conversion do **not** share an environment:

| step | stack |
|---|---|
| `run_<runconfig>.py` (BDSIM) | `/cvmfs/beam-physics.cern.ch/bdsim/...v1.7.7` |
| `convert_edm4hep.py` (EDM4hep) | `/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08` |

`convert.sh` sources key4hep in its own subshell and first clears `PYTHONHOME`,
`PYTHONPATH` and `LD_LIBRARY_PATH` — without that, the BDSIM stack's Python 3.9
leaks into key4hep's 3.13 and the job dies with
`Failed to import encodings module`.

## Where the numbers live

`config/<lattice>/sample_config.py` holds four constants, and nothing else:

```python
BUNCH_INTENSITY = 2.02e11
CHARGE_FRACTION = {"halo": 0.01, "core": 0.99}
NGENERATE       = {"halo": 500000, "core": 1000000}
EVENTS_PER_FILE = 1
```

`ngenerate` comes from here, not from `submit.py` (`--ngenerate` only overrides
it, which is handy for short test runs). The halo tail shape (`xweight`,
`yweight`, `xtail`, `ytail`) stays in `run_halo.py`; masks, solenoid and lattice
options stay in `LCC_v2.py`.

## Adding a run config

To add e.g. `halo1`:

1. add a `"halo1"` key to `CHARGE_FRACTION` and `NGENERATE` in `sample_config.py`
2. copy `run_halo.py` to `run_halo1.py` and point it at `NGENERATE["halo1"]`

Then `--runconfig halo1` works everywhere. Nothing else is hardcoded.

## Normalising the output

Each EDM4hep file carries the normalisation in a `metadata` frame, and every
event carries `nPrimaries`. So:

```python
weight = chargeFraction * bunchIntensity / sum(nPrimaries over all events read)
```

which is correct for any subset of files and any `ngenerate` — no external file
count needed. Counts are **per bunch crossing, one beam**. For a full 1-bunch
halo production `sum(nPrimaries) = 2.02e9` and the weight is exactly 1.0.

## Known limitations

These apply to **both** the core and halo samples and are unchanged in v2 — see
section 7 of `config/LCC_v2/README.md`:

- the beamline stops at QD0AL, so the 7 mm QC1 mask is **not** simulated
- no solenoid, no orbit correctors
- multipoles zeroed (`linear=True`)
- no grazing-incidence X-ray reflection

And two known bugs in the halo beam, left in place on purpose so that v2 isolates
the sampler change — sections 2 and 3 of `config/LCC_v2/README.md`, with the fix
written out as a comment in `generate_4d_distribution`:

- energy spread is 1 MeV instead of 45.6 MeV (`sigma_delta` 2.2e-5, not 1e-3)
- the dispersive orbit `D*delta` is never added to the halo coordinates
