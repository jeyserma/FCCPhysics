# BDSIM

Synchrotron-radiation backgrounds at the FCC-ee IP. BDSIM tracks positrons
through the last ~289 m of the beamline and records the SR photons crossing a
sampler 2.4 m before the IP; those photons are the input to ddsim.

Lattice: LCC_v1 (corresponds to LCC V106.2).

```bash
# how many files are needed for weight 1? (query only, submits nothing)
python submit.py --weightcalc --lattice LCC_v2 --runconfig halo

# dry run (one job, locally, in /tmp/bdsim/<lattice>/<runconfig>/)
python submit.py --dryrun --lattice LCC_v2 --runconfig halo

# production
python submit.py --submit --cms_pool --lattice LCC_v2 --runconfig halo --njobs 2600
```

Output goes to `<storagedir>/<lattice>_<runconfig>/`.

`--weightcalc` reads `sample_config.py`, counts what is already on disk, and
reports the weight that set currently gives plus how many files are still
missing for weight 1:

```
LCC_v2 / halo
  bunch intensity        2.02e+11
  charge fraction        0.01
  particles in one bunch 2.02e+09
  ngenerate per file     500000
  --> files for weight 1 4,040
  ...
  existing files         100
  current weight         40.4
  still needed           3,940
```

`--ngenerate N` overrides `sample_config.py` — useful with `--dryrun` for a
short test run, and with `--weightcalc` to see what a different job size implies.

> `--dryrun` exits 0 even when the job script fails — `dryrun()` does not check
> `subprocess.run`'s return value. Read the log, not the exit code.

## Configs

| | |
|---|---|
| `config/LCC_v1` | the original. Writes HEPEvt text. Kept to reproduce the 2026 campaign. |
| `config/LCC_v2` | **use this.** Fixed halo sampler and energy spread, writes EDM4hep. See `config/LCC_v2/README.md` for the full changelog and the measured job sizing. |

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
NGENERATE       = {"halo": 500000, "core": 500000}
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
