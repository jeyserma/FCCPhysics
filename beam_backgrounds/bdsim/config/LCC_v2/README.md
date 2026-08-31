# LCC_v2

Lattice LCC_v1 (LCC V106.2), same TWISS and aperture files as `config/LCC_v1`.
**The machine model is unchanged** — the generated `input_sequence.gmad`,
`input_components.gmad` and `input_beam.gmad` are byte-identical to those from
`LCC_v1`. What changed is how the beam-halo distribution is generated, plus a
handful of parameters that were silently ignored.

```bash
# dry run
python submit.py --dryrun --lattice LCC_v2 --runconfig halo

# production
python submit.py --submit --cms_pool --lattice LCC_v2 --runconfig halo --njobs 2600
```

Output goes to `<storagedir>/LCC_v2_halo/` and `<storagedir>/LCC_v2_core/`, so
nothing from the v1 campaign is overwritten.

---

## 1. Halo sampling rewritten (`generate_4d_distribution`)

**The physics model is unchanged.** Both versions target an exponential tail in
the Courant–Snyder invariant of each plane, truncated to an annulus:

```
p(Jx) ~ exp(-Jx/scale_x)   on [ (3.5 sigma_x)^2 , (xtail sigma_x)^2 ]
p(Jy) ~ exp(-Jy/scale_y)   on [ (4.0 sigma_y)^2 , (ytail sigma_y)^2 ]
scale_x = (3.5^2 emitx)/xweight        scale_y = (4^2 emity)/yweight
```

Only the numerical method changed.

### v1 — rejection + importance sampling

Scatter candidates uniformly over a 4-D box, discard those outside the annulus,
then draw `ngenerate` of the survivors **with replacement**, with probability
proportional to `exp(...)`. It is a general-purpose recipe — it works for any
density you can evaluate — but it is wasteful, and here it failed in two
compounding ways.

**The candidate pool was too small.** Only 0.36 % of the box lands inside the
annulus, so `init_size = 250 * ngenerate` was chosen to make the number of
survivors come out near `ngenerate`. It came out *under*: **181 129 survivors
for 200 000 requested**. That is why `replace=True` is there — it is
structurally forced, not a stylistic choice. Even with perfectly flat weights,
drawing 200 000 times from a pool of 181 129 yields only ~121 000 distinct
particles.

**The weights then collapsed it further.** The box extends to
`ytail = 513 sigma_y` while the vertical weight scale is only
`sqrt(16/0.06) = 16 sigma_y`, so a candidate at the outer boundary is suppressed
by `exp(-987)`. The median survivor carries weight `exp(-495)`; only 552 of the
181 129 have a weight above `exp(-3)`. The 200 000 draws therefore landed on
1069 distinct points.

**Raising 250 would not have fixed it.** The effective sample size scales
linearly with the pool size, so reaching a genuine ESS of 200 000 would need
`250 -> 706 000`, i.e. 1.4e11 candidates and roughly **7.9 TB** of RAM. This is
not a knob that was set too low; the method is unusable at the target
statistics.

Measured at production settings (`ngenerate = 200000`, `xweight=1.2`,
`yweight=0.06`):

| | v1 | v2 |
|---|---|---|
| candidates thrown | 50 000 000 | — |
| candidates passing the annulus cut | 181 129 (0.36 %) | — |
| **distinct phase-space points in the file** | **1 069** | **200 000** |
| largest multiplicity of a single particle | 9 042 | 1 |
| Kish effective sample size | 70.8 | 200 000 |
| peak RSS of `genGMAD` | **3821 MB** | **746 MB** |
| ... of which module imports | 568 MB | 568 MB |
| ... of which scales with `ngenerate` | 3253 MB | 178 MB (~19 MB sampler, rest `to_csv`) |
| wall time | ~4 s + 50 M-point arrays | ~1 s |

### v2 — direct (inverse-CDF) sampling

Nothing is guessed and checked, because for this density the answer can be
written down. The distribution depends only on the two invariants, and in
action-angle coordinates it factorises completely:

```
p(Jx, Jy, phix, phiy) = p(Jx) . p(Jy) . uniform(phix) . uniform(phiy)
```

That is four independent 1-D problems. The two truncated exponentials invert in
closed form (`truncated_exponential`):

```
J = -a * ln( exp(-lo/a) - u * ( exp(-lo/a) - exp(-hi/a) ) ),   u ~ Uniform(0,1)
```

Every uniform random number maps to exactly one J, guaranteed to lie inside
[lo, hi] and distributed exactly as the truncated exponential — this is the
inverse-transform theorem, not an approximation. The betatron phases are
uniform, and the Courant-Snyder transform maps `(J, phase)` back to `(x, x')`.

One random number per particle, nothing rejected, nothing resampled, every
macroparticle distinct. **There is no `init_size` in v2 because there are no
candidates to reject** — the factor 250 did not become 1, it ceased to exist.

The key step is the change of variables, not the size of any factor: v1 sampled
in the coordinate space (x, x', y, y') where the density cannot be inverted;
v2 samples in the space where it separates.

**No support is lost by dropping the box.** It was only ever a proposal region,
never a physical boundary: the outer ellipse `Jx = (xtail sigma_x)^2` reaches
exactly `xtail sigma_x` in position and `xtail sigma_x'` in angle, so it is
inscribed in the box and touches all four sides. Nothing in the annulus was ever
clipped, and v1 and v2 cover the identical region (see the closure table below).

### Memory

The peak-RSS ratio is just the rejection rate. v1 held 5e7 candidates
(400 MB per float64 array) simultaneously — four coordinate arrays, two
invariants, the weights and the boolean masks, about 3 GB — plus the full-size
temporaries NumPy materialises for every sub-expression of
`gammax*dx**2 + 2*alfx*dx*dxp + betx*dxp**2`, since it does not fuse
expressions. With no rejection, v2's arrays are `ngenerate` long (1.6 MB each),
~19 MB in total; its remaining footprint is the `pd.DataFrame(...).to_csv()`
that writes `inputfile.dat`, unchanged from v1.

### Closure

Verified to reproduce the same distribution:

| | v2 sampled | analytic target |
|---|---|---|
| `<Jx>/emitx` | 22.4476 | 22.4421 |
| `<Jy>/emity` | 282.30 | 282.67 |
| Jx range / emitx | 12.250 … 99.924 | cuts at 12.25 … 100 |

**Expected effect on results: none on the mean.** In the v1 scaling tests the
yield was flat at 0.281 ± 0.004 (ESS≈35) and 0.2855 ± 0.0015 (ESS≈354)
photons/primary — no drift over a factor 10 in effective statistics, so the
v1 campaign mean is sound. What changes is the **precision**: the per-file
photon count scattered ~11 % in v1 (against 0.4 % Poisson, matching 1/sqrt(70.8));
in v2 it should be Poisson. Statistical errors on v1 output must be bootstrapped
over *files*, never over photons — that caveat goes away in v2.

**Side effect:** the RNG changed from the legacy global `np.random.seed` to
`np.random.default_rng(seed)`. The same seed produces a different (correct)
distribution than v1.

## 2. Energy spread fixed

v1:

```python
sampled_E = np.random.normal(45.6, 1.0e-3, size=size)
```

`SIGE = 0.001` in the TWISS header is `dE/E`, so the absolute spread should be
`0.001 x 45.6 GeV = 45.6 MeV`. v1 used **1 MeV**, a factor 45.6 too small
(measured `sigma_rel = 2.2e-5` in `inputfile.dat`). v2 draws
`delta ~ N(0, SIGE)` and sets `E = ENERGY x (1 + delta)`; measured
`sigma_rel = 1.002e-3`, `sigma_E = 45.67 MeV`.

## 3. Dispersive orbit added to the halo (new physics content)

BDSIM applies `distrType="userfile"` coordinates literally — the `dispx`/`dispxp`
in the beam block are only used by `gausstwiss`. So the dispersive offset has to
be added when writing the file, and v1 never did. With `Dx = 0.349 m` at the
start of the beamline and the corrected energy spread this is a **349 µm rms**
horizontal offset on top of an 813 µm betatron sigma — a ~9 % increase of
sigma_x in quadrature. In v1 it would have been 7.7 µm, i.e. nothing, which is
why the two fixes belong together: correcting the energy spread alone would have
had no transverse effect.

Note this still does **not** model an off-momentum halo (a tail *in delta*).
The halo remains a pure betatron-amplitude model with a Gaussian energy core.

## 4. Emittances and energy read from the TWISS header

v1 hardcoded `emitx = 0.7e-9`, `emity = 2.6e-12` and `45.6` GeV inside
`generate_4d_distribution`. They happen to match this lattice
(`@ EX 0.7e-9`, `@ EY 2.6e-12`, `@ ENERGY 45.6`), but would silently disagree on
any other TWISS file. v2 reads `EX`, `EY`, `SIGE`, `ENERGY` from the header
(`read_tfs_header`) and raises if they are missing. The TWISS path is now the
module constant `TWISS_FILE` instead of four copies of the same string literal.

## 5. Arguments that were silently ignored

- **`--withSol`, `--withCorr`** — `run_halo.py`/`run_core.py` parsed them but
  passed the hardcoded `withSol=0, withCorr=0` to `MDIStudy`. Now passed
  through. Defaults are still 0, so **default behaviour is unchanged** (still no
  solenoid), but `--withSol 1` now does something.
- **`--withDip`** — parsed but never passed. Now passed through, *but it still
  has no effect*: `idx_stop` is unconditionally overwritten with `idx(QD0AL)+1`
  further down, so the beamline always ends at QD0AL. A warning comment was
  added at the point of the override. See §7.
- **`XP0`** — the `gausstwiss` branch did
  `text[i].replace("Xp0=0.0", ...)`, but pybdsim writes `Xp0=-0.0`, so the
  substitution never fired and a non-zero `XP0` was silently dropped. Replaced
  by regex substitutions for all four of `X0`, `Y0`, `Xp0`, `Yp0`. The only
  change to the generated core beam block is `Xp0=-0.0` → `Xp0=0`, numerically
  identical, so **the core sample is unaffected**.

## 6. `convert.py`

- **Forward-going photons only.** ~0.11 % of sampled photons have `zp < 0`,
  i.e. they travel back up the beamline away from the IP. They cannot reach the
  detector and are now removed.
- **Correct extrapolation slope.** The `rp > 9 mm` cut extrapolates 8.4 m
  downstream of the sampler (to s = +6.0 m). v1 used `xp` and `yp` directly as
  slopes; they are momentum direction cosines, so the slope is `xp/zp`.
  Identical for `zp ~ 1`, wrong-signed for the backward photons above.
- **Vectorised.** v1 called `result.iloc[i]` ten times per photon inside a
  Python loop (~140 s per production job). Same output format, ~100x faster.

The three physics cuts are otherwise **unchanged**: `E > 2 keV`,
`r < 18 mm` at the sampler, `r > 9 mm` extrapolated to s = +6.0 m.
Note that the third is often quoted as "9 mm at +8.4 m from the IP" — it is
8.4 m from the **sampler**, i.e. +6.0 m from the IP.

## 7. Deliberately NOT changed

These are open modelling decisions, not bugs, and they affect the **core and
halo samples equally**. Several are likely larger than anything fixed above.

- **The beamline ends at QD0AL, 2.4 m before the IP.** `MASK_QC1L` (the 7 mm
  jaw mask), `DRIFT_L0`, `DRIFT_L1`, `DRIFT_SOL` and `DRIFT_R1` are written into
  `input_components.gmad` but never appear in the sequence. Consequence:
  **`maskA[2]`, `maskA[3]` and `roundmask` have no effect**; only `MASK_QC2L`
  (15 x 15 mm, from `maskA[0..1]`) is in the line. This is only correct if the
  downstream ddsim geometry provides the mask.
- **No solenoid** (`withSol=0`). The detector solenoid / anti-solenoid distort
  the vertical orbit through the final doublet, which is the plane the SR fan is
  most sensitive to.
- **Multipoles zeroed.** `MadxTfs2Gmad(..., linear=True)` sets `k2=0` on SDM1L
  (TWISS `K2L = -0.063`) and `knl={0.0,...}` on OCT0L/OCT1L/DEC1L.
- **No grazing-incidence X-ray reflection.** The photons arrive at a median
  1.16 mrad to the axis, right at the total-external-reflection threshold for
  the soft part of the spectrum. BDSIM 1.7.7 supports this
  (`physicsList="xray_reflection ..."`, `xrayAllSurfaceRoughness`); the lines
  are still commented out in `genGMAD`.
- **The halo model parameters** `xweight=1.2`, `yweight=0.06`, `xtail=10`,
  `ytail=513`. Note the tension: `yweight=0.06` puts the vertical tail scale at
  ~16 sigma_y, so a particle at the `ytail=513 sigma_y` boundary is suppressed by
  `exp(-987)`. The outer vertical halo — the part that actually scrapes the
  aperture — carries no weight in this model.
- **The normalisation.** Nothing in the code encodes the 5-minute lifetime or
  the 1 % halo charge fraction; those live only in the sample documentation.
  The weight is
  `w = charge_fraction x 2.02e11 / (ngenerate x n_files)`
  = 0.467 for the halo (1 %) and 15.34 for the core (99 %), giving counts **per
  bunch crossing, one beam**.

## 8. Production sizing (measured)

> **Note:** this section measures the **HEPEvt** path (`convert.py`). If you use
> the EDM4hep converter of section 9 the conversion is ~10x faster, uses ~3x less
> memory and writes ~4.4x less output, which lifts the `ngenerate` ceiling from
> ~865 000 to the runtime limit — see "This moves the job-size limits of section
> 8" at the end of section 9.

Scan of one halo job per point, `ngenerate` = 1e3 … 5e5, LCC_v2, BDSIM v1.7.7,
run on the MIT submit node (idle). The two steps are the ones the condor job
actually runs: `python run_halo.py` (genGMAD + BDSIM) then `python convert.py`.

| ngenerate | gen time | gen RSS | conv time | conv RSS | .root | .hepevt | gamma/primary |
|---|---|---|---|---|---|---|---|
| 1 000 | 11.2 s | 591 MB | 5.8 s | 687 MB | 2.8 MB | 0.04 MB | 0.2510 |
| 5 000 | 21.0 s | 597 MB | 6.0 s | 706 MB | 13.0 MB | 0.24 MB | 0.2834 |
| 20 000 | 57.0 s | 692 MB | 7.8 s | 831 MB | 42.5 MB | 1.00 MB | 0.2890 |
| 50 000 | 128.8 s | 701 MB | 11.6 s | 973 MB | 85.9 MB | 2.43 MB | 0.2825 |
| 100 000 | 248.0 s | 716 MB | 17.4 s | 1157 MB | 158.2 MB | 4.86 MB | 0.2824 |
| 200 000 | 485.3 s | 724 MB | 28.5 s | 1518 MB | 302.6 MB | 9.76 MB | 0.2833 |
| 500 000 | 1203.9 s | 948 MB | 62.3 s | 2600 MB | 736.4 MB | 24.27 MB | 0.2818 |

Least-squares fits over all seven points:

```
gen time [s]   =  8.9  + 0.002389 * N      (max residual 0.5 %)
conv time [s]  =  5.7  + 0.000113 * N      (5.1 %)
gen RSS [MB]   =  632  + 0.000625 * N      (7.0 %)
conv RSS [MB]  =  739  + 0.003770 * N      (8.1 %)   <-- binding constraint
hepevt [MB]    =  0    + 4.85e-5  * N      (48.5 bytes per primary)
root [MB]      =  9    + 0.001458 * N      (worker-node scratch only, never transferred)
```

**Timing caveat:** 2.39 ms/primary is this node. The v1 condor logs give
**6.4–9.0 ms/primary** across grid sites (SPRACE, NCHC, ACNCA, Ultralight), so
plan with **8 ms/primary**, i.e. ~3.3x slower than the numbers above.

Against typical condor limits (4000 MB RAM, few hours, 0.5 GB output transfer):

| limit | caps `ngenerate` at |
|---|---|
| output 0.5 GB | 10 300 000 |
| runtime 3 h at 8 ms/primary | 1 350 000 |
| generation RAM 4000 MB | 5 860 000 |
| **conversion RAM 4000 MB** | **865 000** |

Note this is the opposite of v1, where the *generation* step was the memory hog
(3.8 GB at N = 200 000, see §1). With the sampler fixed, `convert.py` is now the
constraint, because `pybdsim.Data.Load()` reads the whole BDSIM ROOT file into
memory. To go above ~865 000 it would have to stream instead: uproot 4.3.7 is in
the BDSIM stack and `convert.py` needs only ~10 branches of the single QD0AL
sampler plus `SEnd` from the Model tree, so chunked reading would flatten its
memory at ~740 MB and move the cap to runtime (~1.35 M/job).

| ngenerate | grid runtime | conv | peak RAM | output | scratch | jobs for 1 bunch |
|---|---|---|---|---|---|---|
| 200 000 | 0.44 h | 0.5 min | 1493 MB | 9.7 MB | 0.30 GB | 10 100 |
| 300 000 | 0.67 h | 0.7 min | 1870 MB | 14.6 MB | 0.45 GB | 6 733 |
| **500 000** | **1.11 h** | **1.0 min** | **2623 MB** | **24.3 MB** | **0.74 GB** | **4 040** |
| 700 000 | 1.56 h | 1.4 min | 3377 MB | 34.0 MB | 1.03 GB | 2 886 |
| 865 000 | 1.92 h | 1.7 min | 3999 MB | 42.0 MB | 1.27 GB | 2 335 |

### Recommended operating point

`ngenerate = 500000`, `RequestMemory = 3000`.

3000 rather than 4000 puts condor's `1.1 * RequestMemory` hold threshold at
3300 MB, just above the measured 2623 MB — that catches genuine runaways without
holding healthy jobs. 700 000 is not recommended: 3377 MB leaves only 10 % margin
against a 4000 MB ceiling, and grid nodes vary.

### One full bunch of halo

The halo is 1 % of a 2.02e11 bunch, so **1 bunch = 2.02e9 macroparticles**, and
because every v2 macroparticle is an independent particle the weight becomes
exactly **1.0**.

| | |
|---|---|
| jobs at `ngenerate = 500000` | **4 040** |
| total output | **~98 GB** |
| total CPU at 8 ms/primary | ~4 500 CPU-h ≈ 190 CPU-days |

For scale: the v1 campaign had 4.33e9 macroparticles = 2.14 bunches *nominal*,
but only 1.53e6 effective, i.e. 7.6e-4 of a bunch. Reaching a genuine full bunch
therefore takes about half the number of jobs v1 already ran.

### Yield check

Across the seven independent runs the yield is **0.2837 +/- 0.0024
gamma/primary**, against the v1 campaign mean of 0.2843 — agreement to 0.2 %.
This is the high-statistics confirmation that the sampler rewrite moved the
variance and not the physics. The yield is also flat from N = 5000 upwards, so
short jobs are statistically safe in v2; in v1 the N = 1000 ensemble was 2.7x
high because a handful of cloned trajectories dominated it.

## 9. EDM4hep output (`convert_edm4hep.py`)

`convert.py` (HEPEvt text) is superseded by `convert_edm4hep.py`, which writes
EDM4hep directly from the BDSIM ROOT file — consistent with the GuineaPig pairs,
and far easier to merge or overlay later. `convert.py` is kept for reference and
for reproducing v1 output.

### It runs in a different stack

podio and edm4hep are **not** in the BDSIM stack, so the conversion cannot share
a shell with the generation:

| step | stack |
|---|---|
| `run_halo.py` (genGMAD + BDSIM) | `/cvmfs/beam-physics.cern.ch/bdsim/...v1.7.7` |
| `convert_edm4hep.py` | `/cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08` |

Both are on CVMFS, so one condor job can source one and then the other in
separate subshells. `convert.sh` is the wrapper that does this:

```bash
bash convert.sh <seed> [collection]      # sources key4hep, runs the converter
```

To wire it into `submit.py`, replace the conversion line of `make_script` and the
check that follows it, and change `transfer_output_files` from
`output_$(SEED).hepevt` to `output_$(SEED)_edm4hep.root`:

```bash
python convert.py --input "output_${{seed}}.root"        # LCC_v1
bash convert.sh "${{seed}}"                              # LCC_v2

if [ ! -f "output_${{seed}}_edm4hep.root" ]; then ...    # was .hepevt
```

Note the **doubled braces**: `make_script` returns an f-string, so every literal
shell `${...}` must be written `${{...}}` or Python tries to interpolate it and
`make_script()` raises `NameError` before any job is written.

**Stack layering — the reason `convert.sh` unsets three variables.** By the time
the conversion runs, the job script has already sourced the BDSIM stack, which
sets `PYTHONHOME` to its own LCG Python 3.9. Sourcing key4hep on top of that
leaves its Python 3.13 unable to find its stdlib, and the job dies with

```
Fatal Python error: Failed to import encodings module
ModuleNotFoundError: No module named 'encodings'
```

`convert.sh` therefore does `unset PYTHONHOME PYTHONPATH LD_LIBRARY_PATH` before
sourcing key4hep. This only shows up when the converter runs *after* the BDSIM
stack in the same job — testing `convert.sh` from a clean shell does not catch
it. Validate with a full dryrun instead:

```bash
python submit.py --dryrun --lattice LCC_v2 --runconfig halo --ngenerate 2000   # --ngenerate only overrides sample_config.py
```

and **read the log rather than the exit code**: `submit.py`'s `dryrun()` calls
`subprocess.run` without checking the return value, so a failed dryrun still
exits 0. In a real condor job the failure is caught — `set -euo pipefail` aborts
the script and `on_exit_hold` holds the job.

### It does not use pybdsim

The sampler is stored as an unsplit `BDSOutputROOTEventSampler<float>`, but ROOT
reconstructs the class from the **StreamerInfo embedded in the file**, so no BDSIM
dictionaries — and therefore no BDSIM stack — are needed to read it. Two things
that do *not* work, and why `TTree::Draw` is used instead:

- **RDataFrame** sees only the top-level `QD0AL.` column, not its members.
- **uproot** *can* read it (`AsObjects(Model_BDSOutputROOTEventSampler_3c_float_3e_)`)
  but through a pure-Python path: **16 s for 2000 events vs 0.04 s** for
  `TTree::Draw`, i.e. ~400x slower, which would be over an hour at N = 500 000.

The sampler S position comes from `Model.samplerSPosition`, not from looking up
the element in `Model.componentName` — the component names are a
`vector<string>` inside the emulated class and cppyy cannot return them without
the dictionary. Cross-checked against pybdsim's `SEnd`: **286.865599 m** vs
**286.865601 m** (float32 vs double).

### `sample_config.py` and the normalisation metadata

Section 7 notes that in v1 nothing in the data encoded the normalisation — the
1 % / 99 % charge split and the bunch intensity lived only in a text note next to
the samples. In v2 they are defined in one small file and travel with the data.

`sample_config.py` holds four things and nothing else:

```python
BUNCH_INTENSITY = 2.02e11
CHARGE_FRACTION = {"halo": 0.01, "core": 0.99}
NGENERATE       = {"halo": 500000, "core": 1000000}
EVENTS_PER_FILE = 1
```

It is imported by `run_halo.py` / `run_core.py` (for `NGENERATE`) and by
`convert_edm4hep.py` (for the rest). It is deliberately dependency-free —
standard library only — because those two live in different stacks, which is also
why it is not part of `LCC_v2.py`, which cannot be imported at all under key4hep.

Everything else stays where it belongs: the halo tail shape (`xweight`,
`yweight`, `xtail`, `ytail`) in `run_halo.py`, the masks and lattice options in
`LCC_v2.py`.

Note that `BUNCH_INTENSITY` and `CHARGE_FRACTION` are **not simulation inputs** —
BDSIM never sees them. They are applied afterwards to turn simulated counts into
rates per bunch crossing.

The converter writes a `metadata` frame of 18 parameters — the two normalisation
constants, `runConfig`, the weight formula, what was measured from the file
(`nPrimariesTotal`, `nPhotons`, `samplerSPosition_m`) and the skim that was
applied — and **every event frame carries `nPrimaries`**, the number of primary
positrons that event accounts for.

That last point is the important one, because the weight denominator is then a
sum over the events actually read rather than a per-file constant:

```python
w = md.get_parameter("chargeFraction") * md.get_parameter("bunchIntensity") / sum_nPrimaries
```

correct for **any** subset of files, any mix of `ngenerate`, and any
`--events-per-file`. It replaces the v1 rule `weight = fullweight * ntot / n`,
which needed a fixed `ntot` supplied from outside and broke silently when the
file count drifted (as it did: the note said 21 641 files, the directory holds
21 644). For a full 1-bunch halo production `sum(nPrimaries) = 2.02e9` and the
weight is exactly **1.0**.

`--runconfig halo|core` is the converter's only required argument; it selects the
charge fraction and is recorded as a string in the metadata. The job script
passes it automatically:

```bash
bash convert.sh "${{seed}}" {runconfig}      # in submit.py's make_script
```

Note `{runconfig}` has *single* braces: it is an f-string substitution made by
Python when the script is generated, unlike `${{seed}}` which must reach the
shell.

### Splitting a file into several events (`--events-per-file`)

By default one file is one event, matching the old HEPEvt granularity — at
`ngenerate = 500000` that is a single Geant4 event with **~141 000 primary
photons**. For comparison, the v1 files that were already run through
ALFA_VTX108 held ~56 900, so this is 2.5x more per event.

`--events-per-file N` splits the photons over N events instead:

```bash
python convert_edm4hep.py --input output_<seed>.root --events-per-file 100
bash convert.sh <seed> MCParticle 100          # same thing through the wrapper
```

The split is made on **primary positrons**, not on photons: the primaries are
divided into N contiguous blocks and each block's photons become one event. That
keeps every photon radiated by a given positron in the same event, and makes each
output event an independent sample of the halo.

Nothing about the physics or the normalisation changes. The event boundary is
arbitrary in either case — one file is a fixed number of primaries, *not* a bunch
crossing (at `ngenerate = 500000`, one file is 1/4040 of a bunch; a real crossing
is ~5.7e8 photons summed over 4040 files at weight 1.0). The per-photon weight is
untouched, so section 7 still applies unchanged.

Verified: the set of photons written is bit-identical for
`--events-per-file` 1, 10 and 100 — only their distribution over events differs.
Cost is ~141 bytes of podio metadata per extra frame, i.e. +0.25 % on a
141 000-photon file split into 100 events.

Use it if ddsim's per-event memory or wall time becomes awkward, or to avoid
losing a whole file when one ddsim job dies. Leave it at 1 to keep exactly the
v1 behaviour.

### Physics content is unchanged

Same sampler, same skim (photons, forward-going, E > 2 keV, r < 18 mm at the
sampler, r > 9 mm extrapolated to s = +6.0 m), same time reference. One `Frame`
holding one `MCParticleCollection` corresponds exactly to one former HEPEvt file,
so the normalisation bookkeeping of §7 is untouched.

**Units follow the EDM4hep conventions:** momentum GeV, vertex mm, mass GeV, and
**time in ns** — note that HEPEvt used mm/c. A photon at the sampler now has
t = -8.0056 ns instead of -2400 mm/c; both mean "arrives at the IP at t = 0".

The MCParticle collection name defaults to `MCParticle` and is settable with
`--collection`; it must match ddsim's `--edm4hep.mcParticleCollectionName`
(the GuineaPig pairs use `Pairs`).

### Validation against `convert.py`

Same input file, 538 photons out of 4612 sampler hits in both:

| quantity | agreement |
|---|---|
| px, py, pz | **bit-identical** |
| x, y, z | **bit-identical** |
| time | 3.4e-7 relative (the float32 sampler position) |
| PDG / generatorStatus / mass | 22 / 1 / 0 |
| \|p\| vs E | 3e-8 (massless, as required) |

### Size

EDM4hep is **44.9 bytes per photon** against 172 bytes for the HEPEvt text — a
factor **3.8 smaller**. The write loop costs 1.7 us/particle after a one-off ~1 s
cppyy JIT, so the conversion is dominated by reading, not writing.

### Measured at production scale

`ngenerate = 500000` (seed 555000: 1 146 014 sampler hits -> 141 417 photons),
same input file for both converters:

| | `convert.py` -> HEPEvt | `convert_edm4hep.py` -> EDM4hep |
|---|---|---|
| wall time | 62.3 s | **6.3 s** |
| peak RSS | 2600 MB | **849 MB** |
| output size | 24.27 MB | **5.51 MB** |
| bytes per photon | 172 | **39.0** |
| bytes per primary | 48.5 | **11.0** |

### This moves the job-size limits of section 8

The HEPEvt path was capped by `convert.py`'s memory at `ngenerate ~ 865 000`.
With the EDM4hep converter that ceiling is gone, and **runtime becomes the
binding constraint**:

| limit | HEPEvt path | EDM4hep path |
|---|---|---|
| conversion RAM 4000 MB | 865 000 | ~16 000 000 |
| output 0.5 GB | 10 300 000 | ~45 000 000 |
| **runtime 3 h at 8 ms/primary** | 1 350 000 | **1 350 000** |

So the choice is now purely about how long a job you want:

| ngenerate | grid runtime | peak RAM | output | jobs for 1 bunch | total output |
|---|---|---|---|---|---|
| 500 000 | 1.1 h | ~850 MB | 5.5 MB | 4 040 | 22 GB |
| 1 000 000 | 2.2 h | ~950 MB | 11.0 MB | 2 020 | 22 GB |

Both are comfortable. 500 000 stays the safer default — 1 000 000 is 2.2 h at
8 ms/primary and 2.5 h at the slowest site seen in the v1 logs (9 ms), which is
tight against a wall-clock limit for no real gain. Either way the full-bunch
production drops from ~98 GB of HEPEvt to **~22 GB** of EDM4hep.
