# LCC_v2multi — two-stage ("dual sampling") halo production

> **NOT RECOMMENDED — parked. Use `../LCC_v2short/` instead.**
>
> The implementation is *correct*: it reproduces the single-stage photon
> spectrum across the full range including the hard tail (E > 500 keV ratio
> 0.969 ± 0.074 on 9 passes). It is simply **slower than doing nothing**:
> best case 0.80× a single-stage run, because stage 1 costs 0.95 of a full run
> and re-injecting particles costs ~0.5 of a run in geometry-entry alone.
>
> `LCC_v2short` gets the intended saving — **2.13× ± 0.16** — by never running
> the 282 m at all, and its trajectories stay independent so there is no
> statistical ceiling. Replaying here caps at 3.25 bunch crossings of power
> however large K is (section 7).
>
> **Kept for the four BDSIM traps documented below**, all of which cost real
> debugging time and will bite anyone attempting a phase-space handoff again:
> the reference-particle mass silently eating sub-511 keV photons (section 4),
> `samplerDiameter` truncating the handoff plane (section 3), `genGMAD` having
> four option blocks of which only the last is live, and multi-file EDM4hep
> input to ddsim segfaulting.
>
> **If revisiting this for the CORE**, note two things measured on the halo that
> may not carry over. The core's photons are made mostly in the *weak bends*
> (median 3.0 keV, against the bends' Ec = 2.02 keV), not the final focus — so
> (a) `LCC_v2short` is probably **not** valid for the core as it stands, since
> it would drop the bulk of them, and (b) staging would likely be *worse* for
> the core than for the halo, because far more upstream photons cross the split
> plane and each one pays the geometry-entry cost. Both are worth measuring
> rather than assuming.

Same physics as `LCC_v2`. The only difference is *how* the tracking is split:
the beamline is cut at `QF1BL`, the expensive upstream part is run once, and the
final focus is replayed `N_REPLAY` times to buy statistics where the hard
photons are actually made.

Nothing about the lattice, the masks, the halo distribution or the acceptance
skim changes. A run with `N_REPLAY = 1` is physically the same sample as
`LCC_v2`.

---

## 1. Why split there

The tracked region is 289.12 m, but the photons that matter are made in the last
few metres. The weak bends upstream have `B = 1.46 mT`, so `Ec = 2.02 keV`; the
final-focus quads reach `Ec` of tens to hundreds of keV for an off-axis particle,
because they act as dipoles with `B = B'·x`. Measured on the halo: 78 % of the
accepted photons are above 6 keV, i.e. essentially all of them come from the FF.

| split point | stage 2 covers | efficiency ceiling |
|---|---|---|
| `B0BL` | 145.56 m — 50.8 % of the line | 2× |
| **`QF1BL`** | **6.80 m — 2.4 % of the line** | **≈ 42×** |

So one replay of the FF costs 2.4 % of a full primary, and ~40 replays cost one
extra primary.

The split is exact: stage 1 ends after `DRIFT_38` (z = −9.20 m), stage 2 begins
at `QF1BL`. No gap, no overlap. This is the same `idx_start` the pre-existing
`userfile=4` mode used.

---

## 2. How to run it

```bash
source env.sh                                   # BDSIM stack

# stage 1 -- the expensive part, run once
python run_halo_stage1.py --ngenerate 500000 --seed 12345

# bridge -- replicate the phase space at the split plane
python make_stage2_input.py output_12345.root --report

# stage 2 -- the cheap part, replayed N_REPLAY times.
# Always two BDSIM runs (section 4); --species both (the default) does both
# back to back and writes output_<seed>_ep.root and output_<seed>_gamma.root.
python run_halo_stage2.py --seed 54321

# convert BOTH passes into one sample (key4hep, not the BDSIM stack)
python convert_edm4hep.py --input output_54321_ep.root \
                          --input output_54321_gamma.root --runconfig halo
```

`N_REPLAY` lives in `sample_config.py`, not on the command line.

---

## 3. What crosses the split plane

Stage 1 samples **every species**, not just the positrons:

```
e+       4000  (71.4 %)
gamma    1606  (28.6 %)
```

The photons already produced upstream still have to traverse the final focus —
some are stopped by the QC2 mask — so dropping them would silently remove the
soft part of the spectrum. This is why the stage-2 beam file carries a `pdgid`
column. BDSIM 1.7.7 accepts `pdgid` and `weight` as userfile columns; both were
verified to reach the sampler intact before this config was written.

Format, written **one file per species**:

```
pdgid  x[m]  xp[rad]  y[m]  yp[rad]  Ek[GeV]  t[ns]  weight
```

`t` is the absolute stage-1 arrival time, so the time reference at `QD0AL` stays
that of the full beamline — see section 5. Why per species, and why `Ek`: section 4.

Two more things the stage-1 sampler must get right:

* **`samplerDiameter` is 500 mm for stage 1**, not the 36 mm used at `QD0AL`.
  The stage-1 plane is a phase-space handoff, not an acceptance cut, and
  anything it fails to record is silently lost. At 36 mm the recorded photons
  piled up against |x| = 17.98 mm and 22 % of them were missing.
* Backward-going tracks (`zp <= 0`) cannot be injected into a forward line and
  are dropped, with a count reported. In practice this is zero.

---

## 4. Why stage 2 runs twice — the BDSIM species trap

**BDSIM 1.7.7 converts the energy column using the REFERENCE beam particle's
mass, ignoring the per-line `pdgid`.** With `beam, particle="e+"`:

* an `E[GeV]` column makes every photon below 511 keV come out with negative
  kinetic energy, and BDSIM **discards those lines silently** — no warning, the
  event count is just smaller than the file;
* an `Ek[GeV]` column keeps them but corrupts the energies.

Verified directly: a decade scan of photons from 2 keV to 5 MeV injected with
`particle="e+"` lost every one below 511 keV, and the same scan with
`particle="gamma"` round-tripped all of them exactly.

This is why `make_stage2_input.py` writes `stage2_input_ep.dat` and
`stage2_input_gamma.dat`, and stage 2 is run once per species with a matching
reference particle:

`--species both` (the default) runs the two passes back to back, so this is one
command; `--species e+` / `--species gamma` run them individually if you want to
split them across batch jobs.

The two passes share the same primaries — they are two halves of one sample, not
two samples — so they are converted **together** into a single EDM4hep file with
one `nPrimariesTotal`. Converting them separately and adding the files later
would double-count the denominator.

If this ever regresses, the symptom is a stage-2 BDSIM event count lower than
the line count of its input file. Check that first.

---

## 5. Normalisation — the one thing to get right

Uniform replication: every particle is written `N_REPLAY` times with weight 1.0,
and the 1/N is folded into the primary count rather than carried per particle:

```
nPrimariesTotal = NGENERATE * N_REPLAY
```

That choice is deliberate. Per-particle weights would be more CPU-efficient (the
upstream photons drift deterministically and gain nothing from being replayed),
but they would break the single-global-weight assumption that
`convert_edm4hep.py`, the metadata and `sr_weight()` all rely on. With uniform
replication every downstream consumer keeps working unchanged.

Two counts must not be conflated, and `convert_edm4hep.py` keeps them separate:

| | meaning |
|---|---|
| BDSIM events in the stage-2 file | number of *injected particles* — 5606 × N here, split over the two species passes |
| `nPrimariesTotal` | `NGENERATE * N_REPLAY` — the weight denominator |

The converter defaults to the second and prints a note when they differ. Pass
`--n-primaries` to override.

---

## 6. Time reference

Stage 2's own beamline is 6.8 m long, so its `Model` tree reports a small
`samplerSPosition`. Since stage-2 particles are injected carrying their absolute
stage-1 arrival time, the converter must use the **full** beamline value, taken
from `sample_config.SAMPLER_S_FULL_M = 286.865599` rather than read from the
file. Getting this wrong shifts every photon time by ~930 ns.

---

## 7. How to check you actually gained

**Replaying one particle N times does not give N independent samples.** The
copies share their incoming phase space, so the variance goes as

```
(1 + (N-1)·rho) / (M·N)
```

with `rho` the intra-primary correlation. This is structurally the same
pathology as the LCC_v1 sampler bug (1069 distinct particles reused, ESS 70.8,
26× Poisson scatter) — the difference is that a replay re-rolls the downstream
RNG, so it adds *real* but *sub-linear* information.

| rho | effective stats at N=40 | net gain |
|---|---|---|
| 0 | 40 | 20.8× |
| 0.05 | 13.6 | 7.0× |
| 0.2 | 4.5 | 2.4× |
| 0.5 | 2.0 | 1.0× — nothing |

`rho` is likely non-negligible here, because whether a halo particle radiates
hard in QD0A is set by its amplitude, which is fixed upstream. **Measure it
before trusting the gain**: run the same stage-1 dump at `N_REPLAY = 1` and at
`N_REPLAY = 10`, and compare the per-file scatter of the yield against Poisson.
That is the same diagnostic that showed v2 at 1.19× Poisson and v1 at 26×.

Note also the asymmetry the `--report` flag prints: the ~25 % of the dump that is
photons mostly drifts the last 6.8 m deterministically, so those copies are
near-identical and their effective sample size does **not** grow with N, even
though the raw count does. Treat the soft end of the spectrum accordingly.

---

## 8. Validation status

Closure against a single-stage run on an **identical primary list** (same seed,
`cmp`-identical `inputfile.dat`), 4000 primaries, `N_REPLAY = 5`:

| | photons/primary |
|---|---|
| `LCC_v2` single stage | 0.28225 |
| `LCC_v2` production, 4108 files | 0.28223 |
| `LCC_v2multi`, before the species fix | 0.26355  (−6.45 %) |
| **`LCC_v2multi`, after the species fix** | **0.27525  (−2.48 %, 0.6σ)** |

The single-stage path reproduces production to four decimals, which validates the
reference. The −6.45 % was **not** noise: four independent stage-2 seeds over the
same stage-1 dump scattered by only 0.94 %, which is what identified it as a
systematic and led to the two bugs in sections 3 and 4.

After both fixes the closure is 0.6σ. That is consistent, but it is still only a
4000-primary test — **repeat it with enough stage-1 primaries to push the
statistical error below ~1 % before using this for production**, and confirm the
residual −2.5 % does not persist.
