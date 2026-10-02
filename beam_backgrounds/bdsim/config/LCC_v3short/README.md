# LCC_v2short — halo generated at the final focus

Identical physics to `LCC_v2`, but the beamline starts at `QF1BL` (z = −9.20 m)
instead of 289 m upstream, and the halo is generated **there** rather than
tracked to it.

**2.13× faster than `LCC_v2` for the same number of primaries**, with fully
independent trajectories. Exact for detector backgrounds; not for the total
photon flux at the sampler. Both claims are measured below.

---

## 1. Why this is exact, not an approximation

The 282 m it skips does nothing to the positrons but a linear, lossless
transform, so the halo at `QF1BL` is the same distribution you would get by
tracking — and the Courant–Snyder action is invariant under that transport. The
quantum-lifetime tail is a ring-wide property, not something that develops along
the way, so it can be generated at any point from the local optics.

Verified on a 60 000-primary dump taken at the split plane:

| check | result |
|---|---|
| losses in the 282 m | **none** — 60000 / 60000 positrons arrive |
| phase distribution at the plane | **uniform** — χ²/ndf = 1.34 (x), 1.77 (y) |
| vertical tail scale, fitted at the plane | **0.0598** vs generated `yweight` = 0.06 |
| horizontal amplitude range | 3.43–9.96 σ vs generated 3.5–10 σ |

`generate_4d_distribution()` reads its optics from `twiss_file.iloc[idx_start-1]`,
so moving `idx_start` is the whole change: it picks up βx = 3016.6 m,
βy = 5251.6 m at `QF1BL` automatically.

---

## 2. What it drops, and why that is fine

It loses the SR made in the 282 m of weak bends. Those magnets run at
**B = 1.46 mT**, so **Ec = 2.02 keV** and the photons top out near 10 keV.

Measured on the accepted sample: the upstream photons have median **3.1 keV**,
99th percentile **9.4 keV**, and **not one** of 755 is above 300 keV or beyond
2000 µrad. The photons that actually reach the vertex detector have median
**381 keV**. So the dropped component contributes **0 %** of detector hits.

**Use `LCC_v2` instead if you need the total photon flux at the sampler**, where
that soft component is 4.5 % of the count.

---

## 3. Validation against the full line

60 000 short-line primaries vs a 250 000-primary `LCC_v2` reference:

| selection | full line | short line | ratio |
|---|---|---|---|
| all photons | 28.253 % | 26.892 % | 0.952 ± 0.008 |
| **E > 10 keV** | 18.641 % | 18.667 % | **1.001 ± 0.011** |
| E > 50 keV | 7.863 % | 7.870 % | 1.001 ± 0.016 |
| E > 100 keV | 3.832 % | 3.822 % | 0.997 ± 0.023 |
| E > 300 keV | 0.4624 % | 0.4533 % | 0.980 ± 0.066 |
| E > 500 keV | 0.1000 % | 0.0983 % | 0.983 ± 0.142 |
| \|y′\| > 1000 µrad | 0.6936 % | 0.6567 % | 0.947 ± 0.053 |
| \|y′\| > 2000 µrad | 0.0860 % | 0.0883 % | 1.027 ± 0.158 |

The 4.8 % total deficit is exactly the predicted soft component. Above 10 keV
the two agree to 0.1 %.

---

## 4. Cost

Interleaved timing, 60 000 primaries, three rounds (interleaving matters: this
machine runs at load 400+ and identical runs vary by 30 %).

| | time | vs full line |
|---|---|---|
| `LCC_v2` full line | 186 ± 14 s | 1.00 |
| `LCC_v2short` | **87 ± 1 s** | **0.470 ± 0.035** |

**Speedup 2.13× ± 0.16.** Excluding the ~16 s of common setup, the tracking-only
ratio is 0.42.

For the w=1 halo production: the same statistics in ~1930 job-equivalents
instead of 4108, or 2.1× the statistics for the same CPU — which is 1.46× better
precision on detector hit rates, with no ceiling, because every primary is an
independent trajectory.

---

## 5. Compared with the alternatives

| approach | cost | statistics |
|---|---|---|
| `LCC_v2` full line | 1.00 | independent |
| `LCC_v2multi` (stage 1 + K replays), best case K=2 | 1.26 | correlated, **3.25× ceiling** |
| **`LCC_v2short`** | **0.47** | **independent, no ceiling** |

The staged approach was validated as *correct* but is slower than doing nothing,
because stage 1 costs 0.95 of a full run — the 282 m being skipped is nearly the
entire cost — and re-injecting particles costs ~0.5 of a run in geometry entry
alone. See `../LCC_v2multi/README.md`. This config gets the same saving by never
running the 282 m at all.

---

## 6. Time reference

No converter change is needed. The reference particle and the tracked particles
both start at `QF1BL`, so `sampler_s_position` (~6.8 m) and the 2.4 m to the IP
give t0 = 30.7 ns, which is correct for this line. Do **not** apply the
`SAMPLER_S_FULL_M` override used by `LCC_v2multi`; that exists because stage-2
particles there carry absolute times inherited from stage 1.
