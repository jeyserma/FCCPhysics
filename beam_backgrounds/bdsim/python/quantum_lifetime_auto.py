"""
quantum_lifetime_auto.py
========================
Extract the Sands quantum lifetime from a halo beam, taking every lattice
parameter from the run itself instead of hard-coding it.

Same physics and same fit as quantum_lifetime.py; the difference is only where
the numbers come from. That matters because the hard-coded Twiss belongs to the
LCC_v2 generation point (B0DL, betx ~ 945 m); LCC_v2short generates at QF1BL
where betx ~ 3017 m, so the old constants would silently give the wrong action.

Usage
-----
    python quantum_lifetime_auto.py /ceph/submit/data/group/fcc/ee/beam_backgrounds/bdsim/LCC_v2short_halo_standalone_G4SR/tmp/

Takes a run directory, like sr_spectrum.py, and reads GMAD/inputfile.dat from
it. A .dat file or a glob still works.

What is derived, and how
------------------------
beta, alpha   from the beam's own sigma matrix. For any J distribution with
              uniform betatron phase, <x^2> = <J>beta, <xx'> = -<J>alpha,
              <x'^2> = <J>gamma, so sqrt(det sigma) = <J> and beta, alpha follow.
              No need to know which element the beam was generated at.
Dx, Dx'       by regression of x and x' on delta. If the generator applied no
              dispersive orbit -- which LCC_v2 does not, see the comment block in
              generate_4d_distribution -- this comes out ~0 and the correction
              correctly does nothing. The hard-coded version subtracted a
              dispersion that was never put in.
emittances    from the EX / EY header of the TWISS file in GMAD/.
halo bounds   inner from the LCC module (haloNSigmaXInner / YInner), outer from
              the --xtail / --ytail defaults in run_<runconfig>.py. The outer
              bound CANNOT be taken from the data: the tail is exponentially
              suppressed, so the largest sampled amplitude sits far below it
              (ytail = 513 generated, ~51 observed in 60k particles).

Physics
-------
The beam is generated via importance sampling:
    w_i = exp(-wx * Jx/J_inner_x  -  wy * Jy/J_inner_y)
on the annulus  J_inner < J < J_outer  in both planes.

Marginal:  P(t) ∝ exp(-w * t),  t = J/J_inner = (n/n_inner)² ∈ [1, T]
Fit gives w  →  z_Lambert = w * T  →  τ_q = (Tx/2) * exp(z_L) / z_L
"""

import sys
import os
import glob
import re
import numpy as np
import mpmath as mp
from scipy.optimize import curve_fit, brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse

from sample_name import resolve_sample, tried_message

# ── Machine parameters (genuinely fixed) ──────────────────────────────────────

FCC_CIRCUMFERENCE = 90644.816       # m
C_LIGHT           = 299_792_458.0   # m/s
N_TURNS_DAMPING   = 1297
Tx = N_TURNS_DAMPING * FCC_CIRCUMFERENCE / C_LIGHT   # damping time (s)

N_BINS = 60   # histogram bins for the fit


# ── Lattice parameters (filled by derive_parameters) ──────────────────────────
#
# Left as None on purpose: every one of these is a property of the run being
# analysed, not of the script. derive_parameters() sets them before anything
# uses them, and prints what it found.

ALFX = ALFY = BETX = BETY = None
EMITX = EMITY = None
HIN_X = HOUT_X = HIN_Y = HOUT_Y = None
DISPX = DISPXP = None

GAMMAX = GAMMAY = None
J_IN_X = J_OUT_X = J_IN_Y = J_OUT_Y = None
T_X = T_Y = None
SIGMAX = SIGMAY = None

_RAW = {}       # filepath -> loaded array, so the .dat is read only once


def _twiss_from_moments(u, up):
    """
    (beta, alpha, <J>) from the second moments of one plane.

    Works for a halo as well as a gaussian core: the betatron phase is uniform
    either way, so the sigma matrix is <J> * [[beta, -alpha], [-alpha, gamma]]
    and its determinant is <J>^2.
    """
    s11 = np.mean(u * u)
    s12 = np.mean(u * up)
    s22 = np.mean(up * up)
    det = s11 * s22 - s12 * s12
    if det <= 0:
        raise ValueError("non-positive emittance from the beam moments")
    mean_J = np.sqrt(det)
    return s11 / mean_J, -s12 / mean_J, mean_J


def _read_emittances(gmad_dir):
    """EX / EY from the TWISS header of the .tfs sitting in GMAD/."""
    tfs = [f for f in sorted(glob.glob(os.path.join(gmad_dir, "*.tfs")))
           if "aperture" not in os.path.basename(f)]
    if not tfs:
        return None, None, None
    ex = ey = None
    with open(tfs[0]) as fh:
        for line in fh:
            if not line.startswith("@"):
                if line.startswith("*"):
                    break
                continue
            parts = line.split()
            if len(parts) >= 4 and parts[1] == "EX":
                ex = float(parts[-1])
            elif len(parts) >= 4 and parts[1] == "EY":
                ey = float(parts[-1])
    return ex, ey, tfs[0]


def _read_halo_inner(configdir):
    """
    Inner halo bounds from the LCC module (haloNSigmaXInner / YInner).

    Only the INNER bound is read. The outer one is deliberately not taken from
    --xtail / --ytail: those are the edges of the sampling box (ytail = 513),
    not an aperture. See derive_parameters for why that distinction matters.
    """
    inner_x = inner_y = None
    for mod in sorted(glob.glob(os.path.join(configdir, "LCC_*.py"))):
        txt = open(mod).read()
        # Two spellings in use:
        #   LCC_v2*   haloNSigmaXInner, haloNSigmaXOuter = 3.5, self._xtail
        #   LCC_v1p1  haloNSigmaXInner = 3.5
        # The optional ", haloNSigma?Outer" makes one pattern cover both.
        for plane, name in (("X", "inner_x"), ("Y", "inner_y")):
            m = re.search(
                rf"haloNSigma{plane}Inner\s*(?:,\s*haloNSigma{plane}Outer\s*)?=\s*([0-9.]+)",
                txt)
            if m:
                if plane == "X":
                    inner_x = float(m.group(1))
                else:
                    inner_y = float(m.group(1))
        if inner_x is not None:
            break
    return inner_x, inner_y


def derive_parameters(filepaths, configdir, gmad_dir, nout_x, nout_y):
    """Fill the lattice globals from the beam file and the config, and report."""
    global ALFX, ALFY, BETX, BETY, EMITX, EMITY
    global HIN_X, HOUT_X, HIN_Y, HOUT_Y, DISPX, DISPXP
    global GAMMAX, GAMMAY, J_IN_X, J_OUT_X, J_IN_Y, J_OUT_Y, T_X, T_Y, SIGMAX, SIGMAY

    stacked = []
    for fp in filepaths:
        if fp not in _RAW:
            _RAW[fp] = np.loadtxt(fp)
        stacked.append(_RAW[fp])
    data = np.vstack(stacked)
    x, xp, y, yp, E = (data[:, i] for i in range(5))

    # Dispersion by regression on the relative energy deviation, but only
    # applied if it is actually resolved.
    #
    # With sigmaE = 1 MeV at 45.6 GeV the relative spread is 2.2e-5, and the
    # regression's own error is then sigma_x / (sqrt(N) * sigma_delta) ~ 0.22 m
    # -- the same size as the lattice Dx. So on a run like this the fit returns
    # a number that looks like a dispersion and is pure noise, and subtracting
    # it would add a spurious 0.14 %-of-sigma shift rather than remove one.
    # The test below keeps the correction where the energy spread is large
    # enough to constrain it, and drops it where it is not.
    delta = (E - E.mean()) / E.mean()
    sigma_d = float(delta.std())
    n_part = len(x)
    if sigma_d > 0:
        var_d = np.mean(delta * delta)
        dx_fit = float(np.mean(x * delta) / var_d)
        dxp_fit = float(np.mean(xp * delta) / var_d)
        dx_err = float(x.std() / (np.sqrt(n_part) * sigma_d))
        dxp_err = float(xp.std() / (np.sqrt(n_part) * sigma_d))
        significant = abs(dx_fit) > 3.0 * dx_err
    else:
        dx_fit = dxp_fit = dx_err = dxp_err = 0.0
        significant = False

    if significant:
        DISPX, DISPXP = dx_fit, dxp_fit
    else:
        DISPX = DISPXP = 0.0

    x_bet = x - DISPX * delta
    xp_bet = xp - DISPXP * delta

    BETX, ALFX, meanJx = _twiss_from_moments(x_bet, xp_bet)
    BETY, ALFY, meanJy = _twiss_from_moments(y, yp)

    EMITX, EMITY, tfs = _read_emittances(gmad_dir)
    if EMITX is None or EMITY is None:
        sys.exit(f"no EX/EY in a TWISS file under {gmad_dir}; cannot set the "
                 "emittances")

    # Inner bound: a real generation parameter, read from the config.
    # Outer bound: the APERTURE at which the lifetime is quoted, which is a
    # physics choice and not a property of the run. tau_q depends on it
    # exponentially -- for the vertical plane, 50.6 sigma gives 300 s, 60 sigma
    # gives 10600 s, and the 513 sigma sampling box overflows. It is therefore
    # a command-line input, not something derived.
    HIN_X, HIN_Y = _read_halo_inner(configdir)
    HOUT_X, HOUT_Y = nout_x, nout_y
    if HIN_X is None or HIN_Y is None:
        sys.exit(f"could not read haloNSigmaXInner/YInner from {configdir}")

    GAMMAX = (1 + ALFX**2) / BETX
    GAMMAY = (1 + ALFY**2) / BETY
    J_IN_X, J_OUT_X = HIN_X**2 * EMITX, HOUT_X**2 * EMITX
    J_IN_Y, J_OUT_Y = HIN_Y**2 * EMITY, HOUT_Y**2 * EMITY
    T_X = (HOUT_X / HIN_X)**2
    T_Y = (HOUT_Y / HIN_Y)**2
    SIGMAX = np.sqrt(BETX * EMITX)
    SIGMAY = np.sqrt(BETY * EMITY)

    print("Derived from the run:")
    print(f"  beam file(s)   : {len(filepaths)}, {len(data)} particles")
    print(f"  emittances     : ex = {EMITX:.4g}, ey = {EMITY:.4g}   "
          f"({os.path.basename(tfs)})")
    print(f"  twiss x        : betx = {BETX:10.3f} m   alfx = {ALFX:9.4f}   "
          f"sigx = {SIGMAX*1e3:.4f} mm")
    print(f"  twiss y        : bety = {BETY:10.3f} m   alfy = {ALFY:9.4f}   "
          f"sigy = {SIGMAY*1e3:.4f} mm")
    print(f"  energy spread  : sigma_delta = {sigma_d:.3g}")
    if significant:
        print(f"  dispersion     : Dx = {DISPX:.6f} +- {dx_err:.6f} m  "
              f"({abs(dx_fit)/dx_err:.1f} sigma) -- corrected for")
    else:
        print(f"  dispersion     : fit gives {dx_fit:.4f} +- {dx_err:.4f} m "
              f"({abs(dx_fit)/dx_err if dx_err else 0:.1f} sigma) -- NOT resolved,")
        print(f"                   no correction applied. It would have shifted x by "
              f"{abs(dx_fit)*sigma_d*1e3:.4f} mm, "
              f"{100*abs(dx_fit)*sigma_d/x.std():.2f} % of sigma_x.")
    print(f"  halo inner     : x {HIN_X} sigma, y {HIN_Y} sigma   (from the config)")
    print(f"  aperture       : x {HOUT_X} sigma, y {HOUT_Y} sigma   "
          f"(ASSUMED; tau_q depends on it exponentially)")
    print(f"  fit range T    : x {T_X:.1f}, y {T_Y:.1f}")
    nx_obs = np.sqrt((GAMMAX*x_bet**2 + 2*ALFX*x_bet*xp_bet + BETX*xp_bet**2) / EMITX)
    ny_obs = np.sqrt((GAMMAY*y**2 + 2*ALFY*y*yp + BETY*yp**2) / EMITY)
    print(f"  observed n_x   : {nx_obs.min():.2f} - {nx_obs.max():.2f}  "
          f"(inner bound should match {HIN_X})")
    print(f"  observed n_y   : {ny_obs.min():.2f} - {ny_obs.max():.2f}")
    for lbl, obs, hout, hin in (("x", nx_obs, HOUT_X, HIN_X),
                                ("y", ny_obs, HOUT_Y, HIN_Y)):
        frac = ((obs.max() / hin) ** 2 - 1) / ((hout / hin) ** 2 - 1)
        if frac < 0.2:
            print(f"  *** WARNING: in {lbl} the data reaches only {100*frac:.0f} % of "
                  f"the fit range [1, T].")
            print(f"  *** {N_BINS} bins over that range leaves almost everything in "
                  f"bin 1 and the fit degenerates.")
            print(f"  *** Lower --nout-{lbl} towards {obs.max():.0f} sigma, or accept "
                  f"that only the MLE is meaningful.")
    print()


# ── I/O ───────────────────────────────────────────────────────────────────────

def load_file(filepath):
    """Load one .dat file and return (Jx, Jy, nx, ny)."""
    data   = _RAW[filepath] if filepath in _RAW else np.loadtxt(filepath)
    E0     = data[:, 4].mean()
    delta  = (data[:, 4] - E0) / E0
    x_bet  = data[:, 0] - DISPX  * delta
    xp_bet = data[:, 1] - DISPXP * delta
    y_m    = data[:, 2]
    yp_m   = data[:, 3]

    Jx = GAMMAX * x_bet**2 + 2.0 * ALFX * x_bet * xp_bet + BETX * xp_bet**2
    Jy = GAMMAY * y_m**2   + 2.0 * ALFY * y_m   * yp_m   + BETY * yp_m**2

    nx = x_bet / SIGMAX
    ny = y_m   / SIGMAY
    return Jx, Jy, nx, ny


def load_many(filepaths):
    """
    Load and stack particles from multiple .dat files.
    Returns concatenated (Jx, Jy, nx, ny, n_per_file).
    """
    Jx_list, Jy_list, nx_list, ny_list = [], [], [], []
    n_per_file = []
    for fp in filepaths:
        Jx, Jy, nx, ny = load_file(fp)
        Jx_list.append(Jx);  Jy_list.append(Jy)
        nx_list.append(nx);  ny_list.append(ny)
        n_per_file.append(len(Jx))
        print(f"    {fp:40s}  {len(Jx):>8} particles")
    return (np.concatenate(Jx_list), np.concatenate(Jy_list),
            np.concatenate(nx_list), np.concatenate(ny_list),
            n_per_file)


def resolve_paths(args):
    """
    Expand each argument into a list of .dat files.

    Handles three cases:
      1. Plain file path           beam0.dat
      2. Shell glob (pre-expanded) beam0.dat beam1.dat ...
      3. Quoted glob pattern       "path/to/dir/*.dat"
         (shell did NOT expand because of the quotes)

    For case 3 the pattern is split into directory + fnmatch pattern
    and matched via os.listdir, so it works on network/CIFS mounts
    where glob.glob can silently return nothing.
    """
    import fnmatch

    paths = []
    for arg in args:
        # 1. try stdlib glob first (works for local paths)
        expanded = sorted(glob.glob(arg, recursive=True))
        if expanded:
            paths.extend(expanded)
            continue

        # 2. plain existing file
        if os.path.isfile(arg):
            paths.append(arg)
            continue

        # 3. quoted / unresolved pattern — split into dir + fnmatch pattern
        #    e.g. "/eos/user/.../jaw/*.dat" -> dir="/eos/.../jaw", pat="*.dat"
        directory = os.path.dirname(arg) or "."
        pattern   = os.path.basename(arg)

        if os.path.isdir(directory):
            matched = sorted(
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if fnmatch.fnmatch(f, pattern)
                and os.path.isfile(os.path.join(directory, f))
            )
            if matched:
                paths.extend(matched)
            else:
                print(f"  Warning: no files match '{pattern}' in '{directory}' — skipping.")
        else:
            print(f"  Warning: directory not found: '{directory}' — skipping.")

    # deduplicate while preserving order
    seen, unique = set(), []
    for p in paths:
        ap = os.path.abspath(p)
        if ap not in seen:
            seen.add(ap)
            unique.append(p)
    return unique


# ── Physics helpers ───────────────────────────────────────────────────────────

def expected_t(w, T):
    """E[t] for P(t) ∝ exp(-w·t) on [1, T]."""
    A = np.exp(-w);  B = np.exp(-w * T)
    d = A - B
    if abs(d) < 1e-280:
        return 1.0 + 1.0 / w
    return (A - T * B) / d + 1.0 / w


def mle_w(t_arr, T):
    """MLE of w: solve E[t](w) = mean(t)."""
    mean_t = t_arr.mean()
    def eq(w): return expected_t(w, T) - mean_t
    lo = 1e-2
    for _ in range(25):
        if eq(lo) > 0:
            break
        lo /= 5.0
    return brentq(eq, lo, 500.0, xtol=1e-13)


def tq_from_zL(zL):
    return (Tx / 2.0) * np.exp(zL) / zL


def design_weights(tau_q_design=300.0):
    """Generation weights and z_Lambert for the design τ_q."""
    arg = -Tx / (2.0 * tau_q_design)
    zL  = float((-mp.lambertw(arg, -1)).real)
    return (HIN_X / HOUT_X)**2 * zL, (HIN_Y / HOUT_Y)**2 * zL, zL


# ── Fit ───────────────────────────────────────────────────────────────────────

def fit_plane(Jarr, J_in, J_out, T, n_inner, w_guess, label):
    """
    Fit P(t) ∝ exp(-w·t) to the halo action histogram.
    x-axis stored in sigma units  n = n_inner * sqrt(t).
    """
    mask = (Jarr >= J_in) & (Jarr <= J_out)
    t    = Jarr[mask] / J_in
    N    = mask.sum()

    counts, edges = np.histogram(t, bins=N_BINS, range=(1.0, T))
    t_centers = 0.5 * (edges[:-1] + edges[1:])
    n_centers = n_inner * np.sqrt(t_centers)
    errors    = np.sqrt(counts.astype(float));  errors[errors == 0] = 1.0
    good      = counts > 0

    def model(tc, A, w): return A * np.exp(-w * tc)

    popt, pcov = curve_fit(
        model,
        t_centers[good], counts[good].astype(float),
        sigma=errors[good], absolute_sigma=True,
        p0=[counts[good][0] * np.exp(t_centers[good][0] * w_guess), w_guess],
        bounds=([0.0, 0.0], [np.inf, 500.0]),
    )
    A_fit, w_fit = popt
    w_err        = np.sqrt(pcov[1, 1])

    zL     = w_fit * T;  zL_err = w_err * T
    tq     = tq_from_zL(zL)
    dtq    = abs((Tx / 2.0) * np.exp(zL) * (zL - 1.0) / zL**2) * zL_err

    fitted = model(t_centers[good], A_fit, w_fit)
    chi2   = float(np.sum(((counts[good] - fitted) / errors[good])**2))
    ndof   = int(good.sum()) - 2

    w_mle  = mle_w(t, T)
    tq_mle = tq_from_zL(w_mle * T)

    return dict(
        label=label, N=N, t=t,
        counts=counts, t_centers=t_centers, n_centers=n_centers,
        errors=errors, good=good,
        A_fit=A_fit, w_fit=w_fit, w_err=w_err,
        zL=zL, zL_err=zL_err, tq=tq, tq_err=dtq,
        chi2=chi2, ndof=ndof,
        w_mle=w_mle, tq_mle=tq_mle,
        T=T, J_in=J_in, J_out=J_out, n_inner=n_inner,
    )


# ── Output ────────────────────────────────────────────────────────────────────

def print_results(rx, ry, w_x_gen, w_y_gen, zL_design, n_files, n_per_file):
    print()
    print("=" * 65)
    print("  SANDS QUANTUM LIFETIME — FIT RESULTS")
    print("=" * 65)
    print(f"  Files used       : {n_files}  ({sum(n_per_file)} particles total)")
    print(f"  Tx (damping)     = {Tx:.6f} s")
    print(f"  z_L (design)     = {zL_design:.6f}  → τ_q = 300 s = 5 min")
    print()
    for r, w_gen, n_in, n_out in [
        (rx, w_x_gen, HIN_X, HOUT_X),
        (ry, w_y_gen, HIN_Y, HOUT_Y),
    ]:
        lbl = r["label"]
        print(f"  Plane {lbl}  [n_in={n_in}σ, n_out={n_out}σ, T={r['T']:.4f}]:")
        print(f"    N in annulus   : {r['N']}")
        print(f"    w  (fit)       : {r['w_fit']:.6f} ± {r['w_err']:.6f}")
        print(f"    w  (MLE)       : {r['w_mle']:.6f}")
        print(f"    w  (gen)       : {w_gen:.6f}")
        print(f"    z_Lambert      : {r['zL']:.6f} ± {r['zL_err']:.6f}"
              f"  (design {zL_design:.6f})")
        print(f"    τ_q  (fit)     : {r['tq']:.1f} ± {r['tq_err']:.1f} s"
              f"  =  {r['tq']/60:.3f} ± {r['tq_err']/60:.3f} min")
        print(f"    τ_q  (MLE)     : {r['tq_mle']:.1f} s = {r['tq_mle']/60:.3f} min")
        print(f"    τ_q  (design)  : 300.0 s = 5.000 min")
        print(f"    χ²/ndof        : {r['chi2']:.1f} / {r['ndof']}")
        print()
    print("  Note: τ_q ∝ exp(z_L)/z_L — exponentially sensitive to z_L.")
    print(f"  Statistical uncertainty improves as 1/√(N_files).")
    print("=" * 65)


def make_plot(rx, ry, w_x_gen, w_y_gen, zL_design, n_files, outpath):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.patch.set_facecolor("white")

    plane_data = [
        (axes[0], rx, w_x_gen, HIN_X, HOUT_X,
         r"$n_x = x_\beta\,/\,\sigma_x$"),
        (axes[1], ry, w_y_gen, HIN_Y, HOUT_Y,
         r"$n_y = y\,/\,\sigma_y$"),
    ]

    for ax, r, w_gen, n_in, n_out, xlabel in plane_data:
        good = r["good"]

        # data
        ax.errorbar(
            r["n_centers"][good], r["counts"][good], yerr=r["errors"][good],
            fmt="o", ms=4, color="steelblue", ecolor="steelblue",
            elinewidth=1, capsize=2, label="Data (halo)", zorder=3,
        )

        # fit
        n_fine = np.linspace(n_in, n_out * 1.02, 800)
        t_fine = (n_fine / n_in)**2
        ax.plot(
            n_fine, r["A_fit"] * np.exp(-r["w_fit"] * t_fine),
            color="tomato", lw=2.0,
            label=rf"Fit  $w = {r['w_fit']:.4f} \pm {r['w_err']:.4f}$",
        )

        # generation reference
        A_gen = r["A_fit"] * np.exp(-(w_gen - r["w_fit"]) * 1.0)
        ax.plot(
            n_fine, A_gen * np.exp(-w_gen * t_fine),
            color="goldenrod", lw=1.5, ls="--",
            label=rf"Generation  $w_{{\rm gen}} = {w_gen:.4f}$",
        )

        # boundary lines
        ax.axvline(n_in,  color="seagreen",  lw=1.2, ls=":",
                   label=rf"$n_{{\rm in}} = {n_in}\,\sigma$")
        ax.axvline(n_out, color="darkorange", lw=1.2, ls=":",
                   label=rf"$n_{{\rm out}} = {n_out}\,\sigma$")

        ax.set_yscale("log")
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel("Counts / bin", fontsize=11)
        ax.set_title(
            rf"Plane {r['label']}  —  "
            rf"$z_L = {r['zL']:.4f} \pm {r['zL_err']:.4f}$"
            rf"  ($z_{{L,\rm des}} = {zL_design:.4f}$)"
            "\n"
            rf"$\tau_q = {r['tq']:.0f} \pm {r['tq_err']:.0f}$ s"
            rf" $= {r['tq']/60:.2f} \pm {r['tq_err']/60:.2f}$ min"
            rf"   $\chi^2/\rm{{ndof}} = {r['chi2']:.1f}/{r['ndof']}$",
            fontsize=10,
        )
        ax.legend(fontsize=8.5, framealpha=0.85)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.set_facecolor("#f9f9f9")

    n_lbl = f"{n_files} file{'s' if n_files > 1 else ''}"
    fig.suptitle(
        r"Sands Quantum Lifetime  —  $P(t) \propto e^{-w\,t}$,  "
        r"$t = J/J_{\rm inner} = (n/n_{\rm in})^2$,  "
        r"$z_L = w \times T$,  $\tau_q = (T_x/2)\,e^{z_L}/z_L$"
        f"\nFCC-ee Z-pole  |  $T_x = {Tx*1e3:.1f}$ ms  |  "
        f"{n_lbl}  |  Design: 300 s = 5 min",
        fontsize=10, fontweight="bold",
    )
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"\n  Plot saved → {outpath}")


# ── Entry point ───────────────────────────────────────────────────────────────

def resolve_tag(path, configdir=None):
    """
    (tag, runconfig, config directory, sample name) from the run path; the
    parsing lives in sample_name.py, shared with the sr_spectrum scripts.
    """
    res = resolve_sample(path, configdir)
    if res is None:
        sys.exit(tried_message(path) + "\n  or pass --configdir")
    return res


def find_beam(path):
    """
    Accept a run directory, a .dat file, or a glob.

    A directory is resolved to <dir>/GMAD/inputfile.dat, which is what
    genGMAD() writes and what BDSIM was actually fed -- the same convention
    sr_spectrum.py uses for output_*.root.
    """
    if os.path.isdir(path):
        cand = os.path.join(path, "GMAD", "inputfile.dat")
        if not os.path.exists(cand):
            sys.exit(f"no GMAD/inputfile.dat under {path}.\n"
                     "That file is written by genGMAD() for the userfile halo "
                     "(userfile=3); a gaussian-core run does not produce one.")
        return [cand], os.path.join(path, "GMAD")
    hits = resolve_paths([path])
    if not hits:
        sys.exit(f"no beam file matched {path}")
    return hits, os.path.dirname(hits[0])


def main():

    parser = argparse.ArgumentParser(description=__doc__.split("Usage")[0].strip())
    parser.add_argument("input",
                        help="run directory (e.g. /ceph/submit/data/group/fcc/ee/beam_backgrounds/bdsim/LCC_v2short_halo_standalone_G4SR/tmp/), "
                             "or a .dat file / glob")
    parser.add_argument("--configdir", default=None,
                        help="config/<lattice>/ holding the LCC module and "
                             "run_<runconfig>.py (default: inferred from the path)")
    parser.add_argument("--nout-x", type=float, default=9.0,
                        help="horizontal aperture in sigma at which to quote "
                             "tau_q (default 9.0). NOT the --xtail sampling box.")
    parser.add_argument("--nout-y", type=float, default=60.0,
                        help="vertical aperture in sigma at which to quote tau_q "
                             "(default 60.0). NOT the --ytail sampling box, which "
                             "is 513 and would put every particle in one bin.")
    parser.add_argument("-o", "--output",
                        default="/home/submit/jaeyserm/public_html/fccee/bdsim//",
                        help="output base directory")
    args = parser.parse_args()

    runpath, gmad_dir = find_beam(args.input)

    tag, runconfig, configdir, name = resolve_tag(args.input, args.configdir)

    print(f"run dir    : {args.input}")
    print(f"tag        : {tag}   runconfig: {runconfig}   sample: {name}")
    print(f"config     : {configdir}")
    print()

    derive_parameters(runpath, configdir, gmad_dir, args.nout_x, args.nout_y)

    Jx, Jy, nx, ny, n_per_file = load_many(runpath)
    print(f"  Total: {len(Jx)} particles from {len(runpath)} file(s).")

    w_x_gen, w_y_gen, zL_design = design_weights(tau_q_design=300.0)

    print("\nFitting plane x ...")
    rx = fit_plane(Jx, J_IN_X, J_OUT_X, T_X, HIN_X, w_x_gen, "x")

    print("Fitting plane y ...")
    ry = fit_plane(Jy, J_IN_Y, J_OUT_Y, T_Y, HIN_Y, w_y_gen, "y")

    print_results(rx, ry, w_x_gen, w_y_gen, zL_design,
                  len(runpath), n_per_file)

    # <output>/<sample name>/ , matching the production directory naming on ceph
    # (e.g. LCC_v2short_halo, LCC_v2short_halo_standalone_G4SR) so samples made
    # with different BDSIM builds do not overwrite each other's plot.
    outpath = f"{args.output}/{name}/"
    os.makedirs(outpath, exist_ok=True)
    make_plot(rx, ry, w_x_gen, w_y_gen, zL_design, len(runpath),
              f"{outpath}/quantum_lifetime.png")


if __name__ == "__main__":
    main()