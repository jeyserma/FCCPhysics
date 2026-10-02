"""
quantum_lifetime.py
====================
Extract the Sands quantum lifetime from one or more halo beam files.

Usage
-----
    # single file
    python quantum_lifetime.py beam0.dat

    # multiple files (glob or explicit list)
    python quantum_lifetime.py beam*.dat
    python quantum_lifetime.py beam0.dat beam1.dat beam2.dat

All files are stacked into a single particle sample before fitting,
giving better statistical precision on z_Lambert and τ_q.

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
import numpy as np
import mpmath as mp
from scipy.optimize import curve_fit, brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse

# ── Machine / lattice parameters ──────────────────────────────────────────────

FCC_CIRCUMFERENCE = 90644.816       # m
C_LIGHT           = 299_792_458.0   # m/s
N_TURNS_DAMPING   = 1297
Tx = N_TURNS_DAMPING * FCC_CIRCUMFERENCE / C_LIGHT   # damping time (s)

ALFX = 12.135715081527367
ALFY = 12.181892577315056
BETX =  945.3283498584578  # m
BETY =  966.7958864309221   # m

EMITX = 7e-10    # m·rad
EMITY = 2.6e-12   # m·rad

HIN_X  = 3.5;   HOUT_X = 9.0
HIN_Y  = 4.0;   HOUT_Y = 60.0

DISPX  =  0.3489035341023329
DISPXP = -0.0033722486271542612

N_BINS = 60   # histogram bins for the fit

# ── Derived ───────────────────────────────────────────────────────────────────

GAMMAX = (1 + ALFX**2) / BETX
GAMMAY = (1 + ALFY**2) / BETY

J_IN_X  = HIN_X**2  * EMITX
J_OUT_X = HOUT_X**2 * EMITX
J_IN_Y  = HIN_Y**2  * EMITY
J_OUT_Y = HOUT_Y**2 * EMITY

T_X = (HOUT_X / HIN_X)**2
T_Y = (HOUT_Y / HIN_Y)**2

SIGMAX = np.sqrt(BETX * EMITX)
SIGMAY = np.sqrt(BETY * EMITY)


# ── I/O ───────────────────────────────────────────────────────────────────────

def load_file(filepath):
    """Load one .dat file and return (Jx, Jy, nx, ny)."""
    data   = np.loadtxt(filepath)
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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--lattice", type=str, help="Lattice", default="LCC_v1")
    parser.add_argument("--runconfig", type=str, help="Config run file", default="halo")
    parser.add_argument("--suffix", type=str, help="Suffix", default="")
    args = parser.parse_args()

    rundir = f"/tmp/bdsim/{args.lattice}/{args.runconfig}{args.suffix}/GMAD/inputfile.dat"
    runpath = [rundir]


    Jx, Jy, nx, ny, n_per_file = load_many(runpath)
    print(f"  Total: {len(Jx)} particles from {len(runpath)} file(s).")

    w_x_gen, w_y_gen, zL_design = design_weights(tau_q_design=300.0)

    print("\nFitting plane x ...")
    rx = fit_plane(Jx, J_IN_X, J_OUT_X, T_X, HIN_X, w_x_gen, "x")

    print("Fitting plane y ...")
    ry = fit_plane(Jy, J_IN_Y, J_OUT_Y, T_Y, HIN_Y, w_y_gen, "y")

    print_results(rx, ry, w_x_gen, w_y_gen, zL_design,
                  len(runpath), n_per_file)

    outpath = f"/home/submit/jaeyserm/public_html/fccee/bdsim/{args.lattice}/{args.runconfig}{args.suffix}/"
    os.system(f"mkdir -p {outpath}")
    os.system(f"cp /home/submit/jaeyserm/public_html/fccee/bdsim/index.php {outpath}")
    make_plot(rx, ry, w_x_gen, w_y_gen, zL_design, len(runpath), f"{outpath}/quantum_lifetime_fit.png")


if __name__ == "__main__":
    main()