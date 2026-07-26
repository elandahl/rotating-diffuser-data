#!/usr/bin/env python3
"""
Batch-fit all rotating-diffuser runs with the erfc(ln τ) model.

  g2(τ) = A * erfc( (ln τ - B) / C_w ) + D

  τ_50 = exp(B)
  β ≈ 2A   (step height; high plateau ≈ D+2A)

Skip first SKIP_PTS lags (afterpulsing). Fit τ ≤ TAU_FIT_MAX.
Uncertainties from curve_fit covariance, rescaled by reduced χ².
Raw test*.txt / test*.png are not modified.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_g2 import SKIP_PTS, ROOT, load_run  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent / "erfc_fits"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TAU_FIT_MAX = 0.1  # s


def erfc_sigmoid_ln(tau, A, B, Cw, D):
    return A * erfc((np.log(tau) - B) / Cw) + D


def fit_window(tau, g2):
    m = (np.arange(len(tau)) >= SKIP_PTS) & (tau <= TAU_FIT_MAX) & np.isfinite(g2)
    return tau[m], g2[m], m


def fit_erfc(tau_fit, g2_fit, n_boot: int = 400, rng_seed: int = 0):
    x = np.log(tau_fit)
    g_hi = np.percentile(g2_fit, 90)
    g_lo = np.percentile(g2_fit, 10)
    half = 0.5 * (g_hi + g_lo)
    i_half = int(np.argmin(np.abs(g2_fit - half)))
    p0 = [0.5 * (g_hi - g_lo), x[i_half], (x.max() - x.min()) / 8, g_lo]
    bounds = ([0, x.min() - 2, 0.05, 0.5], [2, x.max() + 2, 5.0, 1.5])

    p, cov = curve_fit(
        erfc_sigmoid_ln, tau_fit, g2_fit, p0=p0, bounds=bounds, maxfev=40000
    )
    yhat = erfc_sigmoid_ln(tau_fit, *p)
    resid = g2_fit - yhat
    n = len(g2_fit)
    npar = 4
    dof = max(n - npar, 1)
    chi2 = float(np.sum(resid**2))
    red_chi2 = chi2 / dof
    rms = float(np.sqrt(np.mean(resid**2)))

    A, B, Cw, D = map(float, p)

    # Residual bootstrap: more honest than formal cov when residuals are systematic
    rng = np.random.default_rng(rng_seed)
    boot = []
    for _ in range(n_boot):
        g_b = yhat + rng.choice(resid, size=n, replace=True)
        try:
            pb, _ = curve_fit(
                erfc_sigmoid_ln,
                tau_fit,
                g_b,
                p0=p,
                bounds=bounds,
                maxfev=10000,
            )
            boot.append(pb)
        except Exception:
            continue
    boot = np.asarray(boot)
    if len(boot) >= 20:
        se = np.std(boot, axis=0, ddof=1)
    else:
        # fallback: formal cov × reduced χ²
        cov = cov * red_chi2
        se = np.sqrt(np.diag(cov))

    dA, dB, dCw, dD = map(float, se)
    tau50 = float(np.exp(B))
    dtau50 = tau50 * dB
    beta = 2.0 * A
    dbeta = 2.0 * dA

    return {
        "A": A,
        "dA": dA,
        "B": B,
        "dB": dB,
        "Cw": Cw,
        "dCw": dCw,
        "D": D,
        "dD": dD,
        "tau50_s": tau50,
        "dtau50_s": dtau50,
        "beta": beta,
        "dbeta": dbeta,
        "rms": rms,
        "chi2": chi2,
        "red_chi2": red_chi2,
        "dof": dof,
        "n_fit": n,
        "n_boot": int(len(boot)),
        "yhat": yhat,
        "resid": resid,
        "p": p,
    }


def plot_fit(run, fit, tau_fit, g2_fit) -> Path:
    tau, g2 = run["tau_s"], run["g2"]
    fig, axes = plt.subplots(
        2, 1, figsize=(8.5, 7.0), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
    )
    ax, axr = axes

    skip_m = np.arange(len(tau)) < SKIP_PTS
    fit_m = (np.arange(len(tau)) >= SKIP_PTS) & (tau <= TAU_FIT_MAX)
    late_m = (np.arange(len(tau)) >= SKIP_PTS) & (tau > TAU_FIT_MAX)

    ax.semilogx(tau[skip_m], g2[skip_m], "o", ms=3, color="0.8", label="skipped (afterpulse)")
    ax.semilogx(tau[fit_m], g2[fit_m], "o", ms=4, color="0.15", label="fit data")
    if np.any(late_m):
        ax.semilogx(tau[late_m], g2[late_m], "o", ms=3, color="0.65", label=r"τ > fit max")

    tau_dense = np.logspace(np.log10(tau_fit.min()), np.log10(tau_fit.max()), 500)
    ax.semilogx(
        tau_dense,
        erfc_sigmoid_ln(tau_dense, *fit["p"]),
        "-",
        color="C3",
        lw=2.0,
            label=(
            rf"erfc(ln τ):  τ$_{{50}}$={fit['tau50_s']*1e3:.3f}±{fit['dtau50_s']*1e3:.3f} ms"
            "\n"
            rf"β={fit['beta']:.4f}±{fit['dbeta']:.4f},  "
            rf"$C_w$={fit['Cw']:.3f}±{fit['dCw']:.3f},  "
            rf"D={fit['D']:.4f}±{fit['dD']:.4f}"
        ),
    )
    ax.axhline(1.0, color="0.5", ls="--", lw=0.8)
    ax.axvline(fit["tau50_s"], color="C3", ls=":", lw=1.0, alpha=0.7)

    lo, hi = float(np.min(g2_fit)), float(np.max(g2_fit))
    pad = 0.15 * max(hi - lo, 0.05)
    ax.set_ylim(lo - pad, hi + 1.8 * pad)

    period = run["period_s"]
    period_str = f"{period:.2f} s" if period is not None else "unknown"
    ax.set_title(
        f"{run['stem']}  |  {run['angle_label']}  |  $T_{{\\mathrm{{rot}}}}$={period_str}\n"
        f"$f$={run['f']/1e6:.0f} MHz, skip {SKIP_PTS} lags, "
        rf"fit $\tau\leq${TAU_FIT_MAX:g} s,  RMS={fit['rms']:.4f}"
    )
    ax.set_ylabel(r"$g_2(\tau)$")
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    axr.semilogx(tau_fit, fit["resid"], "o-", ms=3, lw=0.8, color="C3")
    axr.axhline(0, color="0.4", lw=0.8)
    axr.set_xlabel(r"$\tau$ (s)")
    axr.set_ylabel("residual")
    axr.grid(True, which="both", ls=":", alpha=0.5)

    fig.tight_layout()
    out = OUT_DIR / f"{run['stem']}_erfc.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_overview_tau50(rows: list[dict]) -> Path:
    """τ_50 vs rotation period, one panel per angle."""
    by_ang: dict[int, list] = {}
    for r in rows:
        if r["period_s"] is None:
            continue
        by_ang.setdefault(r["angle_deg"], []).append(r)

    angles = sorted(by_ang)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    axes = axes.ravel()
    for ax, ang in zip(axes, angles):
        group = sorted(by_ang[ang], key=lambda r: r["period_s"])
        T = np.array([r["period_s"] for r in group])
        t50 = np.array([r["tau50_ms"] for r in group])
        dt50 = np.array([r["dtau50_ms"] for r in group])
        ax.errorbar(T, t50, yerr=dt50, fmt="o-", ms=6, capsize=3, color="C0")
        # rough guide: τ ∝ T_rot  (since Ω ∝ 1/T)
        if len(T) >= 2 and np.all(t50 > 0):
            # fit t50 = k * T
            k = np.sum(T * t50) / np.sum(T**2)
            TT = np.linspace(T.min() * 0.9, T.max() * 1.05, 50)
            ax.plot(TT, k * TT, "--", color="0.4", lw=1, label=rf"∝ $T_{{\mathrm{{rot}}}}$")
            ax.legend(fontsize=8)
        label = (
            "small angle (nominal 0°)"
            if ang == 0
            else (f"{ang}° (unlabeled→20°)" if ang == 20 else f"{ang}°")
        )
        ax.set_title(label)
        ax.set_xlabel(r"$T_{\mathrm{rot}}$ (s)")
        ax.set_ylabel(r"$\tau_{50}$ (ms)")
        ax.grid(True, ls=":", alpha=0.5)

    fig.suptitle(r"erfc(ln τ) midpoint $\tau_{50}$ vs rotation period", fontsize=12)
    out = OUT_DIR / "overview_tau50_vs_period.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def plot_overview_overlays(rows_full: list[tuple]) -> Path:
    """Overlay data + fit curves by angle (normalized contrast)."""
    by_ang: dict[int, list] = {}
    for run, fit, tau_fit, g2_fit in rows_full:
        by_ang.setdefault(run["angle_deg"], []).append((run, fit, tau_fit, g2_fit))

    angles = sorted(by_ang)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, constrained_layout=True)
    axes = axes.ravel()
    cmap = plt.cm.viridis

    for ax, ang in zip(axes, angles):
        group = sorted(
            by_ang[ang],
            key=lambda t: (t[0]["period_s"] is None, t[0]["period_s"] or 0),
        )
        periods = [t[0]["period_s"] for t in group if t[0]["period_s"] is not None]
        pmin = min(periods) if periods else 1.0
        pmax = max(periods) if periods else 1.0
        for run, fit, tau_fit, g2_fit in group:
            # normalize: (g2 - D) / (2A) so plateau→1, baseline→0
            A, D = fit["A"], fit["D"]
            gN = (g2_fit - D) / (2 * A) if A > 0 else g2_fit - 1
            yN = (fit["yhat"] - D) / (2 * A) if A > 0 else fit["yhat"] - 1
            if run["period_s"] is None:
                color, lab = "0.4", run["stem"]
            else:
                t = np.log(run["period_s"] / pmin) / max(np.log(pmax / pmin), 1e-12)
                color = cmap(t)
                lab = f"T={run['period_s']:.2f}s"
            ax.semilogx(tau_fit, gN, "o", ms=2.5, color=color, alpha=0.55)
            ax.semilogx(tau_fit, yN, "-", lw=1.4, color=color, label=lab)
        ax.axhline(0, color="0.5", ls="--", lw=0.7)
        ax.axhline(1, color="0.5", ls=":", lw=0.7)
        label = (
            "small angle (nominal 0°)"
            if ang == 0
            else (f"{ang}° (unlabeled→20°)" if ang == 20 else f"{ang}°")
        )
        ax.set_title(label)
        ax.set_xlabel(r"$\tau$ (s)")
        ax.set_ylabel(r"$(g_2-D)/(2A)$")
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(fontsize=6.5, loc="best")

    fig.suptitle(
        r"erfc(ln τ) fits — contrast-normalized overlays by angle",
        fontsize=12,
    )
    out = OUT_DIR / "overview_normalized_overlays.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def write_tables(rows: list[dict], plot_names: dict) -> tuple[Path, Path]:
    csv_path = OUT_DIR / "erfc_fit_parameters.csv"
    fieldnames = [
        "stem",
        "angle_deg",
        "angle_label",
        "period_s",
        "rate_Hz",
        "A",
        "dA",
        "B",
        "dB",
        "Cw",
        "dCw",
        "D",
        "dD",
        "tau50_ms",
        "dtau50_ms",
        "beta",
        "dbeta",
        "rms",
        "red_chi2",
        "n_fit",
        "plot",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    md = [
        "# erfc(ln τ) batch fits",
        "",
        r"Model: \(g_2(\tau)=A\,\mathrm{erfc}\big((\ln\tau-B)/C_w\big)+D\)",
        "",
        r"- \(\tau_{50}=e^{B}\) (midpoint on the log-τ sigmoid)",
        r"- \(\beta=2A\) (full step height; high plateau \(\approx D+2A\))",
        f"- Skip first {SKIP_PTS} lags (afterpulsing); fit window τ ≤ {TAU_FIT_MAX:g} s",
        "- 1σ uncertainties from residual bootstrap (400 resamples); "
        "formal OLS errors are much smaller and ignore systematic residual structure",
        "- `test03`: first Octave dump only; no period in filename",
        "",
        "## Parameters",
        "",
        "| File | Angle | T_rot (s) | τ₅₀ (ms) | β | C_w | D | RMS | Plot |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        T = f"{r['period_s']:.2f}" if r["period_s"] is not None else "—"
        md.append(
            f"| `{r['stem']}` | {r['angle_label']} | {T} | "
            f"{r['tau50_ms']:.3f}±{r['dtau50_ms']:.3f} | "
            f"{r['beta']:.4f}±{r['dbeta']:.4f} | "
            f"{r['Cw']:.3f}±{r['dCw']:.3f} | "
            f"{r['D']:.4f}±{r['dD']:.4f} | "
            f"{r['rms']:.4f} | `{r['plot']}` |"
        )
    md += [
        "",
        "## Overviews",
        "",
        "- `erfc_fits/overview_tau50_vs_period.png` — τ₅₀ vs T_rot",
        "- `erfc_fits/overview_normalized_overlays.png` — contrast-normalized data+fits",
        "",
        "Full numeric table: `erfc_fits/erfc_fit_parameters.csv`",
        "",
    ]
    md_path = Path(__file__).resolve().parent / "erfc_fits_summary.md"
    md_path.write_text("\n".join(md))
    return csv_path, md_path


def main():
    files = sorted(ROOT.glob("test*.txt"))
    rows = []
    rows_full = []
    plot_names = {}

    for path in files:
        run = load_run(path)
        tau_fit, g2_fit, _ = fit_window(run["tau_s"], run["g2"])
        if len(tau_fit) < 8:
            print(f"SKIP {run['stem']}: too few points in fit window")
            continue
        try:
            fit = fit_erfc(tau_fit, g2_fit)
        except Exception as e:
            print(f"FAIL {run['stem']}: {e}")
            continue

        out = plot_fit(run, fit, tau_fit, g2_fit)
        plot_names[run["stem"]] = out.name
        row = {
            "stem": run["stem"],
            "angle_deg": run["angle_deg"],
            "angle_label": run["angle_label"],
            "period_s": run["period_s"],
            "rate_Hz": run["rate_hz"],
            "A": fit["A"],
            "dA": fit["dA"],
            "B": fit["B"],
            "dB": fit["dB"],
            "Cw": fit["Cw"],
            "dCw": fit["dCw"],
            "D": fit["D"],
            "dD": fit["dD"],
            "tau50_ms": fit["tau50_s"] * 1e3,
            "dtau50_ms": fit["dtau50_s"] * 1e3,
            "beta": fit["beta"],
            "dbeta": fit["dbeta"],
            "rms": fit["rms"],
            "red_chi2": fit["red_chi2"],
            "n_fit": fit["n_fit"],
            "plot": out.name,
        }
        rows.append(row)
        rows_full.append((run, fit, tau_fit, g2_fit))
        print(
            f"{run['stem']:30s}  τ50={row['tau50_ms']:.3f}±{row['dtau50_ms']:.3f} ms  "
            f"β={row['beta']:.4f}±{row['dbeta']:.4f}  "
            f"Cw={row['Cw']:.3f}±{row['dCw']:.3f}  RMS={row['rms']:.4f}"
        )

    # sort for table: by angle then period
    rows.sort(
        key=lambda r: (
            r["angle_deg"] if r["angle_deg"] is not None else 999,
            r["period_s"] is None,
            r["period_s"] or 0,
            r["stem"],
        )
    )

    ov1 = plot_overview_tau50(rows)
    ov2 = plot_overview_overlays(rows_full)
    csv_path, md_path = write_tables(rows, plot_names)
    print(f"overview τ50 -> {ov1}")
    print(f"overview overlays -> {ov2}")
    print(f"CSV -> {csv_path}")
    print(f"summary -> {md_path}")
    print(f"{len(rows)} fits in {OUT_DIR}")


if __name__ == "__main__":
    main()
