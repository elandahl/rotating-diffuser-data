#!/usr/bin/env python3
"""
Pilot curve fits on a few runs: compare stretched-exp (free α, α=1, α=2),
erfc-sigmoid in ln(τ), and the ideal jinc² speckle model (plain and with a
Gaussian decorrelation envelope). Always display vs log-τ.

Early lags (afterpulsing) are skipped. Raw data files are not modified.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfc, j1

# Reuse loader from plot_g2.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_g2 import SKIP_PTS, ROOT, load_run  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent / "fits"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Pilot files (20° and 15°)
PILOT = [
    "test04_5.70s.txt",
    "test025_7.81s_15deg.txt",
]

# Drop very late noisy baseline beyond this (still on log-τ display of full curve)
TAU_FIT_MAX = 0.1  # s

LAMBDA = 650e-9  # m, diode laser

# X_1/e for [2 J1(X)/X]^2 = 1/e, from the theory note
X_1E = 1.91499


def stretched_exp(tau, beta, tau_c, alpha, C):
    return C + beta * np.exp(-2.0 * (tau / tau_c) ** alpha)


def stretched_exp_fixed_alpha(alpha):
    def f(tau, beta, tau_c, C):
        return C + beta * np.exp(-2.0 * (tau / tau_c) ** alpha)

    return f


def erfc_sigmoid_ln(tau, A, B, Cw, D):
    # y = A * erfc((ln τ - B) / Cw) + D
    x = np.log(tau)
    return A * erfc((x - B) / Cw) + D


def jinc(X):
    """2*J1(X)/X with the X->0 limit of 1."""
    X = np.asarray(X, dtype=float)
    out = np.ones_like(X)
    nz = np.abs(X) > 1e-12
    out[nz] = 2.0 * j1(X[nz]) / X[nz]
    return out


def jinc2(tau, beta, a, C):
    """Ideal centered rotation, short-delay: X = a*tau, a = pi*D*Omega*sin(theta)/lambda."""
    return C + beta * jinc(a * tau) ** 2


def jinc2_damped(tau, beta, a, tau_b, C):
    """jinc^2 times a Gaussian envelope (off-axis / boiling decorrelation)."""
    return C + beta * jinc(a * tau) ** 2 * np.exp(-((tau / tau_b) ** 2))


def rms(y, yhat):
    return float(np.sqrt(np.mean((y - yhat) ** 2)))


def fit_window(tau, g2):
    m = (np.arange(len(tau)) >= SKIP_PTS) & (tau <= TAU_FIT_MAX) & np.isfinite(g2)
    return tau[m], g2[m], m


def fit_all(tau_fit, g2_fit):
    results = {}

    # --- free α ---
    p0 = [g2_fit[0] - 1.0, np.median(tau_fit), 1.5, 1.0]
    bounds = ([0, tau_fit.min() * 0.1, 0.2, 0.5], [2, tau_fit.max(), 5.0, 1.5])
    try:
        p, cov = curve_fit(
            stretched_exp, tau_fit, g2_fit, p0=p0, bounds=bounds, maxfev=20000
        )
        yhat = stretched_exp(tau_fit, *p)
        results["stretch_free"] = {
            "label": rf"stretch free α  (α={p[2]:.2f}, τ$_c$={p[1]*1e3:.2f} ms)",
            "params": {"beta": p[0], "tau_c": p[1], "alpha": p[2], "C": p[3]},
            "yhat": yhat,
            "rms": rms(g2_fit, yhat),
            "fn": lambda t, pp=p: stretched_exp(t, *pp),
        }
    except Exception as e:
        results["stretch_free"] = {"error": str(e)}

    # --- fixed α = 1, 2 ---
    for alpha in (1.0, 2.0):
        key = f"stretch_a{int(alpha)}"
        f = stretched_exp_fixed_alpha(alpha)
        p0 = [g2_fit[0] - 1.0, np.median(tau_fit), 1.0]
        bounds = ([0, tau_fit.min() * 0.1, 0.5], [2, tau_fit.max(), 1.5])
        try:
            p, cov = curve_fit(f, tau_fit, g2_fit, p0=p0, bounds=bounds, maxfev=20000)
            yhat = f(tau_fit, *p)
            results[key] = {
                "label": rf"stretch α={alpha:g}  (τ$_c$={p[1]*1e3:.2f} ms)",
                "params": {"beta": p[0], "tau_c": p[1], "alpha": alpha, "C": p[2]},
                "yhat": yhat,
                "rms": rms(g2_fit, yhat),
                "fn": lambda t, pp=p, ff=f: ff(t, *pp),
            }
        except Exception as e:
            results[key] = {"error": str(e)}

    # --- erfc sigmoid in ln τ ---
    x = np.log(tau_fit)
    # midpoint guess: where g2 halfway from max to min
    g_hi, g_lo = np.percentile(g2_fit, 90), np.percentile(g2_fit, 10)
    half = 0.5 * (g_hi + g_lo)
    i_half = int(np.argmin(np.abs(g2_fit - half)))
    A0 = 0.5 * (g_hi - g_lo)
    B0 = x[i_half]
    Cw0 = (x.max() - x.min()) / 8
    D0 = g_lo
    p0 = [A0, B0, Cw0, D0]
    bounds = (
        [0, x.min() - 2, 0.05, 0.5],
        [2, x.max() + 2, 5.0, 1.5],
    )
    try:
        p, cov = curve_fit(
            erfc_sigmoid_ln, tau_fit, g2_fit, p0=p0, bounds=bounds, maxfev=20000
        )
        yhat = erfc_sigmoid_ln(tau_fit, *p)
        tau50 = float(np.exp(p[1]))
        results["erfc_ln"] = {
            "label": rf"erfc(ln τ)  (τ$_{{50}}$={tau50*1e3:.2f} ms)",
            "params": {"A": p[0], "B": p[1], "C_width": p[2], "D": p[3], "tau50": tau50},
            "yhat": yhat,
            "rms": rms(g2_fit, yhat),
            "fn": lambda t, pp=p: erfc_sigmoid_ln(t, *pp),
        }
    except Exception as e:
        results["erfc_ln"] = {"error": str(e)}

    # --- ideal jinc^2 (short-delay form) ---
    # a sets the decay scale: tau_1/e = X_1E / a
    tau_guess = np.median(tau_fit)
    a0 = X_1E / max(tau_guess, 1e-9)
    p0 = [g2_fit[0] - 1.0, a0, 1.0]
    bounds = ([0, 1e-2, 0.5], [2, 1e9, 1.5])
    try:
        p, cov = curve_fit(
            jinc2, tau_fit, g2_fit, p0=p0, bounds=bounds, maxfev=40000
        )
        yhat = jinc2(tau_fit, *p)
        tau_1e = X_1E / p[1]
        results["jinc2"] = {
            "label": rf"jinc² ideal  (τ$_{{1/e}}$={tau_1e*1e3:.2f} ms)",
            "params": {"beta": p[0], "a": p[1], "C": p[2], "tau_1e": tau_1e},
            "yhat": yhat,
            "rms": rms(g2_fit, yhat),
            "fn": lambda t, pp=p: jinc2(t, *pp),
        }
    except Exception as e:
        results["jinc2"] = {"error": str(e)}

    # --- jinc^2 with Gaussian decorrelation envelope ---
    p0 = [g2_fit[0] - 1.0, a0, tau_guess, 1.0]
    bounds = ([0, 1e-2, 1e-9, 0.5], [2, 1e9, 1e3, 1.5])
    try:
        p, cov = curve_fit(
            jinc2_damped, tau_fit, g2_fit, p0=p0, bounds=bounds, maxfev=40000
        )
        yhat = jinc2_damped(tau_fit, *p)
        tau_1e = X_1E / p[1]
        results["jinc2_damped"] = {
            "label": rf"jinc²×Gauss  (τ$_{{1/e}}$={tau_1e*1e3:.2f} ms, "
            rf"τ$_b$={p[2]*1e3:.2f} ms)",
            "params": {
                "beta": p[0],
                "a": p[1],
                "tau_b": p[2],
                "C": p[3],
                "tau_1e": tau_1e,
            },
            "yhat": yhat,
            "rms": rms(g2_fit, yhat),
            "fn": lambda t, pp=p: jinc2_damped(t, *pp),
        }
    except Exception as e:
        results["jinc2_damped"] = {"error": str(e)}

    return results


def effective_diameter(a, period_s, angle_deg):
    """Invert a = pi*D*Omega*sin(theta)/lambda for the illuminated diameter D."""
    if period_s is None or not angle_deg:
        return None
    omega = 2.0 * np.pi / period_s
    st = np.sin(np.deg2rad(angle_deg))
    if st <= 0:
        return None
    return a * LAMBDA / (np.pi * omega * st)


def plot_pilot(run: dict, results: dict, tau_fit, g2_fit) -> Path:
    tau, g2 = run["tau_s"], run["g2"]
    colors = {
        "stretch_free": "C3",
        "stretch_a1": "C1",
        "stretch_a2": "C2",
        "erfc_ln": "C4",
        "jinc2": "C0",
        "jinc2_damped": "C5",
    }

    fig, axes = plt.subplots(
        2, 1, figsize=(9.5, 8.0), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]}
    )
    ax, axr = axes

    # data: skipped faint, fit-window solid, late-not-fit muted
    skip_m = np.arange(len(tau)) < SKIP_PTS
    fit_m = (np.arange(len(tau)) >= SKIP_PTS) & (tau <= TAU_FIT_MAX)
    late_m = (np.arange(len(tau)) >= SKIP_PTS) & (tau > TAU_FIT_MAX)

    ax.semilogx(tau[skip_m], g2[skip_m], "o", ms=3, color="0.8", label="skipped (afterpulse)")
    ax.semilogx(tau[fit_m], g2[fit_m], "o", ms=4, color="0.15", label="fit data")
    if np.any(late_m):
        ax.semilogx(tau[late_m], g2[late_m], "o", ms=3, color="0.65", label=r"τ > fit max")

    # dense curve for model overlays
    tau_dense = np.logspace(np.log10(tau_fit.min()), np.log10(tau_fit.max()), 400)
    for key, color in colors.items():
        r = results.get(key, {})
        if "fn" not in r:
            continue
        ax.semilogx(
            tau_dense,
            r["fn"](tau_dense),
            "-",
            lw=1.8,
            color=color,
            label=f"{r['label']}  RMS={r['rms']:.4f}",
        )

    ax.axhline(1.0, color="0.5", ls="--", lw=0.8)
    # scale to the fit region; the afterpulse spike would otherwise dominate
    lo, hi = np.min(g2_fit), np.max(g2_fit)
    pad = 0.12 * (hi - lo)
    ax.set_ylim(lo - pad, hi + 2.2 * pad)
    ax.set_ylabel(r"$g_2(\tau)$")
    period = run["period_s"]
    period_str = f"{period:.2f} s" if period is not None else "?"
    ax.set_title(
        f"Pilot fits — {run['stem']}  |  {run['angle_label']}  |  "
        f"$T_{{\\mathrm{{rot}}}}$={period_str}\n"
        f"$f$={run['f']/1e6:.0f} MHz,  skip first {SKIP_PTS} lags "
        rf"(afterpulsing), fit $\tau\leq${TAU_FIT_MAX:g} s"
    )
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(loc="best", fontsize=8)

    # residuals on fit window
    for key, color in colors.items():
        r = results.get(key, {})
        if "yhat" not in r:
            continue
        axr.semilogx(
            tau_fit,
            g2_fit - r["yhat"],
            ".-",
            ms=3,
            lw=0.8,
            color=color,
            label=key,
        )
    axr.axhline(0, color="0.4", lw=0.8)
    axr.set_xlabel(r"$\tau$ (s)")
    axr.set_ylabel("residual")
    axr.grid(True, which="both", ls=":", alpha=0.5)
    axr.legend(fontsize=7, loc="best", ncol=2)

    fig.tight_layout()
    out = OUT_DIR / f"{run['stem']}_pilot_fits.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def summarize(run, results) -> str:
    lines = [f"### `{run['stem']}` — {run['angle_label']}, T_rot={run['period_s']} s", ""]
    lines.append("| Model | contrast | char. time (ms) | shape | baseline | RMS |")
    lines.append("|---|---|---|---|---|---|")
    order = [
        "stretch_free",
        "stretch_a1",
        "stretch_a2",
        "erfc_ln",
        "jinc2",
        "jinc2_damped",
    ]
    for key in order:
        r = results.get(key, {})
        if "error" in r:
            lines.append(f"| {key} | ERROR: {r['error']} | | | | |")
            continue
        p = r["params"]
        if key.startswith("stretch"):
            lines.append(
                f"| {r['label'].split('  (')[0]} | {p['beta']:.4f} | "
                f"τ_c={p['tau_c']*1e3:.3f} | α={p['alpha']:.3f} | {p['C']:.4f} | "
                f"{r['rms']:.5f} |"
            )
        elif key == "erfc_ln":
            lines.append(
                f"| erfc(ln τ) | A={p['A']:.4f} | τ₅₀={p['tau50']*1e3:.3f} | "
                f"C_w={p['C_width']:.3f} | {p['D']:.4f} | {r['rms']:.5f} |"
            )
        elif key == "jinc2":
            D_eff = effective_diameter(p["a"], run["period_s"], run["angle_deg"])
            d_str = f"D_eff={D_eff*1e3:.3f} mm" if D_eff else "—"
            lines.append(
                f"| jinc² ideal | {p['beta']:.4f} | τ_1/e={p['tau_1e']*1e3:.3f} | "
                f"{d_str} | {p['C']:.4f} | {r['rms']:.5f} |"
            )
        else:
            D_eff = effective_diameter(p["a"], run["period_s"], run["angle_deg"])
            d_str = f"D_eff={D_eff*1e3:.3f} mm" if D_eff else "—"
            lines.append(
                f"| jinc²×Gauss | {p['beta']:.4f} | τ_1/e={p['tau_1e']*1e3:.3f} | "
                f"{d_str}, τ_b={p['tau_b']*1e3:.3f} ms | {p['C']:.4f} | {r['rms']:.5f} |"
            )
    lines.append("")
    return "\n".join(lines)


def main():
    md = [
        "# Pilot curve-fit comparison",
        "",
        "Models on the same log-τ axis:",
        "",
        r"1. Stretched exp (free α): \(g_2=C+\beta\exp[-2(\tau/\tau_c)^\alpha]\)",
        r"2. Same with α fixed to 1 and 2",
        r"3. Erfc-sigmoid in \(\ln\tau\): \(g_2=A\,\mathrm{erfc}((\ln\tau-B)/C_w)+D\)",
        r"4. Ideal jinc²: \(g_2=C+\beta[2J_1(a\tau)/(a\tau)]^2\), "
        r"with \(a=\pi D\Omega\sin\theta/\lambda\)",
        r"5. jinc² × Gaussian envelope \(\exp[-(\tau/\tau_b)^2]\) (off-axis / boiling)",
        "",
        r"D_eff is recovered from the fitted \(a\) using the measured rotation period "
        r"and nominal angle, with λ = 650 nm.",
        "",
        f"Skip first {SKIP_PTS} lags (afterpulsing). Fit window: "
        f"τ ≤ {TAU_FIT_MAX:g} s. Clock f = 12 MHz for these runs.",
        "",
        "Equal weight per lag (lags already ~log-spaced).",
        "",
    ]

    for fname in PILOT:
        path = ROOT / fname
        run = load_run(path)
        tau_fit, g2_fit, _ = fit_window(run["tau_s"], run["g2"])
        results = fit_all(tau_fit, g2_fit)
        out = plot_pilot(run, results, tau_fit, g2_fit)
        print(f"{run['stem']}:")
        for k, r in results.items():
            if "error" in r:
                print(f"  {k}: ERROR {r['error']}")
            else:
                print(f"  {k}: RMS={r['rms']:.5f}  params={r['params']}")
        print(f"  -> {out}")
        md.append(summarize(run, results))
        md.append(f"![](fits/{out.name})")
        md.append("")

    readme = OUT_DIR.parent / "pilot_fits.md"
    readme.write_text("\n".join(md))
    print(f"summary -> {readme}")


if __name__ == "__main__":
    main()
