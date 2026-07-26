#!/usr/bin/env python3
"""
Publication figures for the paper/ Overleaf package.

Fig 1: three exemplary erfc(ln τ) fits, stacked horizontally
Fig 2: 1/τ₅₀ vs ω with floated-intercept linear fits only
Fig 3: slopes (floated intercept) vs sinθ; 0° shown but excluded from ∝sinθ guide

Writes PDF figures into paper/figures/ (Overleaf-ready). Intermediate CSV goes
under analysis/scaling/. Raw data untouched.

Run from repo root:
  python3 analysis/make_paper_figures.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_g2 import SKIP_PTS, ROOT, load_run  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
PAPER_FIG = REPO / "paper" / "figures"
PAPER_FIG.mkdir(parents=True, exist_ok=True)
ANALYSIS_SCALING = Path(__file__).resolve().parent / "scaling"
ANALYSIS_SCALING.mkdir(parents=True, exist_ok=True)

CSV = Path(__file__).resolve().parent / "erfc_fits" / "erfc_fit_parameters.csv"
SLOPES_CSV = ANALYSIS_SCALING / "floated_slopes.csv"
TAU_FIT_MAX = 0.1

# Exemplars: mid-speed at 10°, 15°, 20°
EXEMPLARS = [
    ("test010_6.80s_10deg.txt", "10°", "C1"),
    ("test025_7.81s_15deg.txt", "15°", "C2"),
    ("test06_10.74s.txt", "20°", "C3"),
]


def erfc_sigmoid_ln(tau, A, B, Cw, D):
    return A * erfc((np.log(tau) - B) / Cw) + D


def weighted_line(x, y, dy):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    dy = np.maximum(np.asarray(dy, float), 1e-30)
    w = 1.0 / dy**2
    A = np.column_stack([x, np.ones_like(x)])
    W = np.diag(w)
    AtW = A.T @ W
    cov = np.linalg.inv(AtW @ A)
    params = cov @ (AtW @ y)
    m, b = float(params[0]), float(params[1])
    dm, db = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
    return m, b, dm, db


def weighted_line_through_origin(x, y, dy):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    dy = np.maximum(np.asarray(dy, float), 1e-30)
    w = 1.0 / dy**2
    m = float(np.sum(w * x * y) / np.sum(w * x * x))
    chi2 = float(np.sum(w * (y - m * x) ** 2))
    dof = max(len(x) - 1, 1)
    dm = float(np.sqrt((chi2 / dof) / np.sum(w * x * x)))
    return m, dm


def load_csv_rows():
    rows = []
    with CSV.open() as f:
        for r in csv.DictReader(f):
            if not r["period_s"]:
                continue
            T = float(r["period_s"])
            tau50 = float(r["tau50_ms"]) * 1e-3
            dtau50 = float(r["dtau50_ms"]) * 1e-3
            omega = 2 * np.pi / T
            rows.append(
                {
                    "stem": r["stem"],
                    "angle_deg": int(float(r["angle_deg"])),
                    "period_s": T,
                    "omega": omega,
                    "tau50_s": tau50,
                    "dtau50_s": dtau50,
                    "inv_tau": 1.0 / tau50,
                    "d_inv_tau": dtau50 / tau50**2,
                    "sin_theta": float(np.sin(np.deg2rad(float(r["angle_deg"])))),
                    "beta": float(r["beta"]),
                    "dbeta": float(r["dbeta"]),
                    "Cw": float(r["Cw"]),
                    "dCw": float(r["dCw"]),
                    "D": float(r["D"]),
                    "dD": float(r["dD"]),
                    "A": float(r["A"]),
                    "dA": float(r["dA"]),
                    "B": float(r["B"]),
                    "tau50_ms": float(r["tau50_ms"]),
                    "dtau50_ms": float(r["dtau50_ms"]),
                }
            )
    return rows


def fig_exemplary_fits():
    """Three horizontal panels: data + erfc fit + τ50 band."""
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 11,
            "legend.fontsize": 8,
            "figure.dpi": 150,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), sharey=False, constrained_layout=True)

    # reload params from CSV for annotations
    csv_by_stem = {r["stem"]: r for r in load_csv_rows()}

    for ax, (fname, ang_lab, color) in zip(axes, EXEMPLARS):
        run = load_run(ROOT / fname)
        stem = Path(fname).stem
        meta = csv_by_stem[stem]

        tau = run["tau_s"]
        g2 = run["g2"]
        fit_m = (np.arange(len(tau)) >= SKIP_PTS) & (tau <= TAU_FIT_MAX)
        skip_m = np.arange(len(tau)) < SKIP_PTS

        ax.semilogx(tau[skip_m], g2[skip_m], "o", ms=2.5, color="0.75", zorder=1)
        ax.semilogx(tau[fit_m], g2[fit_m], "o", ms=3.2, color="0.15", zorder=2)

        A, B, Cw, D = meta["A"], meta["B"], meta["Cw"], meta["D"]
        # reconstruct from CSV; B = ln(tau50)
        B = np.log(meta["tau50_s"])
        tau_d = np.logspace(np.log10(tau[fit_m].min()), np.log10(tau[fit_m].max()), 400)
        y_d = erfc_sigmoid_ln(tau_d, A, B, Cw, D)
        ax.semilogx(tau_d, y_d, "-", color=color, lw=2.0, zorder=3)

        t50 = meta["tau50_s"]
        dt50 = meta["dtau50_s"]
        ax.axvline(t50, color=color, ls=":", lw=1.2, alpha=0.9)
        ax.axvspan(t50 - dt50, t50 + dt50, color=color, alpha=0.15, lw=0)
        ax.axhline(1.0, color="0.55", ls="--", lw=0.7)

        lo = float(np.min(g2[fit_m]))
        hi = float(np.max(g2[fit_m]))
        pad = 0.18 * (hi - lo)
        ax.set_ylim(lo - pad, hi + 1.6 * pad)
        ax.set_xlim(5e-7, 0.3)
        ax.set_xlabel(r"$\tau$ (s)")
        ax.set_title(
            rf"{ang_lab},  $T_{{\mathrm{{rot}}}}={meta['period_s']:.2f}\,\mathrm{{s}}$"
            "\n"
            rf"$\tau_{{50}}={meta['tau50_ms']:.2f}\pm{meta['dtau50_ms']:.2f}\,\mathrm{{ms}}$"
            "\n"
            rf"$\beta={meta['beta']:.3f}\pm{meta['dbeta']:.3f}$,  "
            rf"$C_w={meta['Cw']:.2f}\pm{meta['dCw']:.2f}$"
        )
        ax.grid(True, which="both", ls=":", alpha=0.45)

    axes[0].set_ylabel(r"$g_2(\tau)$")
    fig.suptitle(
        r"Exemplary erfc$(\ln\tau)$ fits (first 10 lags omitted as afterpulsing)",
        fontsize=12,
    )
    fig.savefig(PAPER_FIG / "fig_exemplary_fits.pdf")
    plt.close(fig)
    print("wrote fig_exemplary_fits")


def fig_inv_tau_vs_omega():
    rows = load_csv_rows()
    by_ang: dict[int, list] = {}
    for r in rows:
        by_ang.setdefault(r["angle_deg"], []).append(r)

    cmap = {0: "C0", 10: "C1", 15: "C2", 20: "C3"}
    labels = {
        0: r"small angle (nom.\ $0^\circ$)",
        10: r"$10^\circ$",
        15: r"$15^\circ$",
        20: r"$20^\circ$",
    }

    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    slope_info = []

    for ang in sorted(by_ang):
        group = sorted(by_ang[ang], key=lambda r: r["omega"])
        w = np.array([r["omega"] for r in group])
        y = np.array([r["inv_tau"] for r in group])
        dy = np.array([r["d_inv_tau"] for r in group])
        color = cmap[ang]
        ax.errorbar(
            w,
            y,
            yerr=dy,
            fmt="o",
            ms=6.5,
            capsize=2.5,
            color=color,
            label=labels[ang],
            zorder=3,
        )
        m, b, dm, db = weighted_line(w, y, dy)
        xx = np.linspace(0, w.max() * 1.06, 80)
        ax.plot(xx, m * xx + b, "-", color=color, lw=1.7, zorder=2)
        slope_info.append((ang, m, dm, b, db))

    ax.set_xlabel(r"$\omega = 2\pi/T_{\mathrm{rot}}$ (rad/s)")
    ax.set_ylabel(r"$1/\tau_{50}$ (s$^{-1}$)")
    ax.set_title(r"$1/\tau_{50}$ versus $\omega$ (floated-intercept linear fits)")
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(PAPER_FIG / "fig_inv_tau50_vs_omega.pdf")
    plt.close(fig)

    # Intermediate slopes for fig 3 (not part of the Overleaf tree)
    with SLOPES_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["angle_deg", "sin_theta", "slope", "dslope", "intercept", "dintercept"])
        for ang, m, dm, b, db in slope_info:
            w.writerow([ang, np.sin(np.deg2rad(ang)), m, dm, b, db])
    print("wrote fig_inv_tau50_vs_omega", slope_info)


def fig_slopes_vs_sintheta():
    rows = list(csv.DictReader(SLOPES_CSV.open()))
    angs = np.array([int(r["angle_deg"]) for r in rows])
    st = np.array([float(r["sin_theta"]) for r in rows])
    m = np.array([float(r["slope"]) for r in rows])
    dm = np.array([float(r["dslope"]) for r in rows])

    cmap = {0: "C0", 10: "C1", 15: "C2", 20: "C3"}
    labels = {
        0: r"small angle (nom.\ $0^\circ$)",
        10: r"$10^\circ$",
        15: r"$15^\circ$",
        20: r"$20^\circ$",
    }

    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    for i, ang in enumerate(angs):
        ax.errorbar(
            st[i],
            m[i],
            yerr=dm[i],
            fmt="o",
            ms=8,
            capsize=3,
            color=cmap[int(ang)],
            label=labels[int(ang)],
            zorder=3,
        )

    mask = angs > 0
    k, dk = weighted_line_through_origin(st[mask], m[mask], dm[mask])
    mf, bf, dmf, dbf = weighted_line(st[mask], m[mask], dm[mask])
    xx = np.linspace(0, max(st.max(), 0.36), 80)
    ax.plot(
        xx,
        k * xx,
        "k--",
        lw=1.4,
        label=rf"$\propto\sin\theta$ (excl.\ $0^\circ$),  $k={k:.0f}\pm{dk:.0f}$",
        zorder=1,
    )

    ax.set_xlabel(r"$\sin\theta$")
    ax.set_ylabel(r"slope $d(1/\tau_{50})/d\omega$  (floated intercept)")
    ax.set_title(r"Angle dependence of the $\omega$-scaling slope")
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_xlim(left=-0.01)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(PAPER_FIG / "fig_slopes_vs_sintheta.pdf")
    plt.close(fig)
    print(
        f"wrote fig_slopes_vs_sintheta  k={k:.3f}±{dk:.3f}  "
        f"free: m={mf:.3f}±{dmf:.3f}, b={bf:.1f}±{dbf:.1f}"
    )


def main():
    fig_exemplary_fits()
    fig_inv_tau_vs_omega()
    fig_slopes_vs_sintheta()


if __name__ == "__main__":
    main()
