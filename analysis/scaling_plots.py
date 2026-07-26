#!/usr/bin/env python3
"""
Scaling analysis from erfc(ln τ) τ₅₀ parameters, following the notes:

  1) 1/τ₅₀ vs ω at each angle (same axes), linear fits → slopes
  2) those slopes vs q (and vs sin θ)

Also produce τ₅₀ vs 1/ω (notes' alternate form).

Uses analysis/erfc_fits/erfc_fit_parameters.csv. Raw data untouched.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
CSV = ROOT / "erfc_fits" / "erfc_fit_parameters.csv"
OUT = ROOT / "scaling"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA = 650e-9  # m


def q_of_angle_deg(theta_deg: float) -> float:
    """Optical scattering vector magnitude (vacuum): q = (4π/λ) sin(θ/2)."""
    return (4.0 * np.pi / LAMBDA) * np.sin(np.deg2rad(theta_deg) / 2.0)


def load_rows():
    rows = []
    with CSV.open() as f:
        for r in csv.DictReader(f):
            if not r["period_s"] or r["period_s"] == "":
                continue  # skip test03 (no period)
            rows.append(
                {
                    "stem": r["stem"],
                    "angle_deg": int(float(r["angle_deg"])),
                    "angle_label": r["angle_label"],
                    "period_s": float(r["period_s"]),
                    "tau50_s": float(r["tau50_ms"]) * 1e-3,
                    "dtau50_s": float(r["dtau50_ms"]) * 1e-3,
                    "beta": float(r["beta"]),
                    "dbeta": float(r["dbeta"]),
                    "Cw": float(r["Cw"]),
                }
            )
    return rows


def weighted_line(x, y, dy):
    """Fit y = m x + b with weights 1/dy^2. Returns m,b,dm,db, yhat."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    dy = np.asarray(dy, float)
    dy = np.maximum(dy, 1e-30)
    w = 1.0 / dy**2
    # design matrix
    A = np.column_stack([x, np.ones_like(x)])
    W = np.diag(w)
    AtW = A.T @ W
    cov = np.linalg.inv(AtW @ A)
    params = cov @ (AtW @ y)
    m, b = float(params[0]), float(params[1])
    dm, db = float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))
    return m, b, dm, db


def weighted_line_through_origin(x, y, dy):
    """Fit y = m x with weights 1/dy^2."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    dy = np.asarray(dy, float)
    dy = np.maximum(dy, 1e-30)
    w = 1.0 / dy**2
    m = float(np.sum(w * x * y) / np.sum(w * x * x))
    # variance: Var(m) = 1 / sum(w x^2)  if model correct; rescale by reduced chi2
    chi2 = float(np.sum(w * (y - m * x) ** 2))
    dof = max(len(x) - 1, 1)
    dm = float(np.sqrt((chi2 / dof) / np.sum(w * x * x)))
    return m, dm


def main():
    rows = load_rows()
    by_ang: dict[int, list] = {}
    for r in rows:
        # derived quantities
        T = r["period_s"]
        omega = 2.0 * np.pi / T  # rad/s
        inv_tau = 1.0 / r["tau50_s"]
        d_inv_tau = r["dtau50_s"] / r["tau50_s"] ** 2
        r.update(
            {
                "omega": omega,
                "f_rot": 1.0 / T,
                "inv_omega": 1.0 / omega,
                "inv_tau": inv_tau,
                "d_inv_tau": d_inv_tau,
                "q": q_of_angle_deg(r["angle_deg"]),
                "sin_theta": np.sin(np.deg2rad(r["angle_deg"])),
            }
        )
        by_ang.setdefault(r["angle_deg"], []).append(r)

    angles = sorted(by_ang)
    # colors
    cmap = {0: "C0", 10: "C1", 15: "C2", 20: "C3"}
    labels = {
        0: "small angle (nom. 0°)",
        10: "10°",
        15: "15°",
        20: "20°",
    }

    slope_rows = []  # per-angle fit results

    # ------------------------------------------------------------------
    # Figure 1: 1/τ₅₀ vs ω  (all angles on one plot)
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.0, 5.8))
    for ang in angles:
        group = sorted(by_ang[ang], key=lambda r: r["omega"])
        w = np.array([r["omega"] for r in group])
        y = np.array([r["inv_tau"] for r in group])
        dy = np.array([r["d_inv_tau"] for r in group])
        color = cmap.get(ang, "k")

        ax.errorbar(
            w,
            y,
            yerr=dy,
            fmt="o",
            ms=7,
            capsize=3,
            color=color,
            label=labels.get(ang, f"{ang}°"),
        )

        # free intercept
        m, b, dm, db = weighted_line(w, y, dy)
        # through origin (notes ideal)
        m0, dm0 = weighted_line_through_origin(w, y, dy)

        xx = np.linspace(0, w.max() * 1.05, 80)
        ax.plot(xx, m * xx + b, "-", color=color, lw=1.5, alpha=0.85)
        ax.plot(xx, m0 * xx, "--", color=color, lw=1.0, alpha=0.55)

        slope_rows.append(
            {
                "angle_deg": ang,
                "label": labels.get(ang, f"{ang}°"),
                "q_m": q_of_angle_deg(ang),
                "sin_theta": float(np.sin(np.deg2rad(ang))),
                "n": len(group),
                "slope_free": m,  # (1/τ) / ω  → dimensionless angle-ish
                "dslope_free": dm,
                "intercept": b,
                "dintercept": db,
                "slope_origin": m0,
                "dslope_origin": dm0,
                # Δθ_c from τ = Δθ_c / ω  ⇒ 1/τ = ω / Δθ_c  ⇒ slope = 1/Δθ_c
                "dtheta_c_rad": 1.0 / m0 if m0 else np.nan,
                "ddtheta_c_rad": (dm0 / m0**2) if m0 else np.nan,
            }
        )

    ax.set_xlabel(r"$\omega = 2\pi / T_{\mathrm{rot}}$ (rad/s)")
    ax.set_ylabel(r"$1/\tau_{50}$ (s$^{-1}$)")
    ax.set_title(
        r"Scaling: $1/\tau_{50}$ vs $\omega$ by angle"
        "\n"
        r"solid = free intercept; dashed = through origin"
    )
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(loc="best", fontsize=9)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    p1 = OUT / "inv_tau50_vs_omega.png"
    fig.savefig(p1, dpi=160)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 2: τ₅₀ vs 1/ω  (notes' alternate)
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8.0, 5.8))
    for ang in angles:
        group = sorted(by_ang[ang], key=lambda r: r["inv_omega"])
        x = np.array([r["inv_omega"] for r in group])
        y = np.array([r["tau50_s"] for r in group])
        dy = np.array([r["dtau50_s"] for r in group])
        color = cmap.get(ang, "k")
        ax.errorbar(x, y, yerr=dy, fmt="o", ms=7, capsize=3, color=color, label=labels[ang])
        m, b, dm, db = weighted_line(x, y, dy)
        xx = np.linspace(0, x.max() * 1.05, 80)
        ax.plot(xx, m * xx + b, "-", color=color, lw=1.5)
        # annotate Δθ_c ≈ m when b≈0 (τ = Δθ_c / ω)
    ax.set_xlabel(r"$1/\omega$ (s/rad)")
    ax.set_ylabel(r"$\tau_{50}$ (s)")
    ax.set_title(r"Alternate form: $\tau_{50}$ vs $1/\omega$ (expect linear through origin)")
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(loc="best", fontsize=9)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    p2 = OUT / "tau50_vs_inv_omega.png"
    fig.savefig(p2, dpi=160)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 3: slopes vs q  and vs sinθ
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8), constrained_layout=True)
    axq, axs = axes

    q = np.array([s["q_m"] for s in slope_rows])
    s0 = np.array([s["slope_origin"] for s in slope_rows])
    ds0 = np.array([s["dslope_origin"] for s in slope_rows])
    sf = np.array([s["slope_free"] for s in slope_rows])
    dsf = np.array([s["dslope_free"] for s in slope_rows])
    st = np.array([s["sin_theta"] for s in slope_rows])
    angs = np.array([s["angle_deg"] for s in slope_rows])

    for i, ang in enumerate(angs):
        c = cmap.get(int(ang), "k")
        axq.errorbar(q[i], s0[i], yerr=ds0[i], fmt="o", ms=8, capsize=3, color=c)
        axq.errorbar(q[i], sf[i], yerr=dsf[i], fmt="s", ms=7, capsize=3, color=c, alpha=0.7)
        axs.errorbar(st[i], s0[i], yerr=ds0[i], fmt="o", ms=8, capsize=3, color=c, label=labels[int(ang)])
        axs.errorbar(st[i], sf[i], yerr=dsf[i], fmt="s", ms=7, capsize=3, color=c, alpha=0.7)

    # guide: slope ∝ sinθ (ideal centered model) — fit excluding 0° for origin-forced slopes
    mask = angs > 0
    if np.count_nonzero(mask) >= 2:
        # s0 = k * sinθ  (through origin in sinθ)
        k, dk = weighted_line_through_origin(st[mask], s0[mask], ds0[mask])
        xx = np.linspace(0, max(st.max(), 0.35), 60)
        axs.plot(xx, k * xx, "k--", lw=1.2, label=rf"∝ $\sin\theta$ (excl. 0°)")
        # also vs q: for small θ, q ≈ (2π/λ)θ ≈ (2π/λ)sinθ roughly related
        # fit s0 = a * q through origin excluding 0
        a_q, da_q = weighted_line_through_origin(q[mask], s0[mask], ds0[mask])
        xq = np.linspace(0, q.max() * 1.05, 60)
        axq.plot(xq, a_q * xq, "k--", lw=1.2, label=rf"∝ $q$ (excl. 0°)")

    axq.set_xlabel(r"$q = (4\pi/\lambda)\sin(\theta/2)$ (m$^{-1}$)")
    axq.set_ylabel(r"slope $d(1/\tau_{50})/d\omega$")
    axq.set_title("Slopes vs scattering vector $q$\n○ through origin  □ free intercept")
    axq.grid(True, ls=":", alpha=0.5)
    axq.legend(fontsize=8)
    axq.set_xlim(left=0)
    axq.set_ylim(bottom=0)

    axs.set_xlabel(r"$\sin\theta$")
    axs.set_ylabel(r"slope $d(1/\tau_{50})/d\omega$")
    axs.set_title(r"Slopes vs $\sin\theta$ (ideal: $\propto\sin\theta$)")
    axs.grid(True, ls=":", alpha=0.5)
    axs.legend(fontsize=8)
    axs.set_xlim(left=0)
    axs.set_ylim(bottom=0)

    p3 = OUT / "slopes_vs_q_and_sintheta.png"
    fig.savefig(p3, dpi=160)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Figure 4: collapse check — G vs ω τ_50  (optional visual from normalized overlays)
    # We'll plot τ50 * ω (should be ~const per angle) as a diagnostic strip
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for ang in angles:
        group = by_ang[ang]
        om = np.array([r["omega"] for r in group])
        t50 = np.array([r["tau50_s"] for r in group])
        dt = np.array([r["dtau50_s"] for r in group])
        prod = om * t50  # ≈ Δθ_c if pure 1/ω
        dprod = om * dt
        ax.errorbar(
            om,
            prod,
            yerr=dprod,
            fmt="o",
            ms=7,
            capsize=3,
            color=cmap[ang],
            label=labels[ang],
        )
        mean = np.average(prod, weights=1 / dprod**2)
        ax.axhline(mean, color=cmap[ang], ls="--", lw=1, alpha=0.6)
    ax.set_xlabel(r"$\omega$ (rad/s)")
    ax.set_ylabel(r"$\omega\,\tau_{50}$ (rad)")
    ax.set_title(
        r"Angular decorrelation scale $\omega\tau_{50}\approx\Delta\theta_c$"
        "\n(should be flat in $\omega$ at fixed angle)"
    )
    ax.grid(True, ls=":", alpha=0.5)
    ax.legend(fontsize=9)
    fig.tight_layout()
    p4 = OUT / "omega_tau50_vs_omega.png"
    fig.savefig(p4, dpi=160)
    plt.close(fig)

    # ------------------------------------------------------------------
    # Write summary tables
    # ------------------------------------------------------------------
    # CSV of slopes
    slope_csv = OUT / "angle_slopes.csv"
    with slope_csv.open("w", newline="") as f:
        fields = [
            "angle_deg",
            "label",
            "q_m",
            "sin_theta",
            "n",
            "slope_origin",
            "dslope_origin",
            "slope_free",
            "dslope_free",
            "intercept",
            "dintercept",
            "dtheta_c_rad",
            "ddtheta_c_rad",
            "dtheta_c_deg",
            "ddtheta_c_deg",
        ]
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in slope_rows:
            s = dict(s)
            s["dtheta_c_deg"] = s["dtheta_c_rad"] * 180 / np.pi
            s["ddtheta_c_deg"] = s["ddtheta_c_rad"] * 180 / np.pi
            w.writerow({k: s.get(k, "") for k in fields})

    md = [
        "# Scaling analysis from erfc τ₅₀",
        "",
        r"Using \(\tau_{50}\) from `erfc(ln τ)` fits. Angular speed "
        r"\(\omega=2\pi/T_{\mathrm{rot}}\). Optical \(q=(4\pi/\lambda)\sin(\theta/2)\) "
        f"with λ = {LAMBDA*1e9:.0f} nm.",
        "",
        r"Notes predict \(\tau_c=\Delta\theta_c(q)/\omega\), so "
        r"\(1/\tau_c\) vs \(\omega\) is linear through the origin with "
        r"slope \(1/\Delta\theta_c\). Ideal centered rotation further gives "
        r"\(\tau^{-1}\propto\Omega\sin\theta\).",
        "",
        "## Plots",
        "",
        f"- `{p1.name}` — \(1/\\tau_{{50}}\) vs \(\\omega\) by angle",
        f"- `{p2.name}` — \(\\tau_{{50}}\) vs \(1/\\omega\)",
        f"- `{p3.name}` — slopes vs \(q\) and vs \(\\sin\\theta\)",
        f"- `{p4.name}` — \(\\omega\\tau_{{50}}\) vs \(\\omega\) (flatness test)",
        "",
        "## Per-angle slopes of \(1/\\tau_{{50}} = m\\,\\omega (+ b)\)",
        "",
        "| Angle | q (µm⁻¹) | sinθ | m (origin) | m (free) | intercept b (s⁻¹) | Δθ_c (deg) |",
        "|---|---|---|---|---|---|---|",
    ]
    for s in slope_rows:
        dth = s["dtheta_c_rad"] * 180 / np.pi
        ddth = s["ddtheta_c_rad"] * 180 / np.pi
        md.append(
            f"| {s['label']} | {s['q_m']*1e-6:.3f} | {s['sin_theta']:.4f} | "
            f"{s['slope_origin']:.3f}±{s['dslope_origin']:.3f} | "
            f"{s['slope_free']:.3f}±{s['dslope_free']:.3f} | "
            f"{s['intercept']:.1f}±{s['dintercept']:.1f} | "
            f"{dth:.3f}±{ddth:.3f} |"
        )
    md += [
        "",
        "Δθ_c = 1/m_origin (characteristic diffuser rotation angle to decorrelate).",
        "",
        f"Numeric slopes: `{slope_csv.name}`",
        "",
    ]
    # fix markdown escaping - use simpler text
    md_path = OUT / "scaling_summary.md"
    # rewrite cleanly without over-escaping
    md_clean = [
        "# Scaling analysis from erfc τ₅₀",
        "",
        "Using τ₅₀ from `erfc(ln τ)` fits. Angular speed ω = 2π/T_rot. "
        f"Optical q = (4π/λ) sin(θ/2) with λ = {LAMBDA*1e9:.0f} nm.",
        "",
        "Notes predict τ_c = Δθ_c(q)/ω, so 1/τ_c vs ω is linear through the origin "
        "with slope 1/Δθ_c. Ideal centered rotation further gives τ⁻¹ ∝ Ω sinθ.",
        "",
        "## Plots",
        "",
        f"- `{p1.name}` — 1/τ₅₀ vs ω by angle (solid = free intercept, dashed = origin)",
        f"- `{p2.name}` — τ₅₀ vs 1/ω",
        f"- `{p3.name}` — slopes vs q and vs sinθ",
        f"- `{p4.name}` — ω·τ₅₀ vs ω (should be flat at fixed angle)",
        "",
        "## Per-angle slopes of 1/τ₅₀ = m·ω (+ b)",
        "",
        "| Angle | q (µm⁻¹) | sinθ | m (origin) | m (free) | intercept b (s⁻¹) | Δθ_c (deg) |",
        "|---|---|---|---|---|---|---|",
    ]
    for s in slope_rows:
        dth = s["dtheta_c_rad"] * 180 / np.pi
        ddth = s["ddtheta_c_rad"] * 180 / np.pi
        md_clean.append(
            f"| {s['label']} | {s['q_m']*1e-6:.3f} | {s['sin_theta']:.4f} | "
            f"{s['slope_origin']:.3f}±{s['dslope_origin']:.3f} | "
            f"{s['slope_free']:.3f}±{s['dslope_free']:.3f} | "
            f"{s['intercept']:.1f}±{s['dintercept']:.1f} | "
            f"{dth:.3f}±{ddth:.3f} |"
        )
    md_clean += [
        "",
        "Δθ_c = 1/m_origin (characteristic diffuser rotation angle to decorrelate).",
        "",
        f"Numeric slopes: `{slope_csv.name}`",
        "",
    ]
    md_path.write_text("\n".join(md_clean))

    print(f"wrote {p1}")
    print(f"wrote {p2}")
    print(f"wrote {p3}")
    print(f"wrote {p4}")
    print(f"wrote {slope_csv}")
    print(f"wrote {md_path}")
    print("\nSlopes:")
    for s in slope_rows:
        print(
            f"  {s['label']:28s}  m0={s['slope_origin']:.3f}±{s['dslope_origin']:.3f}  "
            f"m={s['slope_free']:.3f}±{s['dslope_free']:.3f}  "
            f"b={s['intercept']:.1f}±{s['dintercept']:.1f}  "
            f"Δθ_c={s['dtheta_c_rad']*180/np.pi:.3f}°"
        )


if __name__ == "__main__":
    main()
