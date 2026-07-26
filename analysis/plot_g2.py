#!/usr/bin/env python3
"""
Read FPGA/Octave autocorrelation exports (immutable) and write g2 plots
to analysis/. Does not modify raw test*.txt / test*.png files.

Normalization (from fpga-backup log_with_multi_tau.m):
  shortcorrs = coincidences / total_counts
  g2(tau) = shortcorrs * clock_time / total_counts
          = shortcorrs / <n>,  <n> = total_counts / clock_time

If a file contains a duplicated Octave dump (e.g. test03), only the first
occurrence of each variable is used.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = Path(__file__).resolve().parent / "plots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SKIP_PTS = 10  # match acquisition convention; still plot raw early bins faintly


def parse_octave_text(path: Path) -> dict:
    """Parse Octave -text save; keep first occurrence of each name."""
    text = path.read_text()
    blocks = re.split(r"(?=^# name: )", text, flags=re.M)
    data: dict = {}
    for block in blocks:
        m = re.match(r"# name: (\S+)\n(.*)", block, flags=re.S)
        if not m:
            continue
        name, body = m.group(1), m.group(2)
        if name in data:
            continue  # first dump only
        if "# type: scalar" in body:
            nums = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", body)
            # last numeric line is the value (after header)
            lines = [
                ln.strip()
                for ln in body.splitlines()
                if ln.strip() and not ln.strip().startswith("#")
            ]
            data[name] = float(lines[0]) if lines else float("nan")
        elif "# type: matrix" in body:
            rows = int(re.search(r"# rows: (\d+)", body).group(1))
            cols = int(re.search(r"# columns: (\d+)", body).group(1))
            vals = []
            for ln in body.splitlines():
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                for tok in ln.split():
                    if tok.lower() == "nan":
                        vals.append(np.nan)
                    else:
                        vals.append(float(tok))
            arr = np.array(vals, dtype=float)
            if rows * cols == arr.size:
                arr = arr.reshape(rows, cols)
            data[name] = arr
        elif "# type: sq_string" in body:
            # last non-empty non-# line after length is the string
            lines = [
                ln.rstrip("\n")
                for ln in body.splitlines()
                if ln.strip() and not ln.strip().startswith("#")
            ]
            data[name] = lines[-1] if lines else ""
    return data


def parse_filename(stem: str) -> dict:
    """Extract run id, period (s), angle (deg) from filename stem."""
    info = {"stem": stem, "period_s": None, "angle_deg": None, "angle_label": "?"}
    # test027_17.50s_15deg, test04_5.70s, test03, test021_33.59s_00deg
    m = re.match(
        r"(test\d+)(?:_(\d+(?:\.\d+)?)s)?(?:_(\d+)deg)?$",
        stem,
    )
    if not m:
        return info
    info["run"] = m.group(1)
    if m.group(2):
        info["period_s"] = float(m.group(2))
    if m.group(3) is not None:
        ang = int(m.group(3))
        info["angle_deg"] = ang
        if ang == 0:
            info["angle_label"] = "small angle (nominal 0°)"
        else:
            info["angle_label"] = f"{ang}°"
    else:
        # unlabeled early runs = 20°
        info["angle_deg"] = 20
        info["angle_label"] = "20° (unlabeled)"
    return info


def load_run(path: Path) -> dict:
    raw = parse_octave_text(path)
    meta = parse_filename(path.stem)

    ac_mat = np.asarray(raw.get("AC", raw["shortcorrs"]), dtype=float)
    # Multi-angle leftover: rows = angles; single-angle runs are 1xN
    if ac_mat.ndim == 2:
        ac = ac_mat[0, :]
    else:
        ac = ac_mat.ravel()

    tau_cycles = np.asarray(raw["shortdatax"], dtype=float).ravel()
    f = float(raw["f"])
    clock_time = float(raw.get("clock_time", raw.get("clock_times", np.nan)))
    if np.ndim(clock_time) > 0:
        clock_time = float(np.asarray(clock_time).ravel()[0])
    total_counts = float(raw.get("total_counts", raw.get("ACcounts", np.nan)))
    if np.ndim(total_counts) > 0:
        total_counts = float(np.asarray(total_counts).ravel()[0])
    waittime = float(raw.get("waittime", np.nan))

    n_mean = total_counts / clock_time
    g2 = ac / n_mean
    tau_s = tau_cycles / f
    rate_hz = f * n_mean

    return {
        **meta,
        "path": path,
        "tau_s": tau_s,
        "g2": g2,
        "ac": ac,
        "f": f,
        "clock_time": clock_time,
        "total_counts": total_counts,
        "n_mean": n_mean,
        "rate_hz": rate_hz,
        "waittime": waittime,
        "n_lags": len(tau_s),
    }


def plot_run(run: dict) -> Path:
    tau, g2 = run["tau_s"], run["g2"]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))

    # faint early bins (often skipped)
    if len(tau) > SKIP_PTS:
        ax.semilogx(
            tau[:SKIP_PTS],
            g2[:SKIP_PTS],
            "o",
            ms=3,
            color="0.75",
            label=f"lags 1–{SKIP_PTS} (often skipped)",
        )
        ax.semilogx(
            tau[SKIP_PTS:],
            g2[SKIP_PTS:],
            "o-",
            ms=3.5,
            lw=1.0,
            color="C0",
            label=f"lags {SKIP_PTS+1}–{len(tau)}",
        )
    else:
        ax.semilogx(tau, g2, "o-", ms=3.5, lw=1.0, color="C0")

    ax.axhline(1.0, color="0.5", ls="--", lw=0.8, label="g₂ → 1")
    ax.set_xlabel(r"$\tau$ (s)")
    ax.set_ylabel(r"$g_2(\tau)$")
    period = run["period_s"]
    period_str = f"{period:.2f} s" if period is not None else "unknown"
    ax.set_title(
        f"{run['stem']}  |  {run['angle_label']}  |  T_rot = {period_str}\n"
        f"rate ≈ {run['rate_hz']:.3g} Hz,  N = {run['total_counts']:.3g},  "
        f"T_acq ≈ {run['clock_time']/run['f']:.2f} s"
    )
    ax.grid(True, which="both", ls=":", alpha=0.5)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    out = OUT_DIR / f"{run['stem']}_g2.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_overview(runs: list[dict]) -> Path:
    """One panel per angle, all rotation periods overlaid."""
    by_angle: dict[int, list] = {}
    for r in runs:
        by_angle.setdefault(r["angle_deg"], []).append(r)

    angles = sorted(by_angle.keys())
    fig, axes = plt.subplots(
        2, 2, figsize=(11, 8), sharex=True, sharey=True, constrained_layout=True
    )
    axes = axes.ravel()
    cmap = plt.cm.viridis

    for ax, ang in zip(axes, angles):
        group = sorted(
            by_angle[ang],
            key=lambda r: (r["period_s"] is None, r["period_s"] or 0),
        )
        periods = [r["period_s"] for r in group if r["period_s"] is not None]
        pmin, pmax = (min(periods), max(periods)) if periods else (1.0, 1.0)
        for r in group:
            if r["period_s"] is None:
                color = "0.4"
                lab = r["stem"]
            else:
                # color by log period
                t = np.log(r["period_s"] / pmin) / max(np.log(pmax / pmin), 1e-12)
                color = cmap(t)
                lab = f"T={r['period_s']:.2f}s"
            ax.semilogx(
                r["tau_s"][SKIP_PTS:],
                r["g2"][SKIP_PTS:],
                "-",
                lw=1.2,
                color=color,
                label=lab,
            )
        ax.axhline(1.0, color="0.6", ls="--", lw=0.7)
        label = (
            "small angle (nominal 0°)"
            if ang == 0
            else (f"{ang}° (unlabeled→20°)" if ang == 20 else f"{ang}°")
        )
        ax.set_title(label)
        ax.grid(True, which="both", ls=":", alpha=0.4)
        ax.legend(fontsize=7, loc="best")
        ax.set_xlabel(r"$\tau$ (s)")
        ax.set_ylabel(r"$g_2(\tau)$")

    fig.suptitle(
        r"Rotating diffuser $g_2(\tau)$ by scattering angle "
        f"(skip first {SKIP_PTS} lags)",
        fontsize=12,
    )
    out = OUT_DIR / "overview_by_angle_g2.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def write_index(runs: list[dict], plot_paths: list[Path]) -> Path:
    lines = [
        "# Analysis plots — g₂(τ)",
        "",
        "Generated from immutable `test*.txt` exports. Raw files were not modified.",
        "",
        "**Normalization:** `g2 = AC * clock_time / total_counts` "
        "(AC = coincidences/total_counts from `log_with_multi_tau.m`).",
        "",
        f"**Early lags:** first {SKIP_PTS} points shown faint on individual plots; "
        "overview uses lags after skip.",
        "",
        "**test03:** file contains a duplicated Octave dump; only the first dump was used.",
        "",
        "**Angles:** unlabeled filenames treated as 20°. Nominal `00deg` labeled as "
        "small angle (rough alignment).",
        "",
        "## Runs",
        "",
        "| File | Angle | T_rot (s) | Rate (Hz) | N_counts | T_acq (s) | Plot |",
        "|---|---|---|---|---|---|---|",
    ]
    for r, p in zip(runs, plot_paths):
        period = f"{r['period_s']:.2f}" if r["period_s"] is not None else "—"
        lines.append(
            f"| `{r['stem']}.txt` | {r['angle_label']} | {period} | "
            f"{r['rate_hz']:.3g} | {r['total_counts']:.0f} | "
            f"{r['clock_time']/r['f']:.2f} | `{p.name}` |"
        )
    lines += [
        "",
        "## Overview",
        "",
        "- `plots/overview_by_angle_g2.png`",
        "",
    ]
    out = Path(__file__).resolve().parent / "README.md"
    out.write_text("\n".join(lines))
    return out


def main():
    files = sorted(ROOT.glob("test*.txt"))
    runs = []
    plot_paths = []
    for path in files:
        run = load_run(path)
        runs.append(run)
        plot_paths.append(plot_run(run))
        print(
            f"{run['stem']:30s}  ang={run['angle_label']:28s}  "
            f"T={str(run['period_s']):>6}  rate={run['rate_hz']:.3g} Hz  "
            f"g2[skip]={run['g2'][SKIP_PTS]:.3f}  g2[-1]={run['g2'][-1]:.3f}"
        )

    overview = plot_overview(runs)
    print(f"overview -> {overview}")
    index = write_index(runs, plot_paths)
    print(f"index -> {index}")
    print(f"{len(runs)} plots in {OUT_DIR}")


if __name__ == "__main__":
    main()
