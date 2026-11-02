#!/usr/bin/env python3
"""Create strict-only counterparts to the targeted refinement figures.

The strict targeted-refinement dataset was already calculated and is read-only.
The companion SOC finite differences are rebuilt from strict family data, leaving
missing/out-of-band samples as gaps rather than interpolating across them.
"""
from __future__ import annotations
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results", "strict_policy")
FIGURES = os.path.join(HERE, "figures", "strict_policy")


def main():
    os.makedirs(FIGURES, exist_ok=True)
    with open(os.path.join(RESULTS, "refined_extrema_zoom_strict.json"), encoding="utf-8") as f:
        zoom = json.load(f)
    coarse = zoom["coarse_data"]
    strict_roots = sum(float(r.get("Tc") or 0.0) > 0.0 for r in coarse)

    # Preserve the original diagnostic's four-panel structure but do not imply
    # extrema or finite differences where the strict dataset has no SC samples.
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for j, frac in enumerate(zoom["provenance"]["fractions"]):
        ax = axes[0, j]
        ax.set_xlim(15, 40)
        ax.set_title(f"Strict-policy raw $W/W_N(U)$, $T/T_c={frac:g}$")
        ax.set_xlabel("$U/t$"); ax.set_ylabel("$W/W_N$")
        ax.text(0.5, 0.5, f"No in-band superconducting roots\n{len(coarse) - strict_roots}/{len(coarse)} sampled points",
                ha="center", va="center", transform=ax.transAxes, color="tab:red")
        ax = axes[1, j]
        ax.set_xlim(0.60, 0.90)
        ax.set_title(f"Strict-policy $U$ finite differences, $T/T_c={frac:g}$")
        ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel(r"$\Delta(W/W_N)/\Delta(U/t)$")
        ax.text(0.5, 0.5, "No admissible samples; no differences to plot",
                ha="center", va="center", transform=ax.transAxes, color="tab:red")
    fig.suptitle("Targeted strict-policy scan, $15\\leq U/t\\leq40$, $0.60\\leq V_{SO}/t\\leq0.90$\n"
                 "Out-of-band diagnostic roots excluded; no interpolation", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(FIGURES, "refined_extrema_zoom_U.png"), dpi=180)
    plt.close(fig)

    with open(os.path.join(RESULTS, "strict_family_data.json"), encoding="utf-8") as f:
        family = json.load(f)["uvso_families"]
    selected_u = (30.0, 45.0, 75.0, 95.0, 115.0)
    fracs = (0.3, 0.5)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    cmap = plt.get_cmap("plasma")
    for ax, frac in zip(axes, fracs):
        for i, U in enumerate(selected_u):
            rows = sorted((r for r in family if float(r["U"]) == U), key=lambda r: float(r["VSO"]))
            xs, ys = [], []
            for r in rows:
                idx = (0.3, 0.5, 0.7, 0.9, 0.999).index(frac)
                ratio = r.get("ratio", [None] * 5)[idx]
                xs.append(float(r["VSO"]))
                ys.append(np.nan if ratio is None else float(ratio))
            mids = [(xs[k] + xs[k + 1]) / 2 for k in range(len(xs) - 1)]
            deriv = [(ys[k + 1] - ys[k]) / (xs[k + 1] - xs[k])
                     if np.isfinite(ys[k]) and np.isfinite(ys[k + 1]) else np.nan
                     for k in range(len(xs) - 1)]
            ax.plot(mids, deriv, "-o", ms=2.5, lw=1.1, color=cmap(i / 4), label=f"$U={U:g}t$")
        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(f"Strict-policy SOC finite differences, $T/T_c={frac:g}$")
        ax.set_xlabel("$V_{SO}/t$ (cell midpoint)")
        ax.set_ylabel(r"$\Delta(W/W_N)/\Delta(V_{SO}/t)$")
        ax.legend(fontsize=8, ncol=2)
    fig.suptitle("Strict-policy fixed-U SOC profiles; adjacent-cell differences only, gaps preserved")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(FIGURES, "refined_SOC_finite_differences.png"), dpi=180)
    plt.close(fig)
    print(f"wrote strict refinement figures; targeted strict SC roots={strict_roots}/{len(coarse)}")


if __name__ == "__main__":
    main()
