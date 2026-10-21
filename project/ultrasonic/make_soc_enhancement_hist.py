"""
make_soc_enhancement_hist.py -- histogram view of the spin-orbit response at FIXED
(U, T/T_c): what happens to W/W_N when V_SO is increased, point by point.

Data: results/uvso_families.json (dense grids dU=2.5t, dV_SO panels 0/0.2/0.4/0.6/0.8,
T/T_c = 0.3/0.5/0.7/0.9/0.999; scattering-only, eta=+1 signed, hbar*omega=0.7t).
Parameter sets on the correlated-hopping-dominated root (D0_0 < 0.05) or with Tc = 0 are
EXCLUDED from every panel, so the histograms show the SOC response alone, not the
root-switch dip.

Panels
  (a) WHERE the effect lives: median enhancement R(V_SO=0.8)/R(V_SO=0) per U-bin,
      one bar per T/T_c (log scale).  Saturation region (U <= 30): ~1; decay tail:
      up to ~17.
  (b) WHAT happens as V_SO grows: distribution of the per-point enhancement
      R(V_SO)/R(0) over all (U, T/T_c) points, one step-histogram per V_SO level.
      The whole distribution shifts right with V_SO while keeping its mass at ~1
      (the saturated points).
  (c) PER STEP: for each +0.2 increment of V_SO, the fraction of (U, T/T_c) points
      whose ratio increases, and the median multiplicative factor of those that do.

Output: figures/soc_enhancement_histogram.png
Run:    python make_soc_enhancement_hist.py
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FRACS = (0.3, 0.5, 0.7, 0.9, 0.999)
VSOS = (0.0, 0.2, 0.4, 0.6, 0.8)
FLAG = 0.05  # D0_0 threshold: correlated-hopping-dominated root -> excluded


def load():
    with open(os.path.join(HERE, "results", "uvso_families.json")) as fh:
        d = json.load(fh)
    ok = [r for r in d if r["Tc"] > 0 and r.get("D0_0", 1) >= FLAG]
    tab = {}
    for r in ok:
        tab[(r["U"], r["VSO"])] = r
    return tab


def main():
    tab = load()
    Us = sorted({u for (u, v) in tab})
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.6))

    # ---------------- (a) median enhancement per U-bin, bars per T/T_c
    ax = axes[0]
    edges = [15, 30, 50, 75, 100, 150.01]
    xs = np.arange(len(edges) - 1)
    w = 0.16
    for j, f in enumerate(FRACS):
        heights = []
        for b in range(len(edges) - 1):
            vals = []
            for U in Us:
                if not (edges[b] <= U < edges[b + 1]):
                    continue
                a = tab.get((U, 0.0)); c = tab.get((U, 0.8))
                if a is None or c is None or a["ratio"][j] <= 0:
                    continue
                vals.append(c["ratio"][j] / a["ratio"][j])
            heights.append(np.median(vals) if vals else np.nan)
        ax.bar(xs + (j - 2) * w, heights, width=w, color=cmap(0.05 + 0.9 * j / 4),
               label=f"$T/T_c$={f:g}")
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(["15-30", "30-50", "50-75", "75-100", "100-150"])
    ax.set_xlabel("$U/t$ bin")
    ax.set_ylabel("median $R(V_{SO}=0.8)/R(V_{SO}=0)$")
    ax.set_title("(a) where the SOC effect lives\n(saturated plateau vs decay tail)",
                 fontsize=10)
    ax.axhline(1.0, color="k", lw=0.7, ls=":")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(alpha=0.3, which="both")

    # ---------------- (b) distribution of per-point enhancement, per V_SO
    ax = axes[1]
    Eb = np.array([0.9, 1.05, 1.2, 1.4, 1.7, 2.1, 2.6, 3.5, 5.0, 8.0, 12.0, 20.01])
    cols = plt.get_cmap("plasma")(np.linspace(0.1, 0.85, len(VSOS) - 1))
    for i, v in enumerate(VSOS[1:]):
        vals = []
        for U in Us:
            for j in range(len(FRACS)):
                a = tab.get((U, 0.0)); c = tab.get((U, v))
                if a is None or c is None or a["ratio"][j] <= 0:
                    continue
                vals.append(min(c["ratio"][j] / a["ratio"][j], 20.0))
        cnt, _ = np.histogram(vals, bins=Eb)
        cnt = cnt / cnt.sum()
        ax.step(Eb[:-1], cnt, where="post", color=cols[i], lw=1.8,
                label=f"$V_{{SO}}={v:g}$ vs 0")
    ax.set_xscale("log")
    ax.set_xlabel("per-point enhancement $R(V_{SO})/R(0)$")
    ax.set_ylabel("fraction of $(U, T/T_c)$ points")
    ax.set_title("(b) the enhancement distribution\nshifts right as $V_{SO}$ grows",
                 fontsize=10)
    ax.set_xticks([1, 1.2, 1.4, 1.7, 2.1, 2.6, 3.5, 5, 8, 12, 20])
    ax.set_xticklabels(["1", "1.2", "1.4", "1.7", "2.1", "2.6", "3.5", "5", "8", "12", "20"],
                       fontsize=8)
    ax.axvline(1.0, color="k", lw=0.7, ls=":")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    # ---------------- (c) per +0.2 step: fraction increasing, median factor
    ax = axes[2]
    steps = list(zip(VSOS[:-1], VSOS[1:]))
    frac_inc, med_fac = [], []
    for v0, v1 in steps:
        inc, fac = 0, []
        for U in Us:
            for j in range(len(FRACS)):
                a = tab.get((U, v0)); c = tab.get((U, v1))
                if a is None or c is None or a["ratio"][j] <= 0:
                    continue
                inc += c["ratio"][j] > a["ratio"][j]
                if c["ratio"][j] > a["ratio"][j]:
                    fac.append(c["ratio"][j] / a["ratio"][j])
        tot = sum(1 for U in Us for j in range(len(FRACS))
                  if (U, v0) in tab and (U, v1) in tab and tab[(U, v0)]["ratio"][j] > 0)
        frac_inc.append(inc / tot)
        med_fac.append(np.median(fac))
    xs = np.arange(len(steps))
    ax.bar(xs - 0.2, frac_inc, width=0.4, color="tab:blue", label="fraction of points with $R$ increasing")
    ax.bar(xs + 0.2, med_fac, width=0.4, color="tab:orange", label="median factor of the increasers")
    ax.set_xticks(xs)
    ax.set_xticklabels([f"{a:g}$\\to${b:g}" for a, b in steps], fontsize=9)
    ax.set_xlabel("$V_{SO}/t$ increment")
    ax.set_ylabel("fraction / median factor")
    ax.set_title("(c) each $+0.2$ step: almost always up,\nmedian $\\times$1.05-$\\times$1.25 per step",
                 fontsize=10)
    ax.axhline(1.0, color="k", lw=0.7, ls=":")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    fig.suptitle("Spin-orbit response at FIXED $(U, T/T_c)$: enhancement of $W/W_N$ when "
                 "$V_{SO}$ increases (scattering only, $\\hbar\\omega=0.7t$, $\\Delta_t=4.5t$, "
                 "$n=1.885$; root-switch sets excluded)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = os.path.join(HERE, "figures", "soc_enhancement_histogram.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out)
    # console summary
    for v in VSOS[1:]:
        vals = [tab[(U, v)]["ratio"][j] / tab[(U, 0.0)]["ratio"][j]
                for U in Us for j in range(len(FRACS))
                if (U, 0.0) in tab and (U, v) in tab and tab[(U, 0.0)]["ratio"][j] > 0]
        print(f"V_SO={v:g}: median E={np.median(vals):.3f}  p90={np.percentile(vals,90):.2f}  "
              f"max={max(vals):.2f}  frac>1.05={np.mean(np.array(vals)>1.05):.2f}")


if __name__ == "__main__":
    main()
