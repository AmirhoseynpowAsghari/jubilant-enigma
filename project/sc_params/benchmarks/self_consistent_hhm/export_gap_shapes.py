#!/usr/bin/env python3
"""Export finite-grid HHM gap-shape diagnostics from the self-consistent state CSV.

This is a postprocessor: it never re-solves equations or modifies the base CSV. It exports
per-helicity gap-distribution summaries over the complete Nk x Nk midpoint grid, plus
illustrative high-symmetry-k samples. They are calculated values, not digitized paper data.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "project" / "sc_params"))
from kgrid import build_kgrid  # noqa: E402

SAMPLE_K = (("Gamma", 0.0, 0.0), ("X", math.pi, 0.0),
            ("M", math.pi, math.pi), ("diagonal_mid", math.pi/2, math.pi/2),
            ("axis_mid", math.pi/2, 0.0), ("off_diagonal", math.pi/2, math.pi))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=HERE / "results" / "hhm_self_consistent_grid.csv")
    ap.add_argument("--output", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.input.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"no state rows in {args.input}")

    summary, samples = [], []
    cache = {}
    for r in rows:
        U, dt, vso, nk = float(r["U"]), float(r["Delta_t"]), float(r["V_SO"]), int(r["Nk"])
        key = (U, dt, vso, nk)
        if key not in cache:
            cache[key] = build_kgrid(vso, U, dt, nk)[0]
        eps = cache[key]
        d0, ds = float(r["D0"]), float(r["DS"])
        for h in range(2):
            bar = d0 - ds * eps[:, :, h] / 4.0  # t=1 in the benchmark convention
            summary.append(dict(U=U, Delta_t=dt, n=float(r["n"]), V_SO=vso, Nk=nk,
                                frac=float(r["frac"]), T=float(r["T"]),
                                helicity=("+" if h == 0 else "-"),
                                gapbar_min=float(bar.min()), gapbar_max=float(bar.max()),
                                gapbar_mean=float(bar.mean()), gapbar_mean_abs=float(np.abs(bar).mean()),
                                gapbar_rms=float(np.sqrt(np.mean(bar * bar))),
                                gapbar_negative_fraction=float(np.mean(bar < 0.0))))
        for label, kx, ky in SAMPLE_K:
            eps0 = -2.0 * (math.cos(kx) + math.cos(ky))
            rashba = math.sqrt(math.sin(kx)**2 + math.sin(ky)**2)
            for h, helicity in enumerate(("+", "-")):
                s = 1.0 if h == 0 else -1.0
                energy = eps0 - 2.0 * s * vso * rashba
                samples.append(dict(U=U, Delta_t=dt, n=float(r["n"]), V_SO=vso, Nk=nk,
                                    frac=float(r["frac"]), T=float(r["T"]), k_label=label,
                                    kx=kx, ky=ky, helicity=helicity, epsilon=energy,
                                    gapbar=d0 - ds * energy / 4.0))

    summary_path = args.output / "hhm_gap_shape_grid_summary.csv"
    sample_path = args.output / "hhm_gap_shape_high_symmetry_samples.csv"
    for path, data in ((summary_path, summary), (sample_path, samples)):
        with path.open("w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=list(data[0]), lineterminator="\n")
            wr.writeheader(); wr.writerows(data)
    print(f"Wrote {summary_path} ({len(summary)} rows)")
    print(f"Wrote {sample_path} ({len(samples)} rows)")

    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.4, 4.8))
        groups = sorted({(float(r["V_SO"]), int(r["Nk"])) for r in rows})
        colors = plt.cm.viridis(np.linspace(0.05, 0.95, max(2, len({g[0] for g in groups}))))
        vso_color = {v: colors[i] for i, v in enumerate(sorted({g[0] for g in groups}))}
        styles = {80: "-", 160: "--", 320: ":"}
        for vso, nk in groups:
            subset = [x for x in summary if float(x["V_SO"]) == vso and int(x["Nk"]) == nk
                      and x["helicity"] == "+"]
            subset.sort(key=lambda x: x["frac"])
            ax.fill_between([x["frac"] for x in subset],
                            [x["gapbar_min"] for x in subset],
                            [x["gapbar_max"] for x in subset],
                            color=vso_color[vso], alpha=0.04)
            ax.plot([x["frac"] for x in subset], [x["gapbar_rms"] for x in subset],
                    color=vso_color[vso], linestyle=styles.get(nk, "-"),
                    label=f"VSO={vso:g}, Nk={nk}, helicity + RMS")
        ax.set(xlabel=r"$T/T_c$ (determinant)", ylabel=r"$\bar{\Delta}_{k,+}$ statistics / t",
               title="Calculated gap-shape RMS; shading spans grid min to max")
        ax.grid(alpha=0.25); ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(args.output / "hhm_gap_shape_statistics.png", dpi=180)
        plt.close(fig)
        print(f"Wrote {args.output / 'hhm_gap_shape_statistics.png'}")
    except ImportError as exc:
        print(f"Plot not generated (matplotlib unavailable): {exc}")


if __name__ == "__main__":
    main()
