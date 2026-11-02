#!/usr/bin/env python3
"""Regenerate the targeted low-U/high-SOC scan under the strict mu-band policy.

This is deliberately separate from refine_extrema_zoom.py: it saves strict-policy outputs
under new names and leaves the existing relaxed-policy results and figures untouched.
The strict test uses the analytic continuum Rashba-band edges implemented in sc_params.solvers.

Run from project/ultrasonic:
    PYTHONPATH=../sc_params python3 strict_policy_zoom.py
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import os
import sys
import time
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import refine_extrema_zoom as ref

OUT_DIR = os.path.join(HERE, "results", "strict_policy")
FIG_DIR = os.path.join(HERE, "figures", "strict_policy")


def _set_strict_policy():
    ref.s1sol.MU_WINDOW_POLICY = "strict"


def _evaluate_strict(task):
    _set_strict_policy()
    row = ref._evaluate(task)
    row["policy"] = "strict"
    return row


def _evaluate_strict_high(task):
    _set_strict_policy()
    row = ref._evaluate_high(task)
    row["policy"] = "strict"
    return row


def _strict_figure(records):
    os.makedirs(FIG_DIR, exist_ok=True)
    plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.25, "font.size": 9})
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.6))

    accepted = [r for r in records if float(r.get("Tc") or 0.0) > 0.0]
    rejected = [r for r in records if r.get("no_inband_root")]
    other = [r for r in records if r not in accepted and r not in rejected]
    if accepted:
        axes[0].scatter([r["U"] for r in accepted], [r["VSO"] for r in accepted],
                        c="seagreen", marker="o", s=24, label="strict SC solution found")
    if rejected:
        axes[0].scatter([r["U"] for r in rejected], [r["VSO"] for r in rejected],
                        c="firebrick", marker="x", s=24,
                        label="no in-band T=0 root; relaxed diagnostic retained")
    if other:
        axes[0].scatter([r["U"] for r in other], [r["VSO"] for r in other],
                        c="0.45", marker=".", s=18, label="other/no-Tc outcome")
    axes[0].set_xlim(min(ref.U_GRID) - 0.5, max(ref.U_GRID) + 0.5)
    axes[0].set_ylim(min(ref.VSO_GRID) - 0.01, max(ref.VSO_GRID) + 0.01)
    axes[0].set_xlabel("$U/t$")
    axes[0].set_ylabel("$V_{SO}/t$")
    axes[0].set_title("Strict-policy admissibility map")
    axes[0].legend(fontsize=8, loc="best")

    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, 6))
    ax = axes[1]
    # Overlay both requested reduced temperatures, using distinct line styles.
    for f in ref.FRACS:
        for v, color in zip((0.6, 0.7, 0.75, 0.8, 0.85, 0.9), colors):
            rows = sorted((r for r in records if abs(r["VSO"] - v) < 1e-9),
                          key=lambda r: r["U"])
            vals = [r.get("frac", {}).get(f"{f:g}", {}).get("ratio") for r in rows]
            valid = [(r["U"], y) for r, y in zip(rows, vals) if y is not None]
            if valid:
                ax.plot([x for x, _ in valid], [y for _, y in valid],
                        color=color, lw=1.2, ls="-" if f == ref.FRACS[0] else "--",
                        label=f"$V_{{SO}}={v:g}t$, $T/T_c={f:g}$")
    axes[1].set_xlim(min(ref.U_GRID), max(ref.U_GRID))
    axes[1].set_xlabel("$U/t$")
    axes[1].set_ylabel("$W/W_N$")
    axes[1].set_title("$W/W_N$ only where a strict SC solution exists")
    if not accepted:
        axes[1].text(0.5, 0.5,
                     "No strict in-band superconducting points in this zoom window",
                     ha="center", va="center", transform=axes[1].transAxes,
                     color="firebrick", fontsize=11, wrap=True)
    handles, labels = axes[1].get_legend_handles_labels()
    if handles:
        axes[1].legend(fontsize=7, ncol=2, loc="best")
    fig.suptitle("Strict chemical-potential policy; analytic bare-band edges; "
                 "out-of-band relaxed roots are excluded from physical curves")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = os.path.join(FIG_DIR, "strict_policy_zoom.png")
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main():
    _set_strict_policy()
    tasks = [(u, v) for v in ref.VSO_GRID for u in ref.U_GRID]
    start = time.time()
    with Pool(ref.NPROC, initializer=_set_strict_policy) as pool:
        records = list(pool.imap_unordered(_evaluate_strict, tasks, chunksize=4))
    records.sort(key=lambda r: (r["VSO"], r["U"]))
    if any(r.get("error") for r in records):
        raise RuntimeError(f"coarse strict scan had {sum(bool(r.get('error')) for r in records)} errors")

    summary, cutoffs = ref.summarize(records)
    high_tasks, centers = ref.make_high_tasks(summary)
    high = []
    if high_tasks:
        with Pool(ref.NPROC, initializer=_set_strict_policy) as pool:
            high = list(pool.imap_unordered(_evaluate_strict_high, high_tasks, chunksize=2))
        high.sort(key=lambda r: (r["VSO"], r["U"]))
        if any(r.get("error") for r in high):
            raise RuntimeError(f"high-resolution strict scan had {sum(bool(r.get('error')) for r in high)} errors")
    comparisons, neighborhoods = ref.compare_high(records, high, centers)

    os.makedirs(OUT_DIR, exist_ok=True)
    n_sc = sum(float(r.get("Tc") or 0.0) > 0.0 for r in records)
    n_rejected = sum(bool(r.get("no_inband_root")) for r in records)
    payload = dict(
        provenance=dict(
            stage="targeted strict-policy scan",
            policy="strict",
            band_definition="analytic continuum interval [-4t*sqrt(1+V_SO^2/(2t^2)), +4t*sqrt(1+V_SO^2/(2t^2))]",
            band_edge_tolerance=ref.s1sol.BAND_EDGE_TOL,
            n_fill=ref.N_FILL, Delta_t=ref.DELTA_T, hbar_omega=ref.OMEGA,
            grid=dict(U_min=min(ref.U_GRID), U_max=max(ref.U_GRID), dU=0.5,
                      VSO_min=min(ref.VSO_GRID), VSO_max=max(ref.VSO_GRID), dVSO=0.025,
                      n_U=len(ref.U_GRID), n_VSO=len(ref.VSO_GRID)),
            fractions=list(ref.FRACS),
            coarse_resolution=dict(Nk=ref.NK_STAGE1, Nk_mu=ref.NK_MU,
                                   DOS_Nk=ref.DOS_NK, nell=ref.NELL),
            high_resolution=ref.HIGH,
            n_coarse=len(records), n_strict_sc=len(accepted_count(records)),
            n_no_inband_root=n_rejected, n_high=len(high),
            caveat="Strict root admissibility is applied; it does not rank multiple in-band branches by free energy.",
            utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        ),
        summary=summary,
        screened_cutoff_brackets=cutoffs,
        high_resolution_comparisons=comparisons,
        refined_neighborhood_maxima=neighborhoods,
        coarse_data=records,
        high_resolution_data=high,
    )
    json_path = os.path.join(OUT_DIR, "refined_extrema_zoom_strict.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=1, allow_nan=False)
    csv_path = os.path.join(OUT_DIR, "refined_extrema_zoom_strict_summary.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        fields = ["VSO", "fraction", "n_points", "n_strict_sc", "n_no_inband_root",
                  "U_max", "R_max", "U_min", "R_min"]
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for item in summary:
            ext = item["extrema"]
            writer.writerow(dict(
                VSO=item["VSO"], fraction=item["fraction"], n_points=ext.get("n", 0),
                n_strict_sc=n_sc, n_no_inband_root=n_rejected,
                U_max=(ext.get("global_max") or {}).get("x"),
                R_max=(ext.get("global_max") or {}).get("value"),
                U_min=(ext.get("global_min") or {}).get("x"),
                R_min=(ext.get("global_min") or {}).get("value"),
            ))
    figure_path = _strict_figure(records)
    print(f"Strict policy: {n_sc}/{len(records)} points have an in-band SC solution; "
          f"{n_rejected} have no in-band T=0 root.")
    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {figure_path}")
    print(f"Completed in {time.time()-start:.1f}s")


def accepted_count(records):
    return [r for r in records if float(r.get("Tc") or 0.0) > 0.0]


if __name__ == "__main__":
    main()
