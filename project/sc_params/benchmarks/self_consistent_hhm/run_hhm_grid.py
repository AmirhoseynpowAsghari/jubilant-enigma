#!/usr/bin/env python3
"""Isolated HHM finite-T convergence/observable export (does not touch Fig. 8 outputs).

Example:
  python run_hhm_grid.py --output results --nk 80 160 320 --vso 0 .4 .8

CSV contains the linearized determinant Tc, legacy 1%-gap Tc proxy, the chosen T=0 root,
finite-T coupled states and fitted-gap mu comparison.  Paper values for complete finite-T
curves are not digitized or supplied here; this runner exports calculations for transparent
comparison and does not claim exact curve reproduction.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SC_PARAMS = REPO / "project" / "sc_params"
ULTRA = REPO / "project" / "ultrasonic"
sys.path.insert(0, str(SC_PARAMS))
sys.path.insert(0, str(ULTRA))

from hhm_finite_t import HHMFiniteTemperatureSolver, compare_nonlinear_tc_proxy
from physics_ua import SCState
from solvers import solve_set_core


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=HERE / "results")
    ap.add_argument("--nk", type=int, nargs="+", default=[80, 160, 320])
    ap.add_argument("--vso", type=float, nargs="+", default=[0.0, 0.2, 0.4, 0.6, 0.8])
    ap.add_argument("--fractions", type=float, nargs="+", default=[0.0, 0.25, 0.5, 0.75, 0.9, 0.99])
    ap.add_argument("--U", type=float, default=90.0)
    ap.add_argument("--delta-t", type=float, default=4.5)
    ap.add_argument("--filling", type=float, default=1.875)
    ap.add_argument("--tc-max", type=float, default=8.0)
    ap.add_argument("--tc-scan", type=int, default=401)
    ap.add_argument("--initial-root-index", type=int, default=0,
                    help="select among accepted T=0 roots sorted by decreasing Deff; default 0")
    args = ap.parse_args()
    if any(n < 4 for n in args.nk) or any(not 0.0 <= v < 1.0 for v in args.vso):
        raise ValueError("invalid Nk or SOC")
    if any(not 0.0 <= f < 1.0 for f in args.fractions):
        raise ValueError("fractions must lie in [0,1)")
    args.output.mkdir(parents=True, exist_ok=True)

    fields = ["U", "Delta_t", "n", "V_SO", "Nk", "frac", "T", "Tc_HHM_linearized",
              "Tc_legacy_1pct", "delta_Tc", "relative_delta_Tc", "D0", "DS", "Deff",
              "mu_SC", "mu_N", "mu_fitted_gap", "mu_fit_fraction", "mu_fit_minus_mu_SC", "residual_inf",
              "mu_in_band", "mu_excess", "branch_selection"]
    rows = []
    summary = []
    for vso in args.vso:
        for nk in args.nk:
            print(f"Solving U={args.U:g}, n={args.filling:g}, VSO={vso:g}, Nk={nk}...", flush=True)
            comp = compare_nonlinear_tc_proxy(args.U, args.delta_t, args.filling, vso, Nk=nk)
            solver = HHMFiniteTemperatureSolver(
                args.U, args.delta_t, args.filling, vso, Nk=nk,
                tc_scan_max=args.tc_max, tc_scan_points=args.tc_scan,
                initial_root_index=args.initial_root_index)
            if solver.Tc is None:
                summary.append(dict(U=args.U, Delta_t=args.delta_t, n=args.filling,
                                    V_SO=vso, Nk=nk, status="no_determinant_root"))
                continue
            path = solver.solve_fractions(args.fractions, include_diagnostics=True)
            legacy = comp["legacy_stage1"]
            legacy_state = SCState.from_params(
                vso, float(legacy["Tc"]), float(legacy["D0_0"]), float(legacy["DS_0"]),
                float(legacy["Deff_0"]), n_tgt=args.filling, Nk_mu=max(300, nk))
            tc_info = dict(U=args.U, Delta_t=args.delta_t, n=args.filling, V_SO=vso, Nk=nk,
                           Tc_HHM_linearized=float(solver.Tc),
                           Tc_legacy_1pct=float(legacy["Tc"]),
                           delta_Tc=float(legacy["Tc"] - solver.Tc),
                           initial_root_index=int(args.initial_root_index),
                           relative_delta_Tc=float((legacy["Tc"] - solver.Tc) / solver.Tc),
                           roots=solver.tc_info["roots"],
                           t0_root={k: solver.t0[k] for k in
                                    ("D0", "DS", "Deff", "mu", "residual_inf", "mu_window")},
                           t0_candidates=[{k: r[k] for k in
                                           ("D0", "DS", "Deff", "mu", "residual_inf", "mu_window", "starts")}
                                          for r in solver.t0_roots])
            tc_info["finite_T_diagnostics"] = [
                dict(T=float(state["T"]), frac=float(state["T"] / solver.Tc),
                     selected=dict(D0=float(state["D0"]), DS=float(state["DS"]),
                                  mu=float(state["mu"]), Deff=float(state["Deff"]),
                                  residual_inf=float(state["residual_inf"]),
                                  mu_window=state["mu_window"],
                                  branch_selection=state["branch_selection"]),
                     candidates=state.get("candidates", []))
                for state in path if state.get("requested", False)]
            summary.append(tc_info)
            for state in path:
                frac = float(state["T"] / solver.Tc)
                if not state.get("requested", False):
                    continue
                fit_frac = float(state["T"] / legacy["Tc"])
                fit = legacy_state.at(fit_frac)
                rows.append(dict(
                    U=args.U, Delta_t=args.delta_t, n=args.filling, V_SO=vso, Nk=nk,
                    frac=frac, T=float(state["T"]), Tc_HHM_linearized=float(solver.Tc),
                    Tc_legacy_1pct=float(legacy["Tc"]),
                    delta_Tc=float(legacy["Tc"] - solver.Tc),
                    relative_delta_Tc=float((legacy["Tc"] - solver.Tc) / solver.Tc),
                    D0=float(state["D0"]), DS=float(state["DS"]), Deff=float(state["Deff"]),
                    mu_SC=float(state["mu"]), mu_N=float(state["mu_N"]),
                    mu_fitted_gap=float(fit["mu"]), mu_fit_fraction=fit_frac,
                    mu_fit_minus_mu_SC=float(fit["mu"] - state["mu"]),
                    residual_inf=float(state["residual_inf"]),
                    mu_in_band=bool(state["mu_window"]["inside"]),
                    mu_excess=float(state["mu_window"]["excess"]),
                    branch_selection=state["branch_selection"]))
    csv_path = args.output / "hhm_self_consistent_grid.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)
    json_path = args.output / "hhm_self_consistent_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {csv_path}\nWrote {json_path}")
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        for nk in args.nk:
            subset = [s for s in summary if s.get("Nk") == nk and "Tc_HHM_linearized" in s]
            subset.sort(key=lambda s: s["V_SO"])
            if subset:
                ax.plot([s["V_SO"] for s in subset], [s["Tc_HHM_linearized"] for s in subset],
                        marker="o", label=f"HHM determinant, Nk={nk}")
                ax.plot([s["V_SO"] for s in subset], [s["Tc_legacy_1pct"] for s in subset],
                        marker="x", linestyle="--", label=f"legacy 1% proxy, Nk={nk}")
        ax.set(xlabel=r"$V_{SO}/t$", ylabel=r"$T_c/t$",
               title="Calculated transition temperatures (no digitized paper data)")
        ax.grid(alpha=0.25); ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(args.output / "tc_grid_convergence.png", dpi=180)
        plt.close(fig)

        fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.3), sharex=True)
        observables = (("D0", r"$D_0/t$"), ("DS", r"$D_S/t$"),
                       ("mu_SC_minus_mu_N", r"$(\mu_{SC}-\mu_N)/t$"))
        colors = plt.cm.viridis(np.linspace(0.05, 0.95, max(len(args.vso), 2)))
        styles = {nk: style for nk, style in zip(args.nk, ("-", "--", ":", "-."))}
        for ax, (key, label) in zip(axes, observables):
            for ci, vso in enumerate(args.vso):
                for nk in args.nk:
                    subset = [r for r in rows if r["V_SO"] == vso and r["Nk"] == nk]
                    subset.sort(key=lambda r: r["frac"])
                    if subset:
                        y = ([r["D0"] for r in subset] if key == "D0" else
                             [r["DS"] for r in subset] if key == "DS" else
                             [r["mu_SC"] - r["mu_N"] for r in subset])
                        ax.plot([r["frac"] for r in subset], y, color=colors[ci],
                                linestyle=styles.get(nk, "-"), marker=".",
                                label=f"VSO={vso:g}, Nk={nk}")
            ax.set_xlabel(r"$T/T_c$ (determinant)"); ax.set_ylabel(label); ax.grid(alpha=0.2)
        axes[0].legend(fontsize=6, ncol=2, loc="best")
        fig.suptitle("Self-consistent finite-T observables; lines are calculations, not digitized references")
        fig.tight_layout()
        fig.savefig(args.output / "finite_t_observables.png", dpi=180)
        plt.close(fig)
        print(f"Wrote {args.output / 'tc_grid_convergence.png'}")
        print(f"Wrote {args.output / 'finite_t_observables.png'}")
    except ImportError as exc:
        print(f"Plots not generated (matplotlib unavailable): {exc}")
    for s in summary:
        if "Tc_HHM_linearized" in s:
            print(f"VSO={s['V_SO']:g}, Nk={s['Nk']}: Tc linear={s['Tc_HHM_linearized']:.9f}, "
                  f"1%-proxy={s['Tc_legacy_1pct']:.9f}, delta={s['delta_Tc']:+.3e}, "
                  f"T0 residual={s['t0_root']['residual_inf']:.2e}")


if __name__ == "__main__":
    main()
