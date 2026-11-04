#!/usr/bin/env python3
"""Compare legacy fitted gaps with the optional self-consistent HHM ultrasonic states.

Runs the unchanged absorption.ratio_all kernel for both state models. Outputs are isolated
under this benchmark directory; strict-window results are the primary acceptance diagnostics
and relaxed roots are reported separately, never ranked as thermodynamic solutions.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SC_PARAMS = REPO / "project" / "sc_params"
ULTRA = REPO / "project" / "ultrasonic"
sys.path.insert(0, str(SC_PARAMS))
sys.path.insert(0, str(ULTRA))

from hhm_finite_t import HHMFiniteTemperatureSolver
from solvers import solve_set_core
from absorption import ratio_all
from dos_ua import NormalDOS
from parameters_ua import OMEGA
from physics_ua import SCState


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=HERE / "results")
    ap.add_argument("--U", type=float, nargs="+", default=[30.0, 45.0, 75.0])
    ap.add_argument("--vso", type=float, nargs="+", default=[0.0, 0.4, 0.8])
    ap.add_argument("--fractions", type=float, nargs="+", default=[0.5, 0.9])
    ap.add_argument("--filling", type=float, default=1.885)
    ap.add_argument("--delta-t", type=float, default=4.5)
    ap.add_argument("--nk", type=int, default=80)
    ap.add_argument("--nk-dos", type=int, default=600)
    ap.add_argument("--nell", type=int, default=1400)
    ap.add_argument("--omega", type=float, default=OMEGA)
    ap.add_argument("--tc-max", type=float, default=8.0)
    ap.add_argument("--tc-scan", type=int, default=401)
    ap.add_argument("--initial-root-index", type=int, default=0,
                    help="select among accepted T=0 roots sorted by decreasing Deff")
    args = ap.parse_args()
    if args.nk < 4 or args.nk_dos < 4 or args.nell < 16:
        raise ValueError("Nk and NkDOS must be >=4 and nell >=16")
    if any(not 0.0 <= f < 1.0 for f in args.fractions):
        raise ValueError("fractions must lie in [0,1)")
    args.output.mkdir(parents=True, exist_ok=True)
    rows, summary = [], []
    fields = ["U", "Delta_t", "n", "V_SO", "Nk", "mu_policy", "strict_primary",
              "frac", "legacy_fraction", "T", "Tc_HHM_linearized", "Tc_legacy_1pct", "root_status",
              "mu_SC", "mu_N", "mu_inside", "mu_excess", "D0", "DS", "Deff",
              "residual_inf", "W_WN_HHM", "W_WN_legacy_fit", "difference", "relative_difference",
              "branch_selection"]

    for vso in args.vso:
        print(f"Building normal DOS for VSO={vso:g}...", flush=True)
        dos = NormalDOS(vso, Nk=args.nk_dos)
        for U in args.U:
            print(f"U={U:g}, VSO={vso:g}: solving branches...", flush=True)
            legacy = solve_set_core(U, args.delta_t, args.filling, vso, Nk=args.nk)
            legacy_state = SCState.from_params(
                vso, float(legacy["Tc"]), float(legacy["D0_0"]), float(legacy["DS_0"]),
                float(legacy["Deff_0"]), n_tgt=args.filling, Nk_mu=max(300, args.nk))
            for policy in ("strict", "relaxed"):
                entry = dict(U=U, Delta_t=args.delta_t, n=args.filling, V_SO=vso,
                             Nk=args.nk, mu_policy=policy, strict_primary=(policy == "strict"),
                             initial_root_index=int(args.initial_root_index),
                             Tc_legacy_1pct=float(legacy["Tc"]),
                             Tc_legacy_bracketed=bool(legacy.get("Tc_bracketed", False)))
                try:
                    model = HHMFiniteTemperatureSolver(
                        U, args.delta_t, args.filling, vso, Nk=args.nk, mu_policy=policy,
                        tc_scan_max=args.tc_max, tc_scan_points=args.tc_scan,
                        initial_root_index=args.initial_root_index)
                    entry["Tc_HHM_linearized"] = model.Tc
                    entry["t0_candidates"] = len(model.t0_roots)
                    entry["t0_accepted"] = model.t0 is not None
                    entry["t0_roots"] = [{k: r[k] for k in
                                          ("D0", "DS", "Deff", "mu", "residual_inf", "mu_window")}
                                         for r in model.t0_roots]
                    if model.Tc is None or model.t0 is None:
                        entry["status"] = "no_accepted_nontrivial_T0_root"
                        summary.append(entry)
                        continue
                    states = model.solve_fractions(args.fractions, include_diagnostics=True)
                    requested = [s for s in states if s.get("requested")]
                    entry["finite_T_diagnostics"] = [
                        dict(frac=float(s["T"] / model.Tc), T=float(s["T"]),
                             selected=dict(D0=float(s["D0"]), DS=float(s["DS"]),
                                          mu=float(s["mu"]), Deff=float(s["Deff"]),
                                          residual_inf=float(s["residual_inf"]),
                                          mu_window=s["mu_window"],
                                          branch_selection=s["branch_selection"]),
                             candidates=s.get("candidates", []))
                        for s in requested]
                    for st in requested:
                        frac = float(st["T"] / model.Tc)
                        hhm_state = dict(st, frac=frac)
                        hhm = ratio_all(dos, hhm_state, args.omega,
                                        variant_list=("a_signed",), nell=args.nell)
                        legacy_frac = float(st["T"] / legacy["Tc"])
                        fitted = legacy_state.at(legacy_frac)
                        old = ratio_all(dos, fitted, args.omega,
                                        variant_list=("a_signed",), nell=args.nell)
                        ratio_hhm = float(hhm["a_signed_ratio"])
                        ratio_old = float(old["a_signed_ratio"])
                        rows.append(dict(U=U, Delta_t=args.delta_t, n=args.filling, V_SO=vso,
                                         Nk=args.nk, mu_policy=policy,
                                         strict_primary=(policy == "strict"), frac=frac,
                                         legacy_fraction=legacy_frac, T=float(st["T"]),
                                         Tc_HHM_linearized=float(model.Tc),
                                         Tc_legacy_1pct=float(legacy["Tc"]),
                                         root_status="residual_verified",
                                         mu_SC=float(st["mu"]), mu_N=float(st["mu_N"]),
                                         mu_inside=bool(st["mu_window"]["inside"]),
                                         mu_excess=float(st["mu_window"]["excess"]),
                                         D0=float(st["D0"]), DS=float(st["DS"]),
                                         Deff=float(st["Deff"]),
                                         residual_inf=float(st["residual_inf"]),
                                         W_WN_HHM=ratio_hhm, W_WN_legacy_fit=ratio_old,
                                         difference=ratio_hhm-ratio_old,
                                         relative_difference=(ratio_hhm-ratio_old)/ratio_old
                                         if ratio_old else None,
                                         branch_selection=st["branch_selection"]))
                    entry["status"] = "completed"
                except Exception as exc:
                    entry["status"] = "failed"
                    entry["error"] = f"{type(exc).__name__}: {exc}"
                    print(f"  {policy} policy: {entry['error']}", flush=True)
                summary.append(entry)
    csv_path = args.output / "ultrasonic_hhm_vs_legacy.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)
    json_path = args.output / "ultrasonic_hhm_summary.json"
    json_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Wrote {csv_path}\nWrote {json_path}")
    try:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9.2, 5.6))
        combos = sorted({(r["U"], r["V_SO"]) for r in rows})
        colors = plt.cm.tab20(range(max(len(combos), 1)))
        for ci, (U, vso) in enumerate(combos):
            for policy, style in (("strict", "-"), ("relaxed", "--")):
                subset = [r for r in rows if r["U"] == U and r["V_SO"] == vso
                          and r["mu_policy"] == policy]
                subset.sort(key=lambda r: r["frac"])
                if subset:
                    ax.plot([r["frac"] for r in subset], [r["W_WN_HHM"] for r in subset],
                            marker="o", linestyle=style, color=colors[ci],
                            label=f"U={U:g}, VSO={vso:g}, {policy}, HHM")
                    ax.plot([r["frac"] for r in subset], [r["W_WN_legacy_fit"] for r in subset],
                            marker="s", markerfacecolor="none", linestyle=style,
                            color=colors[ci], label=f"U={U:g}, VSO={vso:g}, {policy}, fit")
        ax.set(xlabel=r"$T/T_c$ (HHM determinant)", ylabel=r"$W/W_N$",
               title="Existing response kernel: self-consistent HHM states")
        ax.grid(alpha=0.25); ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(args.output / "ultrasonic_hhm_vs_legacy.png", dpi=180)
        plt.close(fig)
        print(f"Wrote {args.output / 'ultrasonic_hhm_vs_legacy.png'}")
    except ImportError as exc:
        print(f"Plot not generated (matplotlib unavailable): {exc}")
    print(f"Rows: {len(rows)}; strict-window rows are marked as primary acceptance data.")


if __name__ == "__main__":
    main()
