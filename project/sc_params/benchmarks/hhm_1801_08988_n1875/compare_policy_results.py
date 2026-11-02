#!/usr/bin/env python3
"""Paired relaxed/strict audit for the HHM 2018 Fig. 8 parameter set.

This additive benchmark leaves the production policy and existing results untouched. It tests
both chemical-potential policies on the paper's U=90t, Delta_t=4.5t, n=1.875, V_SO/t
= 0, 0.2, ..., 0.8 set; compares nonlinear 1%-gap Tc estimates with HHM's linearized
Tc determinant; samples the self-consistent gap equations at common fractions of the
linearized Tc; and includes an out-of-band control point where the policies must differ.

Outputs (in this directory): policy_comparison.json, policy_tc_comparison.csv,
policy_gap_functions.csv, and policy_comparison.png.
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import platform
import subprocess
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SC_PARAMS = REPO / "project" / "sc_params"
sys.path.insert(0, str(SC_PARAMS))

import solvers  # noqa: E402
from compare_fig8 import (  # noqa: E402
    DELTA_T, FILLING, LINEARIZED_NK, LOW_T, NONLINEAR_NK, U, VSO_VALUES,
    paper_linearized_tc,
)
from gap_equations import effective_gap, make_residual  # noqa: E402
from kgrid import build_kgrid  # noqa: E402
from parameters_tc import (  # noqa: E402
    BISECTION_TOL, GAP_THRESHOLD_FRAC, GRID_CONV_TOL, MIN_GAP_FOR_SC, T,
)

POLICIES = ("relaxed", "strict")
GAP_FRACTIONS = (0.0, 0.25, 0.5, 0.75, 0.9, 0.99)
RESIDUAL_TOL = solvers.RESID_TOL


def _git_head() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO, text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def _max_residual(residual, state, temp, eps, kernel):
    beta = np.inf if temp == 0.0 else 1.0 / float(temp)
    values = np.asarray(residual([state["D0"], state["DS"], state["mu"]],
                                beta, FILLING, eps, kernel), dtype=float)
    return float(np.max(np.abs(values))) if np.all(np.isfinite(values)) else float("inf")


def _window(mu, vso):
    band = solvers.bare_band_window(np.empty(0), float(vso))
    excess = max(float(mu) - band[1], band[0] - float(mu), 0.0)
    return dict(band=[float(band[0]), float(band[1])], excess=float(excess),
                inside=bool(excess <= solvers.BAND_EDGE_TOL))


def _solve_gap_curve(policy, vso, tc_linear, stage1, nk):
    """Solve and residual-check gaps at common T/Tc(HHM-linearized) for one policy."""
    solvers.MU_WINDOW_POLICY = policy
    eps, kernel = build_kgrid(vso, U, DELTA_T, Nk=nk)
    residual = make_residual(DELTA_T)
    t0 = [float(stage1["D0_0"]), float(stage1["DS_0"]), float(stage1["mu_0"])]
    previous = list(t0)
    band = solvers.bare_band_window(eps, vso) if policy == "strict" else None
    records = []

    for frac in GAP_FRACTIONS:
        temp = float(frac * tc_linear)
        if temp == 0.0:
            candidate = list(t0)
            deff = effective_gap(candidate[0], candidate[1], eps)
        else:
            deff, candidate = solvers._gap_at_T(
                temp, previous, FILLING, eps, kernel, residual,
                t0=t0, full=True, mu_band=band,
            )
            candidate = [float(x) for x in candidate]
        state = dict(D0=float(candidate[0]), DS=float(candidate[1]), mu=float(candidate[2]))
        resid = _max_residual(residual, state, temp, eps, kernel)
        win = _window(state["mu"], vso)
        root_ok = resid <= RESIDUAL_TOL
        policy_ok = (policy == "relaxed") or win["inside"]
        sc_ok = bool(root_ok and policy_ok and float(deff) > MIN_GAP_FOR_SC)
        records.append(dict(
            policy=policy, VSO=float(vso), Nk=int(nk), T_over_Tc_HHM=float(frac),
            T=float(temp), D0=state["D0"] if sc_ok else None,
            DS=state["DS"] if sc_ok else None,
            Deff=float(deff) if sc_ok else None,
            mu=state["mu"] if root_ok else None,
            residual_inf=float(resid), root_verified=bool(root_ok),
            mu_window=win, accepted=bool(sc_ok),
        ))
        if root_ok and policy_ok:
            previous = list(candidate)

    return records


def _solve_paper_point(vso, policy):
    solvers.MU_WINDOW_POLICY = policy
    result = solvers.solve_set(
        U=U, Delta_t=DELTA_T, n_fill=FILLING, V_SO=float(vso),
        Nk=NONLINEAR_NK, check_conv=True,
    )
    result["policy"] = policy
    return result


def _policy_control():
    """A known near-edge point verifies that relaxed and strict really are distinct modes."""
    rows = {}
    for policy in POLICIES:
        solvers.MU_WINDOW_POLICY = policy
        rows[policy] = solvers.solve_set_core(
            U=46.5, Delta_t=4.5, n_fill=1.885, V_SO=0.0, Nk=80,
        )
    relaxed, strict = rows["relaxed"], rows["strict"]
    diag = strict.get("diagnostic_unrestricted") or {}
    expected = (relaxed.get("Tc", 0.0) > 0.0
                and strict.get("Tc", 0.0) == 0.0
                and bool(strict.get("no_inband_root"))
                and not bool(diag.get("mu_inside_band", True)))
    return dict(
        parameter=dict(U=46.5, Delta_t=4.5, n=1.885, VSO=0.0, Nk=80,
                       note="policy-discrimination control; not an HHM Fig. 8 point"),
        relaxed={k: relaxed.get(k) for k in
                 ("Tc", "Tc_bracketed", "D0_0", "DS_0", "Deff_0", "mu_0", "mu_window")},
        strict={k: strict.get(k) for k in
                ("Tc", "Tc_bracketed", "no_inband_root", "mu_0", "diagnostic_unrestricted")},
        expected_policy_difference_observed=bool(expected),
    )


def _plot(tc_rows, gap_rows):
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5))
    vso = [r["VSO"] for r in tc_rows]
    lin = [r["Tc_HHM_linear_Nk320"] for r in tc_rows]
    axes[0, 0].plot(vso, lin, "k-o", ms=4, label="HHM linearized determinant")
    for policy, marker, color in (("relaxed", "s", "tab:blue"),
                                  ("strict", "^", "tab:orange")):
        vals = [r[f"Tc_{policy}_1pct_Nk160"] for r in tc_rows]
        axes[0, 0].plot(vso, vals, marker=marker, color=color, ls="--",
                        ms=4, label=f"{policy}: nonlinear 1% threshold")
    axes[0, 0].set(xlabel="$V_{SO}/t$", ylabel="$T_c/t$",
                   title="Paper criterion vs repository threshold")
    axes[0, 0].legend(fontsize=7)

    colors = plt.get_cmap("viridis")(np.linspace(0.08, 0.92, len(VSO_VALUES)))
    for ax, field, ylabel, title in (
        (axes[0, 1], "D0", "$D_0(T)/D_0(0)$", "On-site gap component"),
        (axes[1, 0], "DS", "$D_s(T)/D_s(0)$", "Correlated-hopping gap component"),
        (axes[1, 1], "Deff", "$\\Delta_{eff}(T)/\\Delta_{eff}(0)$", "Effective gap"),
    ):
        for color, v in zip(colors, VSO_VALUES):
            for policy, ls in (("relaxed", "-"), ("strict", "--")):
                pts = sorted((r for r in gap_rows
                              if r["VSO"] == float(v) and r["policy"] == policy
                              and r["Nk"] == NONLINEAR_NK),
                             key=lambda r: r["T_over_Tc_HHM"])
                base = next((r[field] for r in pts if r["T_over_Tc_HHM"] == 0.0
                             and r["accepted"]), None)
                if base is None or abs(float(base)) < 1e-14:
                    continue
                x, y = [], []
                for r in pts:
                    x.append(r["T_over_Tc_HHM"])
                    val = r[field]
                    y.append(float(val) / float(base) if r["accepted"] and val is not None
                             else np.nan)
                ax.plot(x, y, color=color, ls=ls, marker="o", ms=2.5,
                        label=f"$V_{{SO}}={v:g}$ {policy}")
        ax.set(xlabel="$T/T_c^{HHM}$", ylabel=ylabel, title=title)
        ax.legend(fontsize=6.3, ncol=2)
    fig.suptitle("Paired chemical-potential policies at the HHM Fig. 8 parameter set")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path = HERE / "policy_comparison.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def main():
    original_policy = solvers.MU_WINDOW_POLICY
    tc_rows, gap_rows = [], []
    linearized = {float(v): {
        nk: float(paper_linearized_tc(float(v), nk)) for nk in LINEARIZED_NK
    } for v in VSO_VALUES}

    try:
        for vso in VSO_VALUES:
            paper_tc = linearized[float(vso)][320]
            by_policy, by_policy_80 = {}, {}
            for policy in POLICIES:
                solvers.MU_WINDOW_POLICY = policy
                stage1_80 = solvers.solve_set_core(
                    U=U, Delta_t=DELTA_T, n_fill=FILLING, V_SO=float(vso), Nk=80,
                )
                stage1_160 = _solve_paper_point(vso, policy)
                by_policy_80[policy] = stage1_80
                by_policy[policy] = stage1_160
                gap_rows.extend(_solve_gap_curve(policy, vso, paper_tc, stage1_80, 80))
                gap_rows.extend(_solve_gap_curve(policy, vso, paper_tc, stage1_160, NONLINEAR_NK))

            relaxed, strict = by_policy["relaxed"], by_policy["strict"]
            relaxed_80, strict_80 = by_policy_80["relaxed"], by_policy_80["strict"]
            gaps = {}
            for key in ("Tc", "D0_0", "DS_0", "Deff_0", "mu_0"):
                a, b = relaxed.get(key), strict.get(key)
                gaps[key] = (abs(float(a) - float(b))
                             if a is not None and b is not None else None)
            def get_gap(policy, frac, nk):
                return next(r for r in gap_rows
                            if r["policy"] == policy and r["VSO"] == float(vso)
                            and r["Nk"] == nk and r["T_over_Tc_HHM"] == frac)
            policy_gap_diffs, gap_grid_diffs = {}, {}
            for frac in GAP_FRACTIONS:
                a = get_gap("relaxed", frac, NONLINEAR_NK)
                b = get_gap("strict", frac, NONLINEAR_NK)
                policy_gap_diffs[str(frac)] = {
                    field: (abs(float(a[field]) - float(b[field]))
                            if a["accepted"] and b["accepted"] else None)
                    for field in ("D0", "DS", "Deff", "mu")
                }
                gap_grid_diffs[str(frac)] = {}
                for policy in POLICIES:
                    lo, hi = get_gap(policy, frac, 80), get_gap(policy, frac, NONLINEAR_NK)
                    gap_grid_diffs[str(frac)][policy] = {
                        field: (abs(float(lo[field]) - float(hi[field]))
                                if lo["accepted"] and hi["accepted"] else None)
                        for field in ("D0", "DS", "Deff", "mu")
                    }
            tc_rows.append(dict(
                VSO=float(vso), Tc_HHM_linear_Nk80=linearized[float(vso)][80],
                Tc_HHM_linear_Nk160=linearized[float(vso)][160],
                Tc_HHM_linear_Nk320=paper_tc,
                Tc_relaxed_1pct_Nk80=float(relaxed_80["Tc"]),
                Tc_relaxed_1pct_Nk160=float(relaxed["Tc"]),
                Tc_relaxed_1pct_Nk320=float(relaxed["conv"]["Tc2"]),
                Tc_strict_1pct_Nk80=float(strict_80["Tc"]),
                Tc_strict_1pct_Nk160=float(strict["Tc"]),
                Tc_strict_1pct_Nk320=float(strict["conv"]["Tc2"]),
                nonlinear_Nk80_160_policy_abs_diff={
                    policy: {
                        key: abs(float(by_policy_80[policy][key]) - float(by_policy[policy][key]))
                        for key in ("Tc", "D0_0", "DS_0", "Deff_0", "mu_0")
                    }
                    for policy in POLICIES
                },
                relaxed_Tc_bracketed=bool(relaxed["Tc_bracketed"]),
                strict_Tc_bracketed=bool(strict["Tc_bracketed"]),
                relaxed_mu_window=relaxed.get("mu_window"),
                strict_mu_window=strict.get("mu_window"),
                relaxed_conv=relaxed["conv"], strict_conv=strict["conv"],
                relaxed_vs_strict_stage1_abs_diff=gaps,
                gap_function_abs_differences=policy_gap_diffs,
                gap_function_grid_abs_differences=gap_grid_diffs,
                max_abs_Tc_policy_difference_from_linearized={
                    "relaxed": abs(float(relaxed["Tc"]) - paper_tc),
                    "strict": abs(float(strict["Tc"]) - paper_tc),
                },
            ))

        # Enforce the expected exact-parameter conclusion: the policy constraint is inactive
        # at this published Fig. 8 set, and both calculations match to numerical precision.
        for row in gap_rows:
            if not row["root_verified"] or not row["accepted"]:
                raise AssertionError(f"Rejected/unverified sampled gap point: {row}")
            if row["residual_inf"] > RESIDUAL_TOL:
                raise AssertionError(f"Gap residual exceeds tolerance: {row}")
        for r in tc_rows:
            if not (r["strict_mu_window"] or {}).get("inside", False):
                raise AssertionError(f"Strict Fig. 8 stage-1 root is out of band at VSO={r['VSO']}")
            if not (r["relaxed_mu_window"] or {}).get("inside", False):
                raise AssertionError(f"Relaxed Fig. 8 root is out of band at VSO={r['VSO']}")
            if max(r["relaxed_vs_strict_stage1_abs_diff"].values()) > 1e-7:
                raise AssertionError(f"Policy stage-1 outputs differ at VSO={r['VSO']}: {r}")
            for frac, fields in r["gap_function_abs_differences"].items():
                for field, difference in fields.items():
                    if difference is not None and difference > 1e-7:
                        raise AssertionError(
                            f"Policy gap functions differ at VSO={r['VSO']}, f={frac}, {field}={difference}"
                        )
            if not (r["relaxed_conv"].get("ok") and r["strict_conv"].get("ok")):
                raise AssertionError(f"Nk=80/160 convergence check failed at VSO={r['VSO']}")

        control = _policy_control()
        if not control["expected_policy_difference_observed"]:
            raise AssertionError("The policy-discrimination control did not differ as expected")
    finally:
        solvers.MU_WINDOW_POLICY = original_policy

    provenance = dict(
        stage="paired HHM Fig. 8 chemical-potential policy audit",
        paper="Hutchinson, Hirsch, Marsiglio, arXiv:1801.08988, Fig. 8 parameter set",
        parameters=dict(U=U, Delta_t=DELTA_T, n=FILLING, VSO=list(VSO_VALUES)),
        policy_modes=list(POLICIES), policy_mutation="solvers.MU_WINDOW_POLICY only; restored on exit",
        strict_band="analytic square-lattice Rashba band; edge tolerance 1e-8 t",
        paper_Tc_method="linearized HHM determinant plus normal-state number equation",
        repository_Tc_method="nonlinear gap bisection where Deff(T)=0.01 Deff(0); not the HHM Tc criterion",
        nonlinear_resolution=dict(Nk=NONLINEAR_NK, convergence_Nk2=2 * NONLINEAR_NK,
                                  tolerance=GRID_CONV_TOL),
        gap_resolution=dict(Nk=[80, NONLINEAR_NK]),
        linearized_resolution=list(LINEARIZED_NK),
        gap_sample_fractions_of_linearized_Tc=list(GAP_FRACTIONS),
        residual_tolerance=RESIDUAL_TOL,
        policy_control=control,
        caveat="Matching the Fig. 8 parameter set cannot distinguish strict from relaxed when every accepted root is in-band. It also cannot establish free-energy ordering or exact digitized curve agreement.",
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        git_head=_git_head(), utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    )
    output = dict(provenance=provenance, tc_comparison=tc_rows, gap_functions=gap_rows)
    json_path = HERE / "policy_comparison.json"
    with json_path.open("w", encoding="utf-8") as f:
        json.dump(output, f, indent=1, allow_nan=False)

    tc_csv = HERE / "policy_tc_comparison.csv"
    with tc_csv.open("w", newline="", encoding="utf-8") as f:
        fields = ["VSO", "Tc_HHM_linear_Nk320", "Tc_relaxed_1pct_Nk80",
                  "Tc_relaxed_1pct_Nk160", "Tc_relaxed_1pct_Nk320",
                  "Tc_strict_1pct_Nk80", "Tc_strict_1pct_Nk160",
                  "Tc_strict_1pct_Nk320", "mu_inside", "mu_excess",
                  "max_abs_policy_delta_Tc", "max_abs_policy_delta_D0_0",
                  "max_abs_policy_delta_DS_0", "max_abs_policy_delta_Deff_0",
                  "relaxed_conv_ok", "strict_conv_ok"]
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for r in tc_rows:
            w = r["strict_mu_window"] or {}
            d = r["relaxed_vs_strict_stage1_abs_diff"]
            writer.writerow(dict(
                VSO=r["VSO"], Tc_HHM_linear_Nk320=r["Tc_HHM_linear_Nk320"],
                Tc_relaxed_1pct_Nk80=r["Tc_relaxed_1pct_Nk80"],
                Tc_relaxed_1pct_Nk160=r["Tc_relaxed_1pct_Nk160"],
                Tc_relaxed_1pct_Nk320=r["Tc_relaxed_1pct_Nk320"],
                Tc_strict_1pct_Nk80=r["Tc_strict_1pct_Nk80"],
                Tc_strict_1pct_Nk160=r["Tc_strict_1pct_Nk160"],
                Tc_strict_1pct_Nk320=r["Tc_strict_1pct_Nk320"],
                mu_inside=w.get("inside"), mu_excess=w.get("excess"),
                max_abs_policy_delta_Tc=d["Tc"], max_abs_policy_delta_D0_0=d["D0_0"],
                max_abs_policy_delta_DS_0=d["DS_0"], max_abs_policy_delta_Deff_0=d["Deff_0"],
                relaxed_conv_ok=r["relaxed_conv"].get("ok"),
                strict_conv_ok=r["strict_conv"].get("ok"),
            ))

    gap_csv = HERE / "policy_gap_functions.csv"
    with gap_csv.open("w", newline="", encoding="utf-8") as f:
        fields = ["policy", "VSO", "Nk", "T_over_Tc_HHM", "T", "D0", "DS", "Deff",
                  "mu", "residual_inf", "root_verified", "band_inside", "band_excess", "accepted"]
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for r in gap_rows:
            writer.writerow(dict(
                policy=r["policy"], VSO=r["VSO"], Nk=r["Nk"],
                T_over_Tc_HHM=r["T_over_Tc_HHM"],
                T=r["T"], D0=r["D0"], DS=r["DS"], Deff=r["Deff"], mu=r["mu"],
                residual_inf=r["residual_inf"], root_verified=r["root_verified"],
                band_inside=r["mu_window"]["inside"], band_excess=r["mu_window"]["excess"],
                accepted=r["accepted"],
            ))
    fig_path = _plot(tc_rows, gap_rows)

    print("VSO/t | Tc HHM linear | Tc relaxed(160/320) | Tc strict(160/320) | max policy gap")
    for r in tc_rows:
        diffs = [d for by_field in r["gap_function_abs_differences"].values()
                 for d in by_field.values() if d is not None]
        print(f"{r['VSO']:.1f} | {r['Tc_HHM_linear_Nk320']:.9f} | "
              f"{r['Tc_relaxed_1pct_Nk160']:.9f}/{r['Tc_relaxed_1pct_Nk320']:.9f} | "
              f"{r['Tc_strict_1pct_Nk160']:.9f}/{r['Tc_strict_1pct_Nk320']:.9f} | "
              f"{max(diffs) if diffs else 0.0:.3e}")
    print(f"Control distinguishes policies: {control['expected_policy_difference_observed']}")
    print(f"Wrote {json_path}, {tc_csv}, {gap_csv}, {fig_path}")


if __name__ == "__main__":
    main()
