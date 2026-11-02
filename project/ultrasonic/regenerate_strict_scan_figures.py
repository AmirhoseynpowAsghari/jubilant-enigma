#!/usr/bin/env python3
"""Generate separate strict-policy counterparts for the extrema scan figures.

Relaxed-policy results are read-only. Each coarse (U, V_SO) point is re-solved
under strict bare-band admissibility; stage-2 points are included only when the
finite-temperature superconducting chemical potential is in the analytic band.
The five fixed-U SOC profiles are taken from strict_family_data.json, and their
fixed-absolute-temperature points are independently screened and evaluated.

Outputs are isolated under results/strict_policy and figures/strict_policy.
"""
from __future__ import annotations

import datetime as dt
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import solvers
import scan_extrema_vs_soc as scan
from physics_ua import SCState

OUT_DIR = os.path.join(HERE, "results", "strict_policy")
FIG_DIR = os.path.join(HERE, "figures", "strict_policy")
BAND_TOL = solvers.BAND_EDGE_TOL


def _init_worker():
    solvers.MU_WINDOW_POLICY = "strict"
    scan.s1sol.MU_WINDOW_POLICY = "strict"


def _inside(mu, vso):
    lo, hi = solvers.bare_band_window(np.empty(0), float(vso))
    return mu is not None and lo - BAND_TOL <= float(mu) <= hi + BAND_TOL


def _evaluate(task):
    U, vso = task
    row = dict(U=float(U), VSO=float(vso), policy="strict")
    try:
        sol = solvers.solve_set_core(float(U), scan.DELTA_T, scan.N_FILL, float(vso),
                                     Nk=scan.NK_STAGE1)
    except Exception as exc:
        return row | dict(Tc=0.0, error=f"{type(exc).__name__}: {exc}",
                          mu_window=None, frac={}, fixed_T={})
    tc = float(sol.get("Tc") or 0.0)
    row.update(Tc=tc, mu_0=sol.get("mu_0"), mu_window=sol.get("mu_window"),
               Tc_bracketed=bool(sol.get("Tc_bracketed", tc > 0.0)),
               no_inband_root=bool(sol.get("no_inband_root", False)),
               diagnostic_unrestricted=sol.get("diagnostic_unrestricted"),
               D0_0=sol.get("D0_0"), DS_0=sol.get("DS_0"), Deff_0=sol.get("Deff_0"))
    frac, fixed = {}, {}
    if tc > 0.0:
        sc = SCState(params=dict(VSO=float(vso), Tc=tc, D0_0=float(sol["D0_0"]),
                                 DS_0=float(sol["DS_0"]), Deff_0=float(sol["Deff_0"])),
                     Nk_mu=scan.NK_MU)
        dos = scan.get_dos(vso)
        for f in scan.FRACS:
            state = sc.at(float(f))
            keep = _inside(state["mu"], vso)
            frac[f"{f:g}"] = dict(T=float(state["T"]), mu=float(state["mu"]),
                                  mu_N=float(state["mu_N"]), in_band=bool(keep),
                                  ratio=scan._ratio(dos, state) if keep else None)
        for temp in scan.ABS_T:
            k = f"{temp:.2f}"
            if temp >= tc:
                fixed[k] = dict(T=float(temp), ratio=None, superconducting=False,
                                mu=None, in_band=False)
            else:
                f = float(temp / tc)
                state = sc.at(f)
                keep = _inside(state["mu"], vso)
                fixed[k] = dict(T=float(temp), ratio=scan._ratio(dos, state) if keep else None,
                                superconducting=True, frac=f, mu=float(state["mu"]),
                                mu_N=float(state["mu_N"]), in_band=bool(keep))
    else:
        frac = {f"{f:g}": dict(T=None, ratio=None, mu=None, mu_N=None, in_band=False)
                for f in scan.FRACS}
        fixed = {f"{t:.2f}": dict(T=float(t), ratio=None, superconducting=False,
                                  mu=None, in_band=False) for t in scan.ABS_T}
    return row | dict(frac=frac, fixed_T=fixed)


def _build_strict_soc_data(family_records):
    selected_u = {30.0, 45.0, 75.0, 95.0, 115.0}
    out = []
    for r in family_records:
        U, vso, tc = float(r["U"]), float(r["VSO"]), float(r["Tc"])
        if U not in selected_u:
            continue
        row = dict(U=U, VSO=vso, Tc=tc, policy="strict", frac={}, fixed_T={})
        if tc > 0.0:
            sc = SCState(params=dict(VSO=vso, Tc=tc, D0_0=float(r["D0_0"]),
                                     DS_0=float(r["DS_0"]), Deff_0=float(r["Deff_0"])),
                         Nk_mu=scan.NK_MU)
            dos = scan.get_dos(vso)
            for i, f in enumerate(scan.FRACS):
                mu = r.get("stage2_mu", [None] * len(scan.FRACS))[i]
                keep = bool(r.get("stage2_inside", [False] * len(scan.FRACS))[i])
                row["frac"][f"{f:g}"] = dict(ratio=r["ratio"][i], mu=mu, in_band=keep)
            for temp in scan.ABS_T:
                k = f"{temp:.2f}"
                if temp >= tc:
                    row["fixed_T"][k] = dict(T=float(temp), ratio=None, superconducting=False,
                                            mu=None, in_band=False)
                else:
                    f = float(temp / tc)
                    state = sc.at(f)
                    keep = _inside(state["mu"], vso)
                    row["fixed_T"][k] = dict(T=float(temp),
                                             ratio=scan._ratio(dos, state) if keep else None,
                                             superconducting=True, frac=f, mu=float(state["mu"]),
                                             in_band=bool(keep))
        else:
            row["frac"] = {f"{f:g}": dict(ratio=None, mu=None, in_band=False)
                           for f in scan.FRACS}
            row["fixed_T"] = {f"{t:.2f}": dict(T=float(t), ratio=None,
                                                 superconducting=False, mu=None, in_band=False)
                               for t in scan.ABS_T}
        out.append(row)
    return sorted(out, key=lambda x: (x["U"], x["VSO"]))


def _save_strict_scan_figures(records, summary, soc_data):
    # Match the relaxed analysis layouts while retaining NaN gaps at every rejected sample.
    os.makedirs(FIG_DIR, exist_ok=True)
    plt = scan.plt
    palette = {0.0: "#1f77b4", 0.2: "#ff7f0e", 0.4: "#2ca02c", 0.6: "#d62728",
               0.8: "#9467bd", 0.4082: "#8c564b", 0.1: "#17becf", 0.3: "#bcbd22",
               0.5: "#e377c2", 0.7: "#7f7f7f", 0.9: "#393b79"}
    panels = [("fixed_T_over_Tc", 0.5, "$T/T_c=0.5$"),
              ("fixed_T_over_Tc", 0.9, "$T/T_c=0.9$"),
              ("fixed_absolute_T", 0.1, "$T=0.1t$"),
              ("fixed_absolute_T", 0.2, "$T=0.2t$")]

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    for (ensemble, temp, title), ax in zip(panels, axes.flat):
        for v in (0.0, 0.2, 0.4082, 0.4, 0.6, 0.8, 0.9):
            rows = sorted((r for r in records if r["VSO"] == v), key=lambda r: r["U"])
            if not rows: continue
            key = f"{temp:g}" if ensemble == "fixed_T_over_Tc" else f"{temp:.2f}"
            get = (lambda r, k=key: r["frac"].get(k, {}).get("ratio")) if ensemble == "fixed_T_over_Tc" else (lambda r, k=key: r["fixed_T"].get(k, {}).get("ratio"))
            ys = [np.nan if get(r) is None else float(get(r)) for r in rows]
            xs = [r["U"] for r in rows]
            color = palette[v]
            ax.plot(xs, ys, "-", lw=1.25, color=color, label=f"$V_{{SO}}={v:g}t$")
            scan._mark_extrema(ax, rows, get, color)
        ax.set_title(title); ax.set_xlabel("$U/t$"); ax.set_ylabel("$W/W_N$")
        handles, _ = ax.get_legend_handles_labels()
        if handles: ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Strict-policy global maxima/minima vs $U$; rejected/out-of-band samples are gaps\n"
                 "Absolute-temperature curves include only in-band superconducting states")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(FIG_DIR, "extrema_by_U.png"), dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    for (ensemble, temp, title), ax in zip(panels, axes.flat):
        for idx, U in enumerate((30.0, 45.0, 75.0, 95.0, 115.0)):
            rows = sorted((r for r in soc_data if r["U"] == U), key=lambda r: r["VSO"])
            key = f"{temp:g}" if ensemble == "fixed_T_over_Tc" else f"{temp:.2f}"
            get = (lambda r, k=key: r["frac"].get(k, {}).get("ratio")) if ensemble == "fixed_T_over_Tc" else (lambda r, k=key: r["fixed_T"].get(k, {}).get("ratio"))
            xs = [r["VSO"] for r in rows]
            ys = [np.nan if get(r) is None else float(get(r)) for r in rows]
            color = plt.get_cmap("viridis")(0.1 + 0.8 * idx / 4)
            ax.plot(xs, ys, "-o", ms=2.6, lw=1.2, color=color, label=f"$U={U:g}t$")
            scan._mark_extrema(ax, rows, get, color, x_key="VSO")
        ax.set_title(title); ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("$W/W_N$")
        handles, _ = ax.get_legend_handles_labels()
        if handles: ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Strict-policy $W/W_N$ vs $V_{SO}$ at fixed $U$; rejected samples are gaps")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(FIG_DIR, "soc_scans_by_U.png"), dpi=170); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    for (ensemble, temp, title), ax in zip(panels, axes.flat):
        for which, color, marker in (("global_max", "tab:red", "o"), ("global_min", "tab:blue", "s")):
            rows = sorted((r for r in summary if r["ensemble"] == ensemble and abs(r["value"]-temp) < 1e-10), key=lambda r: r["VSO"])
            xs, ys, endpoints = [], [], []
            for r in rows:
                e = r["extrema"].get(which)
                if e is not None:
                    xs.append(r["VSO"]); ys.append(e["value"]); endpoints.append(e["endpoint"])
            if xs:
                ax.plot(xs, ys, f"-{marker}", color=color, lw=1.2, label="maximum" if which == "global_max" else "minimum")
                for x, y, endpoint in zip(xs, ys, endpoints):
                    if endpoint: ax.scatter([x], [y], marker=marker, s=45, facecolors="none", edgecolors=color, linewidths=1.2)
        ax.set_title(title); ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("extremal $W/W_N$ over $U$")
        handles, _ = ax.get_legend_handles_labels()
        if handles: ax.legend(fontsize=8)
    fig.suptitle("Strict-policy extremal attenuation vs SOC; open markers denote U-boundary extrema")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(FIG_DIR, "extrema_tracks_by_soc.png"), dpi=170); plt.close(fig)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    tasks = [(u, v) for v in scan.VSO_GRID for u in scan.U_GRID]
    t0 = time.time()
    records = []
    print(f"strictly evaluating {len(tasks)} extrema-scan points", flush=True)
    with Pool(scan.NPROC, initializer=_init_worker) as pool:
        for i, row in enumerate(pool.imap_unordered(_evaluate, tasks, chunksize=4), 1):
            records.append(row)
            if i % 100 == 0 or i == len(tasks):
                print(f"{i}/{len(tasks)} points ({time.time()-t0:.0f}s)", flush=True)
    if any(r.get("error") for r in records):
        raise RuntimeError(f"strict extrema scan had {sum(bool(r.get('error')) for r in records)} errors")
    records.sort(key=lambda r: (r["VSO"], r["U"]))
    summary = scan.analyze(records)
    with open(os.path.join(OUT_DIR, "strict_family_data.json"), encoding="utf-8") as f:
        family = json.load(f)["uvso_families"]
    soc_data = _build_strict_soc_data(family)
    soc_profiles = scan.analyze_soc_data(soc_data)
    provenance = dict(
        stage="strict-policy extrema scan",
        policy="strict",
        band_rule="analytic continuum bare-band interval with 1e-8 t edge tolerance",
        grid=dict(U_min=min(scan.U_GRID), U_max=max(scan.U_GRID), dU=2.5,
                  VSO=scan.VSO_GRID, n_U=len(scan.U_GRID), n_VSO=len(scan.VSO_GRID)),
        Delta_t=scan.DELTA_T, n_fill=scan.N_FILL, hbar_omega=scan.OMEGA,
        fractions=list(scan.FRACS), absolute_T=list(scan.ABS_T),
        Nk_stage1=scan.NK_STAGE1, Nk_mu=scan.NK_MU, DOS_Nk=1200, nell=scan.NELL,
        n_records=len(records), n_stage1_sc=sum(r["Tc"] > 0 for r in records),
        n_no_inband_root=sum(bool(r.get("no_inband_root")) for r in records),
        n_stage2_inband=sum(bool(v.get("in_band")) for r in records
                            for v in list(r["frac"].values()) + list(r["fixed_T"].values())),
        caveat="This strict criterion does not compare free energies among in-band roots.",
        utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    )
    output = dict(provenance=provenance, summary=summary, soc_profiles=soc_profiles,
                  soc_data=soc_data, data=records)
    with open(os.path.join(OUT_DIR, "strict_extrema_scan.json"), "w", encoding="utf-8") as f:
        json.dump(output, f, indent=1, allow_nan=False)
    _save_strict_scan_figures(records, summary, soc_data)
    print(f"strict extrema scan complete in {time.time()-t0:.1f}s; "
          f"in-band stage1 roots={provenance['n_stage1_sc']}/{len(records)}", flush=True)


if __name__ == "__main__":
    main()
