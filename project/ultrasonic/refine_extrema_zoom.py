"""Targeted refinement and unsmoothed diagnostics for candidate W/W_N extrema.

This supplements (and never overwrites) scan_extrema_vs_soc.py outputs. It scans
U=15..40 at 0.5t spacing and V_SO=0.60..0.90 at 0.025t spacing for T/Tc=0.3,
0.5 under the relaxed chemical-potential policy. Candidate peaks are then
recomputed with doubled Stage-1, chemical-potential, DOS, and absorption grids.
The fixed-U SOC finite-difference plot uses the existing 0.025-spaced family.
Finite differences are shown without smoothing or interpolation.

Run from project/ultrasonic:
  PYTHONPATH=/tmp/arena-analysis-deps python3 refine_extrema_zoom.py
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import solvers as s1sol
from physics_ua import SCState
from parameters_ua import OMEGA
from dos_ua import NormalDOS
from absorption import S_scattering, S_normal
import scan_extrema_vs_soc as base

N_FILL = 1.885
DELTA_T = 4.5
FRACS = (0.3, 0.5)
U_GRID = [float(round(x, 4)) for x in np.arange(15.0, 40.0001, 0.5)]
VSO_GRID = [float(round(x, 4)) for x in np.arange(0.6, 0.9001, 0.025)]
NK_STAGE1 = 80
NK_MU = 200
DOS_NK = 1200
NELL = 20000
HIGH = dict(NK_STAGE1=160, NK_MU=400, DOS_NK=2400, NELL=40000)
NPROC = 2
_DOS = {}


def get_dos(vso, nk):
    key = (round(float(vso), 4), int(nk))
    if key not in _DOS:
        _DOS[key] = NormalDOS(key[0], Nk=key[1])
    return _DOS[key]


def ratio(dos, state, nell):
    sn = S_normal(dos, state["T"], OMEGA, state["mu_N"])
    sc, _ = S_scattering(dos, state["T"], OMEGA, state["D0"], state["DS"],
                         state["mu"], "a_signed", nell=nell)
    return float(sc / sn)


def _evaluate(task):
    U, vso = task
    row = dict(U=float(U), VSO=float(vso), policy=s1sol.MU_WINDOW_POLICY)
    try:
        sol = s1sol.solve_set_core(float(U), DELTA_T, N_FILL, float(vso), Nk=NK_STAGE1)
    except Exception as exc:
        return row | dict(error=f"{type(exc).__name__}: {exc}", Tc=0.0,
                          mu_window=None, frac={})
    tc = float(sol.get("Tc") or 0.0)
    row.update(Tc=tc, Tc_bracketed=bool(sol.get("Tc_bracketed", tc > 0)),
               mu_0=sol.get("mu_0"), mu_window=sol.get("mu_window"),
               D0_0=sol.get("D0_0"), DS_0=sol.get("DS_0"),
               Deff_0=sol.get("Deff_0"))
    if tc <= 0.0:
        return row | dict(frac={f"{f:g}": dict(T=None, ratio=None) for f in FRACS})
    sc = SCState(params=dict(VSO=float(vso), Tc=tc, D0_0=float(sol["D0_0"]),
                             DS_0=float(sol["DS_0"]), Deff_0=float(sol["Deff_0"])),
                 Nk_mu=NK_MU)
    dos = get_dos(vso, DOS_NK)
    values = {}
    for f in FRACS:
        state = sc.at(float(f))
        values[f"{f:g}"] = dict(T=float(state["T"]),
                                ratio=ratio(dos, state, NELL),
                                mu=float(state["mu"]))
    return row | dict(frac=values)


def _evaluate_high(task):
    U, vso = task
    row = dict(U=float(U), VSO=float(vso), policy=s1sol.MU_WINDOW_POLICY)
    try:
        sol = s1sol.solve_set_core(float(U), DELTA_T, N_FILL, float(vso),
                                   Nk=HIGH["NK_STAGE1"])
    except Exception as exc:
        return row | dict(error=f"{type(exc).__name__}: {exc}", Tc=0.0, frac={})
    tc = float(sol.get("Tc") or 0.0)
    row.update(Tc=tc, Tc_bracketed=bool(sol.get("Tc_bracketed", tc > 0)),
               mu_0=sol.get("mu_0"), mu_window=sol.get("mu_window"),
               D0_0=sol.get("D0_0"), DS_0=sol.get("DS_0"),
               Deff_0=sol.get("Deff_0"))
    if tc <= 0.0:
        return row | dict(frac={f"{f:g}": dict(T=None, ratio=None) for f in FRACS})
    sc = SCState(params=dict(VSO=float(vso), Tc=tc, D0_0=float(sol["D0_0"]),
                             DS_0=float(sol["DS_0"]), Deff_0=float(sol["Deff_0"])),
                 Nk_mu=HIGH["NK_MU"])
    dos = get_dos(vso, HIGH["DOS_NK"])
    values = {}
    for f in FRACS:
        state = sc.at(float(f))
        values[f"{f:g}"] = dict(T=float(state["T"]),
                                ratio=ratio(dos, state, HIGH["NELL"]),
                                mu=float(state["mu"]))
    return row | dict(frac=values)


def _prominence(y, i, kind):
    if kind == "max":
        return float(y[i] - max(min(y[:i + 1]), min(y[i:])))
    return float(min(max(y[:i + 1]), max(y[i:])) - y[i])


def summarize(records):
    out = []
    for vso in VSO_GRID:
        rows = [r for r in records if r["VSO"] == vso]
        for f in FRACS:
            k = f"{f:g}"
            e = base.extrema_for_profile(rows,
                    lambda r, k=k: r["frac"].get(k, {}).get("ratio"), x_key="U")
            out.append(dict(VSO=vso, fraction=f, extrema=e,
                            maximum_mu_inside=window_flag(rows, e.get("global_max"), k),
                            minimum_mu_inside=window_flag(rows, e.get("global_min"), k)))
    cutoffs = []
    for f in FRACS:
        hits = [r["VSO"] for r in out if r["fraction"] == f and r["extrema"]["robust_maxima"]]
        first = min(hits) if hits else None
        previous = max((v for v in VSO_GRID if first is not None and v < first), default=None)
        cutoffs.append(dict(fraction=f, first_screened_VSO=first,
                            preceding_grid_VSO=previous,
                            bracket=([previous, first] if previous is not None else None)))
    return out, cutoffs


def window_flag(rows, extremum, fraction_key):
    if extremum is None:
        return None
    row = next((r for r in rows if r["U"] == extremum["x"]), None)
    return (row.get("mu_window") or {}).get("inside") if row else None


def make_high_tasks(summary):
    tasks = set()
    centers = []
    for s in summary:
        e = s["extrema"].get("global_max")
        if e is None:
            continue
        center = float(e["x"])
        centers.append(dict(VSO=s["VSO"], fraction=s["fraction"], center_U=center))
        for u in (center - 0.5, center, center + 0.5):
            if min(U_GRID) <= u <= max(U_GRID):
                tasks.add((float(round(u, 4)), float(s["VSO"])))
    return sorted(tasks), centers


def compare_high(coarse, high, centers):
    coarse_index = {(r["U"], r["VSO"]): r for r in coarse}
    high_index = {(r["U"], r["VSO"]): r for r in high}
    comparisons = []
    for h in high:
        c = coarse_index[(h["U"], h["VSO"])]
        q = dict(U=h["U"], VSO=h["VSO"],
                 Tc_delta=(h["Tc"] - c["Tc"]),
                 D0_0_delta=(h.get("D0_0") - c.get("D0_0")
                             if h.get("D0_0") is not None and c.get("D0_0") is not None else None),
                 DS_0_delta=(h.get("DS_0") - c.get("DS_0")
                             if h.get("DS_0") is not None and c.get("DS_0") is not None else None),
                 mu_0_delta=(h.get("mu_0") - c.get("mu_0")
                             if h.get("mu_0") is not None and c.get("mu_0") is not None else None),
                 ratio={})
        for f in FRACS:
            k = f"{f:g}"
            a, b = h["frac"].get(k, {}).get("ratio"), c["frac"].get(k, {}).get("ratio")
            q["ratio"][k] = dict(coarse=b, refined=a,
                                  abs_delta=(abs(a-b) if a is not None and b is not None else None))
        comparisons.append(q)
    neighborhoods = []
    for c in centers:
        f, v, u0 = c["fraction"], c["VSO"], c["center_U"]
        k = f"{f:g}"
        candidates = [high_index[(float(round(u, 4)), v)]
                      for u in (u0 - 0.5, u0, u0 + 0.5)
                      if (float(round(u, 4)), v) in high_index
                      and high_index[(float(round(u, 4)), v)]["frac"].get(k, {}).get("ratio") is not None]
        if candidates:
            best = max(candidates, key=lambda r: r["frac"][k]["ratio"])
            neighborhoods.append(dict(VSO=v, fraction=f, coarse_center_U=u0,
                                      refined_neighborhood_U=[r["U"] for r in candidates],
                                      refined_local_max_U=best["U"],
                                      refined_local_max_R=best["frac"][k]["ratio"]))
    return comparisons, neighborhoods


def _profile_matrix(records, fraction):
    mat = []
    for v in VSO_GRID:
        rows = sorted((r for r in records if r["VSO"] == v), key=lambda r: r["U"])
        mat.append([r["frac"][f"{fraction:g}"].get("ratio") for r in rows])
    return np.asarray(mat, dtype=float)


def make_figures(records, soc_family):
    plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.25, "font.size": 9})
    selected_vso = (0.6, 0.7, 0.75, 0.8, 0.85, 0.9)
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.95, len(selected_vso)))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for col, f in enumerate(FRACS):
        ax = axes[0, col]
        for v, color in zip(selected_vso, colors):
            rows = sorted((r for r in records if r["VSO"] == v), key=lambda r: r["U"])
            vals = [r["frac"][f"{f:g}"].get("ratio") for r in rows]
            valid = [(r["U"], y) for r, y in zip(rows, vals) if y is not None]
            ax.plot([x for x, _ in valid], [y for _, y in valid], "-", lw=1.3,
                    color=color, label=f"$V_{{SO}}={v:g}t$")
            if valid:
                x, y = max(valid, key=lambda z: z[1])
                ax.scatter([x], [y], marker="o", s=24, color=color, edgecolors="black", zorder=5)
        ax.set_xlim(min(U_GRID), max(U_GRID))
        ax.set_title(f"Raw $W/W_N(U)$, $T/T_c={f:g}$")
        ax.set_xlabel("$U/t$"); ax.set_ylabel("$W/W_N$")
        ax.legend(fontsize=7, ncol=2)

        ax = axes[1, col]
        matrix = _profile_matrix(records, f)
        dmat = np.diff(matrix, axis=1) / 0.5
        uma = (np.asarray(U_GRID[:-1]) + np.asarray(U_GRID[1:])) / 2
        lim = float(np.nanmax(np.abs(dmat)))
        mesh = ax.pcolormesh(VSO_GRID, uma, dmat.T, shading="auto", cmap="coolwarm",
                             norm=TwoSlopeNorm(vmin=-lim, vcenter=0.0, vmax=lim))
        fig.colorbar(mesh, ax=ax, label=r"finite difference $\Delta(W/W_N)/\Delta(U/t)$")
        ax.set_title(f"Unsmoothed $U$-finite-difference map, $T/T_c={f:g}$")
        ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("$U/t$ (cell midpoint)")
    fig.suptitle("Targeted relaxed-policy scan, $15\u2264U/t\u226440$, $0.60\u2264V_{SO}/t\u22640.90$\n"
                 "Raw curves and first differences; no smoothing/interpolation. Markers show sampled maxima in the zoom window")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(HERE, "figures", "refined_extrema_zoom_U.png"), dpi=180)
    plt.close(fig)

    # Existing detailed fixed-U family: plot adjacent-grid finite differences, not a fit.
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    selected_u = (30.0, 45.0, 75.0, 95.0, 115.0)
    cmap = plt.get_cmap("plasma")
    for ax, f in zip(axes, FRACS):
        for i, U in enumerate(selected_u):
            rows = sorted((r for r in soc_family if r["U"] == U), key=lambda r: r["VSO"])
            vs = np.asarray([r["VSO"] for r in rows], dtype=float)
            ys = np.asarray([r["ratio"][FRACS.index(f)] for r in rows], dtype=float)
            xmid = 0.5 * (vs[:-1] + vs[1:])
            dydv = np.diff(ys) / np.diff(vs)
            ax.plot(xmid, dydv, "-o", ms=2.5, lw=1.1, color=cmap(i / 4), label=f"$U={U:g}t$")
        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(f"$\u0394(W/W_N)/\u0394(V_{{SO}}/t)$, $T/T_c={f:g}$")
        ax.set_xlabel("$V_{SO}/t$ (cell midpoint)"); ax.set_ylabel("finite difference")
        ax.legend(fontsize=8, ncol=2)
    fig.suptitle("Fixed-U SOC profiles: adjacent 0.025t finite differences, unsmoothed")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(HERE, "figures", "refined_SOC_finite_differences.png"), dpi=180)
    plt.close(fig)


def main():
    if s1sol.MU_WINDOW_POLICY != "relaxed":
        raise RuntimeError(f"expected relaxed policy, received {s1sol.MU_WINDOW_POLICY!r}")
    tasks = [(u, v) for v in VSO_GRID for u in U_GRID]
    t0 = time.time()
    records = []
    with Pool(NPROC) as pool:
        for i, row in enumerate(pool.imap_unordered(_evaluate, tasks, chunksize=4), 1):
            records.append(row)
            if i % 100 == 0 or i == len(tasks):
                print(f"coarse {i}/{len(tasks)} ({time.time()-t0:.0f}s)", flush=True)
    records.sort(key=lambda r: (r["VSO"], r["U"]))
    if any(r.get("error") for r in records):
        raise RuntimeError(f"coarse scan had {sum(bool(r.get('error')) for r in records)} errors")
    summary, cutoffs = summarize(records)

    high_tasks, centers = make_high_tasks(summary)
    high = []
    with Pool(NPROC) as pool:
        for i, row in enumerate(pool.imap_unordered(_evaluate_high, high_tasks, chunksize=2), 1):
            high.append(row)
            if i % 20 == 0 or i == len(high_tasks):
                print(f"refined {i}/{len(high_tasks)} ({time.time()-t0:.0f}s)", flush=True)
    high.sort(key=lambda r: (r["VSO"], r["U"]))
    if any(r.get("error") for r in high):
        raise RuntimeError(f"high-resolution checks had {sum(bool(r.get('error')) for r in high)} errors")
    comparisons, neighborhoods = compare_high(records, high, centers)

    family_path = os.path.join(HERE, "results", "uvso_families.json")
    with open(family_path) as f:
        family = json.load(f)
    fam = {(float(r["U"]), float(r["VSO"])): r for r in family}
    cross_diffs = []
    for row in records:
        old = fam.get((row["U"], row["VSO"]))
        if old is None:
            continue
        for frac in FRACS:
            cross_diffs.append(abs(row["frac"][f"{frac:g}"]["ratio"] - old["ratio"][FRACS.index(frac)]))

    # Raw finite-difference data for selected fixed-U SOC profiles.
    soc_diff_data = []
    for U in (30.0, 45.0, 75.0, 95.0, 115.0):
        rows = sorted((r for r in family if r["U"] == U), key=lambda r: r["VSO"])
        vs = [float(r["VSO"]) for r in rows]
        for f in FRACS:
            values = [float(r["ratio"][FRACS.index(f)]) for r in rows]
            soc_diff_data.append(dict(U=U, fraction=f,
                                      VSO_midpoints=[0.5 * (a+b) for a,b in zip(vs[:-1],vs[1:])],
                                      dR_dVSO=[(b-a)/(vb-va) for a,b,va,vb in zip(values[:-1],values[1:],vs[:-1],vs[1:])]))

    provenance = dict(stage="targeted-extrema-refinement", policy=s1sol.MU_WINDOW_POLICY,
                      grid=dict(U_min=min(U_GRID), U_max=max(U_GRID), dU=0.5,
                                VSO_min=min(VSO_GRID), VSO_max=max(VSO_GRID), dVSO=0.025,
                                n_U=len(U_GRID), n_VSO=len(VSO_GRID)),
                      fractions=list(FRACS), N_fill=N_FILL, Delta_t=DELTA_T, hbar_omega=OMEGA,
                      coarse_resolution=dict(Nk=NK_STAGE1, Nk_mu=NK_MU, DOS_Nk=DOS_NK, nell=NELL),
                      refined_resolution=HIGH, n_coarse=len(records), n_high=len(high),
                      n_high_comparisons=len(comparisons), n_out_of_band=sum(
                          (r.get("mu_window") or {}).get("inside") is False for r in records),
                      reference_crosscheck=dict(reference="uvso_families.json", n_values=len(cross_diffs),
                                                max_abs_difference=max(cross_diffs) if cross_diffs else None),
                      caveat="Combined resolution comparison, not an isolated error budget; relaxed roots are tagged, not stability-certified.",
                      utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    result = dict(provenance=provenance, summary=summary, screened_cutoff_brackets=cutoffs,
                  high_resolution_comparisons=comparisons,
                  refined_neighborhood_maxima=neighborhoods,
                  soc_finite_differences=soc_diff_data,
                  coarse_data=records, high_resolution_data=high)
    outpath = os.path.join(HERE, "results", "refined_extrema_zoom.json")
    with open(outpath, "w") as f:
        json.dump(result, f, indent=1, allow_nan=False)
    csvpath = os.path.join(HERE, "results", "refined_extrema_zoom_summary.csv")
    with open(csvpath, "w", newline="") as f:
        fields = ["fraction", "VSO", "n", "U_max", "R_max", "max_mu_inside", "U_min", "R_min",
                  "min_mu_inside", "robust_maxima", "robust_minima"]
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for s in summary:
            e = s["extrema"]
            mx, mn = e.get("global_max") or {}, e.get("global_min") or {}
            w.writerow(dict(fraction=s["fraction"], VSO=s["VSO"], n=e["n"],
                            U_max=mx.get("x"), R_max=mx.get("value"), max_mu_inside=s["maximum_mu_inside"],
                            U_min=mn.get("x"), R_min=mn.get("value"), min_mu_inside=s["minimum_mu_inside"],
                            robust_maxima=json.dumps(e.get("robust_maxima", []), separators=(",", ":")),
                            robust_minima=json.dumps(e.get("robust_minima", []), separators=(",", ":"))))
    make_figures(records, family)
    print(f"Done in {time.time()-t0:.0f}s; coarse={len(records)}, high={len(high)}, "
          f"reference max |dR|={provenance['reference_crosscheck']['max_abs_difference']}")


if __name__ == "__main__":
    main()
