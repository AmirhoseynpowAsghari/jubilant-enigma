"""Scan W/W_N extrema versus U and V_SO under the relaxed D1 policy.

This is a targeted follow-up to the U- and V_SO-family figures. It evaluates
both fixed reduced temperatures (T/Tc = 0.3, 0.5, 0.7, 0.9, 0.999) and fixed
physical temperatures (T = 0.1t, 0.2t), so these two ensembles are not mixed.
Out-of-band roots are retained and tagged; this script makes no phase-stability
claim about them.

Grid: U=15..150 in 2.5t steps; V_SO=0..0.9 in 0.1t steps, plus 0.4082t.
The grid provides cutoff brackets, not a continuum proof. Global extrema include
endpoints. Interior extrema are also reported with their prominence; only those
above the explicit screening threshold are called robust candidates.

Outputs: results/extrema_vs_soc_scan.json, results/extrema_vs_soc_summary.csv,
         figures/extrema_by_U.png, figures/soc_scans_by_U.png,
         figures/extrema_tracks_by_soc.png
Run: PYTHONPATH=../sc_params python3 scan_extrema_vs_soc.py
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import solvers as s1sol
from physics_ua import SCState
from parameters_ua import OMEGA
from dos_ua import NormalDOS
from absorption import S_scattering, S_normal

N_FILL = 1.885
DELTA_T = 4.5
NELL = 20000
NK_STAGE1 = 80
NK_MU = 200
FRACS = (0.3, 0.5, 0.7, 0.9, 0.999)
ABS_T = (0.1, 0.2)
U_GRID = [float(round(x, 4)) for x in np.arange(15.0, 150.0001, 2.5)]
VSO_GRID = sorted(set([float(round(x, 4)) for x in np.arange(0.0, 0.9001, 0.1)] + [0.4082]))
NPROC = 2
# Screening rule for local extrema. Raw local extrema and prominence are kept too.
PROM_ABS = 5e-4
PROM_FRAC = 0.02
_DOS = {}


def get_dos(vso):
    key = round(float(vso), 4)
    if key not in _DOS:
        _DOS[key] = NormalDOS(key, Nk=1200)
    return _DOS[key]


def _ratio(dos, state):
    sn = S_normal(dos, state["T"], OMEGA, state["mu_N"])
    sc, _ = S_scattering(dos, state["T"], OMEGA, state["D0"], state["DS"],
                         state["mu"], "a_signed", nell=NELL)
    return float(sc / sn)


def evaluate(task):
    U, vso = task
    base = dict(U=float(U), VSO=float(vso), policy=s1sol.MU_WINDOW_POLICY)
    try:
        s = s1sol.solve_set(float(U), DELTA_T, N_FILL, float(vso), Nk=NK_STAGE1)
    except Exception as exc:
        return base | dict(Tc=0.0, error=f"{type(exc).__name__}: {exc}",
                           mu_window=None, frac={}, fixed_T={})
    tc = float(s.get("Tc") or 0.0)
    base.update(Tc=tc, mu_0=s.get("mu_0"), mu_window=s.get("mu_window"),
                Tc_bracketed=bool(s.get("Tc_bracketed", tc > 0.0)),
                D0_0=s.get("D0_0"), DS_0=s.get("DS_0"), Deff_0=s.get("Deff_0"))
    if tc <= 0.0:
        return base | dict(frac={}, fixed_T={})
    sc = SCState(params=dict(VSO=float(vso), Tc=tc, D0_0=float(s["D0_0"]),
                             DS_0=float(s["DS_0"]), Deff_0=float(s["Deff_0"])),
                 Nk_mu=NK_MU)
    dos = get_dos(vso)
    fr = {}
    for f in FRACS:
        state = sc.at(float(f))
        fr[f"{f:g}"] = dict(T=float(state["T"]), ratio=_ratio(dos, state),
                            mu=float(state["mu"]), mu_N=float(state["mu_N"]))
    abs_t = {}
    for temp in ABS_T:
        key = f"{temp:.2f}"
        if temp >= tc:
            abs_t[key] = dict(T=float(temp), ratio=None, superconducting=False,
                              frac=None, mu=None, mu_N=None)
        else:
            f = float(temp / tc)
            state = sc.at(f)
            abs_t[key] = dict(T=float(temp), ratio=_ratio(dos, state),
                              superconducting=True, frac=f,
                              mu=float(state["mu"]), mu_N=float(state["mu_N"]))
    return base | dict(frac=fr, fixed_T=abs_t)


def _prominence(y, i, kind):
    if kind == "max":
        return float(y[i] - max(min(y[:i + 1]), min(y[i:])))
    return float(min(max(y[:i + 1]), max(y[i:])) - y[i])


def extrema_for_profile(points, value_getter, x_key="U"):
    """Global extrema plus raw/prominence-screened interior extrema for ordered points."""
    pts = sorted([p for p in points if value_getter(p) is not None], key=lambda p: p[x_key])
    if not pts:
        return dict(n=0, global_max=None, global_min=None,
                    local_maxima=[], local_minima=[], robust_maxima=[], robust_minima=[])
    x = [float(p[x_key]) for p in pts]
    y = [float(value_getter(p)) for p in pts]
    span = max(y) - min(y)
    threshold = max(PROM_ABS, PROM_FRAC * span)
    imax, imin = max(range(len(y)), key=y.__getitem__), min(range(len(y)), key=y.__getitem__)
    global_max = dict(x=x[imax], value=y[imax], endpoint=imax in (0, len(y) - 1),
                      below_1e_14=(y[imax] <= 1e-14))
    global_min = dict(x=x[imin], value=y[imin], endpoint=imin in (0, len(y) - 1),
                      below_1e_14=(y[imin] <= 1e-14))
    maxima, minima = [], []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1]:
            maxima.append(dict(x=x[i], value=y[i], prominence=_prominence(y, i, "max")))
        if y[i] < y[i - 1] and y[i] <= y[i + 1]:
            minima.append(dict(x=x[i], value=y[i], prominence=_prominence(y, i, "min")))
    robust_max = [p for p in maxima if p["prominence"] >= threshold]
    robust_min = [p for p in minima if p["prominence"] >= threshold]
    return dict(n=len(pts), x_min=x[0], x_max=x[-1], y_min=min(y), y_max=max(y),
                prominence_threshold=threshold, global_max=global_max, global_min=global_min,
                local_maxima=maxima, local_minima=minima,
                robust_maxima=robust_max, robust_minima=robust_min)


def contiguous_intervals(points, predicate, x_key="U"):
    points = sorted(points, key=lambda p: p[x_key])
    out, start, last = [], None, None
    for p in points:
        if predicate(p):
            if start is None:
                start = float(p[x_key])
            last = float(p[x_key])
        elif start is not None:
            out.append([start, last])
            start = last = None
    if start is not None:
        out.append([start, last])
    return out


def cutoff_intervals(records, temp, x_key="U"):
    """Report Tc>T intervals and linearly interpolated Tc=T boundary brackets."""
    pts = sorted(records, key=lambda p: p[x_key])
    out, start, last = [], None, None
    for i, p in enumerate(pts):
        active = float(p.get("Tc") or 0.0) > temp
        if active and start is None:
            start = i
        if not active and start is not None:
            end = i - 1
            out.append(_cutoff_interval(pts, start, end, temp, x_key))
            start = None
    if start is not None:
        out.append(_cutoff_interval(pts, start, len(pts) - 1, temp, x_key))
    return out


def _cutoff_interval(pts, i0, i1, temp, x_key):
    def cross(i_active, i_other):
        a, b = pts[i_active], pts[i_other]
        ta, tb = float(a.get("Tc") or 0.0), float(b.get("Tc") or 0.0)
        xa, xb = float(a[x_key]), float(b[x_key])
        if ta == tb:
            estimate = 0.5 * (xa + xb)
        else:
            estimate = xa + (temp - ta) * (xb - xa) / (tb - ta)
        return dict(estimate=estimate, bracket=sorted([xa, xb]), method="linear Tc crossing")
    left = cross(i0, i0 - 1) if i0 > 0 else dict(estimate=None, bracket=None, method="scan boundary")
    right = cross(i1, i1 + 1) if i1 < len(pts) - 1 else dict(estimate=None, bracket=None, method="scan boundary")
    return dict(sampled_U=[float(pts[i0][x_key]), float(pts[i1][x_key])],
                lower_cutoff=left, upper_cutoff=right)


def analyze(records):
    summary = []
    for v in VSO_GRID:
        rows = [r for r in records if r["VSO"] == v]
        for f in FRACS:
            k = f"{f:g}"
            e = extrema_for_profile(rows, lambda r, k=k: r["frac"].get(k, {}).get("ratio"))
            summary.append(dict(ensemble="fixed_T_over_Tc", value=float(f), VSO=v,
                                extrema=e, robust_peak_count=len(e["robust_maxima"]),
                                robust_valley_count=len(e["robust_minima"])))
        for temp in ABS_T:
            k = f"{temp:.2f}"
            e = extrema_for_profile(rows, lambda r, k=k: r["fixed_T"].get(k, {}).get("ratio"))
            valid_intervals = contiguous_intervals(
                rows, lambda r, k=k: bool(r["fixed_T"].get(k, {}).get("superconducting")))
            summary.append(dict(ensemble="fixed_absolute_T", value=float(temp), VSO=v,
                                extrema=e, superconducting_U_intervals=valid_intervals,
                                Tc_crossing_intervals=cutoff_intervals(rows, temp),
                                robust_peak_count=len(e["robust_maxima"]),
                                robust_valley_count=len(e["robust_minima"])))
    return summary


def build_soc_data(family_records):
    """Use the companion 0.025-spaced V_SO family for detailed fixed-U profiles."""
    selected_u = {30.0, 45.0, 75.0, 95.0, 115.0}
    data = []
    for r in family_records:
        if float(r["U"]) not in selected_u:
            continue
        row = {"U": float(r["U"]), "VSO": float(r["VSO"]), "Tc": float(r["Tc"]),
               "frac": {f"{f:g}": {"ratio": float(r["ratio"][i])}
                        for i, f in enumerate(FRACS)}, "fixed_T": {}}
        sc = SCState(params=dict(VSO=float(r["VSO"]), Tc=float(r["Tc"]),
                                 D0_0=float(r["D0_0"]), DS_0=float(r["DS_0"]),
                                 Deff_0=float(r["Deff_0"])), Nk_mu=NK_MU)
        dos = get_dos(r["VSO"])
        for temp in ABS_T:
            k = f"{temp:.2f}"
            if temp >= row["Tc"]:
                row["fixed_T"][k] = dict(T=float(temp), ratio=None, superconducting=False)
            else:
                state = sc.at(float(temp / row["Tc"]))
                row["fixed_T"][k] = dict(T=float(temp), ratio=_ratio(dos, state),
                                         superconducting=True, frac=float(temp / row["Tc"]))
        data.append(row)
    data.sort(key=lambda r: (r["U"], r["VSO"]))
    return data


def analyze_soc_data(data):
    profiles = []
    for U in (30.0, 45.0, 75.0, 95.0, 115.0):
        rows = [r for r in data if r["U"] == U]
        for f in FRACS:
            k = f"{f:g}"
            e = extrema_for_profile(rows, lambda r, k=k: r["frac"].get(k, {}).get("ratio"), x_key="VSO")
            profiles.append(dict(ensemble="fixed_T_over_Tc", U=U, value=float(f), extrema=e))
        for temp in ABS_T:
            k = f"{temp:.2f}"
            e = extrema_for_profile(rows, lambda r, k=k: r["fixed_T"].get(k, {}).get("ratio"), x_key="VSO")
            profiles.append(dict(ensemble="fixed_absolute_T", U=U, value=float(temp), extrema=e))
    return profiles


def _mark_extrema(ax, profile, getter, color, x_key="U"):
    e = extrema_for_profile(profile, getter, x_key=x_key)
    for name, marker in (("global_max", "o"), ("global_min", "s")):
        q = e[name]
        if q:
            ax.scatter([q["x"]], [q["value"]], marker=marker, s=42,
                       facecolors=(color if not q["endpoint"] else "none"),
                       edgecolors=color, linewidths=1.2, zorder=5)


def make_figures(records, summary, soc_data):
    plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.25, "font.size": 9})
    palette = {0.0: "#1f77b4", 0.2: "#ff7f0e", 0.4: "#2ca02c", 0.6: "#d62728",
               0.8: "#9467bd", 0.4082: "#8c564b", 0.1: "#17becf", 0.3: "#bcbd22",
               0.5: "#e377c2", 0.7: "#7f7f7f", 0.9: "#393b79"}
    selected_vso = (0.0, 0.2, 0.4082, 0.4, 0.6, 0.8, 0.9)
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    panels = [("fixed_T_over_Tc", 0.5, axes[0, 0], "$T/T_c=0.5$"),
              ("fixed_T_over_Tc", 0.9, axes[0, 1], "$T/T_c=0.9$"),
              ("fixed_absolute_T", 0.1, axes[1, 0], "$T=0.1t$"),
              ("fixed_absolute_T", 0.2, axes[1, 1], "$T=0.2t$")]
    for ensemble, temp, ax, title in panels:
        for v in selected_vso:
            rows = sorted((r for r in records if r["VSO"] == v), key=lambda r: r["U"])
            if not rows:
                continue
            if ensemble == "fixed_T_over_Tc":
                k = f"{temp:g}"
                ys = [r["frac"].get(k, {}).get("ratio") for r in rows]
            else:
                k = f"{temp:.2f}"
                ys = [r["fixed_T"].get(k, {}).get("ratio") for r in rows]
            xs = [r["U"] for r in rows]
            valid = [(x, y) for x, y in zip(xs, ys) if y is not None]
            if not valid:
                continue
            color = palette[v]
            ax.plot([x for x, _ in valid], [y for _, y in valid], "-", lw=1.25,
                    color=color, label=f"$V_{{SO}}={v:g}t$")
            _mark_extrema(ax, rows, lambda r, k=k, ensemble=ensemble, temp=temp:
                          r["frac"].get(k, {}).get("ratio") if ensemble == "fixed_T_over_Tc"
                          else r["fixed_T"].get(k, {}).get("ratio"), color)
        ax.set_title(title)
        ax.set_xlabel("$U/t$")
        ax.set_ylabel("$W/W_N$")
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Global maxima (circles) and minima (squares) vs U; open markers are valid-profile endpoints\n"
                 "Some endpoints are T=Tc cutoffs; values <=1e-14 are flagged, not clipped\n"
                 "Relaxed policy: out-of-band chemical potentials are plotted but tagged in the JSON")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(HERE, "figures", "extrema_by_U.png"), dpi=170)
    plt.close(fig)

    # W/W_N vs SOC at the five U values in the existing companion-family figure.
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    panels = [("fixed_T_over_Tc", 0.5, axes[0, 0], "$T/T_c=0.5$"),
              ("fixed_T_over_Tc", 0.9, axes[0, 1], "$T/T_c=0.9$"),
              ("fixed_absolute_T", 0.1, axes[1, 0], "$T=0.1t$"),
              ("fixed_absolute_T", 0.2, axes[1, 1], "$T=0.2t$")]
    cmap = plt.get_cmap("viridis")
    for ensemble, temp, ax, title in panels:
        for idx, U in enumerate((30.0, 45.0, 75.0, 95.0, 115.0)):
            rows = sorted((r for r in soc_data if r["U"] == U), key=lambda r: r["VSO"])
            k = f"{temp:g}" if ensemble == "fixed_T_over_Tc" else f"{temp:.2f}"
            ys = [(r["frac"].get(k, {}).get("ratio") if ensemble == "fixed_T_over_Tc"
                  else r["fixed_T"].get(k, {}).get("ratio")) for r in rows]
            valid = [(r["VSO"], y) for r, y in zip(rows, ys) if y is not None]
            if not valid:
                continue
            color = cmap(0.1 + 0.8 * idx / 4)
            ax.plot([x for x, _ in valid], [y for _, y in valid], "-o", ms=2.6, lw=1.2,
                    color=color, label=f"$U={U:g}t$")
            _mark_extrema(ax, rows, lambda r, k=k, ensemble=ensemble:
                          r["frac"].get(k, {}).get("ratio") if ensemble == "fixed_T_over_Tc"
                          else r["fixed_T"].get(k, {}).get("ratio"), color, x_key="VSO")
        ax.set_title(title)
        ax.set_xlabel("$V_{SO}/t$")
        ax.set_ylabel("$W/W_N$")
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle("W/W_N vs V_SO at fixed U (0.025t spacing): global maxima/minima marked\n"
                 "Open symbols are valid-profile endpoints; absolute-T curves include only T<Tc")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(HERE, "figures", "soc_scans_by_U.png"), dpi=170)
    plt.close(fig)

    # Extrema values and U positions as SOC changes for each ensemble/temperature.
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    panels = [("fixed_T_over_Tc", 0.5, axes[0, 0], "$T/T_c=0.5$"),
              ("fixed_T_over_Tc", 0.9, axes[0, 1], "$T/T_c=0.9$"),
              ("fixed_absolute_T", 0.1, axes[1, 0], "$T=0.1t$"),
              ("fixed_absolute_T", 0.2, axes[1, 1], "$T=0.2t$")]
    for ensemble, temp, ax, title in panels:
        for which, color, marker in (("global_max", "tab:red", "o"),
                                     ("global_min", "tab:blue", "s")):
            rows = [r for r in summary if r["ensemble"] == ensemble and abs(r["value"] - temp) < 1e-10]
            rows.sort(key=lambda r: r["VSO"])
            xs, ys, labels = [], [], []
            for r in rows:
                e = r["extrema"].get(which)
                if e is not None:
                    xs.append(r["VSO"]); ys.append(e["value"]); labels.append(e["endpoint"])
            if xs:
                ax.plot(xs, ys, f"-{marker}", color=color, lw=1.2,
                        label="maximum" if which == "global_max" else "minimum")
                for x, y, endpoint in zip(xs, ys, labels):
                    if endpoint:
                        ax.scatter([x], [y], marker=marker, s=45, facecolors="none",
                                   edgecolors=color, linewidths=1.2, zorder=5)
        ax.set_title(title)
        ax.set_xlabel("$V_{SO}/t$")
        ax.set_ylabel("extremal $W/W_N$ over U")
        ax.legend(fontsize=8)
    fig.suptitle("Extremal attenuation values vs SOC; open markers denote U-boundary extrema")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(HERE, "figures", "extrema_tracks_by_soc.png"), dpi=170)
    plt.close(fig)


def main():
    if s1sol.MU_WINDOW_POLICY != "relaxed":
        raise RuntimeError(f"expected relaxed D1 default, got {s1sol.MU_WINDOW_POLICY!r}")
    tasks = [(u, v) for v in VSO_GRID for u in U_GRID]
    t0 = time.time()
    records = []
    with Pool(NPROC) as pool:
        for i, r in enumerate(pool.imap_unordered(evaluate, tasks, chunksize=4), 1):
            records.append(r)
            if i % 100 == 0 or i == len(tasks):
                print(f"{i}/{len(tasks)} points ({time.time()-t0:.0f}s)", flush=True)
    records.sort(key=lambda r: (r["VSO"], r["U"]))
    # Load the denser companion family both for exact-grid verification and for
    # fixed-U SOC curves (0.025 spacing rather than the coarse 0.1 scan grid).
    old_path = os.path.join(HERE, "results", "uvso_families.json")
    with open(old_path) as fh:
        old = json.load(fh)
    summary = analyze(records)
    soc_data = build_soc_data(old)
    soc_profiles = analyze_soc_data(soc_data)
    old_index = {(r["U"], r["VSO"]): r for r in old}
    diffs = []
    fixed_t_diffs = []
    soc_index = {(r["U"], r["VSO"]): r for r in soc_data}
    for r in records:
        old_r = old_index.get((r["U"], r["VSO"]))
        if old_r is None:
            continue
        for f in FRACS:
            k = f"{f:g}"
            if k in r["frac"]:
                diffs.append(abs(r["frac"][k]["ratio"] - old_r["ratio"][FRACS.index(f)]))
        soc_r = soc_index.get((r["U"], r["VSO"]))
        if soc_r is not None:
            for temp in ABS_T:
                k = f"{temp:.2f}"
                new_val = r["fixed_T"].get(k, {}).get("ratio")
                dense_val = soc_r["fixed_T"].get(k, {}).get("ratio")
                if new_val is not None and dense_val is not None:
                    fixed_t_diffs.append(abs(new_val - dense_val))
    provenance = dict(stage="extrema-vs-SOC", policy=s1sol.MU_WINDOW_POLICY,
                      grid=dict(U_min=min(U_GRID), U_max=max(U_GRID), dU=2.5, n_U=len(U_GRID),
                                VSO=VSO_GRID, n_VSO=len(VSO_GRID)),
                      Delta_t=DELTA_T, n=N_FILL, hbar_omega=OMEGA, fracs=list(FRACS),
                      absolute_T=list(ABS_T), Nk_stage1=NK_STAGE1, Nk_mu=NK_MU,
                      DOS_Nk=1200, nell=NELL, small_value_flag_threshold=1e-14,
                      prominence_rule=dict(abs=PROM_ABS, frac_of_profile_range=PROM_FRAC),
                      n_records=len(records), n_out_of_band=sum(
                          1 for r in records if (r.get("mu_window") or {}).get("inside") is False),
                      crosscheck=dict(reference="uvso_families.json", n_values=len(diffs),
                                      max_abs_difference=max(diffs) if diffs else None,
                                      fixed_T_n_values=len(fixed_t_diffs),
                                      fixed_T_max_abs_difference=max(fixed_t_diffs) if fixed_t_diffs else None),
                      utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"))
    out = dict(provenance=provenance, summary=summary, soc_profiles=soc_profiles,
               soc_data=soc_data, data=records)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "extrema_vs_soc_scan.json"), "w") as fh:
        json.dump(out, fh, indent=1, allow_nan=False)
    with open(os.path.join(HERE, "results", "extrema_vs_soc_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["ensemble", "value", "VSO", "n", "U_max", "R_max",
                                            "max_endpoint", "U_min", "R_min", "min_endpoint",
                                            "robust_maxima", "robust_minima", "SC_U_intervals",
                                            "Tc_crossing_U_brackets"], lineterminator="\n")
        w.writeheader()
        for s in summary:
            e = s["extrema"]
            w.writerow(dict(ensemble=s["ensemble"], value=s["value"], VSO=s["VSO"], n=e["n"],
                            U_max=e.get("global_max", {}).get("x"),
                            R_max=e.get("global_max", {}).get("value"),
                            max_endpoint=e.get("global_max", {}).get("endpoint"),
                            U_min=e.get("global_min", {}).get("x"),
                            R_min=e.get("global_min", {}).get("value"),
                            min_endpoint=e.get("global_min", {}).get("endpoint"),
                            robust_maxima=json.dumps(e.get("robust_maxima", []), separators=(",", ":")),
                            robust_minima=json.dumps(e.get("robust_minima", []), separators=(",", ":")),
                            SC_U_intervals=json.dumps(s.get("superconducting_U_intervals", [])),
                            Tc_crossing_U_brackets=json.dumps(s.get("Tc_crossing_intervals", []),
                                                              separators=(",", ":"))))
    make_figures(records, summary, soc_data)
    print(f"Done in {time.time()-t0:.0f}s; {len(records)} points; "
          f"max cross-check |ΔR|={provenance['crosscheck']['max_abs_difference']}")


if __name__ == "__main__":
    main()
