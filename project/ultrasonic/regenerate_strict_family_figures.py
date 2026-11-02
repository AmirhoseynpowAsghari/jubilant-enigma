#!/usr/bin/env python3
"""Regenerate the W/W_N parameter-family figures with strict bare-band admissibility.

Existing relaxed-policy inputs are read-only. Stage-1 roots are re-solved for every unique
(U, V_SO) point under the strict policy. If the strict solution matches the archived relaxed
root, archived W/W_N samples are reused; changed strict branches are re-integrated. Stage-2
chemical potentials are checked against the analytic band interval and points outside it are
omitted (not plotted as zero).

New strict outputs go only under results/strict_policy and figures/strict_policy.
Run from project/ultrasonic:
    PYTHONPATH=../sc_params python3 regenerate_strict_family_figures.py
"""
from __future__ import annotations

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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import solvers
from absorption import S_normal, S_scattering
from dos_ua import NormalDOS
from parameters_ua import OMEGA
from physics_ua import SCState
import make_uvso_families as uv
import make_delta_families as delta

OUT_DIR = os.path.join(HERE, "results", "strict_policy")
FIG_DIR = os.path.join(HERE, "figures", "strict_policy")
N_FILL = 1.885
DELTA_T = 4.5
NK_STAGE1 = 80
NK_MU = 200
NELL = 20000
BAND_TOL = solvers.BAND_EDGE_TOL


def _strict_worker_init():
    solvers.MU_WINDOW_POLICY = "strict"


def _solve_strict_stage1(task):
    _strict_worker_init()
    U, vso = task
    try:
        s = solvers.solve_set_core(float(U), DELTA_T, N_FILL, float(vso), Nk=NK_STAGE1)
    except Exception as exc:
        return dict(U=float(U), VSO=float(vso), Tc=0.0, error=f"{type(exc).__name__}: {exc}")
    return dict(U=float(U), VSO=float(vso), Tc=float(s.get("Tc") or 0.0),
                Tc_bracketed=bool(s.get("Tc_bracketed", False)),
                D0_0=float(s.get("D0_0") or 0.0), DS_0=float(s.get("DS_0") or 0.0),
                Deff_0=float(s.get("Deff_0") or 0.0), mu_0=s.get("mu_0"),
                mu_window=s.get("mu_window"), no_inband_root=bool(s.get("no_inband_root", False)),
                diagnostic_unrestricted=s.get("diagnostic_unrestricted"), policy="strict")


def _same_root(strict_row, relaxed_row):
    if not relaxed_row or float(strict_row.get("Tc") or 0.0) <= 0.0:
        return False
    keys = ("Tc", "D0_0", "DS_0", "Deff_0", "mu_0")
    for key in keys:
        a, b = strict_row.get(key), relaxed_row.get(key)
        if a is None or b is None:
            return False
        if abs(float(a) - float(b)) > 2e-7 * max(1.0, abs(float(a)), abs(float(b))):
            return False
    return True


def _inside(mu, band):
    return mu is not None and band[0] - BAND_TOL <= float(mu) <= band[1] + BAND_TOL


def _stage2_curve(stage1, fractions, old_curve=None, old_x=None, old_y=None):
    """Return ratios and stage-2 mu flags for a strict stage-1 state.

    Reuses archived ratios only when the strict and relaxed stage-1 solutions agree. In either
    case, the finite-temperature SC chemical potential is screened by the same strict band rule.
    """
    vso = float(stage1["VSO"])
    band = solvers.bare_band_window(np.empty(0), vso)
    scstate = SCState(params=dict(VSO=vso, Tc=float(stage1["Tc"]),
                                 D0_0=float(stage1["D0_0"]), DS_0=float(stage1["DS_0"]),
                                 Deff_0=float(stage1["Deff_0"])),
                      Nk_mu=NK_MU, n_tgt=N_FILL)
    same = old_curve is not None
    ratios, mus, accepted = [], [], []
    dos = None
    for i, frac in enumerate(fractions):
        state = scstate.at(float(frac))
        keep = _inside(state["mu"], band)
        mus.append(float(state["mu"]))
        accepted.append(bool(keep))
        if not keep:
            ratios.append(None)
            continue
        if same and old_curve is not None and i < len(old_curve):
            ratios.append(float(old_curve[i]))
            continue
        if dos is None:
            dos = NormalDOS(vso, Nk=1200)
        sn = S_normal(dos, state["T"], OMEGA, state["mu_N"])
        ss, _ = S_scattering(dos, state["T"], OMEGA, state["D0"], state["DS"],
                             state["mu"], "a_signed", nell=NELL)
        ratios.append(float(ss / sn))
    return ratios, mus, accepted


def _strict_uv_row(stage1, relaxed):
    if stage1["Tc"] <= 0.0:
        return stage1 | dict(ratio=[None] * len(uv.FRACS), stage2_mu=[], stage2_inside=[])
    old_ratio = relaxed.get("ratio") if _same_root(stage1, relaxed) else None
    ratios, mus, inside = _stage2_curve(stage1, uv.FRACS, old_curve=old_ratio)
    return stage1 | dict(ratio=ratios, stage2_mu=mus, stage2_inside=inside,
                          reused_relaxed_ratio=old_ratio is not None)


def _strict_delta_row(stage1, relaxed):
    if stage1["Tc"] <= 0.0:
        n = len(delta.FR_GRID)
        return stage1 | dict(frac=list(delta.FR_GRID), ratio=[None] * n,
                             D0=[None] * n, Emin=[None] * n, T=[None] * n,
                             deffn=[None] * n, stage2_mu=[], stage2_inside=[])
    same = _same_root(stage1, relaxed)
    old_ratio = relaxed.get("ratio") if same else None
    if same:
        ratios, mus, inside = _stage2_curve(stage1, delta.FR_GRID, old_curve=old_ratio)
        result = dict(stage1)
        for key in ("D0", "Emin", "T", "deffn"):
            result[key] = list(relaxed[key])
        result.update(frac=list(delta.FR_GRID), ratio=ratios, stage2_mu=mus,
                      stage2_inside=inside, reused_relaxed_ratio=True)
        return result

    band = solvers.bare_band_window(np.empty(0), float(stage1["VSO"]))
    scstate = SCState(params=dict(VSO=float(stage1["VSO"]), Tc=float(stage1["Tc"]),
                                 D0_0=float(stage1["D0_0"]), DS_0=float(stage1["DS_0"]),
                                 Deff_0=float(stage1["Deff_0"])),
                      Nk_mu=NK_MU, n_tgt=N_FILL)
    dos = NormalDOS(float(stage1["VSO"]), Nk=1200)
    eps, _ = uv.build_kgrid(float(stage1["VSO"]), float(stage1["U"]), DELTA_T, Nk=NK_STAGE1)
    ratios, d0s, emins, temps, deffns, mus, inside = [], [], [], [], [], [], []
    for frac in delta.FR_GRID:
        state = scstate.at(float(frac))
        keep = _inside(state["mu"], band)
        mus.append(float(state["mu"]))
        inside.append(bool(keep))
        if not keep:
            ratios.append(None); d0s.append(None); emins.append(None); temps.append(float(state["T"])); deffns.append(None)
            continue
        sn = S_normal(dos, state["T"], OMEGA, state["mu_N"])
        ss, _ = S_scattering(dos, state["T"], OMEGA, state["D0"], state["DS"],
                             state["mu"], "a_signed", nell=NELL)
        p = delta.qp_constants(state["D0"], state["DS"], state["mu"])
        ratios.append(float(ss / sn)); d0s.append(float(state["D0"]))
        emins.append(float(p["Emin"])); temps.append(float(state["T"]))
        deffns.append(float(np.mean(np.abs(state["D0"] - state["DS"] * eps / 4.0))) /
                       float(stage1["Deff_0"]))
    return stage1 | dict(frac=list(delta.FR_GRID), ratio=ratios, D0=d0s, Emin=emins,
                         T=temps, deffn=deffns, stage2_mu=mus, stage2_inside=inside,
                         reused_relaxed_ratio=False)


def _band_status(row):
    if float(row.get("Tc") or 0.0) <= 0.0 or row.get("mu_0") is None:
        return False
    band = solvers.bare_band_window(np.empty(0), float(row["VSO"]))
    return _inside(row["mu_0"], band)


def _plot_uvso(data, panels, panel_kind, x_kind, filename, title):
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), squeeze=False)
    fracs = uv.FRACS
    for i, pv in enumerate(panels):
        ax = axes[i // 3, i % 3]
        for j, frac in enumerate(fracs):
            col = cmap(0.05 + 0.9 * j / (len(fracs) - 1))
            rows = sorted((r for r in data if r[panel_kind] == pv), key=lambda r: r[x_kind])
            xs = [float(r[x_kind]) for r in rows]
            ys = [r["ratio"][j] if r.get("Tc", 0.0) > 0 and r.get("ratio") else np.nan for r in rows]
            if xs:
                ax.plot(xs, ys, "-", lw=1.4, color=col, label=f"$T/T_c$={frac:g}")
                small = [(x, y) for x, y, r in zip(xs, ys, rows)
                         if np.isfinite(y) and float(r.get("D0_0", 1.0)) < 0.05]
                if small:
                    ax.plot([p[0] for p in small], [p[1] for p in small], "s", ms=2.8,
                            mfc="none", mec="tab:red", mew=0.8, ls="none")
        kl = "$U/t$" if x_kind == "U" else "$V_{SO}/t$"
        pl = "$V_{SO}$" if panel_kind == "VSO" else "$U$"
        ax.set_xlabel(kl); ax.set_ylabel("$W/W_N$"); ax.set_title(f"{pl} = {pv:g}", fontsize=10)
        ax.grid(alpha=0.3)
        subset = [r for r in data if r[panel_kind] == pv]
        n_excluded = sum(float(r.get("Tc") or 0.0) <= 0.0 for r in subset)
        if n_excluded:
            ax.text(0.02, 0.98, f"no strict SC root: {n_excluded}/{len(subset)} sampled points",
                    transform=ax.transAxes, va="top", fontsize=7, color="tab:red")
    leg = axes[1, 2]
    leg.axis("off")
    handles = [leg.plot([], [], "-", color=cmap(0.05 + 0.9 * j / (len(fracs) - 1)),
                        label=f"$T/T_c$={f:g}")[0] for j, f in enumerate(fracs)]
    handles.append(leg.plot([], [], "s", ms=4, mfc="none", mec="tab:red", ls="none",
                            label=r"$D_0(0)<0.05t$")[0])
    leg.legend(handles=handles, loc="center", title="Strict-policy curves")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(FIG_DIR, filename)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out, flush=True)


def _plot_delta(data, rows, row_kind, curve_kind, filename, title):
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(len(rows), 2, figsize=(11.6, 3.1 * len(rows)), sharex="col", squeeze=False)
    cvals = sorted({r[curve_kind] for r in data})
    for i, rv in enumerate(rows):
        for j, (xkey, xlabel) in enumerate((("frac", "$T/T_c$"), ("deffn", "$\\Delta_{eff}(T)/\\Delta_{eff}(0)$"))):
            ax = axes[i, j]
            for r in data:
                if r[row_kind] != rv or r.get("Tc", 0.0) <= 0.0:
                    continue
                ci = cvals.index(r[curve_kind])
                color = cmap(0.05 + 0.9 * ci / max(len(cvals) - 1, 1))
                xx = np.asarray(r[xkey], dtype=float)
                yy = np.asarray([np.nan if y is None else y for y in r["ratio"]], dtype=float)
                ax.plot(xx, yy, "-o", ms=2.2, lw=1.2, color=color,
                        label=f"$U={r[curve_kind]:g}$" if curve_kind == "U" else f"$V_{{SO}}={r[curve_kind]:g}$")
            ax.set_xlabel(xlabel); ax.set_ylabel("$W/W_N$" if j == 0 else "")
            if xkey == "deffn": ax.set_xlim(0.0, 1.02)
            ax.grid(alpha=0.3)
            if i == 0 and j == 0: ax.legend(fontsize=7.2, ncol=2, loc="upper left")
        axes[i, 0].set_title(("$V_{SO}$" if row_kind == "VSO" else "$U$") + f" = {rv:g}", loc="left", fontsize=10)
        rows_here = [r for r in data if r[row_kind] == rv]
        no_sc = sum(float(r.get("Tc") or 0.0) <= 0.0 for r in rows_here)
        if no_sc:
            axes[i, 1].text(0.98, 0.03, f"no strict SC root at {no_sc}/{len(rows_here)} curve sets",
                            transform=axes[i, 1].transAxes, ha="right", va="bottom",
                            fontsize=7, color="tab:red")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = os.path.join(FIG_DIR, filename)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print("wrote", out, flush=True)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(FIG_DIR, exist_ok=True)
    relaxed_uv = json.load(open(os.path.join(HERE, "results", "uvso_families.json"), encoding="utf-8"))
    relaxed_delta = json.load(open(os.path.join(HERE, "results", "delta_families.json"), encoding="utf-8"))
    tasks = sorted({(float(r["U"]), float(r["VSO"])) for r in relaxed_uv + relaxed_delta})
    print(f"strictly re-solving {len(tasks)} unique family parameter points", flush=True)
    t0 = time.time()
    with Pool(2, initializer=_strict_worker_init) as pool:
        strict_roots = list(pool.imap_unordered(_solve_strict_stage1, tasks, chunksize=8))
    if any(r.get("error") for r in strict_roots):
        raise RuntimeError(f"strict stage-1 scan failed at {sum(bool(r.get('error')) for r in strict_roots)} points")
    roots = {(r["U"], r["VSO"]): r for r in strict_roots}
    old_uv = {(float(r["U"]), float(r["VSO"])): r for r in relaxed_uv}
    old_delta = {(float(r["U"]), float(r["VSO"])): r for r in relaxed_delta}
    strict_uv, strict_delta = [], []
    for key in sorted(old_uv):
        strict_uv.append(_strict_uv_row(roots[key], old_uv[key]))
    for key in sorted(old_delta):
        strict_delta.append(_strict_delta_row(roots[key], old_delta[key]))

    provenance = dict(
        stage="strict-policy family figure regeneration",
        policy="strict",
        band_rule="analytic Rashba continuum band interval with 1e-8 t tolerance at T=0 and in Tc solver",
        stage1_Nk=NK_STAGE1, stage2_Nk_mu=NK_MU, DOS_Nk=1200, NELL=NELL,
        n_unique_stage1=len(strict_roots),
        n_strict_sc=sum(r["Tc"] > 0.0 for r in strict_roots),
        n_no_inband_root=sum(bool(r.get("no_inband_root")) for r in strict_roots),
        n_uv_points=len(strict_uv), n_delta_points=len(strict_delta),
        n_reused_uv=sum(bool(r.get("reused_relaxed_ratio")) for r in strict_uv),
        n_reused_delta=sum(bool(r.get("reused_relaxed_ratio")) for r in strict_delta),
        interpretation="Relaxed archives are read-only. Every unique stage-1 point is re-solved strictly; ratios are reused only if the strict root matches the archived relaxed root. Stage-2 points with out-of-band SC mu are omitted, not set to zero.",
        caveat="Strict admissibility does not compare free energies between multiple in-band roots.",
        utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    )
    output = dict(provenance=provenance, strict_roots=strict_roots,
                  uvso_families=strict_uv, delta_families=strict_delta)
    json_path = os.path.join(OUT_DIR, "strict_family_data.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=1, allow_nan=False)

    _plot_uvso(strict_uv, list(uv.FAM_A_VSO), "VSO", "U",
               "WWN_vs_U_families_VSO.png",
               "$W/W_N$ vs $U$: strict analytic bare-band policy; $\\Delta_t=4.5t$, $n=1.885$")
    _plot_uvso(strict_uv, list(uv.FAM_B_U), "U", "VSO",
               "WWN_vs_VSO_families_U.png",
               "$W/W_N$ vs $V_{SO}$: strict analytic bare-band policy; $\\Delta_t=4.5t$, $n=1.885$")
    _plot_delta(strict_delta, list(delta.FAM1_VSO), "VSO", "U",
                "WWN_vs_TTc_and_Delta_Ufamilies.png",
                "$W/W_N$ vs $T/T_c$ and normalized gap: strict bare-band policy; curves by $U$")
    _plot_delta(strict_delta, list(delta.FAM2_U), "U", "VSO",
                "WWN_vs_TTc_and_Delta_VSOfamilies.png",
                "$W/W_N$ vs $T/T_c$ and normalized gap: strict bare-band policy; curves by $V_{SO}$")
    summary = dict(provenance=provenance,
                   stage1_points=[dict(U=r["U"], VSO=r["VSO"], Tc=r["Tc"],
                                       mu_0=r.get("mu_0"), mu_window=r.get("mu_window"),
                                       no_inband_root=r.get("no_inband_root", False),
                                       diagnostic_unrestricted=r.get("diagnostic_unrestricted"))
                                  for r in strict_roots])
    with open(os.path.join(OUT_DIR, "strict_family_provenance.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1, allow_nan=False)
    print(f"strict family regeneration finished in {time.time()-t0:.1f}s; "
          f"stage1 valid={provenance['n_strict_sc']}/{len(strict_roots)}")
    print("Wrote", json_path, flush=True)


if __name__ == "__main__":
    main()
