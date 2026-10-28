"""
make_delta_families.py -- W/W_N vs T/T_c AND vs the gap Delta^0(T), for two requested
family grids (Delta_t = 4.5 t, n = 1.885, hbar*omega = 0.7 t, eta=+1 signed, scattering
only -- the standing conventions):

  Family 1: V_SO = 0.0, 0.4, 0.8   x  U = 15, 35, 55, 75, 95, 105, 115
  Family 2: U    = 35, 75, 115     x  V_SO = 0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9

Layout: one figure per family, 3 rows (the outer parameter) x 2 columns
(left: x = T/T_c, right: x = Delta^0(T)/t, the self-consistent on-site gap amplitude at
that temperature -- Delta is monotone in T so the right column is the same physics
re-parameterised by the gap, showing how directly the ratio tracks Delta).

Sets without a residual-verified root are reported and skipped. Under the adopted relaxed
window, out-of-band roots are retained and flagged with open gray triangles; red-square notes
have the separate legacy criterion D0(0) < 0.05t. Out-of-band retention does not establish
thermodynamic stability or a Fermi-surface interpretation.

Outputs: results/delta_families.json,
         figures/WWN_vs_TTc_and_Delta_Ufamilies.png,
         figures/WWN_vs_TTc_and_Delta_VSOfamilies.png
Run: python make_delta_families.py    (~6-8 min on 2 cores)
"""
from __future__ import annotations

import os
import sys
import json
import time

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
from absorption import S_scattering, S_normal, qp_constants
from gap_equations import effective_gap
from kgrid import build_kgrid
from param_scans_peaks import get_dos

FR_GRID = np.concatenate([np.linspace(0.05, 0.95, 28), [0.99, 0.999]])
NELL = 20000

FAM1_VSO = (0.0, 0.4, 0.8)
FAM1_U = (15.0, 35.0, 55.0, 75.0, 95.0, 105.0, 115.0)
FAM2_U = (35.0, 75.0, 115.0)
FAM2_VSO = (0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9)


def one_set(task):
    U, VSO = task
    try:
        s = s1sol.solve_set(U, 4.5, 1.885, VSO, Nk=80)
    except Exception:
        return dict(U=U, VSO=VSO, Tc=0.0, policy=s1sol.MU_WINDOW_POLICY, mu_window=None)
    if s["Tc"] <= 1e-9:
        return dict(U=U, VSO=VSO, Tc=0.0, policy=s1sol.MU_WINDOW_POLICY, mu_window=s.get("mu_window"))
    sc = SCState(params=dict(VSO=float(s["V_SO"]), Tc=float(s["Tc"]), D0_0=float(s["D0_0"]),
                             DS_0=float(s["DS_0"]), Deff_0=float(s["Deff_0"])), Nk_mu=200)
    dos = get_dos(VSO)
    eps_grid, _ = build_kgrid(VSO, U, 4.5, 80)      # paper's Delta_eff definition
    fr, d0, em, tt, dn = [], [], [], [], []
    for frac in FR_GRID:
        st = sc.at(frac)
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                            "a_signed", nell=NELL)
        fr.append(S / SN)
        d0.append(st["D0"])
        em.append(P["Emin"])
        tt.append(st["T"])
        dn.append(effective_gap(st["D0"], st["DS"], eps_grid) / s["Deff_0"])
    return dict(U=U, VSO=VSO, Tc=s["Tc"], D0_0=s["D0_0"], DS_0=s["DS_0"],
                Deff_0=s["Deff_0"], mu_0=s["mu_0"], mu_window=s.get("mu_window"),
                policy=s.get("policy", s1sol.MU_WINDOW_POLICY), frac=list(FR_GRID), ratio=fr,
                D0=d0, Emin=em, T=tt, deffn=dn)


def run(tasks):
    from multiprocessing import Pool
    with Pool(2) as p:
        return p.map(one_set, tasks, chunksize=1)


def family_figure(data, rows, row_kind, curve_kind, path, title):
    """rows: values of row_kind ('VSO' or 'U'); curves: values of curve_kind."""
    cols = plt.get_cmap("tab10").colors
    cmap2 = plt.get_cmap("viridis")
    fig, axes = plt.subplots(len(rows), 2, figsize=(11.6, 3.1 * len(rows)), sharex="col",
                             squeeze=False)
    cvals = sorted({r[curve_kind] for r in data if r["Tc"] > 0})
    for i, rv in enumerate(rows):
        for j, (xk, xl) in enumerate((("frac", "$T/T_c$"),
                                      ("deffn", "$\\Delta_{eff}(T)/\\Delta_{eff}(0)$"))):
            ax = axes[i, j]
            for r in data:
                if r[row_kind] != rv or r["Tc"] <= 0:
                    continue
                ci = cvals.index(r[curve_kind])
                col = cmap2(0.05 + 0.9 * ci / max(len(cvals) - 1, 1))
                lbl = (f"$U$={r[curve_kind]:g}" if curve_kind == "U"
                       else f"$V_{{SO}}$={r[curve_kind]:g}")
                ax.plot(r[xk], r["ratio"], "-o", ms=2.4, lw=1.4, color=col, label=lbl)
                mu_out = (r.get("mu_window") or {}).get("inside") is False
                if mu_out:
                    ax.plot(r[xk], r["ratio"], "^", ms=3.0, mfc="none", mec="0.25",
                            mew=0.65, ls="none",
                            label="mu outside band (relaxed)" if i == 0 and j == 0 else None)
                if r.get("D0_0", 1.0) < 0.05:
                    ax.plot(r[xk], r["ratio"], "s", ms=2.8, mfc="none", mec="tab:red",
                            mew=0.7, ls="none",
                            label="D0(0)<0.05t" if i == 0 and j == 0 else None)
            ax.set_xlabel(xl)
            ax.set_ylabel("$W/W_N$" if j == 0 else "")
            if xk == "deffn":
                ax.set_xlim(0.0, 1.02)
            ax.grid(alpha=0.3)
            if i == 0 and j == 0:
                ax.legend(fontsize=7.5, ncol=2, loc="upper left")
        kind_lbl = "$V_{SO}$" if row_kind == "VSO" else "$U$"
        axes[i, 0].set_title(f"{kind_lbl} = {rv:g}", fontsize=10, loc="left")
        # flag non-SC / dominated roots inside the row
        notes = []
        for r in data:
            if r[row_kind] == rv:
                if r["Tc"] <= 0:
                    notes.append(f"U={r['U']:g}: no SC")
                elif r.get("D0_0", 1) < 0.05:
                    notes.append(f"U={r['U']:g}: $D_0(0)<0.05t$")
                if r.get("Tc", 0) > 0 and (r.get("mu_window") or {}).get("inside") is False:
                    notes.append(f"U={r['U']:g}: mu outside band")
        if notes and row_kind == "VSO":
            axes[i, 1].text(0.98, 0.02, "; ".join(notes), transform=axes[i, 1].transAxes,
                            ha="right", va="bottom", fontsize=7, color="tab:red")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path)


if __name__ == "__main__":
    os.chdir(HERE)
    t0 = time.time()
    tasks = sorted({(U, v) for v in FAM1_VSO for U in FAM1_U} |
                   {(U, v) for v in FAM2_VSO for U in FAM2_U})
    print(f"{len(tasks)} parameter sets", flush=True)
    data = [r for r in run(tasks) if r is not None]
    print(f"done in {time.time()-t0:.0f}s", flush=True)
    with open("results/delta_families.json", "w") as fh:
        json.dump(data, fh, indent=1)
    for r in data:
        if r["Tc"] <= 0 or r.get("D0_0", 1) < 0.05:
            print(f"  note: U={r['U']:g} V_SO={r['VSO']:g} -> Tc={r['Tc']:.4f} "
                  f"D0_0={r.get('D0_0', 0):.4f} DS_0={r.get('DS_0', 0):.4f}", flush=True)
    family_figure([r for r in data if r["VSO"] in FAM1_VSO and r["U"] in FAM1_U],
                  list(FAM1_VSO), "VSO", "U",
                  "figures/WWN_vs_TTc_and_Delta_Ufamilies.png",
                  "$W/W_N$ vs $T/T_c$ and normalized gap; relaxed $\\mu$ window (out-of-band roots marked); "
                  "$\\Delta_t=4.5t$, $n=1.885$, $\\hbar\\omega=0.7t$; rows: $V_{SO}$, curves: $U$")
    family_figure([r for r in data if r["U"] in FAM2_U and r["VSO"] in FAM2_VSO],
                  list(FAM2_U), "U", "VSO",
                  "figures/WWN_vs_TTc_and_Delta_VSOfamilies.png",
                  "$W/W_N$ vs $T/T_c$ and normalized gap; relaxed $\\mu$ window (out-of-band roots marked); "
                  "$\\Delta_t=4.5t$, $n=1.885$, $\\hbar\\omega=0.7t$; rows: $U$, curves: $V_{SO}$")
