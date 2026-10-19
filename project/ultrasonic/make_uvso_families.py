"""
make_uvso_families.py -- requested companion scans at Delta_t = 4.5 t, n = 1.885,
hbar*omega = 0.7 t, eta=+1 signed, scattering only:

  Fig A: W/W_N vs U      (U = 15...150) for V_SO = 0.0, 0.2, 0.4, 0.6, 0.8,
         one panel per V_SO, curves = T/T_c = 0.3, 0.5, 0.7, 0.9, 0.999
  Fig B: W/W_N vs V_SO   (V_SO = 0...0.9) for U = 30, 45, 75, 95, 115,
         one panel per U,   curves = the same T/T_c values

Sets with no SC solution (Tc = 0) are omitted from the curves and noted; sets on the
correlated-hopping-dominated root (Delta^0_0 << Delta^S_0) are plotted but flagged.

Outputs: results/uvso_families.json,
         figures/WWN_vs_U_families_VSO.png, figures/WWN_vs_VSO_families_U.png
Run: python make_uvso_families.py   (~4-6 min on 2 cores)
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
from absorption import S_scattering, S_normal
from param_scans_peaks import get_dos

FRACS = (0.3, 0.5, 0.7, 0.9, 0.999)
NELL = 20000

UGRID = [15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 75.0, 85.0, 95.0, 105.0, 115.0, 125.0,
         135.0, 150.0]
FAM_A_VSO = (0.0, 0.2, 0.4, 0.6, 0.8)
FAM_B_U = (30.0, 45.0, 75.0, 95.0, 115.0)
VSOGRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


def one_set(task):
    U, VSO = task
    try:
        s = s1sol.solve_set(U, 4.5, 1.885, VSO, Nk=80)
    except Exception:
        return dict(U=U, VSO=VSO, Tc=0.0)
    if s["Tc"] <= 1e-9:
        return dict(U=U, VSO=VSO, Tc=0.0)
    sc = SCState(params=dict(VSO=float(s["V_SO"]), Tc=float(s["Tc"]), D0_0=float(s["D0_0"]),
                             DS_0=float(s["DS_0"]), Deff_0=float(s["Deff_0"])), Nk_mu=200)
    dos = get_dos(VSO)
    rr = []
    for frac in FRACS:
        st = sc.at(frac)
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                            "a_signed", nell=NELL)
        rr.append(S / SN)
    return dict(U=U, VSO=VSO, Tc=s["Tc"], D0_0=s["D0_0"], DS_0=s["DS_0"],
                Deff_0=s["Deff_0"], ratio=rr)


def family_figure(data, panels, panel_kind, x_kind, path, title):
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), squeeze=False)
    for i, pv in enumerate(panels):
        ax = axes[i // 3, i % 3]
        for j, frac in enumerate(FRACS):
            xx, yy = [], []
            for r in data:
                if r[panel_kind] != pv or r["Tc"] <= 0:
                    continue
                xx.append(r[x_kind])
                yy.append(r["ratio"][j])
            if not xx:
                continue
            o = np.argsort(xx)
            ax.plot(np.array(xx)[o], np.array(yy)[o], "-o", ms=3, lw=1.5,
                    color=cmap(0.05 + 0.9 * j / (len(FRACS) - 1)),
                    label=f"$T/T_c$={frac:g}")
        kl = "$U/t$" if x_kind == "U" else "$V_{SO}/t$"
        pl = "$V_{SO}$" if panel_kind == "VSO" else "$U$"
        ax.set_xlabel(kl)
        ax.set_ylabel("$W/W_N$")
        ax.set_title(f"{pl} = {pv:g}", fontsize=10)
        ax.grid(alpha=0.3)
        notes = [f"U={r['U']:g}: no SC" for r in data
                 if r[panel_kind] == pv and r["Tc"] <= 0]
        notes += [f"U={r['U']:g}: $\\Delta^0\\ll\\Delta^S$" for r in data
                  if r[panel_kind] == pv and r["Tc"] > 0 and r.get("D0_0", 1) < 0.05]
        if notes:
            ax.text(0.02, 0.98, "; ".join(notes), transform=ax.transAxes, va="top",
                    fontsize=7, color="tab:red")
    axl = axes[1, 2]
    axl.axis("off")
    hs = [axl.plot([], [], "-o", ms=3, lw=1.5, color=cmap(0.05 + 0.9 * j / (len(FRACS) - 1)),
                   label=f"$T/T_c$={f:g}")[0] for j, f in enumerate(FRACS)]
    axl.legend(handles=hs, fontsize=10, loc="center", title="curves", title_fontsize=10)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path)


if __name__ == "__main__":
    os.chdir(HERE)
    t0 = time.time()
    tasks = sorted({(U, v) for v in FAM_A_VSO for U in UGRID} |
                   {(U, v) for v in VSOGRID for U in FAM_B_U})
    print(f"{len(tasks)} parameter sets", flush=True)
    from multiprocessing import Pool
    with Pool(2) as p:
        data = [r for r in p.map(one_set, tasks, chunksize=2) if r is not None]
    print(f"done in {time.time()-t0:.0f}s", flush=True)
    with open("results/uvso_families.json", "w") as fh:
        json.dump(data, fh, indent=1)
    for r in data:
        if r["Tc"] <= 0 or r.get("D0_0", 1) < 0.05:
            print(f"  note: U={r['U']:g} V_SO={r['VSO']:g} -> Tc={r['Tc']:.4f} "
                  f"D0_0={r.get('D0_0', 0):.4f} DS_0={r.get('DS_0', 0):.4f}", flush=True)
    family_figure([r for r in data if r["VSO"] in FAM_A_VSO and r["U"] in UGRID],
                  list(FAM_A_VSO), "VSO", "U", "figures/WWN_vs_U_families_VSO.png",
                  "$W/W_N$ vs $U$ for each $T/T_c$; panels: $V_{SO}$; "
                  "$\\Delta_t=4.5t$, $n=1.885$, $\\hbar\\omega=0.7t$")
    family_figure([r for r in data if r["U"] in FAM_B_U and r["VSO"] in VSOGRID],
                  list(FAM_B_U), "U", "VSO", "figures/WWN_vs_VSO_families_U.png",
                  "$W/W_N$ vs $V_{SO}$ for each $T/T_c$; panels: $U$; "
                  "$\\Delta_t=4.5t$, $n=1.885$, $\\hbar\\omega=0.7t$")
