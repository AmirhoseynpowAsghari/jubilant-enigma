"""
make_uvso_families.py -- requested companion scans at Delta_t = 4.5 t, n = 1.885,
hbar*omega = 0.7 t, eta=+1 signed, scattering only:

  Fig A: W/W_N vs U      (U = 15...150, step 0.5 -> 271 points per V_SO) for V_SO = 0.0,
         0.2, 0.4, 0.6, 0.8; curves = T/T_c = 0.3, 0.5, 0.7, 0.9, 0.999
  Fig B: W/W_N vs V_SO   (V_SO = 0...0.9, step 0.025 -> 37 points) for U = 30, 45, 75,
         95, 115; curves = the same T/T_c values

The relaxed chemical-potential policy includes residual-verified roots even when mu is outside
 the bare band. Such points are explicitly flagged by open gray triangles; they are algebraic
 mean-field solutions, not automatically stable Fermi-surface states. Small red squares retain
 their separate legacy meaning: D0(0) < 0.05t.

Sets without a residual-verified superconducting root (Tc = 0) are omitted and noted.

Outputs: results/uvso_families.json, results/uvso_families_PROVENANCE.json,
         figures/WWN_vs_U_families_VSO.png, figures/WWN_vs_VSO_families_U.png
Run: python make_uvso_families.py
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
from parameters_ua import OMEGA, NBINS, SMOOTH
from absorption import S_scattering, S_normal
from dos_ua import NormalDOS

# Bounded DOS cache (LRU, 3 entries): the dense V_SO grid has 42 distinct values and the
# machine has 2 GB RAM -- an unbounded cache OOM-kills the pool workers.  Construction
# parameters are IDENTICAL to param_scans_peaks.get_dos, so all numbers are unchanged.
_DOS, _DOS_ORDER = {}, []


def get_dos(VSO):
    key = round(float(VSO), 4)
    if key not in _DOS:
        if len(_DOS_ORDER) >= 3:
            del _DOS[_DOS_ORDER.pop(0)]
        _DOS[key] = NormalDOS(key, Nk=1200)
        _DOS_ORDER.append(key)
    return _DOS[key]

FRACS = (0.3, 0.5, 0.7, 0.9, 0.999)
NELL = 20000

UGRID = [float(round(x, 3)) for x in np.arange(15.0, 150.01, 0.5)]   # 271 points (relaxed window)
FAM_A_VSO = (0.0, 0.2, 0.4, 0.6, 0.8)
FAM_B_U = (30.0, 45.0, 75.0, 95.0, 115.0)
VSOGRID = [float(round(x, 4)) for x in np.arange(0.0, 0.9001, 0.025)]  # 37 points


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
                Deff_0=s["Deff_0"], mu_0=s.get("mu_0"),
                mu_window=s.get("mu_window"), policy=s.get("policy", s1sol.MU_WINDOW_POLICY),
                ratio=rr)


def family_figure(data, panels, panel_kind, x_kind, path, title):
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.4), squeeze=False)
    for i, pv in enumerate(panels):
        ax = axes[i // 3, i % 3]
        for j, frac in enumerate(FRACS):
            xx, yy, flag, mu_out = [], [], [], []
            for r in data:
                if r[panel_kind] != pv or r["Tc"] <= 0:
                    continue
                xx.append(r[x_kind])
                yy.append(r["ratio"][j])
                flag.append(r.get("D0_0", 1) < 0.05)
                mw = r.get("mu_window") or {}
                mu_out.append(mw.get("inside") is False)
            if not xx:
                continue
            o = np.argsort(xx)
            xx, yy = np.array(xx)[o], np.array(yy)[o]
            flag, mu_out = np.array(flag)[o], np.array(mu_out)[o]
            col = cmap(0.05 + 0.9 * j / (len(FRACS) - 1))
            ax.plot(xx, yy, "-", lw=1.4, color=col, label=f"$T/T_c$={frac:g}")
            if flag.any():
                ax.plot(xx[flag], yy[flag], "s", ms=2.6, mfc="none", mec="tab:red",
                        mew=0.8, ls="none")
            if mu_out.any():
                ax.plot(xx[mu_out], yy[mu_out], "^", ms=2.8, mfc="none", mec="0.25",
                        mew=0.65, ls="none")
        kl = "$U/t$" if x_kind == "U" else "$V_{SO}/t$"
        pl = "$V_{SO}$" if panel_kind == "VSO" else "$U$"
        ax.set_xlabel(kl)
        ax.set_ylabel("$W/W_N$")
        ax.set_title(f"{pl} = {pv:g}", fontsize=10)
        ax.grid(alpha=0.3)
        def _ranges(vals):
            vals = sorted(vals)
            out, lo, prev = [], vals[0], vals[0]
            for v in vals[1:]:
                if v - prev > (3.0 if x_kind == "U" else 0.03):
                    out.append((lo, prev)); lo = v
                prev = v
            out.append((lo, prev))
            key = "U" if x_kind == "U" else "V_{SO}"
            return ", ".join((f"${key}$={a:g}" if abs(a - b) < 1e-9 else
                              f"${key}$={a:g}–{b:g}") for a, b in out)
        notes = []
        nosc = [r[x_kind] for r in data if r[panel_kind] == pv and r["Tc"] <= 0]
        flagged = [r[x_kind] for r in data
                   if r[panel_kind] == pv and r["Tc"] > 0 and r.get("D0_0", 1) < 0.05]
        if nosc:
            notes.append("no SC: " + _ranges(nosc))
        if flagged:
            notes.append(r"$\Delta^0\ll\Delta^S$ (red squares): " + _ranges(flagged))
        if notes:
            ax.text(0.02, 0.98, "; ".join(notes), transform=ax.transAxes, va="top",
                    fontsize=7, color="tab:red")
    axl = axes[1, 2]
    axl.axis("off")
    hs = [axl.plot([], [], "-", lw=1.4, color=cmap(0.05 + 0.9 * j / (len(FRACS) - 1)),
                   label=f"$T/T_c$={f:g}")[0] for j, f in enumerate(FRACS)]
    hs.append(axl.plot([], [], "s", ms=4, mfc="none", mec="tab:red", mew=0.9, ls="none",
                       label=r"$D_0(0)<0.05t$")[0])
    hs.append(axl.plot([], [], "^", ms=4, mfc="none", mec="0.25", mew=0.8, ls="none",
                       label=r"$\mu_0$ outside band (relaxed)")[0])
    axl.legend(handles=hs, fontsize=10, loc="center", title="curves", title_fontsize=10)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path)


if __name__ == "__main__":
    os.chdir(HERE)
    t0 = time.time()
    # V_SO-major task order + large chunks: each worker only ever holds the DOS objects
    # of one or two adjacent V_SO values (bounded cache above).
    tasks = sorted({(U, v) for v in FAM_A_VSO for U in UGRID} |
                   {(U, v) for v in VSOGRID for U in FAM_B_U}, key=lambda tv: (tv[1], tv[0]))
    print(f"{len(tasks)} parameter sets", flush=True)
    from multiprocessing import Pool
    data = []
    with Pool(2) as p:
        for i, r in enumerate(p.imap(one_set, tasks, chunksize=16), 1):
            if r is not None:
                data.append(r)
            if i % 25 == 0:
                print(f"  {i}/{len(tasks)} sets, {time.time()-t0:.0f}s", flush=True)
    print(f"done in {time.time()-t0:.0f}s", flush=True)
    with open("results/uvso_families.json", "w") as fh:
        json.dump(data, fh, indent=1, allow_nan=False)
    import datetime, subprocess
    try:
        git_head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE,
                                           text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        git_head = "unknown"
    provenance = dict(
        stage="U-family", policy=s1sol.MU_WINDOW_POLICY, grid="midpoint",
        U_min=min(UGRID), U_max=max(UGRID), dU=0.5, n_U=len(UGRID),
        stage1_Nk=80, stage2_Nk_mu=200, DOS_Nk=1200, DOS_nbins=NBINS,
        DOS_smooth=SMOOTH, nell=NELL, fracs=list(FRACS), n_records=len(data),
        git_head=git_head,
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"))
    with open("results/uvso_families_PROVENANCE.json", "w") as fh:
        json.dump(provenance, fh, indent=1, allow_nan=False)
    for r in data:
        if r["Tc"] <= 0 or r.get("D0_0", 1) < 0.05:
            print(f"  note: U={r['U']:g} V_SO={r['VSO']:g} -> Tc={r['Tc']:.4f} "
                  f"D0_0={r.get('D0_0', 0):.4f} DS_0={r.get('DS_0', 0):.4f}", flush=True)
    family_figure([r for r in data if r["VSO"] in FAM_A_VSO and r["U"] in UGRID],
                  list(FAM_A_VSO), "VSO", "U", "figures/WWN_vs_U_families_VSO.png",
                  "$W/W_N$ vs $U$ for each $T/T_c$; panels: $V_{SO}$; relaxed $\\mu$ window "
                  "(out-of-band roots marked); $\\Delta_t=4.5t$, $n=1.885$, $\\hbar\\omega=0.7t$")
    family_figure([r for r in data if r["U"] in FAM_B_U and r["VSO"] in VSOGRID],
                  list(FAM_B_U), "U", "VSO", "figures/WWN_vs_VSO_families_U.png",
                  "$W/W_N$ vs $V_{SO}$ for each $T/T_c$; panels: $U$; relaxed $\\mu$ window "
                  "(out-of-band roots marked); $\\Delta_t=4.5t$, $n=1.885$, $\\hbar\\omega=0.7t$")
