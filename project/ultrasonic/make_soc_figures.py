"""
make_soc_figures.py -- the spin-orbit (V_SO) effect IN THE SCATTERING CHANNEL
(eta=+1 signed coherence factor, scattering only, hbar*omega = 0.7 t), visualised so that
it cannot be missed, with a mechanism decomposition.

Why V_SO must show up here (all four routes are internal to the scattering channel):
  (1) band splitting eps_ks = eps_k - 2 s V_SO |sin k|  ->  the DOS N(eps) = SUM_s N_s(eps)
      splits its van-Hove structure and shifts weight around the Fermi level;
  (2) the self-consistent coupling ratio D0_0/Tc FALLS strongly with V_SO (0.906 -> 0.500
      for V_SO = 0 -> 0.8 at U=95), so at fixed T/Tc the suppression factor f(Delta/T)
      weakens and W/W_N rises;  Tc itself roughly doubles (0.286 -> 0.649), which is the
      dominant route at FIXED physical T;
  (3) the branch structure eps_sigma(E) = eps* + sigma*sqrt(E^2-Emin^2)/kappa moves
      (eps*, kappa, E_min all depend on the self-consistent DS, mu), redistributing the
      (sigma, sigma') branch weights of the scattering sum;
  (4) the gap anisotropy Dbar(eps) sampled on the two helicity sheets changes with the
      splitting, i.e. the coherence factor is evaluated at different (xi, Dbar) pairs.

Figures:
  figures/soc_effect_scattering.png  (2x2)
     (a) W/W_N vs T/Tc for V_SO = 0, 0.2, 0.4082, 0.6, 0.8   (U=95, Delta_t=4.5)
     (b) W/W_N vs V_SO at fixed T/Tc = 0.5, 0.9, 0.999 (U=95 and U=125)
     (c) relative SOC enhancement (R(V_SO)-R(0))/R(0) vs T/Tc
     (d) W/W_N vs V_SO at FIXED T = 0.1t and 0.2t (scattering only)
  figures/soc_mechanism.png  (2x2)
     (a) DOS N(eps) for three V_SO (van-Hove splitting)
     (b) Tc and D0_0/Tc vs V_SO  (routes 2)
     (c) (sigma,sigma') branch-weight fractions of S_sc at T/Tc=0.5 vs V_SO (route 3)
     (d) E_min(T/Tc=0.5) and mu vs V_SO (routes 3-4)
Run: python make_soc_figures.py   (~3 min)
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
from dos_ua import NormalDOS
from physics_ua import SCState
from parameters_ua import OMEGA
from absorption import S_scattering, S_normal, qp_constants
from param_scans_peaks import get_dos

VSO_LIST = (0.0, 0.2, 0.4082, 0.6, 0.8)
FR_GRID = np.concatenate([np.linspace(0.05, 0.95, 28), [0.99, 0.999]])
TS = (0.10, 0.20)
NELL = 20000


def one_set(task):
    U, VSO, Dt = task
    s = s1sol.solve_set(U, Dt, 1.885, VSO, Nk=80)
    if s["Tc"] <= 0:
        return None
    sc = SCState(params=dict(VSO=float(s["V_SO"]), Tc=float(s["Tc"]), D0_0=float(s["D0_0"]),
                             DS_0=float(s["DS_0"]), Deff_0=float(s["Deff_0"])), Nk_mu=200)
    dos = get_dos(VSO)
    fr_r, fr_w, fr_em = [], [], []
    for frac in FR_GRID:
        st = sc.at(frac)
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S, parts = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                                "a_signed", nell=NELL)
        tot = sum(parts.values())
        fr_r.append(S / SN)
        fr_w.append({f"{k[0]:+d}{k[1]:+d}": v / tot for k, v in parts.items()} if tot > 0 else {})
        fr_em.append(P["Emin"])
    tr = {}
    for T in TS:
        st = sc.at(T / s["Tc"])
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                            "a_signed", nell=NELL)
        tr[f"{T:.2f}"] = S / SN
    st5 = sc.at(0.5)
    return dict(U=U, VSO=VSO, Dt=Dt, Tc=s["Tc"], D0_0=s["D0_0"], DS_0=s["DS_0"],
                mu_0=s["mu_0"], ratio_vs_frac=fr_r, weights_vs_frac=fr_w,
                emin_vs_frac=fr_em, ratio_vs_T=tr,
                mu5=float(st5["mu"]), emin5=float(qp_constants(st5["D0"], st5["DS"],
                                                               st5["mu"])["Emin"]))


if __name__ == "__main__":
    os.chdir(HERE)
    t0 = time.time()
    tasks = [(U, v, 4.5) for U in (95.0, 125.0) for v in VSO_LIST]
    from multiprocessing import Pool
    with Pool(2) as p:
        res = [r for r in p.map(one_set, tasks, chunksize=1) if r]
    print(f"solved {len(res)} sets in {time.time()-t0:.0f}s", flush=True)
    with open("results/soc_scan.json", "w") as fh:
        json.dump(res, fh, indent=1)
    R = {(r["U"], round(r["VSO"], 4)): r for r in res}

    # ---------------------------------------------------------------- figure 1
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.2))
    cols = plt.get_cmap("tab10").colors
    ax = axes[0, 0]
    for i, v in enumerate(VSO_LIST):
        r = R[(95.0, round(v, 4))]
        ax.plot(FR_GRID, r["ratio_vs_frac"], "-o", ms=2.5, lw=1.5, color=cols[i],
                label=f"$V_{{SO}}$={v:g}")
    ax.set_xlabel("$T/T_c$"); ax.set_ylabel("$W/W_N$ (scattering, $\\eta=+1$ signed)")
    ax.set_title("(a) $U=95t,\\ \\Delta_t=4.5t,\\ \\hbar\\omega=0.7t$: SOC lifts every curve",
                 fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[0, 1]
    for i, (U, mk) in enumerate(((95.0, "o"), (125.0, "s"))):
        for j, frac in enumerate((0.5, 0.9, 0.999)):
            ii = int(np.argmin(abs(FR_GRID - frac)))
            yy = [R[(U, round(v, 4))]["ratio_vs_frac"][ii] for v in VSO_LIST]
            ax.plot(VSO_LIST, yy, "-" + mk, ms=4, lw=1.4, color=cols[j + 1],
                    label=f"$U$={U:g}, $T/T_c$={frac}")
    ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("$W/W_N$")
    ax.set_title("(b) fixed $T/T_c$: monotone SOC enhancement, up to $\\times$2-5", fontsize=9)
    ax.legend(fontsize=7.5); ax.grid(alpha=0.3)
    ax = axes[1, 0]
    for i, v in enumerate(VSO_LIST[1:]):
        r0 = R[(95.0, 0.0)]["ratio_vs_frac"]
        rv = R[(95.0, round(v, 4))]["ratio_vs_frac"]
        ax.plot(FR_GRID, (np.array(rv) - np.array(r0)) / np.maximum(np.array(r0), 1e-12),
                "-o", ms=2.5, lw=1.5, color=cols[i + 1], label=f"$V_{{SO}}$={v:g}")
    ax.set_xlabel("$T/T_c$"); ax.set_ylabel("$(R(V_{SO})-R(0))/R(0)$")
    ax.set_title("(c) relative SOC enhancement of the scattering ratio", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[1, 1]
    for j, T in enumerate(TS):
        for i, U in enumerate((95.0, 125.0)):
            yy = [R[(U, round(v, 4))]["ratio_vs_T"][f"{T:.2f}"] for v in VSO_LIST]
            ax.plot(VSO_LIST, yy, "-" + ("o" if U == 95 else "s"), ms=4, lw=1.4,
                    color=cols[j], ls=("-" if U == 95 else "--"),
                    label=f"$U$={U:g}, $T$={T:g}t")
    ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("$W/W_N$")
    ax.set_title("(d) fixed physical $T$: SOC effect survives (via $T_c$, dome)", fontsize=9)
    ax.legend(fontsize=7.5); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("figures/soc_effect_scattering.png", dpi=150)
    plt.close(fig)
    print("wrote figures/soc_effect_scattering.png")

    # ---------------------------------------------------------------- figure 2
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.0))
    ax = axes[0, 0]
    eps = np.linspace(-4.0, 4.0, 800)
    for i, v in enumerate((0.0, 0.4082, 0.8)):
        dos = NormalDOS(v, Nk=1400, nbins=6000, smooth=0.02)
        ax.plot(eps, dos.N(eps), lw=1.4, color=cols[i], label=f"$V_{{SO}}$={v:g}")
    ax.set_xlabel("$\\varepsilon/t$"); ax.set_ylabel("$N(\\varepsilon)$")
    ax.set_title("(a) route 1: van-Hove splitting of the DOS", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[0, 1]
    vv = [v for v in VSO_LIST]
    tc = [R[(95.0, round(v, 4))]["Tc"] for v in vv]
    ratio0 = [R[(95.0, round(v, 4))]["D0_0"] / R[(95.0, round(v, 4))]["Tc"] for v in vv]
    ax.plot(vv, tc, "-o", color=cols[0], label="$T_c/t$")
    ax.plot(vv, ratio0, "-s", color=cols[1], label="$\\Delta^0_0/T_c$")
    ax.set_xlabel("$V_{SO}/t$")
    ax.set_title("(b) route 2: $T_c$ up, coupling ratio down $\\Rightarrow$ less suppression",
                 fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[1, 0]
    keys = ["-1-1", "-1+1", "+1-1", "+1+1"]
    labs = ["$(-,-)$", "$(-,+)$", "$(+,-)$", "$(+,+)$"]
    i5 = int(np.argmin(abs(FR_GRID - 0.5)))
    bot = np.zeros(len(vv))
    for ki, k in enumerate(keys):
        w = np.array([R[(95.0, round(v, 4))]["weights_vs_frac"][i5].get(k, 0.0) for v in vv])
        ax.bar([f"{v:g}" for v in vv], w, bottom=bot, color=cols[ki], label=labs[ki], width=0.6)
        bot = bot + w
    ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("fraction of $S_{sc}$")
    ax.set_title("(c) route 3: branch-weight redistribution at $T/T_c=0.5$", fontsize=9)
    ax.legend(fontsize=7.5, ncol=4); ax.grid(alpha=0.3)
    ax = axes[1, 1]
    em = [R[(95.0, round(v, 4))]["emin5"] for v in vv]
    mu = [R[(95.0, round(v, 4))]["mu5"] for v in vv]
    ax.plot(vv, em, "-o", color=cols[0], label="$E_{min}(0.5T_c)/t$")
    ax.plot(vv, mu, "-s", color=cols[1], label="$\\mu(0.5T_c)/t$")
    ax.set_xlabel("$V_{SO}/t$")
    ax.set_title("(d) routes 3-4: gap minimum and chemical potential move", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("figures/soc_mechanism.png", dpi=150)
    plt.close(fig)
    print("wrote figures/soc_mechanism.png")

    # ---------------------------------------------------------------- numbers
    print("\nV_SO effect in the scattering channel (U=95, Delta_t=4.5, hbar*omega=0.7t):")
    print(f"{'V_SO':>6} {'Tc/t':>7} {'D0_0/Tc':>8} {'R(0.5)':>8} {'R(0.9)':>8} {'R(0.999)':>9} "
          f"{'R(0.1t)':>8} {'R(0.2t)':>8}")
    for v in VSO_LIST:
        r = R[(95.0, round(v, 4))]
        i5 = int(np.argmin(abs(FR_GRID - 0.5))); i9 = int(np.argmin(abs(FR_GRID - 0.9)))
        print(f"{v:6g} {r['Tc']:7.4f} {r['D0_0']/r['Tc']:8.3f} {r['ratio_vs_frac'][i5]:8.5f} "
              f"{r['ratio_vs_frac'][i9]:8.5f} {r['ratio_vs_frac'][-1]:9.5f} "
              f"{r['ratio_vs_T']['0.10']:8.5f} {r['ratio_vs_T']['0.20']:8.5f}")
