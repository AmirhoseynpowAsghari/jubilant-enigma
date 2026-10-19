"""
make_soc_explicit.py -- how EXPLICIT is the spin-orbit effect in our scattering kernel?

The kernel C_I(xi,xi',Dbar,Dbar'; E,E') carries no V_SO: SOC enters only through
(i) band splitting -> DOS / phase space / Fermi volumes, (ii) self-consistent scales
(Tc, D0, DS, mu), (iii) E_min and thresholds.  In addition, the inherited vertex
approximation (paper, Walker, Sardwal homogeneous limit) takes inter-helicity matrix
elements EQUAL to intra-helicity ones, while a strict q->0 density/strain vertex is
helicity-diagonal (Rashba spinor overlap -> delta_ss').

This script quantifies both statements with k-space lattice sums that track helicity s:

  (a) inter-helicity fraction of S_sc vs V_SO  (how much weight rides on channels a
      spinor-conserving vertex would suppress);
  (b) R_A (isotropic vertex, our standard) vs R_B (helicity-diagonal vertex, strict
      q->0) vs V_SO at T/Tc = 0.5, 0.9  -- the vertex systematic band of the SOC trend;
  (c) signed difference curves dR = R(V_SO) - R(0) vs T/T_c (the explicit effect size);
  (d) R_B/R_A - 1 vs V_SO.
Normal-state references use the SAME vertex rule in each variant (S_N_A: full DOS
convolution; S_N_B: helicity-diagonal convolution), so each ratio is internally
consistent; at V_SO=0 the two variants coincide exactly (sanity check).

Output: figures/soc_explicit_vertex.png  (copied to SOCfigures/35_...),
        results/soc_explicit.json
Run: python make_soc_explicit.py   (~2-4 min)
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
from parameters_ua import OMEGA, T_HOP
from absorption import qp_constants, S_normal
from param_scans_peaks import get_dos

FRACS = (0.5, 0.9)
NBIN = 2500
WS = (0.05, 0.035, 0.025)
VSO_LIST = (0.0, 0.2, 0.4082, 0.6, 0.8)


def hel_states(VSO, Nk, D0, DS, mu, t=T_HOP):
    """Particle-branch quasiparticle states on the quarter zone, helicity-resolved."""
    kx = np.linspace(0, np.pi, Nk)
    KX, KY = np.meshgrid(kx, kx)
    w = (np.pi / Nk) ** 2 / (np.pi ** 2 / 4)          # quarter-zone weight -> 1
    cs = np.cos(KX) + np.cos(KY)
    g = np.sqrt(np.sin(KX) ** 2 + np.sin(KY) ** 2)
    out = []
    for s in (+1, -1):
        eps = -2 * t * cs - 2 * s * VSO * g
        xi = eps - mu
        Db = D0 - (DS / (4 * t)) * eps
        E = np.hypot(xi, Db)
        out.append(dict(s=s, eps=eps.ravel(), xi=xi.ravel(), Db=Db.ravel(),
                        E=E.ravel(), w=w))
    return out


def sums(states, T, omega, mu, nbin=NBIN):
    """Per-(s,s') scattering sums and per-variant normal references, Gaussian delta."""
    Eall = np.concatenate([st["E"] for st in states])
    Emin, Emax = Eall.min() * 0.995, Eall.max() * 1.002
    de = (Emax - Emin) / nbin
    Eb = Emin + (np.arange(nbin) + 0.5) * de
    f = 1.0 / (np.exp(Eb / T) + 1.0)
    therm = f[:, None] - f[None, :]
    blocks = {}
    for st in states:
        idx = ((st["E"] - Emin) / de).astype(int)
        ok = (idx >= 0) & (idx < nbin)
        idx = idx[ok]
        cnt = np.bincount(idx, weights=np.full(ok.sum(), st["w"]), minlength=nbin)
        xi = np.bincount(idx, weights=st["w"] * st["xi"][ok], minlength=nbin) / np.maximum(cnt, 1e-300)
        Db = np.bincount(idx, weights=st["w"] * st["Db"][ok], minlength=nbin) / np.maximum(cnt, 1e-300)
        epsb = np.bincount(idx, weights=st["w"] * st["eps"][ok], minlength=nbin) / np.maximum(cnt, 1e-300)
        blocks[st["s"]] = dict(cnt=cnt, xi=xi, Db=Db, eps=epsb)
    epsall = np.concatenate([st["eps"] for st in states])
    lo, hi = epsall.min() * 0.999, epsall.max() * 1.001
    nen = 4000
    den = (hi - lo) / nen
    en = lo + (np.arange(nen) + 0.5) * den
    fn = 1.0 / (np.exp((en - mu) / T) + 1.0)
    tn = fn[:, None] - fn[None, :]
    cntN = {}
    for st in states:
        ie = ((st["eps"] - lo) / den).astype(int)
        ok = (ie >= 0) & (ie < nen)
        cntN[st["s"]] = np.bincount(ie[ok], weights=st["w"][ok] if np.ndim(st["w"]) else np.full(ok.sum(), st["w"]), minlength=nen)
    cntN_tot = cntN[+1] + cntN[-1]
    res = {}
    for w in WS:
        dlt = np.exp(-(((Eb[None, :] - Eb[:, None] - omega) / w) ** 2)) / (w * np.sqrt(np.pi))
        base = therm * dlt
        ss = {}
        for sa in (+1, -1):
            for sb in (+1, -1):
                A, B = blocks[sa], blocks[sb]
                EE = np.maximum(Eb[:, None] * Eb[None, :], 1e-300)
                C = 0.5 * (1 + (A["xi"][:, None] * B["xi"][None, :]
                                - A["Db"][:, None] * B["Db"][None, :]) / EE)
                ss[(sa, sb)] = float((A["cnt"][:, None] * B["cnt"][None, :] * C * base).sum())
        # normal references with the SAME width and weights: A = all pairs, B = diag in s
        dn = np.exp(-(((en[None, :] - en[:, None] - omega) / w) ** 2)) / (w * np.sqrt(np.pi))
        bn = tn * dn
        SA = float((cntN_tot[:, None] * cntN_tot[None, :] * bn).sum())
        SB = sum(float((cntN[sq][:, None] * cntN[sq][None, :] * bn).sum()) for sq in cntN)
        ss["SA"], ss["SB"] = SA, SB
        res[w] = ss
    return res


def rich(vals):
    w = np.array(WS)
    A = np.vstack([np.ones_like(w), w ** 2]).T
    c, *_ = np.linalg.lstsq(A, np.array(vals), rcond=None)
    return c[0]


if __name__ == "__main__":
    os.chdir(HERE)
    t0 = time.time()
    out = {}
    for v in VSO_LIST:
        s = s1sol.solve_set(95.0, 4.5, 1.885, v, Nk=80)
        sc = SCState(params=dict(VSO=float(s["V_SO"]), Tc=float(s["Tc"]),
                                 D0_0=float(s["D0_0"]), DS_0=float(s["DS_0"]),
                                 Deff_0=float(s["Deff_0"])), Nk_mu=200)
        rec = {}
        for frac in FRACS:
            st = sc.at(frac)
            states = hel_states(v, 240, st["D0"], st["DS"], st["mu"])
            res = sums(states, st["T"], OMEGA, st["mu"])
            A_tot = [sum(res[w][k] for k in ((1, 1), (1, -1), (-1, 1), (-1, -1))) for w in WS]
            B_dia = [res[w][(1, 1)] + res[w][(-1, -1)] for w in WS]
            SA = [res[w]["SA"] for w in WS]
            SB = [res[w]["SB"] for w in WS]
            RA = rich(A_tot) / rich(SA)
            RB = rich(B_dia) / rich(SB)
            inter = 1.0 - rich(B_dia) / rich(A_tot)
            rec[str(frac)] = dict(RA=float(RA), RB=float(RB), inter=float(inter))
            print(f"V_SO={v:6.4g} T/Tc={frac}: R_A={rec[str(frac)]['RA']:.5f} "
                  f"R_B={rec[str(frac)]['RB']:.5f} inter-helicity fraction={inter:.4f}",
                  flush=True)
        out[f"{v:.4f}"] = rec
    with open("results/soc_explicit.json", "w") as fh:
        json.dump(out, fh, indent=1)

    # existing E-space curves for the difference panel
    soc = {(r["U"], round(r["VSO"], 4)): r for r in json.load(open("results/soc_scan.json"))}
    FRG = np.concatenate([np.linspace(0.05, 0.95, 28), [0.99, 0.999]])

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.2))
    cols = plt.get_cmap("tab10").colors
    ax = axes[0, 0]
    for j, frac in enumerate(FRACS):
        ax.plot(VSO_LIST, [out[f"{v:.4f}"][str(frac)]["inter"] for v in VSO_LIST],
                "-o", color=cols[j], label=f"$T/T_c$={frac}")
    ax.set_xlabel("$V_{SO}/t$")
    ax.set_ylabel("inter-helicity fraction of $S_{sc}$")
    ax.set_title("(a) weight on channels a spinor-conserving q$\\to$0 vertex suppresses",
                 fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[0, 1]
    for j, frac in enumerate(FRACS):
        ax.plot(VSO_LIST, [out[f"{v:.4f}"][str(frac)]["RA"] for v in VSO_LIST],
                "-o", color=cols[j], label=f"$R_A$ (isotropic), $T/T_c$={frac}")
        ax.plot(VSO_LIST, [out[f"{v:.4f}"][str(frac)]["RB"] for v in VSO_LIST],
                "--s", color=cols[j], label=f"$R_B$ (helicity-diag), $T/T_c$={frac}")
    ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("$W/W_N$")
    ax.set_title("(b) vertex approximation vs strict q$\\to$0: same SOC trend", fontsize=9)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    ax = axes[1, 0]
    r0 = soc[(95.0, 0.0)]["ratio_vs_frac"]
    for i, v in enumerate((0.2, 0.4082, 0.6, 0.8)):
        rv = soc[(95.0, round(v, 4))]["ratio_vs_frac"]
        ax.plot(FRG, np.array(rv) - np.array(r0), "-o", ms=2.5, color=cols[i + 1],
                label=f"$V_{{SO}}$={v:g}")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("$T/T_c$"); ax.set_ylabel("$R(V_{SO})-R(0)$")
    ax.set_title("(c) signed SOC effect vs $T/T_c$ (explicit, fixed $T/T_c$ ensemble)",
                 fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[1, 1]
    for j, frac in enumerate(FRACS):
        ax.plot(VSO_LIST, [out[f"{v:.4f}"][str(frac)]["RB"] / out[f"{v:.4f}"][str(frac)]["RA"] - 1
                           for v in VSO_LIST], "-o", color=cols[j], label=f"$T/T_c$={frac}")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("$V_{SO}/t$"); ax.set_ylabel("$R_B/R_A - 1$")
    ax.set_title("(d) size of the vertex systematic vs the SOC trend", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.suptitle("Explicit vs indirect spin-orbit content of the scattering kernel "
                 "($U=95t$, $\\Delta_t=4.5t$, $\\hbar\\omega=0.7t$)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig("figures/soc_explicit_vertex.png", dpi=150)
    plt.close(fig)
    print("wrote figures/soc_explicit_vertex.png", f"({time.time()-t0:.0f}s)")
