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
    """All (helicity h, branch sigma) quasiparticle states on the quarter zone.

    F-14 fix: the k-grid is the MIDPOINT rule (kgrid_ua.kgrid) with uniform weight 1/Nk^2
    (the old endpoint linspace double-counted boundary nodes).  Every state carries
      h     = Rashba helicity (+1/-1), which ONLY the helicity-diagonal vertex R_B uses;
      sigma = branch label of the E -> eps inversion, sigma = sign(eps - eps*), which the
              coherence factor needs (bruteforce_ua docstring: binning by helicity mixes the
              two branches inside one bin and makes the coherence factor wrong).
    """
    from kgrid_ua import kgrid
    e0, R = kgrid(Nk)
    eps_base = e0                              # e0 = -2 t (cos kx + cos ky)
    b = DS / (4 * t)
    estar = (mu + b * D0) / (1.0 + b * b)
    out = []
    for h in (+1, -1):
        eps = (eps_base - 2 * h * VSO * R).ravel()
        xi = eps - mu
        Db = D0 - b * eps
        E = np.hypot(xi, Db)
        sig = np.where(eps >= estar, +1, -1)
        for s_ in (+1, -1):
            m = sig == s_
            out.append(dict(h=h, s=s_, eps=eps[m], xi=xi[m], Db=Db[m], E=E[m],
                            w=np.full(m.sum(), 1.0 / (Nk * Nk))))
    return out


def sums(states, T, omega, mu, nbin=NBIN):
    """Per-(sigma,sigma') scattering sums, helicity-resolved for the vertex, and normal refs.

    Blocks are the (h, sigma) groups of states.  Coherence C is evaluated from the bin-mean
    (xi, Dbar) of each (h, sigma) block, so every block contains only ONE branch.
    Vertex A (isotropic): all (h,sigma) x (h',sigma') pairs.
    Vertex B (helicity-diagonal): only pairs with h = h'.
    """
    Eall = np.concatenate([st["E"] for st in states])
    Emin, Emax = Eall.min() * 0.995, Eall.max() * 1.002
    de = (Emax - Emin) / nbin
    Eb = Emin + (np.arange(nbin) + 0.5) * de
    f = 1.0 / (np.exp(Eb / T) + 1.0)
    therm = f[:, None] - f[None, :]
    blocks = {}
    for st in states:
        key = (st["h"], st["s"])
        idx = ((st["E"] - Emin) / de).astype(int)
        ok = (idx >= 0) & (idx < nbin)
        idx = idx[ok]
        wv = st["w"][ok] if np.ndim(st["w"]) else np.full(ok.sum(), st["w"])
        cnt = np.bincount(idx, weights=wv, minlength=nbin)
        den = np.maximum(cnt, 1e-300)
        xi = np.bincount(idx, weights=wv * st["xi"][ok], minlength=nbin) / den
        Db = np.bincount(idx, weights=wv * st["Db"][ok], minlength=nbin) / den
        blocks[key] = dict(cnt=cnt, xi=xi, Db=Db)
    epsall = np.concatenate([st["eps"] for st in states])
    lo, hi = epsall.min() * 0.999, epsall.max() * 1.001
    nen = 4000
    den_e = (hi - lo) / nen
    en = lo + (np.arange(nen) + 0.5) * den_e
    fn = 1.0 / (np.exp((en - mu) / T) + 1.0)
    tn = fn[:, None] - fn[None, :]
    cntN = {}
    for st in states:
        ie = ((st["eps"] - lo) / den_e).astype(int)
        ok = (ie >= 0) & (ie < nen)
        wv = st["w"][ok] if np.ndim(st["w"]) else np.full(ok.sum(), st["w"])
        key = st["h"]
        cntN[key] = cntN.get(key, np.zeros(nen)) + np.bincount(ie[ok], weights=wv, minlength=nen)
    cntN_tot = cntN[+1] + cntN[-1]
    res = {}
    ROW = 400                                   # row chunk: bounds peak memory (F-11 pattern)
    keys = list(blocks.keys())
    for w in WS:
        ss = {}                                 # keyed by (h,s),(h',s')
        acc = {(kA, kB): 0.0 for kA in keys for kB in keys}
        for i0 in range(0, nbin, ROW):
            i1 = min(i0 + ROW, nbin)
            Ebi = Eb[i0:i1]
            therm_c = therm[i0:i1]
            dlt = np.exp(-(((Eb[None, :] - Ebi[:, None] - omega) / w) ** 2)) / (w * np.sqrt(np.pi))
            base = therm_c * dlt
            EE = np.maximum(Ebi[:, None] * Eb[None, :], 1e-300)
            for kA in keys:
                A = blocks[kA]
                for kB in keys:
                    B = blocks[kB]
                    C = 0.5 * (1 + (A["xi"][i0:i1, None] * B["xi"][None, :]
                                    - A["Db"][i0:i1, None] * B["Db"][None, :]) / EE)
                    acc[(kA, kB)] += float((A["cnt"][i0:i1, None] * B["cnt"][None, :] * C * base).sum())
        ss.update(acc)
        # normal references, chunked the same way
        SA = 0.0; SB = {h: 0.0 for h in (+1, -1)}
        for i0 in range(0, nen, ROW):
            i1 = min(i0 + ROW, nen)
            dn = np.exp(-(((en[None, :] - en[i0:i1, None] - omega) / w) ** 2)) / (w * np.sqrt(np.pi))
            bn = tn[i0:i1] * dn
            SA += float((cntN_tot[i0:i1, None] * cntN_tot[None, :] * bn).sum())
            for h in (+1, -1):
                SB[h] += float((cntN[h][i0:i1, None] * cntN[h][None, :] * bn).sum())
        ss["SA"] = SA
        ss["SB"] = SB[+1] + SB[-1]
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
            keysA = [(a, b) for a in ((1, 1), (1, -1), (-1, 1), (-1, -1))
                     for b in ((1, 1), (1, -1), (-1, 1), (-1, -1))]
            A_tot = [sum(res[w][k] for k in keysA) for w in WS]
            B_dia = [sum(res[w][(a, b)] for a in ((1, 1), (1, -1), (-1, 1), (-1, -1))
                         for b in ((1, 1), (1, -1), (-1, 1), (-1, -1)) if a[0] == b[0])
                     for w in WS]
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
