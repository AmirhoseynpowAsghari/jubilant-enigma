"""
calibrate_smooth.py -- calibration of the DOS Gaussian smoothing SMOOTH (audit F-03).

The normal DOS is a histogram of the quarter-zone midpoint grid, Gaussian-smoothed with width
sigma (SMOOTH).  Two error sources compete:
  (i)  smoothing BIAS: the smoothed DOS differs from the sharp one by O(sigma^2 N'') -> the
       extrapolation R(sigma) = R0 + c*sigma^2 (sigma -> 0) estimates the unsmoothed value;
  (ii) histogram NOISE: at small sigma the k-grid histogram is noisy; detected by comparing two
       k-grid resolutions Nk_DOS at the same sigma.

DECISION RULE (fixed before the calibration was run, to avoid choosing a value by its result):
    sigma* = argmin over the candidate set of  E(sigma) = |R(sigma) - R0| / R0
                                              + |R_{Nk2}(sigma) - R_{Nk1}(sigma)| / R_{Nk1}(sigma)
evaluated for the master state (U=95, Dt=4.5, n=1.885, V_SO=0.4082, T/Tc = 0.5 and 0.9),
with R0 the sigma->0 extrapolation (least-squares fit of R0 + c*sigma^2 to the three smallest
sigma).  The adopted SMOOTH is the candidate minimising E, and the reported systematic is
|R(sigma*) - R0|/R0 (per T/Tc, maximised).

Output: results/calibrate_smooth.json
Run:    python calibrate_smooth.py     (~5 min)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import solvers as s1sol
from physics_ua import SCState
from dos_ua import NormalDOS
from absorption import S_scattering, S_normal
from parameters_ua import OMEGA, NBINS, SMOOTH

SIGMAS = (0.0025, 0.005, 0.01, 0.02, 0.04)
BINS_SWEEP = (6000, 12000, 24000, 48000)   # histogram resolution test (F-03 follow-up)
NKS = (1200, 2000)
FRACS = (0.5, 0.9)
NELL = 12000
VSO = 0.4082


def ratio(dos, st):
    SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
    S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], "a_signed", nell=NELL)
    return S / SN


def main():
    t0 = time.time()
    r = s1sol.solve_set_core(95.0, 4.5, 1.885, VSO, Nk=80)
    sc = SCState(params=dict(VSO=VSO, Tc=r["Tc"], D0_0=r["D0_0"], DS_0=r["DS_0"],
                             Deff_0=r["Deff_0"]), n_tgt=1.885, Nk_mu=200)
    states = {f: sc.at(f) for f in FRACS}
    table = {}
    for nk in NKS:
        for sg in SIGMAS:
            dos = NormalDOS(VSO, Nk=nk, nbins=NBINS, smooth=sg)
            for f in FRACS:
                table[(nk, sg, f)] = ratio(dos, states[f])
            print(f"Nk={nk} sigma={sg}: " + "  ".join(f"R({f})={table[(nk, sg, f)]:.6f}" for f in FRACS),
                  f"({time.time()-t0:.0f}s)", flush=True)
    out = dict(VSO=VSO, Tc=r["Tc"], D0_0=r["D0_0"], DS_0=r["DS_0"], fracs=list(FRACS),
               sigmas=list(SIGMAS), Nk=list(NKS), NELL=NELL, rule=(
                   "argmin_sigma [ |R-R0|/R0 + |R_Nk2-R_Nk1|/R_Nk1 ], R0 = least-squares "
                   "fit of R0 + c*sigma^2 on the three smallest sigma"),
               table={f"{nk}|{sg}|{f}": v for (nk, sg, f), v in table.items()})
    # ---- histogram-resolution test at sigma = 0.005 (discretisation must be converged) ----
    binsw = {}
    for nb in BINS_SWEEP:
        dos = NormalDOS(VSO, Nk=NKS[0], nbins=nb, smooth=SMOOTH)
        binsw[str(nb)] = [ratio(dos, states[f]) for f in FRACS]
        print(f"bins={nb} sigma={SMOOTH}: " + "  ".join(f"R({f})={v:.7f}" for f, v in zip(FRACS, binsw[str(nb)])), flush=True)
    out["bins_sweep_sigma0p005"] = binsw
    d_bin = [abs(binsw[str(BINS_SWEEP[i+1])][j] - binsw[str(BINS_SWEEP[i])][j]) / binsw[str(BINS_SWEEP[i])][j]
             for i in range(len(BINS_SWEEP)-1) for j in range(len(FRACS))]
    out["bins_change_last_doubling_max"] = float(max(d_bin[-len(FRACS):]))
    E = {}
    R0 = {}
    for f in FRACS:
        small = SIGMAS[:3]
        A = np.vstack([np.ones(3), np.array(small) ** 2]).T
        y = np.array([table[(NKS[0], s, f)] for s in small])
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        R0[f] = float(coef[0])
    out["R0_extrapolated"] = {str(f): R0[f] for f in FRACS}
    for sg in SIGMAS:
        e_bias = max(abs(table[(NKS[0], sg, f)] - R0[f]) / R0[f] for f in FRACS)
        e_noise = max(abs(table[(NKS[1], sg, f)] - table[(NKS[0], sg, f)]) / table[(NKS[0], sg, f)]
                      for f in FRACS)
        E[sg] = e_bias + e_noise
        out.setdefault("per_sigma", {})[str(sg)] = dict(bias=e_bias, noise=e_noise, E=E[sg])
    adopted = min(E, key=E.get)
    sys_at = max(abs(table[(NKS[0], adopted, f)] - R0[f]) / R0[f] for f in FRACS)
    out["SMOOTH_adopted"] = float(adopted)
    out["SMOOTH_in_code"] = float(SMOOTH)
    out["NBINS_in_code"] = int(NBINS)
    # sigma-convergence status: relative change per halving of sigma at the two smallest sigma
    h = []
    for f in FRACS:
        a = table[(NKS[0], 0.0025, f)]; b = table[(NKS[0], 0.005, f)]; c = table[(NKS[0], 0.01, f)]
        h.append(max(abs(b - a) / a, abs(c - b) / b))
    out["sigma_halving_change_max"] = float(max(h))
    out["sigma_limit_converged"] = bool(max(h) < 0.005)   # threshold stated in CORRECTION_LOG.md
    out["systematic_at_adopted"] = float(sys_at)
    out["systematic_at_0p02_old"] = float(max(abs(table[(NKS[0], 0.02, f)] - R0[f]) / R0[f]
                                              for f in FRACS))
    with open(os.path.join(HERE, "results", "calibrate_smooth.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("R0 (sigma->0):", {f: round(R0[f], 6) for f in FRACS})
    print("E(sigma):", {s: round(v, 6) for s, v in E.items()})
    print(f"SMOOTH in code={SMOOTH} NBINS={NBINS}; sigma-limit converged={out['sigma_limit_converged']} "
          f"(max change per sigma halving {out['sigma_halving_change_max']:.2%})")
    print(f"ADOPTED (pre-declared rule) SMOOTH={adopted}  systematic={sys_at:.2e}  "
          f"(old SMOOTH=0.02 systematic={out['systematic_at_0p02_old']:.2e})")


if __name__ == "__main__":
    main()
