"""
scanning_tc.py -- production stage-1 runs: group sets, master V_SO scan, U scans.
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
Verbatim logic of the repository's stage-1 solver (branch master, commit cda0d56),
parameterised over (U, Delta_t, n, V_SO) so that every parameter set of arXiv:2607.02812
can be solved with the code that produced the paper's own stage-1 numbers.
"""
from __future__ import annotations

import json
import time
import numpy as np
from scipy.optimize import fsolve, brentq

from parameters_tc import SETS, MASTER, SCAN_U, VSO_GRID, Nk_DEF
from solvers import solve_set, sweep_T
import io_utils_tc as io


def run_groups():
    out = []
    for s in SETS:
        r = solve_set(U=s["U"], Delta_t=s["Delta_t"], n_fill=s["n"], V_SO=s["V_SO"])
        r.update(set=s)
        print(f"  {s['group']} {s['label']:8s} Tc={r['Tc']:.4f} Deff0={r['Deff_0']:.4f}", flush=True)
        out.append(r)
    io.save("groups_stage1.json", out)
    return out


def run_master_scan(npts=50):
    V = list(np.linspace(0.0, 0.8, npts))
    rec = dict(V_SO=V, Tc=[], D0_T0=[], DS_T0=[], Deff_T0=[])
    for v in V:
        r = solve_set(U=MASTER["U"], Delta_t=MASTER["Delta_t"], n_fill=MASTER["n"], V_SO=v)
        rec["Tc"].append(r["Tc"]); rec["D0_T0"].append(r["D0_0"])
        rec["DS_T0"].append(r["DS_0"]); rec["Deff_T0"].append(r["Deff_0"])
    io.save("tc_vs_vso.json", dict(meta=MASTER, results=rec))
    return rec


def run_u_scans():
    out = {}
    for U in SCAN_U:
        rec = dict(V_SO=VSO_GRID, Tc=[], Deff_0=[], D0_0=[], DS_0=[])
        for v in VSO_GRID:
            r = solve_set(U=U, Delta_t=MASTER["Delta_t"], n_fill=MASTER["n"], V_SO=v)
            rec["Tc"].append(r["Tc"]); rec["Deff_0"].append(r["Deff_0"])
            rec["D0_0"].append(r["D0_0"]); rec["DS_0"].append(r["DS_0"])
            print(f"  U={U:6.0f} VSO={v:.3f} Tc={r['Tc']:.4f}", flush=True)
        out[str(int(U))] = rec
    io.save("u_scans_stage1.json", out)
    return out
