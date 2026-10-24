"""
scanning_tc.py -- production stage-1 runs: group sets, master V_SO scan, U scans.
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
Stage-1 logic of the repository (branch master, commit cda0d56), parameterised over
(U, Delta_t, n, V_SO).  Since the 2026-10-24 audit campaign every solve runs on the
MIDPOINT k-grid (F-02) with the Nk -> 2Nk convergence check enabled (recorded per set in
the JSONs), and sets whose Tc cannot be bracketed are flagged with Tc=null instead of a
fabricated value (F-08).
"""
from __future__ import annotations

import json
import time
import numpy as np
from scipy.optimize import fsolve, brentq

from parameters_tc import (SETS, MASTER, SCAN_U, VSO_GRID, VSO_EXACT, Nk_DEF,
                           GRID_CONV_TOL)
from solvers import solve_set, sweep_T
import io_utils_tc as io


def _fmt_tc(r):
    return "FLAGGED (not bracketed)" if r["Tc"] is None else f"{r['Tc']:.4f}"


def run_groups():
    out = []
    for s in SETS:
        r = solve_set(U=s["U"], Delta_t=s["Delta_t"], n_fill=s["n"], V_SO=s["V_SO"],
                      check_conv=True)
        r.update(set=s)
        print(f"  {s['group']} {s['label']:8s} Tc={_fmt_tc(r)} Deff0={r['Deff_0']:.4f} "
              f"conv_ok={r['conv']['ok']} dTc={r['conv']['dTc']}", flush=True)
        out.append(r)
    io.save("groups_stage1.json", out)
    io.save("PROVENANCE.json", io.provenance(dict(
        stage="1", grid="midpoint", Nk_DEF=Nk_DEF, check_conv=True,
        GRID_CONV_TOL=GRID_CONV_TOL,
        conv_all_ok=bool(all(r["conv"]["ok"] for r in out)))))
    return out


def run_master_scan(npts=50):
    # F-09 support: the exact V_SO values requested downstream (VSO_EXACT, incl. the
    # master 0.4082) are merged into the scan grid; sc_params_at no longer snaps.
    V = sorted(set(np.linspace(0.0, 0.8, npts).tolist()) | set(VSO_EXACT))
    rec = dict(V_SO=V, Tc=[], D0_T0=[], DS_T0=[], Deff_T0=[], conv_ok=[], dTc=[])
    for v in V:
        r = solve_set(U=MASTER["U"], Delta_t=MASTER["Delta_t"], n_fill=MASTER["n"],
                      V_SO=v, check_conv=True)
        rec["Tc"].append(r["Tc"]); rec["D0_T0"].append(r["D0_0"])
        rec["DS_T0"].append(r["DS_0"]); rec["Deff_T0"].append(r["Deff_0"])
        rec["conv_ok"].append(bool(r["conv"]["ok"]))
        rec["dTc"].append(r["conv"]["dTc"])
        print(f"  VSO={v:.6f} Tc={_fmt_tc(r)} conv_ok={r['conv']['ok']}", flush=True)
    io.save("tc_vs_vso.json", dict(
        meta=dict(MASTER, grid="midpoint", Nk=Nk_DEF, check_conv=True,
                  npts=len(V), vso_exact=list(VSO_EXACT)),
        results=rec))
    return rec


def run_u_scans():
    out = {}
    for U in SCAN_U:
        rec = dict(V_SO=VSO_GRID, Tc=[], Deff_0=[], D0_0=[], DS_0=[],
                   Tc_bracketed=[], conv_ok=[])
        for v in VSO_GRID:
            r = solve_set(U=U, Delta_t=MASTER["Delta_t"], n_fill=MASTER["n"], V_SO=v,
                          check_conv=True)
            rec["Tc"].append(r["Tc"]); rec["Deff_0"].append(r["Deff_0"])
            rec["D0_0"].append(r["D0_0"]); rec["DS_0"].append(r["DS_0"])
            rec["Tc_bracketed"].append(bool(r["Tc_bracketed"]))
            rec["conv_ok"].append(bool(r["conv"]["ok"]))
            print(f"  U={U:6.0f} VSO={v:.3f} Tc={_fmt_tc(r)}", flush=True)
        out[str(int(U))] = rec
    io.save("u_scans_stage1.json", out)
    return out
