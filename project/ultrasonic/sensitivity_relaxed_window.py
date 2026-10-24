"""
sensitivity_relaxed_window.py -- what the chemical-potential WINDOW DECISION (D1) changes.

The default policy (parameters_tc.MU_WINDOW_POLICY = "strict") rejects superconducting roots whose
chemical potential lies outside the band.  This script evaluates the SAME sets under the
"relaxed" policy and records what they would give.  It writes a sensitivity record only; the
main results always use the strict default.  Nothing here changes a shipped number.

Sets: every point of the companion family grid (uvso_families.json) with Tc = 0 under the strict
policy, plus group C (Dt = 4.5) and U = -75 (V_SO = 0 ... 0.8).
Run:  python3 sensitivity_relaxed_window.py     (~15 min on 2 cores)
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import solvers as S
import make_uvso_families as MU
from kgrid import build_kgrid
from gap_equations import make_residual

OUT = os.path.join(HERE, "results", "sensitivity", "relaxed_window_sensitivity.json")


def one(task):
    """Relaxed-policy evaluation of one (U, V_SO, Delta_t) set."""
    U, V, Dt = task
    S.MU_WINDOW_POLICY = "relaxed"                    # module-level switch read by _band()
    try:
        r = S.solve_set_core(U, Dt, 1.885, V, Nk=80)
    finally:
        S.MU_WINDOW_POLICY = "strict"
    rec = dict(U=float(U), VSO=float(V), Dt=float(Dt), Tc_relaxed=float(r["Tc"]),
               D0_0=float(r["D0_0"]), DS_0=float(r["DS_0"]), Deff_0=float(r["Deff_0"]),
               mu_0=float(r["mu_0"]), mu_window=r.get("mu_window"))
    if not (r["Tc"] > 0):
        return rec
    # ratios with the production evaluation (same code as the families), relaxed state
    S.MU_WINDOW_POLICY = "relaxed"
    try:
        o = MU.one_set((float(U), float(V))) if Dt == 4.5 else None
    finally:
        S.MU_WINDOW_POLICY = "strict"
    if o is not None and o.get("Tc", 0) > 0:
        rec["ratio_FRACS"] = dict(zip([str(f) for f in MU.FRACS], o["ratio"]))
    return rec


def main():
    t0 = time.time()
    fam = json.load(open(os.path.join(HERE, "results", "uvso_families.json")))
    tasks = sorted({(float(r["U"]), float(r["VSO"]), 4.5) for r in fam if r["Tc"] <= 0})
    extra = [(75.0, 0.8, 4.5), (-75.0, 0.0, 4.5), (-75.0, 0.4, 4.5), (-75.0, 0.8, 4.5)]
    tasks = sorted(set(tasks) | set(extra))
    print(f"{len(tasks)} sets to evaluate under the relaxed policy", flush=True)
    with Pool(2) as p:
        recs = p.map(one, tasks, chunksize=2)
    n_sc = sum(1 for r in recs if r["Tc_relaxed"] > 0)
    out = dict(policy="relaxed", note="sensitivity only; main results use policy=strict",
               n_sets=len(recs), n_superconducting_under_relaxed=n_sc, sets=recs)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1, allow_nan=False)
    print(f"relaxed policy: {n_sc}/{len(recs)} sets superconducting; wrote {OUT} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
