"""
param_scans_peaks.py -- systematic (U, V_SO, Delta_t) scans of W/W_N looking for peaks,
local maxima, resonances and other non-monotonic features.

PHYSICAL MECHANISM EXPECTED TO PRODUCE PEAKS (not put in by hand):
the effective gap of the HHM model is energy-dependent,

    Dbar(eps) = D0 - (Delta_t / 4t) * eps        (b = DS/4t after self-consistency),

so it can vanish INSIDE the band.  The quasiparticle gap is

    E_min = |D0 - b*mu| / kappa ,     kappa = sqrt(1+b^2),

i.e. E_min = 0 exactly when Dbar(mu) = 0: a node of the gap on the Fermi shell.
Across that point the system goes gapped -> gapless -> gapped, and W/W_N (which is
exponentially suppressed by f(E_min) at fixed T/Tc) must rise towards its normal-state
value and fall again: a genuine peak.  The ridge D0(T) = b(T)*mu(T) runs through the
(U, V_SO, Delta_t) space, so the peak can appear in any of the three scans.
A second marked (not plotted as in-range answer) feature is the pair-creation threshold
2*E_min(T) = hbar*omega: beyond it the physical attenuation acquires the S_pair channel
(REPORT Appendix A); the scattering-only ratio plotted here is continuous across it but
the *total* physical answer changes character, so the crossing is annotated.

Protocol (professor's guidelines):
  * broad, physically meaningful ranges (registry of arXiv:2607.02812: U>0 attractive
    convention at n = 1.885; U = -75 excluded as the unphysical strong-coupling root, E7);
  * several fixed values of each companion parameter;
  * automatic local-maximum detection with prominence threshold;
  * refinement scans (13 then 25 points) around every candidate + parabola vertex;
  * artifact checks at each confirmed peak: (i) doubled quadrature resolution nell,
    (ii) INDEPENDENT k-space lattice evaluation (bruteforce_ua, T7 machinery);
  * monotonic families are reported as monotonic -- no peaks are invented.

Outputs: results/param_scans_coarse.json, results/param_scans_peaks.json,
         figures/param_scan_U.png, param_scan_VSO.png, param_scan_Dt.png,
         figures/param_scan_diagnostics.png
Run:  python param_scans_peaks.py        (~10-15 min on 2 cores)
"""
from __future__ import annotations

import os
import sys
import json
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))
import solvers as s1sol

from dos_ua import NormalDOS
from physics_ua import SCState
from parameters_ua import OMEGA, T_HOP
from absorption import S_scattering, S_normal, qp_constants
import bruteforce_ua as bf

FRACS = (0.5, 0.9)
NELL = 20000
N_TARGET = 1.885

_DOS = {}


def get_dos(VSO):
    key = round(float(VSO), 4)
    if key not in _DOS:
        _DOS[key] = NormalDOS(key, Nk=1200, nbins=6000, smooth=0.02)
    return _DOS[key]


def evaluate(task):
    """One (U, V_SO, Delta_t) point: stage-1 solve + stage-2 ratio at both FRACS."""
    U, VSO, Dt = task
    try:
        s = s1sol.solve_set(U, Dt, N_TARGET, VSO, Nk=80)
    except Exception:
        return None
    if not (s["Tc"] > 1e-6 and 0.0 < s["D0_0"] < 10.0 and 0.0 <= s["DS_0"] < 20.0):
        return None
    sc = SCState(params=dict(VSO=float(s["V_SO"]), Tc=float(s["Tc"]), D0_0=float(s["D0_0"]),
                             DS_0=float(s["DS_0"]), Deff_0=float(s["Deff_0"])), Nk_mu=200)
    dos = get_dos(VSO)
    fr = []
    for frac in FRACS:
        st = sc.at(frac)
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                            "a_signed", nell=NELL)
        fr.append(dict(frac=frac, T=float(st["T"]), ratio=float(S / SN),
                       Emin=float(P["Emin"]), mu=float(st["mu"]),
                       D0=float(st["D0"]), DS=float(st["DS"]),
                       pair_open=bool(2.0 * P["Emin"] < OMEGA)))
    return dict(U=float(U), VSO=float(VSO), Dt=float(Dt), Tc=float(s["Tc"]),
                D0_0=float(s["D0_0"]), DS_0=float(s["DS_0"]), Deff_0=float(s["Deff_0"]),
                mu_0=float(s["mu_0"]), fracs=fr)


# --------------------------------------------------------------------------- scans
def coarse_tasks():
    t = []
    for VSO in (0.0, 0.4082, 0.8):                 # P1: vs U
        for Dt in (3.5, 4.5, 5.5):
            for U in np.arange(50.0, 200.1, 10.0):
                t.append((float(U), VSO, Dt))
    for U in (75.0, 95.0, 125.0, 150.0):           # P2: vs V_SO
        for Dt in (3.5, 4.5):
            for VSO in list(np.arange(0.0, 1.01, 0.1)) + [0.4082]:
                t.append((U, float(VSO), Dt))
    for (U, VSO) in ((95, 0.0), (95, 0.4082), (95, 0.8),
                     (125, 0.0), (125, 0.4082), (150, 0.4082)):   # P3: vs Delta_t
        for Dt in np.arange(0.5, 6.01, 0.25):
            t.append((float(U), float(VSO), float(Dt)))
    seen, out = set(), []
    for x in t:
        k = (round(x[0], 4), round(x[1], 4), round(x[2], 4))
        if k not in seen:
            seen.add(k)
            out.append(x)
    return out


def run_pool(tasks, nproc=2):
    from multiprocessing import Pool
    with Pool(nproc) as p:
        return [r for r in p.map(evaluate, tasks, chunksize=4) if r is not None]


# --------------------------------------------------------------------------- peaks
def find_local_maxima(x, y, prom_rel=0.10, prom_abs=0.015):
    """Interior strict local maxima with prominence; returns list of indices."""
    idx = []
    span = np.max(y) - np.min(y)
    prom = max(prom_abs, prom_rel * span)
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] >= y[i + 1]:
            left = np.min(y[:i + 1])
            right = np.min(y[i:])
            if y[i] - max(left, right) >= prom:
                idx.append(i)
    return idx


def is_monotone(y, tol=1e-9):
    d = np.diff(y)
    return bool(np.all(d >= -tol) or np.all(d <= tol))


def refine_peak(task_of_x, lo, hi, n=13):
    """Dense rescan of [lo, hi]."""
    xs = np.linspace(lo, hi, n)
    res = run_pool([task_of_x(v) for v in xs])
    return [r for r in res if r is not None and "fracs" in r]


def parabola_vertex(xv, yv, i):
    """Vertex of the parabola through the three points around index i."""
    xa, xb, xc = xv[i - 1], xv[i], xv[i + 1]
    ya, yb, yc = yv[i - 1], yv[i], yv[i + 1]
    den = ya * (xb - xc) + yb * (xc - xa) + yc * (xa - xb)
    if abs(den) < 1e-14:
        return float(xb)
    num = ya * (xb * xb - xc * xc) + yb * (xc * xc - xa * xa) + yc * (xa * xa - xb * xb)
    return float(num / (2.0 * den))


# --------------------------------------------------------------------------- main
if __name__ == "__main__":
    os.chdir(HERE)
    t_start = time.time()
    tasks = coarse_tasks()
    print(f"coarse scan: {len(tasks)} parameter sets x {len(FRACS)} temperatures", flush=True)
    coarse = run_pool(tasks)
    print(f"coarse done: {len(coarse)} valid points in {time.time()-t_start:.0f}s", flush=True)
    with open("results/param_scans_coarse.json", "w") as fh:
        json.dump(coarse, fh, indent=1)

    # organise into curves
    def curve_key(r, which):
        if which == "U":
            return ("U", round(r["VSO"], 4), round(r["Dt"], 4))
        if which == "VSO":
            return ("VSO", round(r["U"], 4), round(r["Dt"], 4))
        return ("Dt", round(r["U"], 4), round(r["VSO"], 4))

    def xval(r, which):
        return r["U"] if which == "U" else (r["VSO"] if which == "VSO" else r["Dt"])

    peaks, report = [], {}
    for which in ("U", "VSO", "Dt"):
        fams = {}
        for r in coarse:
            fams.setdefault(curve_key(r, which), []).append(r)
        report[which] = {}
        for key, rs in sorted(fams.items()):
            rs.sort(key=lambda r: xval(r, which))
            xv = np.array([xval(r, which) for r in rs])
            for fi, frac in enumerate(FRACS):
                yv = np.array([r["fracs"][fi]["ratio"] for r in rs])
                em = np.array([r["fracs"][fi]["Emin"] for r in rs])
                entry = dict(monotone=is_monotone(yv), maxima=[], emin_min=float(em.min()),
                             emin_arg=float(xv[int(np.argmin(em))]))
                for i in find_local_maxima(xv, yv):
                    entry["maxima"].append(dict(i=int(i), x=float(xv[i]), y=float(yv[i]),
                                                emin=float(em[i])))
                report[which][f"{key} frac={frac}"] = entry
                for m in entry["maxima"]:
                    peaks.append(dict(scan=which, key=list(key), frac=frac,
                                      x=m["x"], y=m["y"], emin=m["emin"],
                                      x_left=float(xv[m["i"] - 1]), x_right=float(xv[m["i"] + 1])))
    print(f"candidate maxima: {len(peaks)}", flush=True)

    # --------------------------------------------------------------- refinement
    def task_of(which, key):
        if which == "U":
            _, VSO, Dt = key
            return lambda v: dict(evaluate((float(v), VSO, Dt)) or {}, _x=float(v))
        if which == "VSO":
            _, U, Dt = key
            return lambda v: dict(evaluate((U, float(v), Dt)) or {}, _x=float(v))
        _, U, VSO = key
        return lambda v: dict(evaluate((U, VSO, float(v))) or {}, _x=float(v))

    refined = []
    for p in peaks:
        f = task_of(p["scan"], tuple(p["key"]))
        fi = FRACS.index(p["frac"])
        res = refine_peak(f, p["x_left"], p["x_right"], n=13)
        if not res:
            continue
        xv = np.array([r["_x"] for r in res])
        yv = np.array([r["fracs"][fi]["ratio"] for r in res])
        em = np.array([r["fracs"][fi]["Emin"] for r in res])
        i = int(np.argmax(yv))
        j = int(np.argmin(em))
        vert = parabola_vertex(xv, yv, i) if 0 < i < len(yv) - 1 else float(xv[i])
        # second, denser refinement on the half window (resolution check)
        lo, hi = xv[max(i - 1, 0)], xv[min(i + 1, len(xv) - 1)]
        res2 = refine_peak(f, lo, hi, n=25)
        xv2 = np.array([r["_x"] for r in res2]) if res2 else xv
        yv2 = np.array([r["fracs"][fi]["ratio"] for r in res2]) if res2 else yv
        i2 = int(np.argmax(yv2)) if len(yv2) else i
        p2 = dict(p, x_ref=float(xv[i]), y_ref=float(yv[i]), vertex=float(vert),
                  x_ref2=float(xv2[i2]), y_ref2=float(yv2[i2]),
                  emin_ref=float(em[j]), emin_arg_ref=float(xv[j]),
                  shift=abs(float(xv2[i2]) - float(vert)))
        refined.append(p2)
        print(f"  peak {p['scan']} {p['key']} frac={p['frac']}: coarse x={p['x']:.3f} "
              f"y={p['y']:.4f} -> refined x={p2['x_ref']:.3f} y={p2['y_ref']:.4f} "
              f"(vertex {vert:.3f}, dense {p2['x_ref2']:.3f}, E_min,min={p2['emin_ref']:.4f} "
              f"@ {p2['emin_arg_ref']:.3f})", flush=True)

    # --------------------------------------------------------------- artifact checks
    for p in refined:
        # rebuild the parameter triple from the scan type
        if p["scan"] == "U":
            U, VSO, Dt = p["x_ref"], p["key"][1], p["key"][2]
        elif p["scan"] == "VSO":
            U, VSO, Dt = p["key"][1], p["x_ref"], p["key"][2]
        else:
            U, VSO, Dt = p["key"][1], p["key"][2], p["x_ref"]
        r = evaluate((float(U), float(VSO), float(Dt)))
        if r is None:
            continue
        fi = FRACS.index(p["frac"])
        st_ratio = r["fracs"][fi]["ratio"]
        # (i) doubled quadrature resolution
        sc = SCState(params=dict(VSO=r["VSO"], Tc=r["Tc"], D0_0=r["D0_0"], DS_0=r["DS_0"],
                                 Deff_0=r["Deff_0"]), Nk_mu=200)
        st = sc.at(p["frac"])
        dos = get_dos(VSO)
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S_hi, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                               "a_signed", nell=2 * NELL)
        r_hi = S_hi / SN
        # (ii) independent k-space lattice evaluation
        lst = bf.lattice_states(VSO, 320, st["D0"], st["DS"], st["mu"])
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        Eg = (P["Emin"] * 0.995, float(lst["E"].max()) * 1.002)
        ws = np.array([0.08, 0.05, 0.03])
        vals = np.array([bf.S_sc_bf(lst, st["T"], OMEGA, "a_signed", Eg, 3000, w) for w in ws])
        A = np.vstack([np.ones_like(ws), ws ** 2]).T
        coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
        r_bf = float(coef[0]) / SN
        p["check_nell"] = float(r_hi)
        p["check_kpace"] = float(r_bf)
        p["check_dev_nell"] = abs(r_hi - st_ratio) / st_ratio
        p["check_dev_kpace"] = abs(r_bf - st_ratio) / st_ratio
        print(f"  artifact check @({U:.2f},{VSO:.3f},{Dt:.2f}) frac={p['frac']}: "
              f"E-space={st_ratio:.5f}  nell x2={r_hi:.5f} (dev {p['check_dev_nell']:.1e})  "
              f"k-space={r_bf:.5f} (dev {p['check_dev_kpace']:.1e})", flush=True)

    with open("results/param_scans_peaks.json", "w") as fh:
        json.dump(dict(report=report, peaks=refined), fh, indent=1)
    print(f"TOTAL {time.time()-t_start:.0f}s", flush=True)
