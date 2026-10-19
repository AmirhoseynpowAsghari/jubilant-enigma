"""
param_scan_features.py -- focused (U, V_SO, Delta_t) scans for NON-MONOTONIC FEATURES.

Follow-up to param_scans_peaks.py, whose coarse fixed-T/Tc scan of the primary channel
(eta=+1, signed, scattering only) found NO local maxima: every family is monotone or
sigmoidal, because the E_min -> 0 ridge (D0 = b*mu) coincides with the COLLAPSE of the
superconducting dome (Tc -> 0), so at fixed T/Tc the ratio dies smoothly instead of
peaking.  Genuine features can still appear in two physical channels, scanned here:

  (R1) HEISENBERG-SLICHTER / coherence resonance, Case-II channel (eta = -1 kernel
       C_II = 1/2[1 + (xi xi' + Dbar Dbar')/(EE')]): at FIXED temperature the ratio is
       enhanced when 2*Delta(T) ~ hbar*omega -- a true resonance in any parameter that
       moves Delta(T) across hbar*omega/2 (U, Delta_t, V_SO all do).
  (R2) PAIR-CREATION THRESHOLD (WITHDRAWN as quantitative claim, REPORT App. A2): the
       combination S_sc + S_pair was examined as an OUTLOOK; its >1 features are withdrawn.
       The scattering-only curve remains the quantitative result; the pair-open region is
       shaded for orientation only.

Ensembles per parameter set:
   frac = 0.5, 0.9      : a_signed (primary; reused from the coarse scan where possible)
   T = 0.10 t, 0.20 t   : a_signed, a_signed+pair (total), case2_signed

Families:
   F_U   : U = 50..200 step 7.5,  (V_SO, Delta_t) in {(0,4.5), (0.4082,4.5)}
   F_Dt  : Delta_t = 0.5..6 step 0.25, (U,V_SO) in {(95,0),(95,0.4082),(125,0),(125,0.4082)}
   F_VSO : V_SO = 0..1 step 0.1 (+0.4082), (U,Delta_t) in {(95,4.5),(125,4.5)}

Peak protocol as before: prominence filter -> 13-pt and 25-pt refinement -> parabola
vertex -> artifact checks (doubled nell; independent k-space lattice sum, which now also
supports the case2_signed kernel).

Outputs: results/param_features.json, figures/param_scan_{U,VSO,Dt}.png (2x2),
         figures/param_scan_diagnostics.png
Run:  python param_scan_features.py        (~6-10 min on 2 cores)
"""
from __future__ import annotations

import os
import sys
import json
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))

import solvers as s1sol
from dos_ua import NormalDOS
from physics_ua import SCState
from parameters_ua import OMEGA
from absorption import S_scattering, S_normal, S_paircreation, qp_constants
import bruteforce_ua as bf
from param_scans_peaks import (get_dos, find_local_maxima, is_monotone, parabola_vertex,
                               FRACS)


def run_pool_local(tl, nproc=2):
    """Pool over THIS module's evaluate (param_scans_peaks.run_pool would use its own)."""
    from multiprocessing import Pool
    with Pool(nproc) as p:
        return [r for r in p.map(evaluate, tl, chunksize=4) if r is not None]

TS = (0.10, 0.20)          # fixed physical temperatures (in t)
NELL = 20000


def evaluate(task):
    U, VSO, Dt = task
    try:
        s = s1sol.solve_set(U, Dt, 1.885, VSO, Nk=80)
    except Exception:
        return None
    if not (s["Tc"] > 1e-9 and 0.0 < s["D0_0"] < 10.0 and 0.0 <= s["DS_0"] < 20.0):
        return None
    sc = SCState(params=dict(VSO=float(s["V_SO"]), Tc=float(s["Tc"]), D0_0=float(s["D0_0"]),
                             DS_0=float(s["DS_0"]), Deff_0=float(s["Deff_0"])), Nk_mu=200)
    dos = get_dos(VSO)
    out = dict(U=float(U), VSO=float(VSO), Dt=float(Dt), Tc=float(s["Tc"]),
               D0_0=float(s["D0_0"]), DS_0=float(s["DS_0"]), mu_0=float(s["mu_0"]),
               Deff_0=float(s["Deff_0"]), fracs={}, temps={})
    for frac in FRACS:
        st = sc.at(frac)
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                            "a_signed", nell=NELL)
        out["fracs"][str(frac)] = dict(T=float(st["T"]), a=float(S / SN),
                                       Emin=float(P["Emin"]))
    for T in TS:
        frac = T / s["Tc"] if s["Tc"] > 0 else 1e9
        st = sc.at(frac)                    # frac>1 -> normal state (phi=0), T exact
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                            "a_signed", nell=NELL)
        S2, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                             "case2_signed", nell=NELL)
        Sp, _ = S_paircreation(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                               "a_signed", nell=NELL)
        out["temps"][f"{T:.2f}"] = dict(
            T=float(st["T"]), sc=float(st["T"] < s["Tc"]),
            a=float(S / SN), tot=float((S + Sp) / SN), c2=float(S2 / SN),
            pair=float(Sp / SN), Emin=float(P["Emin"]), D0=float(st["D0"]),
            DS=float(st["DS"]), mu=float(st["mu"]),
            pair_open=bool(2.0 * P["Emin"] < OMEGA))
    return out


def tasks():
    t = []
    for VSO in (0.0, 0.4082):
        for U in np.arange(50.0, 200.1, 7.5):
            t.append((float(U), VSO, 4.5))
    for (U, VSO) in ((95, 0.0), (95, 0.4082), (125, 0.0), (125, 0.4082)):
        for Dt in np.arange(0.5, 6.01, 0.25):
            t.append((float(U), float(VSO), float(Dt)))
    for U in (95.0, 125.0):
        for VSO in list(np.arange(0.0, 1.01, 0.1)) + [0.4082]:
            t.append((U, float(VSO), 4.5))
    seen, out = set(), []
    for x in t:
        k = tuple(round(v, 4) for v in x)
        if k not in seen:
            seen.add(k)
            out.append(x)
    return out


# --------------------------------------------------------------------------- curves
CHAN = [("fracs", "0.5", "a", "a_signed, T/Tc=0.5"),
        ("fracs", "0.9", "a", "a_signed, T/Tc=0.9"),
        ("temps", "0.10", "a", "a_signed, T=0.1t"),
        ("temps", "0.10", "tot", "total (sc+pair), T=0.1t"),
        ("temps", "0.10", "c2", "Case-II, T=0.1t"),
        ("temps", "0.20", "a", "a_signed, T=0.2t"),
        ("temps", "0.20", "tot", "total (sc+pair), T=0.2t"),
        ("temps", "0.20", "c2", "Case-II, T=0.2t")]


def curves(data, which):
    """{family_key: {chan_label: (x, y, emin, pair_open)}}"""
    fams = {}
    for r in data:
        if which == "U":
            key = (round(r["VSO"], 4), round(r["Dt"], 4))
            x = r["U"]
        elif which == "VSO":
            key = (round(r["U"], 4), round(r["Dt"], 4))
            x = r["VSO"]
        else:
            key = (round(r["U"], 4), round(r["VSO"], 4))
            x = r["Dt"]
        f = fams.setdefault(key, {c[3]: ([], [], [], []) for c in CHAN})
        for grp, gk, ch, lab in CHAN:
            d = r[grp][gk]
            f[lab][0].append(x)
            f[lab][1].append(d[ch])
            f[lab][2].append(d.get("Emin", np.nan))
            f[lab][3].append(bool(d.get("pair_open", False)))
    for key in fams:
        for lab in fams[key]:
            order = np.argsort(fams[key][lab][0])
            fams[key][lab] = [np.asarray(v)[order] for v in fams[key][lab]]
    return fams


if __name__ == "__main__":
    os.chdir(HERE)
    t0 = time.time()
    tl = tasks()
    data = None
    if os.path.exists("results/param_features.json"):
        with open("results/param_features.json") as fh:
            cand = json.load(fh)
        if cand and "temps" in cand[0]:
            data = cand
            print(f"loaded cached scan: {len(data)} points", flush=True)
    if data is None:
        print(f"feature scan: {len(tl)} sets", flush=True)
        data = run_pool_local(tl)
        print(f"done {len(data)} valid in {time.time()-t0:.0f}s", flush=True)
        with open("results/param_features.json", "w") as fh:
            json.dump(data, fh, indent=1)

    summary = {}
    peaks = []
    for which in ("U", "VSO", "Dt"):
        fams = curves(data, which)
        summary[which] = {}
        for key, chans in sorted(fams.items()):
            for lab, (x, y, em, po) in chans.items():
                mono = is_monotone(y)
                mx = find_local_maxima(x, y, prom_rel=0.06, prom_abs=0.01)
                summary[which][f"{key} | {lab}"] = dict(
                    monotone=mono, maxima=[dict(x=float(x[i]), y=float(y[i])) for i in mx])
                for i in mx:
                    peaks.append(dict(scan=which, key=list(key), chan=lab,
                                      x=float(x[i]), y=float(y[i]),
                                      xl=float(x[i - 1]), xr=float(x[i + 1])))
    with open("results/param_features_summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    n_mono = sum(1 for w in summary.values() for e in w.values() if e["monotone"])
    n_tot = sum(len(w) for w in summary.values())
    print(f"curves: {n_tot}, monotone: {n_mono}, candidate maxima: {len(peaks)}", flush=True)
    for p in peaks:
        print(f"  cand {p['scan']} {p['key']} {p['chan']}: x={p['x']:.3f} y={p['y']:.4f}",
              flush=True)

    # ------------------------------------------------------------- refinement
    def ychan(r, chan):
        for grp, gk, ch, lab in CHAN:
            if lab == chan:
                return r[grp][gk][ch]
        raise KeyError(chan)

    def xof(r, which):
        return r["U"] if which == "U" else (r["VSO"] if which == "VSO" else r["Dt"])

    def triples(which, key, xs):
        out = []
        for v in xs:
            if which == "U":
                out.append((float(v), key[0], key[1]))
            elif which == "VSO":
                out.append((key[0], float(v), key[1]))
            else:
                out.append((key[0], key[1], float(v)))
        return out

    refined = []
    for p in peaks:
        key = tuple(p["key"])
        res = run_pool_local(triples(p["scan"], key, np.linspace(p["xl"], p["xr"], 13)))
        if not res:
            continue
        xv = np.array([xof(r, p["scan"]) for r in res])
        yv = np.array([ychan(r, p["chan"]) for r in res])
        i = int(np.argmax(yv))
        vert = parabola_vertex(xv, yv, i) if 0 < i < len(yv) - 1 else float(xv[i])
        lo, hi = xv[max(i - 1, 0)], xv[min(i + 1, len(xv) - 1)]
        res2 = run_pool_local(triples(p["scan"], key, np.linspace(lo, hi, 25)))
        xv2 = np.array([xof(r, p["scan"]) for r in res2]) if res2 else xv
        yv2 = np.array([ychan(r, p["chan"]) for r in res2]) if res2 else yv
        i2 = int(np.argmax(yv2)) if len(yv2) else i
        p2 = dict(p, x_ref=float(xv[i]), y_ref=float(yv[i]), vertex=float(vert),
                  x_ref2=float(xv2[i2]), y_ref2=float(yv2[i2]),
                  shift=abs(float(xv2[i2]) - float(vert)))
        refined.append(p2)
        print(f"  peak {p['scan']} {p['key']} {p['chan']}: x={p2['x_ref']:.3f} "
              f"y={p2['y_ref']:.4f} vertex={vert:.3f} dense={p2['x_ref2']:.3f} "
              f"(shift {p2['shift']:.3f})", flush=True)

    # ------------------------------------------------------------- artifact checks
    for p in refined:
        if p["scan"] == "U":
            U, VSO, Dt = p["x_ref"], p["key"][0], p["key"][1]
        elif p["scan"] == "VSO":
            U, VSO, Dt = p["key"][0], p["x_ref"], p["key"][1]
        else:
            U, VSO, Dt = p["key"][0], p["key"][1], p["x_ref"]
        r = evaluate((float(U), float(VSO), float(Dt)))
        if r is None:
            continue
        var = "case2_signed" if p["chan"].startswith("Case-II") else "a_signed"
        tot = p["chan"].startswith("total")
        T = 0.10 if "0.1t" in p["chan"] else (0.20 if "0.2t" in p["chan"] else None)
        if T is None:
            continue
        sc = SCState(params=dict(VSO=r["VSO"], Tc=r["Tc"], D0_0=r["D0_0"], DS_0=r["DS_0"],
                                 Deff_0=r["Deff_0"]), Nk_mu=200)
        st = sc.at(T / r["Tc"])
        dos = get_dos(VSO)
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        S_hi, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], var,
                               nell=2 * NELL)
        if tot:
            Sp_hi, _ = S_paircreation(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"],
                                      var, nell=2 * NELL)
            S_hi += Sp_hi
        lst = bf.lattice_states(VSO, 320, st["D0"], st["DS"], st["mu"])
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        Eg = (P["Emin"] * 0.995, float(lst["E"].max()) * 1.002)
        Eg_pair = (max(P["Emin"] * 0.995, 1e-9), OMEGA * 1.002)
        ws = np.array([0.08, 0.05, 0.03])
        vals = np.array([bf.S_sc_bf(lst, st["T"], OMEGA, var, Eg, 3000, w) for w in ws])
        A = np.vstack([np.ones_like(ws), ws ** 2]).T
        coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
        k_sc = float(coef[0])
        if tot:   # pair leg: fixed small width (Richardson invalid, see S_pair_bf doc)
            k_sc += bf.S_pair_bf(lst, st["T"], OMEGA, var, Eg_pair, 3000, 0.005)
        p["check_nell"] = float(S_hi / SN)
        p["check_kpace"] = float(k_sc / SN)
        p["dev_nell"] = abs(p["check_nell"] - p["y_ref"]) / max(p["y_ref"], 1e-12)
        p["dev_kpace"] = abs(p["check_kpace"] - p["y_ref"]) / max(p["y_ref"], 1e-12)
        print(f"  check @({U:.1f},{VSO:.3f},{Dt:.2f}) {p['chan']}: E={p['y_ref']:.5f} "
              f"nell2={p['check_nell']:.5f} k={p['check_kpace']:.5f} "
              f"devs {p['dev_nell']:.1e}/{p['dev_kpace']:.1e}", flush=True)

    with open("results/param_features_peaks.json", "w") as fh:
        json.dump(refined, fh, indent=1)
    print(f"TOTAL {time.time()-t0:.0f}s", flush=True)
