"""
scanning_ua.py -- production stage-2 runs: T sweeps and V_SO scans of the corrected W/W_N.
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
"""
from __future__ import annotations

import json
import os
import time
import numpy as np

from parameters_ua import OMEGA, N_TARGET
from dos_ua import NormalDOS
from physics_ua import SCState
from absorption import ratio_all
import io_utils_ua as io

FRACS = list(np.concatenate([np.linspace(0.02, 0.98, 49), [0.99, 0.995, 0.999]]))
_dos = {}


def dos_for(v):
    if v not in _dos:
        _dos[v] = NormalDOS(v)
    return _dos[v]


def sc_from_stage1(p, n_tgt=N_TARGET):
    return SCState.from_params(p["VSO"] if "VSO" in p else p["V_SO"], p["Tc"], p["D0_0"],
                               p["DS_0"], p["Deff_0"], n_tgt=n_tgt)


def temp_sweep(sc, fracs=FRACS, variants=("a_signed", "b_abs", "c_paper", "d_std"),
               with_pair=False):
    """corrected W/W_N(T/Tc) for one SC state; optionally the pair-creation channel too."""
    from absorption import S_paircreation
    dos = dos_for(sc.VSO)
    from physics_ua import phi_gap
    out = dict(V_SO=sc.VSO, Tc=sc.Tc, Deff_0=sc.Deff_0, D0_0=sc.D0_0, DS_0=sc.DS_0,
               hbar_omega=OMEGA, n_target=sc.n_tgt, fracs=[], ratios={}, Emin=[],
               two_Emin=[], Deff=[])
    for v in variants:
        out["ratios"][v] = []
    out["ratios"]["pair"] = []
    for f in fracs:
        st = sc.at(f)
        r = ratio_all(dos, st, OMEGA, variant_list=list(variants))
        out["fracs"].append(f)
        for v in variants:
            out["ratios"][v].append(r[v + "_ratio"])
        if with_pair:
            Sp, _ = S_paircreation(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], "a_signed")
            out["ratios"]["pair"].append(Sp / r["S_N"])
        else:
            out["ratios"]["pair"].append(None)
        out["Emin"].append(r["Emin"])
        out["two_Emin"].append(2.0 * r["Emin"])
        out["Deff"].append(phi_gap(st["T"], sc.Tc) * sc.Deff_0)
        out["ratios"].setdefault("total", []).append(
            (r["a_signed_ratio"] + (Sp / r["S_N"] if with_pair else 0.0)) if with_pair
            else r["a_signed_ratio"])
    return out


def run_master(with_pair=True):
    """the repository's master set (U=95t, Dt=4.5t, n=1.885) at V_SO = 0, 0.4082, 0.8."""
    import sc_loader
    man = {}
    for v in (0.0, 0.4082, 0.8):
        p = sc_loader.sc_params_at(v)
        sw = temp_sweep(sc_from_stage1(p), with_pair=with_pair)
        io.save(f"temp_sweep_VSO_{v:.4f}.json", sw)
        man[f"{v:.4f}"] = f"temp_sweep_VSO_{v:.4f}.json"
        print(f"  V_SO={v:.4f}  lim W/W_N = {sw['ratios']['a_signed'][-1]:.5f}", flush=True)
    io.save("temp_sweep_manifest.json", man)


def run_groups(groups_stage1, n_of_set):
    """corrected sweeps for every group set of the paper (stage-1 input from sc_params)."""
    out = []
    for p in groups_stage1:
        s = p["set"]
        sc = sc_from_stage1(p, n_tgt=n_of_set(s))
        sw = temp_sweep(sc)
        sw["set"] = s
        sw["stage1"] = {k: p[k] for k in ("Tc", "D0_0", "DS_0", "Deff_0")}
        out.append(sw)
        print(f"  {s['group']} {s['label']:8s} lim = {sw['ratios']['a_signed'][-1]:.5f}", flush=True)
    io.save("groups_WWN.json", out)
    return out


def run_vso_scans(scans_stage1, fracs=(0.3, 0.5, 0.7, 0.9)):
    """W/W_N vs V_SO at fixed T/Tc for U = 95, 125, -75 (stage-1 from sc_params)."""
    out = {}
    for key, rec in scans_stage1.items():
        entry = dict(V_SO=rec["V_SO"], ratio={str(f): [] for f in fracs})
        for i, v in enumerate(rec["V_SO"]):
            p = dict(VSO=v, Tc=rec["Tc"][i], D0_0=rec["D0_0"][i], DS_0=rec["DS_0"][i],
                     Deff_0=rec["Deff_0"][i])
            if p["Tc"] <= 0:
                for f in fracs:
                    entry["ratio"][str(f)].append(None)
                continue
            sc = sc_from_stage1(p)
            for f in fracs:
                st = sc.at(f)
                r = ratio_all(dos_for(v), st, OMEGA, variant_list=["a_signed"])
                entry["ratio"][str(f)].append(r["a_signed_ratio"])
        out[key] = entry
        print(f"  U={key} scan done", flush=True)
    io.save("vso_scan.json", out)
    return out
