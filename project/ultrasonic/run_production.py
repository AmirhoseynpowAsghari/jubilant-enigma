"""
run_production.py -- recompute ALL production numbers with the corrected energy-space
measure (Levy-1963 validated; see levy1963_check.py and the absorption.py docstring).
Writes results/production_raw.json; export_results.py fans it out to every consumer.
Run:  python run_production.py        (about 10 min)
"""
import json
import os
import time
import numpy as np

import sc_loader
import io_utils_ua as io
from parameters_ua import OMEGA
from absorption import ratio_all, S_paircreation, S_normal
from scanning_ua import dos_for, FRACS, sc_from_stage1

FOUR = ["a_signed", "b_abs", "c_paper", "d_std"]
SCAN_FRACS = (0.3, 0.5, 0.7, 0.9)


def closed_R(T):
    E = np.linspace(1e-12, 400.0 * max(T, OMEGA), 400001)
    f1 = 1.0 / (np.exp(E / T) + 1.0)
    f2 = 1.0 / (np.exp((E + OMEGA) / T) + 1.0)
    return float(2.0 * np.trapezoid(f1 - f2, E) / OMEGA)


def state_rows(sc, fracs, variants):
    dos = dos_for(sc.VSO)
    rows = []
    for f in fracs:
        st = sc.at(f)
        r = ratio_all(dos, st, OMEGA, variant_list=variants, nell=6000)
        rn_mu = S_normal(dos, st["T"], OMEGA, st["mu"])
        row = dict(r)
        for v in variants:
            row[v + "_parts"] = {f"{a},{b}": w for (a, b), w in row[v + "_parts"].items()}
        row["phi"] = st["Deff"] / sc.Deff_0
        row["S_N_muSC"] = rn_mu
        row["R_delta0"] = closed_R(st["T"])
        for v in variants:
            Sp, _ = S_paircreation(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], v, nell=6000)
            row[v + "_pair"] = Sp
            row[v + "_pair_ratio"] = Sp / r["S_N"]
            row[v + "_total_ratio"] = (r[v] + Sp) / r["S_N"]
            row[v + "_ratio_muSC"] = r[v] / rn_mu
        rows.append(row)
        print(f"    frac={f:.3f}  a={row['a_signed_ratio']:.6f}", flush=True)
    return rows


def main():
    t0 = time.time()
    raw = {"master": {}, "groups": [], "scans": {}}
    print("== master sweeps ==", flush=True)
    for v in (0.0, 0.4082, 0.8):
        p = sc_loader.sc_params_at(v)
        sc = sc_from_stage1(p)
        raw["master"][f"{v:.4f}"] = dict(meta=p, rows=state_rows(sc, FRACS, FOUR))
        io.save("production_raw_partial.json", raw)
    print("== group sets ==", flush=True)
    gp = io.load("../sc_params/results/groups_stage1.json") if False else json.load(
        open(os.path.join(os.path.dirname(__file__), "..", "sc_params", "results",
                          "groups_stage1.json")))
    N = {"A": 0.8, "B": 1.885, "C": 1.885, "D": 1.885}
    for p in gp:
        s = p["set"]
        sc = sc_from_stage1(p, n_tgt=N[s["group"]])
        rows = state_rows(sc, FRACS, ["a_signed"])
        io.save("production_raw_partial.json", raw)
        raw["groups"].append(dict(set=s, stage1={k: p[k] for k in
                                 ("U", "Delta_t", "n", "V_SO", "Tc", "D0_0", "DS_0", "Deff_0", "mu_0")},
                                  rows=rows))
    print("== V_SO scans ==", flush=True)
    us = json.load(open(os.path.join(os.path.dirname(__file__), "..", "sc_params", "results",
                                     "u_scans_stage1.json")))
    for key, rec in us.items():
        entry = dict(V_SO=[], Tc=[], Deff_0=[], D0_0=[], DS_0=[], rows=[])
        for i, v in enumerate(rec["V_SO"]):
            p = dict(VSO=v, Tc=rec["Tc"][i], D0_0=rec["D0_0"][i], DS_0=rec["DS_0"][i],
                     Deff_0=rec["Deff_0"][i])
            entry["V_SO"].append(v)
            entry["Tc"].append(rec["Tc"][i]); entry["Deff_0"].append(rec["Deff_0"][i])
            entry["D0_0"].append(rec["D0_0"][i]); entry["DS_0"].append(rec["DS_0"][i])
            if rec["Tc"][i] <= 0:
                entry["rows"].append(None)
                continue
            sc = sc_from_stage1(p)
            entry["rows"].append(state_rows(sc, SCAN_FRACS, ["a_signed"]))
            print(f"  U={key} V_SO={v:.3f} done", flush=True)
        raw["scans"][key] = entry
    io.save("production_raw.json", raw)
    print(f"production_raw.json written in {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
