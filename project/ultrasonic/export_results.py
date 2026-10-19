"""
export_results.py -- fan out results/production_raw.json to every consumer, all paths
relative to this repository:
  ultrasonic/results/results.json + results.csv        (consumer schema of the wwn_eta1 pipeline)
  ultrasonic/results/groups.json + scans.json          (inputs of paper_correction figure builders)
  ultrasonic/results/*                                 (this repository)
  ultrasonic_U/results/*                               (this repository)
Run after run_production.py:  python export_results.py
"""
import json
import os

import numpy as np

import io_utils_ua as io
from parameters_ua import OMEGA, N_TARGET, GAP_FIT_P, GAP_FIT_Q
from scanning_ua import dos_for

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = io.load("production_raw.json")
VARS = ["a_signed", "b_abs", "c_paper", "d_std"]


def pkey(k):
    a, b = k if isinstance(k, tuple) else tuple(int(x) for x in k.split(","))
    return f"{a:+d}{b:+d}"


def row_out(r):
    o = dict(frac=r["frac"], T=r["T"], phi=r["phi"], D0=r["D0"], DS=r["DS"], Deff=r["Deff"],
             mu=r["mu"], mu_N=r["mu_N"], Emin=r["Emin"], two_Emin=r["2Emin"], kappa=r["kappa"],
             estar=r["estar"], S_N=r["S_N"], S_N_muSC=r["S_N_muSC"],
             R_delta0=r["R_delta0"], R_delta0_exact=r["R_delta0"])
    for v in VARS:
        if v in r:
            o[v] = r[v]
            o[v + "_ratio"] = r[v + "_ratio"]
            o[v + "_ratio_muSC"] = r[v + "_ratio_muSC"]
            o[v + "_parts"] = {pkey(k): val for k, val in r[v + "_parts"].items()}
            o[v + "_pair"] = r[v + "_pair"]
            o[v + "_pair_ratio"] = r[v + "_pair_ratio"]
            o[v + "_total_ratio"] = r[v + "_total_ratio"]
    return o


# ------------------------------------------------------------------ wwn_eta1/results.json
wj = {}
for v, rec in RAW["master"].items():
    p = rec["meta"]
    wj[v] = dict(VSO=p["VSO"], Tc=p["Tc"], D0_0=p["D0_0"], DS_0=p["DS_0"], Deff_0=p["Deff_0"],
                 Wband=dos_for(p["VSO"]).Wband, meta=p, hbar_omega=OMEGA, n_target=N_TARGET,
                 gap_fit_p=GAP_FIT_P, gap_fit_q=GAP_FIT_Q,
                 rows=[row_out(r) for r in rec["rows"]])
json.dump(wj, open(os.path.join(HERE, "results", "results.json"), "w"), indent=1)

# ------------------------------------------------------------------ csv / W_WN_table.csv
hdr = ["VSO", "T_over_Tc", "T", "Delta_eff", "E_min", "2E_min", "mu", "mu_N", "R_delta0_exact",
       "W/WN_(a)_signed", "W/WN_(b)_abs", "W/WN_(c)_paper", "W/WN_(d)_textbook",
       "pair_(a)", "total_(a)", "part_++", "part_+-", "part_-+", "part_--"]
lines = [",".join(hdr)]
for v, rec in RAW["master"].items():
    for r in rec["rows"]:
        pa = {pkey(k): v for k, v in r["a_signed_parts"].items()}
        key = lambda a, b: f"{a:+d}{b:+d}"
        lines.append(",".join(str(x) for x in [
            float(v), r["frac"], r["T"], r["Deff"], r["Emin"], r["2Emin"], r["mu"], r["mu_N"],
            r["R_delta0"], r["a_signed_ratio"], r["b_abs_ratio"], r["c_paper_ratio"],
            r["d_std_ratio"], r["a_signed_pair_ratio"], r["a_signed_total_ratio"],
            pa[key(1, 1)], pa[key(-1, 1)], pa[key(1, -1)], pa[key(-1, -1)]]))
csv_txt = "\n".join(lines) + "\n"
open(os.path.join(HERE, "results", "results.csv"), "w").write(csv_txt)
open(f"{HERE}/results/W_WN_table.csv", "w").write(csv_txt)

# ------------------------------------------------------------------ repository temp sweeps
man = {}
for v, rec in RAW["master"].items():
    p = rec["meta"]
    rows = rec["rows"]
    out = dict(V_SO=p["VSO"], Tc=p["Tc"], Deff_0=p["Deff_0"], D0_0=p["D0_0"], DS_0=p["DS_0"],
               hbar_omega=OMEGA, n_target=N_TARGET,
               fracs=[r["frac"] for r in rows],
               ratios={var: [r[var + "_ratio"] for r in rows] for var in VARS} | {
                   "pair": [r["a_signed_pair_ratio"] for r in rows],
                   "total": [r["a_signed_total_ratio"] for r in rows]},
               Emin=[r["Emin"] for r in rows], two_Emin=[r["2Emin"] for r in rows],
               Deff=[r["Deff"] for r in rows])
    name = f"temp_sweep_VSO_{float(v):.4f}.json"
    json.dump(out, open(f"{HERE}/results/{name}", "w"), indent=1)
    man[f"{float(v):.4f}"] = name
json.dump(man, open(f"{HERE}/results/temp_sweep_manifest.json", "w"), indent=1)

# ------------------------------------------------------------------ groups
gw = []
for g in RAW["groups"]:
    rows = g["rows"]
    sw = dict(set=g["set"], stage1=g["stage1"], fracs=[r["frac"] for r in rows],
              a_ratio=[r["a_signed_ratio"] for r in rows],
              ratios=dict(a_signed=[r["a_signed_ratio"] for r in rows]),
              Emin=[r["Emin"] for r in rows],
              aux=dict(T=[r["T"] for r in rows], Deff=[r["Deff"] for r in rows],
                       mu=[r["mu"] for r in rows], mu_N=[r["mu_N"] for r in rows]))
    gw.append(sw)
    if g["set"]["group"] in ("A", "B"):
        json.dump(sw, open(f"{HERE}/../ultrasonic_U/results/"
                           f"temp_sweep_U_{g['set']['U']:.2f}.json", "w"), indent=1)
json.dump(gw, open(f"{HERE}/results/groups_WWN.json", "w"), indent=1)
json.dump(gw, open(f"{HERE}/../ultrasonic_U/results/groups_WWN_AB.json", "w"), indent=1)
json.dump(gw, open(os.path.join(HERE, "results", "groups.json"), "w"), indent=1)

# ------------------------------------------------------------------ scans
vs, ps = {}, {}
for key, rec in RAW["scans"].items():
    fr = [f for f in (0.3, 0.5, 0.7, 0.9)]
    ratio = {str(f): [] for f in fr}
    for i, rows in enumerate(rec["rows"]):
        for f in fr:
            ratio[str(f)].append(None if rows is None else
                                 [r for r in rows if abs(r["frac"] - f) < 1e-9][0]["a_signed_ratio"])
    vs[key] = dict(V_SO=rec["V_SO"], ratio=ratio)
    ps[key] = dict(U=float(key), Delta_t=4.5, n=1.885, V_SO=rec["V_SO"], Tc=rec["Tc"],
                   Deff_0=rec["Deff_0"], D0_0=rec["D0_0"], DS_0=rec["DS_0"], ratio=ratio)
json.dump(vs, open(f"{HERE}/results/vso_scan.json", "w"), indent=1)
json.dump(ps, open(os.path.join(HERE, "results", "scans.json"), "w"), indent=1)
print("export done:", len(wj), "master sweeps,", len(gw), "groups,", len(vs), "scans")
