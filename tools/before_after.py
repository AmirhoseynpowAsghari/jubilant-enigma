"""
before_after.py -- builds the before/after comparison of the 2026-10-24 correction campaign
directly from the archived pre-correction outputs (archive/pre_f02_baseline/) and the current
regenerated outputs.  Writes results/BEFORE_AFTER.json and CORRECTION_TABLES.md.

Nothing here is typed by hand: every number in CORRECTION_TABLES.md is read from a file.
Run:  python3 tools/before_after.py
"""
import json, os, math
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OLD = os.path.join(ROOT, "archive", "pre_f02_baseline")
NEW_S1 = os.path.join(ROOT, "project", "sc_params", "results")
NEW_S2 = os.path.join(ROOT, "project", "ultrasonic", "results")
NEW_U = os.path.join(ROOT, "project", "ultrasonic_U", "results")


def J(p):
    with open(p) as fh:
        return json.load(fh)


def nearest(vlist, v):
    return min(range(len(vlist)), key=lambda i: abs(vlist[i] - v))


def pct(new, old):
    if old in (0, None) or new is None:
        return None
    return 100.0 * (new - old) / abs(old)


def main():
    out = {"stage1_master": [], "stage1_groups": [], "production_master": [],
           "groups_limits": [], "soc_explicit": [], "param_peaks": {}, "files": {}}
    # ---------------- stage 1, master V_SO points
    old = J(os.path.join(OLD, "sc_params_results", "tc_vs_vso.json"))["results"]
    new = J(os.path.join(NEW_S1, "tc_vs_vso.json"))["results"]
    for v in (0.0, 0.4082, 0.8):
        i_o = nearest(old["V_SO"], v); i_n = nearest(new["V_SO"], v)
        row = dict(V_SO=v)
        for key_o, key_n, lab in (("Tc", "Tc", "Tc"), ("D0_T0", "D0_T0", "D0_0"),
                                  ("DS_T0", "DS_T0", "DS_0"), ("Deff_T0", "Deff_T0", "Deff_0")):
            a, b = old[key_o][i_o], new[key_n][i_n]
            row[lab] = dict(old=a, new=b, pct=pct(b, a))
        out["stage1_master"].append(row)
    # ---------------- stage 1, group sets
    og = J(os.path.join(OLD, "sc_params_results", "groups_stage1.json"))
    ng = J(os.path.join(NEW_S1, "groups_stage1.json"))
    for o, n in zip(og, ng):
        s = n["set"]
        out["stage1_groups"].append(dict(label=f"{s['group']} {s['label']}",
                                         Tc_old=o["Tc"], Tc_new=n["Tc"],
                                         pct=pct(n["Tc"], o["Tc"]),
                                         no_inband_root=bool(n.get("no_inband_root", False))))
    # ---------------- production, master V_SO sweeps
    for v in ("0.0000", "0.4082", "0.8000"):
        o = J(os.path.join(OLD, "ultrasonic_results", f"temp_sweep_VSO_{v}.json"))
        n = J(os.path.join(NEW_S2, f"temp_sweep_VSO_{v}.json"))
        row = dict(V_SO=float(v))
        for var in ("a_signed", "b_abs", "c_paper", "d_std"):
            row[var] = {}
            for fr in (0.5, 0.9, 0.999):
                io_ = nearest(o["fracs"], fr); in_ = nearest(n["fracs"], fr)
                a = o["ratios"][var][io_]; b = n["ratios"][var][in_]
                row[var][str(fr)] = dict(old=a, new=b, pct=pct(b, a))
        out["production_master"].append(row)
    # ---------------- group limits (a_signed at T/Tc=0.999)
    og = J(os.path.join(OLD, "ultrasonic_results", "groups_WWN.json"))
    ng = J(os.path.join(NEW_S2, "groups_WWN.json"))
    # the 13 sets are in the same order in both files; labels differ only in spacing
    assert len(og) == len(ng) == 13 and all(a["set"]["group"] == b["set"]["group"] for a, b in zip(og, ng))
    for o, n in zip(og, ng):
        lab = f"{o['set']['group']} {o['set']['label']}"
        a_old = o["a_ratio"][-1] if isinstance(o.get("a_ratio"), list) else None
        a_new = (n["a_ratio"][-1] if (n and isinstance(n.get("a_ratio"), list) and n["a_ratio"]
                                      and n["a_ratio"][-1] is not None) else None)
        out["groups_limits"].append(dict(label=lab, old=a_old, new=a_new, pct=pct(a_new, a_old),
                                         no_superconductivity=(a_new is None)))
    # ---------------- SOC explicit decomposition
    try:
        so = J(os.path.join(OLD, "ultrasonic_results", "soc_explicit.json"))
        sn = J(os.path.join(NEW_S2, "soc_explicit.json"))
        for v in ("0.0000", "0.2000", "0.4082", "0.6000", "0.8000"):
            for fr in ("0.5", "0.9"):
                if v in so and v in sn:
                    out["soc_explicit"].append(dict(V_SO=float(v), frac=float(fr),
                        RA_old=so[v][fr]["RA"], RA_new=sn[v][fr]["RA"],
                        RB_old=so[v][fr]["RB"], RB_new=sn[v][fr]["RB"],
                        inter_old=so[v][fr]["inter"], inter_new=sn[v][fr]["inter"]))
    except (FileNotFoundError, KeyError) as e:
        out["soc_explicit_error"] = str(e)
    # ---------------- parameter peaks (counts and locations)
    for name in ("param_features_peaks.json", "param_scans_peaks.json"):
        po = os.path.join(OLD, "ultrasonic_results", name)
        pn = os.path.join(NEW_S2, name)
        if os.path.exists(pn):
            n = J(pn)
            n_pk = len(n["peaks"]) if isinstance(n, dict) and "peaks" in n else len(n)
            o_pk = None
            if os.path.exists(po):
                o = J(po)
                o_pk = len(o["peaks"]) if isinstance(o, dict) and "peaks" in o else len(o)
            out["param_peaks"][name] = dict(old=o_pk, new=n_pk)
    # ---------------- lists of files compared
    out["files"] = dict(old_dir="archive/pre_f02_baseline", new_dirs=["project/sc_params/results",
                        "project/ultrasonic/results", "project/ultrasonic_U/results"])
    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    with open(os.path.join(ROOT, "results", "BEFORE_AFTER.json"), "w") as fh:
        json.dump(out, fh, indent=1, allow_nan=False, default=float)

    # ---------------- Markdown tables
    L = ["# Before / after: numerical consequences of the correction campaign", "",
         "Generated by `tools/before_after.py` from `archive/pre_f02_baseline/` (before) and the",
         "regenerated outputs (after).  Percentages are (new - old)/|old|.", ""]
    L += ["## Stage 1: Tc, gap and chemical potential on the master set (U=95t, Dt=4.5t, n=1.885)", "",
          "| V_SO | Tc old | Tc new | d%Tc | D0(0) old | D0(0) new | Deff(0) old | Deff(0) new | d%Deff |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in out["stage1_master"]:
        L.append(f"| {r['V_SO']:.4f} | {r['Tc']['old']:.5f} | {r['Tc']['new']:.5f} | {r['Tc']['pct']:+.2f} | "
                 f"{r['D0_0']['old']:.5f} | {r['D0_0']['new']:.5f} | {r['Deff_0']['old']:.5f} | "
                 f"{r['Deff_0']['new']:.5f} | {r['Deff_0']['pct']:+.2f} |")
    L += ["", "## Stage 1: the 13 parameter sets of arXiv:2607.02812", "",
          "| set | Tc old | Tc new | d%Tc | flag |", "|---|---|---|---|---|"]
    for r in out["stage1_groups"]:
        p = "" if r["pct"] is None else f"{r['pct']:+.2f}"
        flag = "no superconducting state in the window (mu outside band)" if r["no_inband_root"] else ""
        L.append(f"| {r['label']} | {r['Tc_old']:.5f} | {r['Tc_new']:.5f} | {p} | {flag} |")
    L += ["", "## Production W/W_N on the master set (scattering channel, all four coherence variants)", "",
          "| V_SO | variant | T/Tc=0.5 old | new | T/Tc=0.9 old | new | T/Tc=0.999 old | new | d%(0.999) |",
          "|---|---|---|---|---|---|---|---|---|"]
    for r in out["production_master"]:
        for var in ("a_signed", "b_abs", "c_paper", "d_std"):
            c = r[var]
            L.append(f"| {r['V_SO']:.4f} | {var} | {c['0.5']['old']:.5f} | {c['0.5']['new']:.5f} | "
                     f"{c['0.9']['old']:.5f} | {c['0.9']['new']:.5f} | {c['0.999']['old']:.5f} | "
                     f"{c['0.999']['new']:.5f} | {c['0.999']['pct']:+.2f} |")
    L += ["", "## W/W_N limits T/Tc -> 0.999 for the 13 parameter sets (variant a)", "",
          "| set | old | new | d% | note |", "|---|---|---|---|---|"]
    for r in out["groups_limits"]:
        old_s = "-" if r["old"] is None else f"{r['old']:.5f}"
        new_s = "-" if r["new"] is None else f"{r['new']:.5f}"
        p = "" if r["pct"] is None else f"{r['pct']:+.2f}"
        note = "no SC state (flagged)" if r["no_superconductivity"] else ""
        L.append(f"| {r['label']} | {old_s} | {new_s} | {p} | {note} |")
    if out["soc_explicit"]:
        L += ["", "## SOC vertex decomposition (F-14: branch-consistent binning)", "",
              "| V_SO | T/Tc | R_A old | R_A new | R_B old | R_B new | inter-helicity old | new |",
              "|---|---|---|---|---|---|---|---|"]
        for r in out["soc_explicit"]:
            L.append(f"| {r['V_SO']:.4f} | {r['frac']} | {r['RA_old']:.5f} | {r['RA_new']:.5f} | "
                     f"{r['RB_old']:.5f} | {r['RB_new']:.5f} | {r['inter_old']:.4f} | {r['inter_new']:.4f} |")
    if out["param_peaks"]:
        L += ["", "## Parameter-peak study: number of refined peak records", "",
              "| file | old | new |", "|---|---|---|"]
        for k, v in out["param_peaks"].items():
            L.append(f"| {k} | {v['old']} | {v['new']} |")
    with open(os.path.join(ROOT, "CORRECTION_TABLES.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote results/BEFORE_AFTER.json and CORRECTION_TABLES.md")


if __name__ == "__main__":
    main()
