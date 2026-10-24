"""
addendum_numbers.py -- every number printed in the corrected addendum, read from the regenerated
results (NOTHING typed by hand).  Used by build_corrected_pdf.py and build_fa_pdf.py.

Route A (standing constraint): the pair channel is an OUTLOOK only.  No pair-channel number is
produced here; the builders must not print any.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
S1 = os.path.join(ROOT, "project", "sc_params", "results")
S2 = os.path.join(ROOT, "project", "ultrasonic", "results")
BA = os.path.join(ROOT, "results", "BEFORE_AFTER.json")


def _j(p):
    with open(p) as fh:
        return json.load(fh)


def _nearest(xs, v):
    return min(range(len(xs)), key=lambda i: abs(xs[i] - v))


def numbers():
    n = {}
    # ---- stage 1 groups
    g1 = _j(os.path.join(S1, "groups_stage1.json"))
    gw = _j(os.path.join(S2, "groups_WWN.json"))
    rows = []
    for a, b in zip(g1, gw):
        s = a["set"]
        lim = b["a_ratio"][-1] if b.get("a_ratio") and b["a_ratio"][-1] is not None else None
        mn = min([x for x in (b.get("a_ratio") or []) if x is not None], default=None)
        rows.append(dict(group=s["group"], label=s["label"], Tc=a["Tc"], Deff=a["Deff_0"],
                         lim=lim, min=mn, nosc=bool(lim is None),
                         mu_window=(a.get("mu_window") or a.get("diagnostic_unrestricted") or {})))
    n["groups"] = rows
    # ---- master V_SO sweeps
    mv = []
    m = _j(os.path.join(S1, "tc_vs_vso.json"))["results"]
    for v in ("0.0000", "0.4082", "0.8000"):
        t = _j(os.path.join(S2, f"temp_sweep_VSO_{v}.json"))
        vv = float(v)
        i = _nearest(m["V_SO"], vv)
        ia = _nearest(t["fracs"], 0.999)
        i7 = _nearest(t["fracs"], 0.7)
        mv.append(dict(V_SO=vv, Tc=m["Tc"][i], Deff=m["Deff_T0"][i], D0=m["D0_T0"][i],
                       DS=m["DS_T0"][i], lim=t["ratios"]["a_signed"][ia],
                       R07=t["ratios"]["a_signed"][i7],
                       d=t["ratios"]["d_std"][ia], c=t["ratios"]["c_paper"][ia],
                       b=t["ratios"]["b_abs"][ia]))
    n["master"] = mv
    # ---- V_SO scans at T/Tc = 0.9 (U = 95, 125)
    vs = _j(os.path.join(S2, "vso_scan.json"))
    scans = {}
    for k, rec in vs.items():
        vals = [x for x in rec["ratio"]["0.9"] if x is not None]
        scans[k] = dict(lo=min(vals) if vals else None, hi=max(vals) if vals else None,
                        n_sc=len(vals), n=len(rec["ratio"]["0.9"]))
    n["scans"] = scans
    # ---- before/after
    b = _j(BA)
    n["ba_master"] = b["production_master"]
    n["ba_stage1"] = b["stage1_master"]
    # ---- families / window sensitivity
    fam = _j(os.path.join(S2, "uvso_families.json"))
    n["fam_total"] = len(fam)
    n["fam_sc"] = sum(1 for r in fam if r["Tc"] > 0)
    return n


def fmt(x, p=4):
    return "n/a" if x is None else f"{x:.{p}f}"


if __name__ == "__main__":
    import pprint
    pprint.pprint({k: v for k, v in numbers().items() if k in ("master", "scans", "fam_total", "fam_sc")})
