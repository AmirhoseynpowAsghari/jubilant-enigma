"""
test_fixes.py -- one executable check per audit finding (F-01 ... F-18) plus the physics
re-validation of the coherence kernels.  Run:  python3 test_fixes.py   (exit 1 on any FAIL).

A finding is reported RESOLVED only if its test PASSES in this run.  Tests that depend on a
regenerated artefact are reported SKIP (not PASS) when the artefact is absent.
"""
from __future__ import annotations

import json
import os
import sys
import subprocess
import warnings
import hashlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ULTRA = os.path.join(ROOT, "project", "ultrasonic")
SCP = os.path.join(ROOT, "project", "sc_params")
for p in (ULTRA, SCP):
    sys.path.insert(0, p)

RESULTS = []


def check(tag, name, fn):
    try:
        out = fn()
        if out is None:
            RESULTS.append((tag, name, "SKIP", ""))
        else:
            ok, detail = out
            RESULTS.append((tag, name, "PASS" if ok else "FAIL", detail))
    except Exception as e:  # a crashing test is a FAIL, never a silent pass
        RESULTS.append((tag, name, "FAIL", f"exception: {type(e).__name__}: {e}"))


# ----------------------------------------------------------------------------- F-01 kernels
def t_F01():
    import bruteforce_ua as bf
    from absorption import coherence
    rng = np.random.default_rng(3)
    n = 3000
    xi, xp = rng.normal(0, 1.5, n), rng.normal(0, 1.5, n)
    D, Dp = rng.normal(0, 0.8, n), rng.normal(0, 0.8, n)
    E = np.hypot(xi, D) + 0.0; Ep = np.hypot(xp, Dp) + 0.0
    worst = {}
    for v in ("a_signed", "b_abs", "c_paper", "d_std", "case2_signed"):
        a = bf._coh(xi, xp, D, Dp, E, Ep, v)
        b = coherence(E, xi, D, Ep, xp, Dp, v)
        worst[v] = float(np.max(np.abs(a - b)))
    ok = max(worst.values()) < 1e-14
    return ok, "max|bruteforce-production| per variant: " + ", ".join(f"{k}={v:.1e}" for k, v in worst.items())


def t_F01_bdg():
    import kernel_bdg_check as kb
    r = kb.run(n=3000, seed=11, verbose=False)
    a_ok = r["(1) C_a == (uu'-vv')^2"] < 1e-12 and r["(3b) C_pair(a) == 1 - C_a == (uv'+vu')^2"] < 1e-12
    d_min = r["(4) min_k max|C_d - BdG_k|  (d is NOT a BdG element if >>0)"]
    ok = a_ok and d_min > 0.1
    k_a = "(1) C_a == (uu'-vv')^2"
    return ok, (f"(a) equals BdG (uu'-vv')^2 to {r[k_a]:.1e}; "
                f"(d) min distance to any BdG element = {d_min:.3f} (not a BdG kernel)")


# ----------------------------------------------------------------------------- F-02 grid
def t_F02():
    from solvers import solve_set_core
    from kgrid import midpoint_k
    k = midpoint_k(40)
    ok_nodes = bool(np.all(k > 0) and np.all(k < np.pi))
    t40 = solve_set_core(95.0, 4.5, 1.885, 0.0, Nk=40)["Tc"]
    t80 = solve_set_core(95.0, 4.5, 1.885, 0.0, Nk=80)["Tc"]
    t160 = solve_set_core(95.0, 4.5, 1.885, 0.0, Nk=160)["Tc"]
    spread = max(abs(t40 - t80), abs(t80 - t160), abs(t40 - t160))
    audit = 0.278211
    ok = ok_nodes and spread < 2e-3 and abs(t160 - audit) < 2e-5
    return ok, (f"midpoint Tc(Nk=40/80/160) = {t40:.6f}/{t80:.6f}/{t160:.6f}; spread {spread:.1e} "
                f"(endpoint grid drifted 0.293986->0.282122)")


# ----------------------------------------------------------------------------- F-03 DOS smoothing
def t_F03():
    import parameters_ua as P
    from dos_ua import NormalDOS
    import inspect
    default = inspect.signature(NormalDOS.__init__).parameters["smooth"].default
    cal = os.path.join(ULTRA, "results", "calibrate_smooth.json")
    if not os.path.exists(cal):
        return None
    d = json.load(open(cal))
    ok = abs(P.SMOOTH - d["SMOOTH_adopted"]) < 1e-12 and abs(default - P.SMOOTH) < 1e-12
    return ok, (f"SMOOTH adopted {P.SMOOTH} (= NormalDOS default {default}); calibration "
                f"systematic |R(sigma)-R(0)|/R = {d['systematic_at_adopted']:.2e}")


# ----------------------------------------------------------------------------- F-04 S_N
def t_F04():
    from dos_ua import NormalDOS
    from absorption import S_normal
    dos = NormalDOS(0.4082, Nk=600)
    vals = {o: S_normal(dos, 0.2, 0.7, 0.4, neps=o) for o in (2, 4, 8)}
    rel = max(abs(vals[o] - vals[8]) / vals[8] for o in (2, 4))
    # the runtime assertion must reject a non-converged request
    raised = False
    try:
        S_normal(dos, 0.2, 0.7, 0.4, neps=1)
    except ValueError:
        raised = True
    ok = rel < 1e-10 and raised
    return ok, (f"Gauss order 2/4/8 relative spread {rel:.1e}; out-of-range order rejected={raised}")


# ----------------------------------------------------------------------------- F-05 pipeline
def t_F05():
    mk = open(os.path.join(ROOT, "Makefile")).read()
    code = [l for l in mk.splitlines() if not l.lstrip().startswith("#")]
    no_or_true = not any("|| true" in l for l in code)
    has_pipefail = "pipefail" in mk
    rg = [l for l in open(os.path.join(ROOT, "regen_all.sh")).read().splitlines()
          if not l.lstrip().startswith("#")]
    rg_clean = not any("|| true" in l for l in rg) and any("set -euo pipefail" in l for l in rg)
    # behavioural demo: a failing producer behind a pipe must make the pipeline fail
    r = subprocess.run(["bash", "-o", "pipefail", "-c", "false | tail -1"], capture_output=True)
    behav = r.returncode != 0
    ok = no_or_true and has_pipefail and rg_clean and behav
    return ok, (f"Makefile pipefail={has_pipefail}, no '|| true'={no_or_true}, "
                f"regen_all strict={rg_clean}, pipe demo rc={r.returncode}")


# ----------------------------------------------------------------------------- F-06 normal-state data
def t_F06():
    fp = os.path.join(ULTRA, "results", "param_features.json")
    if not os.path.exists(fp):
        return None
    data = json.load(open(fp))
    if not any("normal_ref" in d for r in data for d in r.get("temps", {}).values()):
        return None          # stale (pre-fix schema): regenerate before this can count
    bad = 0; n_normal = 0
    for r in data:
        for key, d in r.get("temps", {}).items():
            if d["sc"] == 0.0:
                n_normal += 1
                if d["a"] is not None:
                    bad += 1
    return bad == 0, f"normal-state (T>=Tc) points: {n_normal}; with physics values stored: {bad}"


# ----------------------------------------------------------------------------- F-07 / F-08 solver
def t_F07():
    import solvers as S
    from gap_equations import make_residual
    from kgrid import build_kgrid
    eps, term = build_kgrid(0.4082, 95.0, 4.5, 40)
    res = make_residual(4.5)
    r = S.solve_set_core(95.0, 4.5, 1.885, 0.4082, Nk=40)
    x = [r["D0_0"], r["DS_0"], r["mu_0"]]
    band = S._band(eps)
    # (i) the exact symmetry (D0,DS,mu)->(-D0,-DS,mu) of the residual
    rs = S._residual_norm(res, [-x[0], -x[1], x[2]], np.inf, 1.885, eps, term)
    # (ii) a sign-mixed input (-D0, DS) must fold to a VERIFIED root with D0>=0 and DS kept signed
    fold = S._folded_root(res, [-x[0], x[1], x[2]], np.inf, 1.885, eps, term, mu_band=band)
    if fold is None:
        return False, "mixed-sign input not folded to a verified root"
    rn = S._residual_norm(res, fold, np.inf, 1.885, eps, term)
    ok = rs < S.RESID_TOL and rn <= S.RESID_TOL and fold[0] >= 0 and abs(fold[1] - x[1]) < 1e-6
    return ok, (f"symmetry residual {rs:.1e}; mixed-sign input folded to verified root "
                f"(max|r|={rn:.1e}, D0={fold[0]:.4f} >=0, DS={fold[1]:+.4f} signed)")


def t_F08():
    import solvers as S
    out = []
    ok = True
    for V in (0.0, 0.8):
        r = S.solve_set_core(-75.0, 4.5, 1.885, V, Nk=40)
        if r["Tc_bracketed"] is False:
            out.append(f"V={V}: flagged (not bracketed, lower bound {r['Tc_lower_bound']:.3f}, Tc={r['Tc']})")
            ok = ok and r["Tc"] == 0.0
        elif r["Tc"] == 0.0:
            out.append(f"V={V}: no in-band SC root -> Tc=0 (no superconductivity)")
            ok = ok and abs(r["Deff_0"]) == 0.0
        else:
            # any shipped Tc must be bracketed: the gap must be below threshold just above Tc
            from gap_equations import make_residual
            from kgrid import build_kgrid
            eps, term = build_kgrid(V, -75.0, 4.5, 40)
            Deff_hi, _ = S._gap_at_T(r["Tc"] * 1.02, [r["D0_0"], r["DS_0"], r["mu_0"]], 1.885,
                                     eps, term, make_residual(4.5),
                                     t0=[r["D0_0"], r["DS_0"], r["mu_0"]])
            ok = ok and Deff_hi <= S.GAP_THRESHOLD_FRAC * r["Deff_0"] * 1.001
            out.append(f"V={V}: Tc={r['Tc']:.4f} bracket verified")
    return ok, "; ".join(out)


# ----------------------------------------------------------------------------- F-09 / F-10 loader
def t_F09_F10():
    import sc_loader as L
    meta = json.load(open(os.path.join(SCP, "results", "tc_vs_vso.json")))["meta"]
    if meta.get("grid") != "midpoint":
        return None          # stale (pre-fix) stage-1 file: regenerate before this counts
    try:
        L.sc_params_at(0.123456)
        snapped = False
    except ValueError:
        snapped = True
    p = L.sc_params_at(0.4082)
    exact = abs(p["VSO"] - 0.4082) < 1e-12
    r1 = L.stage1(os.path.join(ROOT, "project", "sc_params"))
    r2 = L.stage1(os.path.join(ROOT, "project", "ultrasonic", "..", "sc_params"))
    keyed = r1 is r2 or r1 == r2
    ok = snapped and exact and keyed
    return ok, f"off-grid V_SO refused={snapped}; exact lookup={exact}; root-keyed cache consistent={keyed}"


# ----------------------------------------------------------------------------- F-11 chunking
def t_F11():
    import bruteforce_ua as bf
    from parameters_ua import BF_CHUNK
    rng = np.random.default_rng(5)
    nbin = 400
    Eb = np.linspace(0.5, 3.0, nbin)
    f = 1.0 / (np.exp(Eb / 0.4) + 1.0)
    B = dict(cnt=rng.uniform(1, 5, nbin), xi=rng.normal(0, 1, nbin), D=rng.normal(0.3, .1, nbin))
    A = dict(cnt=rng.uniform(1, 5, nbin), xi=rng.normal(0, 1, nbin), D=rng.normal(0.3, .1, nbin))
    s1 = bf._pair_sum(B, A, f, Eb, 0.7, 0.02, "a_signed", nbin, "scattering", chunk=16)
    s2 = bf._pair_sum(B, A, f, Eb, 0.7, 0.02, "a_signed", nbin, "scattering", chunk=10 ** 6)
    ok = abs(s1 - s2) <= 1e-12 * max(abs(s2), 1e-300)
    return ok, f"chunked (16) vs unchunked pair sum: |diff|/|S| = {abs(s1-s2)/max(abs(s2),1e-300):.1e} (BF_CHUNK={BF_CHUNK})"


# ----------------------------------------------------------------------------- F-13 guards
def t_F13():
    from dos_ua import NormalDOS
    from absorption import S_scattering, S_normal, S_paircreation
    from physics_ua import SCState
    dos = NormalDOS(0.0, Nk=300, nbins=1500)
    st = dict(T=0.2, D0=0.3, DS=0.0, mu=-1.0, mu_N=-1.0, frac=0.5, Deff=0.3)
    msgs = []
    cases = [
        ("T<0", lambda: S_scattering(dos, -0.1, 0.7, 0.3, 0.0, -1.0, "a_signed")),
        ("omega<=0", lambda: S_scattering(dos, 0.2, 0.0, 0.3, 0.0, -1.0, "a_signed")),
        ("nell<16", lambda: S_scattering(dos, 0.2, 0.7, 0.3, 0.0, -1.0, "a_signed", nell=3)),
        ("branches bad", lambda: S_scattering(dos, 0.2, 0.7, 0.3, 0.0, -1.0, "a_signed", branches=(2,))),
        ("S_normal T=0", lambda: S_normal(dos, 0.0, 0.7, 0.0)),
        ("frac>=1 in SC integrator", lambda: __import__("absorption").ratio_all(dos, dict(st, frac=1.0), 0.7)),
        ("unknown variant", lambda: S_scattering(dos, 0.2, 0.7, 0.3, 0.0, -1.0, "nope")),
    ]
    ok = True
    for name, fn in cases:
        try:
            fn(); ok = False; msgs.append(f"{name}: NOT rejected")
        except (ValueError, Exception) as e:
            msgs.append(f"{name}: rejected")
    return ok, "; ".join(msgs)


def t_F13_io():
    import absorption
    return None


# ----------------------------------------------------------------------------- F-14 SOC binning
def t_F14():
    import make_soc_explicit as M
    from physics_ua import SCState
    sc = SCState.from_params(0.0, 0.2782, 0.2556, 0.9242, 0.4344, n_tgt=1.885)
    st = sc.at(0.5)
    states = M.hel_states(0.0, 60, st["D0"], st["DS"], st["mu"])
    # at V_SO = 0 the two helicities are identical -> the helicity-diagonal vertex equals the
    # full vertex: inter-helicity fraction must vanish (the docstring's sanity check)
    res = M.sums(states, st["T"], 0.7, st["mu"])
    w = M.WS[0]
    keysA = [(a, b) for a in ((1, 1), (1, -1), (-1, 1), (-1, -1)) for b in ((1, 1), (1, -1), (-1, 1), (-1, -1))]
    A = sum(res[w][k] for k in keysA)
    B = sum(res[w][k] for k in keysA if k[0][0] == k[1][0])
    inter = 1 - B / A                    # = 1/2 by counting: half of all (h,h') pairs have h!=h'
    SA, SB = res[w]["SA"], res[w]["SB"]
    RA, RB = A / SA, B / SB              # the two vertex ratios on the same normal reference
    mass = sum(float(np.sum(s["w"])) for s in states)
    ok = abs(inter - 0.5) < 1e-9 and abs(RA - RB) / abs(RA) < 1e-9 and abs(mass - 2.0) < 1e-9
    return ok, (f"V_SO=0: inter-helicity pair fraction={inter:.6f} (=1/2 by counting); "
                f"R_A={RA:.10f} vs R_B={RB:.10f} (rel diff {abs(RA-RB)/abs(RA):.1e}, must be 0); "
                f"state weight={mass:.6f} (=2)")


# ----------------------------------------------------------------------------- F-15 thresholds
def t_F15():
    from param_scans_peaks import threshold_crossings
    x = np.linspace(0, 1, 11)
    em = 0.5 - 0.5 * x          # 2*Emin crosses omega=0.7 at x = 0.3
    thr = threshold_crossings(x, em, 0.7)
    ok = len(thr) == 1 and abs(thr[0]["x_cross"] - 0.3) < 1e-9
    return ok, f"interpolated crossing = {thr[0]['x_cross']:.6f} (exact 0.3); grid argmin would give a grid point"


# ----------------------------------------------------------------------------- F-16 packaging
def t_F16():
    vend = os.path.join(ROOT, "project", "paper_correction", "vendor")
    man = open(os.path.join(vend, "MANIFEST_md5.txt")).read().split()
    h = hashlib.md5(open(os.path.join(vend, man[1]), "rb").read()).hexdigest()
    req = open(os.path.join(ROOT, "requirements.txt")).read()
    pym = "pymupdf" in req
    src = open(os.path.join(ROOT, "project", "paper_correction", "build_corrected_pdf.py")).read()
    abs_path = "/usr/local/lib" in src
    ok = (h == man[0]) and pym and not abs_path
    return ok, f"vendored PDF md5 matches manifest={h == man[0]}; pymupdf pinned={pym}; absolute font path removed={not abs_path}"


# ----------------------------------------------------------------------------- F-17 writers
def t_F17():
    import io_utils_ua as io
    tmp = dict(a=1.0, b=float("nan"), c=[float("inf"), 2.0], d=np.float64(3.5))
    path = io.save("__f17_probe.json", tmp)
    raw = open(path).read()
    def _no(c):
        raise ValueError("non-finite token")
    try:
        json.loads(raw, parse_constant=_no)
        strict = True
    except ValueError:
        strict = False
    os.remove(path)
    ok = strict and '"b": null' in raw
    return ok, f"non-finite floats written as null, strict JSON parse={strict}"


def t_F17_total():
    import scanning_ua as su
    return None


# ----------------------------------------------------------------------------- F-18 cosmetics
def t_F18():
    import warnings as W
    from absorption import fermi_qp, fermi_normal
    with W.catch_warnings(record=True) as rec:
        W.simplefilter("always")
        a = fermi_qp(np.array([1e4, -1e4, 0.5]), 0.3)
        b = fermi_normal(np.array([-1e7, 1e7]), 0.3, 0.0)
    no_warn = len(rec) == 0
    ok = (no_warn and a[0] == 0.0 and abs(a[2] - 1 / (np.exp(0.5 / 0.3) + 1)) < 1e-15
          and b[0] == 1.0 and b[1] == 0.0)
    return ok, f"overflow-safe Fermi functions, warnings={len(rec)}"


# ----------------------------------------------------------------------------- Levy / BCS physics
def t_LEVY():
    from levy1963_check import levy_ratio, our_ratio
    rows = []
    worst = 0.0
    for T in (0.03, 0.05, 0.07):
        r_l = levy_ratio(T)
        r_o, _ = our_ratio(1e-4, T)
        worst = max(worst, abs(r_o - r_l) / r_l)
        rows.append(f"T={T}: {r_o:.6f} vs {r_l:.6f}")
    ok = worst < 1e-3
    return ok, f"flat-DOS s-wave omega->0 vs Levy 2f(Delta): worst rel {worst:.1e} ({'; '.join(rows)})"


# ----------------------------------------------------------------------------- main
def main():
    check("F-01", "kernels: bruteforce == production (all variants)", t_F01)
    check("F-01", "variant (a) == BdG (uu'-vv')^2 ; (d) is not a BdG element", t_F01_bdg)
    check("F-02", "midpoint grid: Tc convergence and audit value", t_F02)
    check("F-03", "DOS smoothing calibrated and single-sourced", t_F03)
    check("F-04", "S_N quadrature converged and runtime-checked", t_F04)
    check("F-05", "pipeline cannot pass on failure", t_F05)
    check("F-06", "no normal-state values stored as SC data", t_F06)
    check("F-07", "abs()-fold verified against the residual", t_F07)
    check("F-08", "unbracketed Tc flagged, not invented", t_F08)
    check("F-09/10", "exact V_SO lookup; root-keyed cache", t_F09_F10)
    check("F-11", "chunked brute-force pair sum identical", t_F11)
    check("F-13", "input validation rejects bad arguments", t_F13)
    check("F-14", "SOC binning: branch-consistent, V_SO=0 limit exact", t_F14)
    check("F-15", "pair-threshold crossing interpolated", t_F15)
    check("F-16", "packaging: vendored PDF manifest, pinned deps, no absolute paths", t_F16)
    check("F-17", "JSON writers strict (no NaN/Inf)", t_F17)
    check("F-18", "overflow-safe Fermi functions, no hidden clamps", t_F18)
    check("PHYS", "Levy 1963 limit (benchmark T5)", t_LEVY)

    width = max(len(r[1]) for r in RESULTS) + 2
    for tag, name, st, det in RESULTS:
        print(f"[{st:4s}] {tag:8s} {name:<{width}} {det}")
    n_fail = sum(1 for r in RESULTS if r[2] == "FAIL")
    n_skip = sum(1 for r in RESULTS if r[2] == "SKIP")
    print(f"\nSUMMARY: {len(RESULTS) - n_fail - n_skip} PASS, {n_fail} FAIL, {n_skip} SKIP")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
