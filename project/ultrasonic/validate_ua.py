"""
validate.py -- every check that the recalculation is trustworthy.  Run:  python3 validate.py

T0  quadrature            : the ell = sqrt(E^2 - Emin^2) substitution used by S_scattering is
                            identical to a high-resolution integral over E with the explicit
                            Jacobian E/(kappa*ell).
T1  Delta -> 0            : closed form.  At Delta = 0 the signed eta = +1 coherence factor
                            becomes a BRANCH PROJECTOR, and the scattering channel keeps only
                            part of the normal-state thermal window:
                                W_sc/W_N |_{Delta=0} = R(hw,T) = (2/hw) int_0^inf [f(E)-f(E+hw)] dE
                            R -> 1 for hw << T, R -> 2 T ln2 / hw for hw >> T.
T2  textbook BCS limit    : flat DOS, CONSTANT gap, mu = 0, hbar*omega -> 0, eta = +1
                                W_sc/W_N -> 2 f(Delta) = 2/(exp(Delta/T)+1)      (exact)
T3  lattice cross-check   : energy space (branch inversion + Jacobian) versus an explicit sum
                            over all (k, helicity) x (k', helicity') pairs on the lattice, with
                            sigma = sign(eps - eps*) as the branch label.
T4  Cauchy-Schwarz bound  : for the SIGNED factor |xi xi' - Dbar Dbar'| <= E E', so 0 <= C <= 1
                            and W/W_N <= 1 necessarily.  Variants (b),(c) violate C <= 1.
"""
from __future__ import annotations

import numpy as np

from dos_ua import NormalDOS
from physics_ua import SCState, solve_mu
from parameters_ua import OMEGA, T_HOP
from absorption import (S_scattering, S_normal, qp_constants, fermi_qp, fermi_normal,
                          coherence, branch, E_of_ell, pair_range, VARIANTS)
import bruteforce_ua as bf


class FlatDOS:
    """N(eps) = 2/(2 Wc) on [-Wc, Wc]; integral = 2.  For T2."""

    def __init__(self, Wc=12.0, nb=12000):
        self.Wband = Wc
        self.grid = np.linspace(-Wc, Wc, nb)
        self.de = self.grid[1] - self.grid[0]
        self.dens = np.full(nb, 2.0 / (2 * Wc))
        self.eps_i, self.w_i = self.grid, self.dens * self.de

    def N(self, eps):
        eps = np.asarray(eps, float)
        return np.where(np.abs(eps) <= self.Wband, np.interp(eps, self.grid, self.dens), 0.0)


def R_delta0(omega, T, npts=2000001):
    E = np.linspace(1e-12, 400.0 * max(T, omega), npts)
    return 2.0 * np.trapezoid(fermi_qp(E, T) - fermi_qp(E + omega, T), E) / omega


# ---------------------------------------------------------------------------------------------
def T0_quadrature():
    print("=" * 92)
    print("T0  the ell-substitution equals the explicit E-integral with the Jacobian E/(kappa*ell)")
    print("=" * 92)
    Delta, T, Wc = 0.6, 0.3, 12.0
    dos = FlatDOS(Wc=Wc, nb=12000)
    worst = 0.0
    for om in (0.05, 0.01, 0.002):
        S, p = S_scattering(dos, T, om, Delta, 0.0, 0.0, "a_signed", nell=400000)
        P = qp_constants(Delta, 0.0, 0.0)
        lo, hi = pair_range(dos, P, om, +1, +1)
        E = np.linspace(lo, hi, 4000001)
        dE = E[1] - E[0]
        l = np.sqrt(np.maximum(E ** 2 - P["Emin"] ** 2, 0.0))
        Ep = E + om
        lp = np.sqrt(np.maximum(Ep ** 2 - P["Emin"] ** 2, 0.0))
        C = coherence(E, l, np.full_like(E, Delta), Ep, lp, np.full_like(E, Delta), "a_signed")
        th = fermi_qp(E, T) - fermi_qp(Ep, T)
        jac = (np.where(l > 0, E / np.maximum(l, 1e-300), 0.0) *
               np.where(lp > 0, Ep / np.maximum(lp, 1e-300), 0.0))
        J = np.sum(jac * C * th * (Ep <= hi + 1e-12)) * dE
        ref = p[(1, 1)] / dos.N(0.0) ** 2
        worst = max(worst, abs(J - ref) / ref)
        print(f"   hbar*omega={om:<7} channel (+1,+1): E-space integral={J:.8e}  "
              f"module={ref:.8e}  rel.dev={abs(J-ref)/ref:.2e}")
    print(f"   --> worst relative deviation {worst:.2e}")
    return worst


def T1_delta_to_zero():
    print("=" * 92)
    print("T1  Delta -> 0.  Closed form R(hw,T) = (2/hw) int_0^inf [f(E)-f(E+hw)] dE")
    print("    (flat DOS).  For the real lattice DOS the same projector acts, but N(eps) is")
    print("    not flat over the thermal window, so the number differs from R.")
    print("=" * 92)
    dos = FlatDOS(Wc=12.0, nb=12000)
    worst = 0.0
    print(f"   {'T':>7} {'hw':>6} | {'(a) numeric':>12} {'R closed':>10} | {'(b)':>10} {'(c)':>10}")
    for T in (0.05, 0.2, 0.5, 1.0):
        for om in (0.7, 0.2):
            SN = S_normal(dos, T, om, 0.0)
            r = {}
            for v, lab in VARIANTS:
                S, _ = S_scattering(dos, T, om, 0.0, 0.0, 0.0, v, nell=200000)
                r[v] = S / SN
            worst = max(worst, abs(r["a_signed"] - R_delta0(om, T)))
            print(f"   {T:7.3f} {om:6.2f} | {r['a_signed']:12.6f} {R_delta0(om,T):10.6f} | "
                  f"{r['b_abs']:10.6f} {r['c_paper']:10.6f}")
    print(f"   --> worst deviation of (a) from the closed form: {worst:.2e}")
    print("   NOTE: (b) and (c) are exactly 2 x (a) at Delta = 0 and EXCEED 1 whenever hw < ~2T.")
    return worst


def T2_bcs_limit():
    print("=" * 92)
    print("T2  textbook BCS limit: flat DOS, constant gap, mu=0, eta=+1, hbar*omega -> 0")
    print("    expected  W/W_N -> 2 f(Delta) = 2/(exp(Delta/T)+1)")
    print("=" * 92)
    worst = 0.0
    print(f"   {'Delta':>6} {'T':>6} | {'hw=0.02':>9} {'hw=0.005':>9} {'hw=0.001':>9} "
          f"{'hw=0.0002':>10} | {'2f(Delta)':>10} {'rel.dev':>9}")
    for Delta in (0.3, 0.6, 1.0):
        for T in (0.15, 0.3, 0.6):
            dos = FlatDOS(Wc=12.0, nb=12000)
            want = 2.0 / (np.exp(Delta / T) + 1.0)
            row = []
            for om in (0.02, 0.005, 0.001, 0.0002):
                S, _ = S_scattering(dos, T, om, Delta, 0.0, 0.0, "a_signed", nell=400000)
                SN = S_normal(dos, T, om, 0.0)
                row.append(S / SN)
            worst = max(worst, abs(row[-1] - want) / want)
            print(f"   {Delta:6.2f} {T:6.2f} | " + "".join(f"{x:9.5f}" for x in row[:3])
                  + f"{row[3]:10.5f} | {want:10.6f} {abs(row[-1]-want)/want:9.2e}")
    print(f"   --> worst relative deviation in the hbar*omega -> 0 limit: {worst:.2e}")
    return worst


def T3_lattice(VSO=0.4082, fracs=(0.5, 0.7, 0.9, 0.99, 0.999), Nk=320, nbin=3000):
    print("=" * 92)
    print(f"T3  energy space vs explicit lattice sum over (k,helicity) x (k',helicity),  V_SO={VSO}")
    print("    sigma = sign(eps - eps*);  delta(E'-E-hw) Gaussian, Richardson in its width")
    print("=" * 92)
    dos = NormalDOS(VSO, Nk=1200)
    sc = SCState(VSO, Nk_mu=200)
    out = {}
    for frac in fracs:
        st = sc.at(frac)
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        lst = bf.lattice_states(VSO, Nk, st["D0"], st["DS"], st["mu"])
        Eg = (P["Emin"] * 0.995, float(lst["E"].max()) * 1.002)
        print(f"   T/Tc={frac:<8} E_min={P['Emin']:.5f}")
        for v, lab in VARIANTS:
            S, _ = S_scattering(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], v, nell=40000)
            r_es = S / SN
            ws = np.array([0.08, 0.05, 0.03])
            vals = np.array([bf.S_sc_bf(lst, st["T"], OMEGA, v, Eg, nbin, w) for w in ws])
            A = np.vstack([np.ones_like(ws), ws ** 2]).T
            coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
            r_bf = coef[0] / SN
            dev = abs(r_es - r_bf) / max(r_es, 1e-12)
            out[(frac, v)] = (r_es, r_bf, dev)
            print(f"       {lab:<18} energy space={r_es:.6f}  lattice={r_bf:.6f}  rel.dev={dev:.2e}")
    return out


def T4_bound():
    print("=" * 92)
    print("T4  (i) 0<=C_a<=1 pointwise on a STRUCTURED grid (gap edge, mixed xi signs);")
    print("    (ii) ratio-level bound at Delta=0: (a),(d) = R <= 1, (b),(c) = 2R (can exceed 1).")
    print("    NOTE: random sampling can never show C>1 because |xi xi'| <= EE' identically;")
    print("    the bound that physics breaks is the RATIO bound, violated by (b),(c) at Delta=0.")
    print("=" * 92)
    D = 0.35
    E = np.concatenate([D + np.logspace(-4, 0.7, 800), np.linspace(D, 6, 800)])
    worst = 0.0
    for sa in (+1, -1):
        for sb in (+1, -1):
            xi = sa * np.sqrt(np.maximum(E ** 2 - D ** 2, 0.0))
            xip = sb * np.sqrt(np.maximum((E + 0.2) ** 2 - D ** 2, 0.0))
            C = coherence(E, xi, np.full_like(E, D), E + 0.2, xip, np.full_like(E + 0.2, D),
                          "a_signed")
            worst = max(worst, float(max(np.max(C) - 1.0, -np.min(C))))
    print(f"   (i) signed C on structured grid: bound residual {worst:.2e}  ->  "
          f"{'ok' if worst < 1e-12 else 'VIOLATION'}")
    dos = FlatDOS(Wc=12.0, nb=12000)
    for v in ("a_signed", "b_abs", "c_paper", "d_std"):
        S, _ = S_scattering(dos, 1.0, 0.2, 0.0, 0.0, 0.0, v, nell=40000)
        SN = S_normal(dos, 1.0, 0.2, 0.0)
        print(f"   (ii) Delta=0, hw=0.2, T=1: {v:9s} W/W_N = {S/SN:.6f}"
              + ("   <= 1 ok" if S / SN <= 1 + 1e-9 else "   > 1  (bound broken by this variant)"))
    return worst


if __name__ == "__main__":
    import sys
    t0 = T0_quadrature()
    print()
    t1 = T1_delta_to_zero()
    print()
    t2 = T2_bcs_limit()
    print()
    t3 = T3_lattice()
    print()
    T4_bound()
    print()
    print("=" * 92)
    byv = {}
    for k, v in t3.items():
        byv.setdefault(k[1] if isinstance(k, tuple) else k, []).append(v[2])
    rng = "  ".join(f"{kk}:{min(vv):.1e}-{max(vv):.1e}" for kk, vv in sorted(byv.items()))
    print(f"SUMMARY  T0 {t0:.1e}   T1 {t1:.1e}   T2 {t2:.1e}   T3 per variant {rng}")
    print("=" * 92)
    # ------------------------------------------------------------------ pass/fail (F-05)
    # Thresholds fixed BEFORE the recomputation; a failing check makes the script exit 1 so
    # that `make test` (pipefail) cannot report success on a broken validation.
    TOL = dict(T0=1e-3, T1=1e-4, T2=1e-4, T3_a=3e-2)
    t3a = max(v[2] for k, v in t3.items() if (k[1] if isinstance(k, tuple) else k) == "a_signed")
    checks = dict(T0=t0 < TOL["T0"], T1=t1 < TOL["T1"], T2=t2 < TOL["T2"], T3_a=t3a < TOL["T3_a"])
    print("VALIDATE_UA " + " ".join(f"{k}={'PASS' if ok else 'FAIL'}" for k, ok in checks.items())
          + "   (T3 is reported for d_std but not gated: that variant is not a BdG matrix element,"
          + " see kernel_bdg_check.py)")
    sys.exit(0 if all(checks.values()) else 1)
