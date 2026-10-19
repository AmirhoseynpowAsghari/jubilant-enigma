"""
dualspace_check.py -- T7: dual-formulation cross-validation (k-space vs energy space),
mirroring the benchmarking strategy of

    A. Sardwal, A. Kreisel and S. Mukherjee,
    "Ultrasonic attenuation in inhomogeneous superconductors",
    arXiv:2608.19079 (2026), Sec. III and Fig. 3,

who validate their real-space BdG attenuation formula against (i) a fully momentum-space
formulation and (ii) the analytic Walker-Smith-Samokhin treatment [PRB 65, 014517 (2001)],
concluding "good agreement between the various approaches".

Sardwal et al. start from the phonon-number equation of motion + linear response and obtain
(their Eqs. 15-17)

    alpha_s(q,lambda) = -(l_q^2 / 2 hbar v_s N) SUM_{m,n} |F_mn(q)|^2 Im chi_mn(omega),
    chi_mn(omega)     = [f(E_m) - f(E_n)] / (hbar*omega + E_m - E_n + i eta),
    F_mn(q)           = SUM_ij sigma  dt_ij phi*_ij(q) (u^m*_i u^n_j - v^n_i v^m*_j).

For a homogeneous superconductor the coherence factor of Eq. (16) reduces to

    (u u' - v v')^2 = 1/2 [ 1 + (xi xi' - Delta Delta') / (E E') ]  =  our  C_I  (signed, eta=+1),

and their PP-channel processes III/IV (Appendix A, Eqs. A6-A9) carry the weight
(1 - f_m - f_n) -- our pair-creation channel S_pair with C_pair = 1 - C_I.  Their statement
that the pair channel is "negligible for homogeneous superconductors ... across the
transition" is the published counterpart of our S_pair = 0 inside the 2*Delta > hbar*omega
window.  So the same golden-rule structure, the same signed coherence factor and the same
channel decomposition appear in an independent 2026 calculation.

Our triangle (this script):
  (i)   k-space  : explicit double lattice sum over ALL (k,helicity) x (k',helicity') pairs
                   with a Gaussian delta and Richardson extrapolation in its width
                   (bruteforce_ua.py -- no analytic inversion, no Jacobian).
  (ii)  E-space  : analytic branch inversion eps_sigma(E) + Jacobian + midpoint quadrature
                   (absorption.py S_scattering -- the production pipeline).
  (iii) analytic : Levy, PR 131, 1497 (1963), Eq. (8):  alpha_s/alpha_n = 2 f(Delta(T))
                   (s-wave, flat DOS, omega -> 0) -- leg (ii) reproduces it to 3e-6 (T5).

T7a  HHM production model (V_SO = 0.4082 and 0, hbar*omega = 0.7 t):
     (i) vs (ii) as a function of the k-grid size Nk -- the user-visible statement
     "solved in both spaces, same results", with the convergence made explicit.
T7b  s-wave triangle at omega -> 0 with a low-density chemical potential mu = -3.6 t
     (smooth band DOS far from the van-Hove point; Sardwal et al. likewise keep the DOS
     smooth at mu by choosing t' = 0.3, n = 0.8):
     (i) vs (ii) [same band DOS: pure implementation cross-check]
     (ii) vs (iii) [band DOS vs flat DOS: the physical band-structure correction to Levy]
     (ii-flat) vs (iii) [T5 identity, reproduced here as the fourth column].

Run:  python dualspace_check.py      (~5 min)
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dos_ua import NormalDOS
from physics_ua import SCState
from parameters_ua import OMEGA, T_HOP
from absorption import S_scattering, S_normal, qp_constants
import bruteforce_ua as bf
from levy1963_check import gap_BCS, levy_ratio, our_ratio, TC

MU_TRI = -3.6 * T_HOP        # low-density mu for the s-wave triangle (smooth DOS region)
OMEGA_TRI = 1e-4             # omega -> 0 for the Levy leg


# --------------------------------------------------------------------------- helpers
def _rich(vals, ws):
    """Richardson extrapolation w -> 0 of S(w) = S0 + c w^2."""
    ws = np.asarray(ws, float)
    A = np.vstack([np.ones_like(ws), ws ** 2]).T
    coef, *_ = np.linalg.lstsq(A, np.asarray(vals, float), rcond=None)
    return float(coef[0])


def lattice_ratio(VSO, Nk, T, omega, D0, DS, mu, SN, nbin=3000,
                  ws=(0.08, 0.05, 0.03), variant="a_signed", Eg=None, diag=False):
    """Leg (i): k-space lattice sum, Richardson-extrapolated, divided by the SAME S_N."""
    lst = bf.lattice_states(VSO, Nk, D0, DS, mu)
    if Eg is None:
        P = qp_constants(D0, DS, mu)
        Eg = (P["Emin"] * 0.995, float(lst["E"].max()) * 1.002)
    ws = np.asarray(ws, float)
    vals = np.array([bf.S_sc_bf(lst, T, omega, variant, Eg, nbin, w) for w in ws])
    r = _rich(vals, ws) / SN
    if diag:
        # consistency of the w^2 extrapolation: all widths vs the smallest widths
        if len(ws) >= 3:
            r2 = _rich(vals[1:], ws[1:]) / SN
            return r, abs(r - r2) / max(abs(r), 1e-300)
        return r, float("nan")
    return r


def tri_window(T, D):
    """Thermal window of the lower-state energy for the omega->0 triangle leg.

    At hbar*omega -> 0 the thermal factor f(E)-f(E+omega) = omega*(-f'(E)) confines the
    integrand to E in [Delta, Delta + ~12 T] (weight beyond: e^-12).  Windowing is what
    makes the small Richardson widths affordable (de << w on the FULL range would need
    nbin ~ 1e5 -> a 1e10 dense matrix).
    """
    lo = 0.995 * D if D > 0 else 1e-6
    hi = (D + min(12.0 * T, 0.6)) * 1.001 if D > 0 else max(0.6, 12.0 * T)
    return (lo, hi)


def tri_ws(T, window):
    """Adaptive Gaussian widths for the omega->0 leg.

    The w^2 Richardson scheme is only valid while the smearing is small on the scale on
    which the integrand varies, i.e. w << T (the thermal window) -- with the fixed
    production widths w = 0.03..0.08 the T7b leg FAILED at low T (dev ~ 16 % at
    T = 0.03 t, tracking w/T), which is exactly the diagnostic.  Fix: w1 = 0.2 T,
    three geometrically spaced widths, and de = w1/16 so that even the smallest width
    covers >= 4 bins (drop to two widths if the nbin cap forces a coarser de).
    """
    w1 = min(0.06, 0.2 * T)
    ws = np.array([w1, 0.5 * w1, 0.25 * w1])
    de_target = w1 / 16.0
    nbin = int(min(max(np.ceil((window[1] - window[0]) / de_target), 800), 4000))
    de = (window[1] - window[0]) / nbin
    keep = ws / de >= 3.5
    ws = ws[keep]
    return ws, nbin


def energy_ratio(dos, T, omega, D0, DS, mu, SN, nell=40000, variant="a_signed"):
    """Leg (ii): production energy-space pipeline, divided by the SAME S_N."""
    S, _ = S_scattering(dos, T, omega, D0, DS, mu, variant, nell=nell)
    return S / SN


# --------------------------------------------------------------------------- T7a
def T7a(VSOs=(0.4082, 0.0), fracs=(0.5, 0.9, 0.999), Nks=(160, 240, 320, 452)):
    print("=" * 96)
    print("T7a  k-space lattice sum vs energy-space pipeline, HHM production model,")
    print(f"     hbar*omega = {OMEGA} t, variant a_signed; Richardson in the Gaussian delta width")
    print("=" * 96)
    out = {}
    for VSO in VSOs:
        dos = NormalDOS(VSO, Nk=1200, nbins=6000, smooth=0.02)
        sc = SCState(VSO, Nk_mu=200)
        print(f"  V_SO = {VSO}   (Tc = {sc.Tc:.4f} t)")
        print(f"    {'T/Tc':>6} " + "".join(f"{'Nk=' + str(n):>34}" for n in Nks))
        print(f"    {'':>6} " + "".join(f"{'E-space':>11}{'k-space':>11}{'dev':>12}" for n in Nks))
        for frac in fracs:
            st = sc.at(frac)
            SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
            r_es = energy_ratio(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], SN)
            row = f"    {frac:6.3f} "
            for Nk in Nks:
                r_bf = lattice_ratio(VSO, Nk, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], SN)
                dev = abs(r_es - r_bf) / max(abs(r_es), 1e-12)
                out[(VSO, frac, Nk)] = (r_es, r_bf, dev)
                row += f"{r_es:11.6f}{r_bf:11.6f}{dev:12.2e}"
            print(row)
    worst = max(v[2] for v in out.values())
    ok = worst < 0.03
    print(f"  worst rel. dev. over the whole table: {worst:.2e}   ->  {'PASS' if ok else 'FAIL'}")
    return out, ok


# --------------------------------------------------------------------------- T7b
def _tri_state(T):
    D = gap_BCS(T)
    dos = NormalDOS(0.0, Nk=1200, nbins=6000, smooth=0.02)
    return D, dos


def T7b(fracs=(0.3, 0.5, 0.7, 0.9), Nk=320):
    print("=" * 96)
    print("T7b  s-wave triangle at hbar*omega -> 0:  (i) k-space lattice  vs  (ii) E-space (band DOS)")
    print(f"     vs  (iii) Levy 2f(Delta)  vs  (ii-flat) E-space (flat DOS).   mu = {MU_TRI} t (low density,")
    print("     smooth DOS), V_SO = 0, Delta^s = 0, Delta(T) from the BCS gap equation, Tc = 0.1 t.")
    print("     mu_N = mu (shift negligible: Delta << distance to the band edge).")
    print("     k-space leg: thermal window [Delta, Delta+min(12T, 0.6)] with ADAPTIVE Gaussian")
    print("     widths w1 = 0.2T (the fixed production widths w ~ 0.03-0.08 >> T fail at omega->0;")
    print("     Rich.cons = consistency of the w^2 extrapolation, all widths vs two smallest.")
    print("=" * 96)
    dos = NormalDOS(0.0, Nk=1200, nbins=6000, smooth=0.02)
    print(f"    {'T/Tc':>6} {'2f(D) Levy':>11} {'E-sp flat':>11} {'E-sp band':>11} {'k-space':>11}"
          f" {'dev(i,ii)':>10} {'dev(ii,iii)':>12} {'Rich.cons':>10}")
    out = {}
    worst_ab = 0.0
    for frac in fracs:
        T = frac * TC
        D = gap_BCS(T)
        SN = S_normal(dos, T, OMEGA_TRI, MU_TRI)
        r_lev = levy_ratio(T)
        r_flat, _ = our_ratio(OMEGA_TRI, T)
        r_es = energy_ratio(dos, T, OMEGA_TRI, D, 0.0, MU_TRI, SN)
        win = tri_window(T, D)
        ws, nb = tri_ws(T, win)
        r_bf, cons = lattice_ratio(0.0, Nk, T, OMEGA_TRI, D, 0.0, MU_TRI, SN,
                                   nbin=nb, ws=ws, Eg=win, diag=True)
        dev_ab = abs(r_es - r_bf) / max(r_es, 1e-12)
        dev_bc = abs(r_es - r_lev) / max(r_lev, 1e-12)
        worst_ab = max(worst_ab, dev_ab)
        out[frac] = (r_lev, r_flat, r_es, r_bf, dev_ab, dev_bc)
        print(f"    {frac:6.2f} {r_lev:11.6f} {r_flat:11.6f} {r_es:11.6f} {r_bf:11.6f}"
              f" {dev_ab:10.2e} {dev_bc:12.2e} {cons:10.2e}")
    ok = worst_ab < 0.03
    print(f"  worst (i)-(ii) rel. dev.: {worst_ab:.2e}  ->  {'PASS' if ok else 'FAIL'}")
    print("  note: dev(ii,iii) is the PHYSICAL band-DOS correction to Levy's free-electron limit")
    print("        (his derivation assumes a flat DOS at the Fermi level); the flat-DOS column")
    print("        reproduces Levy to T5 accuracy, as it must.")
    return out, ok


# --------------------------------------------------------------------------- figure
def make_figure(t7a, path="figures/dual_space_k_vs_energy.png"):
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.3))

    # (a) s-wave triangle sweep
    ax = axes[0]
    fr = np.linspace(0.05, 0.99, 30)
    dos = NormalDOS(0.0, Nk=1200, nbins=6000, smooth=0.02)
    rA, rB, rC = [], [], []
    for f in fr:
        T = f * TC
        D = gap_BCS(T)
        SN = S_normal(dos, T, OMEGA_TRI, MU_TRI)
        rB.append(energy_ratio(dos, T, OMEGA_TRI, D, 0.0, MU_TRI, SN, nell=20000))
        win = tri_window(T, D)
        ws, nb = tri_ws(T, win)
        rA.append(lattice_ratio(0.0, 240, T, OMEGA_TRI, D, 0.0, MU_TRI, SN,
                                nbin=nb, ws=ws, Eg=win))
        rC.append(levy_ratio(T))
    ax.plot(fr, rC, "k--", lw=1.4, label="(iii) Levy  $2f(\\Delta)$  [analytic]")
    ax.plot(fr, rB, "-", color="tab:blue", lw=1.6, label="(ii) E-space (band DOS)")
    ax.plot(fr, rA, "o", color="tab:red", ms=3.6, label="(i) k-space lattice, $N_k=240$")
    ax.set_xlabel("$T/T_c$"); ax.set_ylabel("$\\alpha_s/\\alpha_n$")
    ax.set_title("(a) s-wave limit, $\\hbar\\omega\\to0$, $\\mu=-3.6t$\nthree complementary formulations")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.3)

    # (b) HHM production model: E-space vs k-space at Nk=320
    ax = axes[1]
    VSO = 0.4082
    fracs = (0.3, 0.5, 0.7, 0.9, 0.99, 0.999)
    dos = NormalDOS(VSO, Nk=1200, nbins=6000, smooth=0.02)
    sc = SCState(VSO, Nk_mu=200)
    res, ksp = [], []
    for f in fracs:
        st = sc.at(f)
        SN = S_normal(dos, st["T"], OMEGA, st["mu_N"])
        res.append(energy_ratio(dos, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], SN))
        ksp.append(lattice_ratio(VSO, 320, st["T"], OMEGA, st["D0"], st["DS"], st["mu"], SN))
    ax.plot(fracs, res, "-", color="tab:blue", lw=1.6, label="(ii) E-space")
    ax.plot(fracs, ksp, "o", color="tab:red", ms=5, label="(i) k-space lattice, $N_k=320$")
    ax.set_xlabel("$T/T_c$"); ax.set_ylabel("$W/W_N$  (a_signed)")
    ax.set_title("(b) HHM model, $V_{SO}=0.4082$, $\\hbar\\omega=0.7t$\nboth spaces, same result")
    ax.legend(fontsize=9); ax.grid(alpha=0.3)

    # (c) convergence of the k-space leg with Nk
    ax = axes[2]
    Nks = sorted({k[2] for k in t7a})
    for frac, c, m in ((0.5, "tab:green", "s"), (0.999, "tab:purple", "^")):
        dev = [t7a[(VSO, frac, n)][2] for n in Nks]
        ax.plot(Nks, dev, m + "-", color=c, ms=5, label=f"$T/T_c={frac}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$N_k$  (quarter-zone grid)"); ax.set_ylabel("rel. dev.  (i) vs (ii)")
    ax.set_title("(c) convergence, $V_{SO}=0.4082$, $\\hbar\\omega=0.7t$\nfloor = NormalDOS construction error")
    ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  figure -> {path}")
    return path


if __name__ == "__main__":
    import os, sys
    here = os.path.dirname(os.path.abspath(__file__))
    os.chdir(here)
    t7a, ok_a = T7a()
    t7b, ok_b = T7b()
    make_figure(t7a)
    print("=" * 96)
    print(f"T7  dual-space validation:  T7a {'PASS' if ok_a else 'FAIL'}   T7b {'PASS' if ok_b else 'FAIL'}")
    sys.exit(0 if (ok_a and ok_b) else 1)
