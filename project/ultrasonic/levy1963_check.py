"""
levy1963_check.py -- reduction of our model to plain s-wave BCS and comparison with
  M. Levy, "Ultrasonic Attenuation in Superconductors for ql < 1",
  Phys. Rev. 131, 1497 (1963), Eq. (8):   alpha_s/alpha_n = 2 f_0(eps_0(T)).

Levy's derivation (Boltzmann equation with collision drag, ql < 1, omega*tau -> 0) and the
golden-rule derivation used in this project must coincide in the limit
    (i)   s-wave:      Delta^s = 0  (b = 0 -> kappa = 1, eps* = mu, E_min = Delta),
    (ii)  V_SO = 0,    flat DOS N(eps) = const, mu = 0  (free-electron limit of Levy),
    (iii) omega -> 0   (Levy's acoustic frequency is negligible on the gap scale).
In that limit our branch-resolved integrand collapses identically:

    nu_sigma(E)^2 C_I(sigma,sigma; E,E) = N(eps_sigma)^2 (E/ell)^2 * xi^2/E^2 = N(eps_sigma)^2

because xi_sigma^2 = ell^2 = E^2 - Delta^2, so

    W/W_N  ->  2 INT_Delta^inf (-df/dE) dE / INT (-df_N/dE) dE  =  2 f(Delta(T)),

Levy's Eq. (8).  Run:  python levy1963_check.py
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from absorption import S_scattering, S_normal, S_paircreation
from validate_ua import FlatDOS

LAM_CUT = 30.0            # flat band half-width Lambda / t
TC = 0.1                  # Tc / t (arbitrary; ratios are universal in T/Tc)


# --------------------------------------------------------------------------- BCS s-wave input
def _lam_for_Tc():
    global _LAM
    if _LAM is not None:
        return _LAM
    x = np.linspace(1e-9, LAM_CUT, 200001)
    f = lambda lam: 1.0 - lam * np.trapezoid(np.tanh(x / (2 * TC)) / x, x)
    _LAM = brentq(f, 0.01, 5.0)
    return _LAM


_LAM = None


def gap_BCS(T: float) -> float:
    """solve 1 = lam * int_0^Lambda tanh(E/2T)/E_deps  with E = sqrt(eps^2+D^2)."""
    global _LAM
    if _LAM is None:
        _LAM = _lam_for_Tc()
    if T >= TC:
        return 0.0
    eps = np.linspace(1e-9, LAM_CUT, 60001)

    def g(D):
        E = np.hypot(eps, D)
        return 1.0 - _LAM * np.trapezoid(np.tanh(E / (2 * T)) / E, eps)
    if g(1e-12) > 0:
        return 0.0
    return brentq(g, 1e-12, 4.0 * TC)


def gap0():
    _lam_for_Tc()
    eps = np.linspace(1e-9, LAM_CUT, 60001)
    return brentq(lambda D: 1.0 - _LAM * np.trapezoid(1.0 / np.hypot(eps, D), eps), 1e-6, 4.0 * TC)


def levy_ratio(T: float) -> float:
    """Levy 1963 Eq. (8): 2 f_0(eps_0(T))."""
    D = gap_BCS(T)
    if D <= 0:
        return 1.0
    return 2.0 / (np.exp(D / T) + 1.0)


# --------------------------------------------------------------------------- our pipeline
def our_ratio(omega: float, T: float, with_pair: bool = False):
    dos = FlatDOS(Wc=LAM_CUT)
    D = gap_BCS(T)
    S, _ = S_scattering(dos, T, omega, D, 0.0, 0.0, "a_signed", nell=20000)
    SN = S_normal(dos, T, omega, 0.0)
    r = S / SN
    if with_pair and D > 0 and omega > 2 * D:
        Sp, _ = S_paircreation(dos, T, omega, D, 0.0, 0.0, "a_signed", nell=20000)
        return (S + Sp) / SN, r
    return r, r


def run_table():
    print("=" * 96)
    print("T5  s-wave reduction vs Levy 1963, PR 131, 1497, Eq. (8):  alpha_s/alpha_n = 2 f(eps_0(T))")
    print("    model reduced: DS=0 (b=0), V_SO=0, flat DOS, mu=0;  Delta(T) from the BCS gap eq.")
    print("=" * 96)
    D0 = gap0()
    print(f"    weak-coupling check: 2*Delta(0)/Tc = {2*D0/TC:.3f}  (BCS 3.53; Levy's Fig.1 uses 3.5)")
    omegas = (0.02, 0.005, 0.001, 2e-4, 5e-5)
    fracs = (0.3, 0.5, 0.7, 0.9)
    print(f"    {'T/Tc':>6} {'Levy 2f(D)':>11} " + " ".join(f"{'w=%.0e' % w:>11}" for w in omegas))
    worst = 0.0
    for fr in fracs:
        T = fr * TC
        lev = levy_ratio(T)
        vals = [our_ratio(w, T)[0] for w in omegas]
        worst = max(worst, abs(vals[-1] - lev) / lev)
        print(f"    {fr:6.2f} {lev:11.6f} " + " ".join(f"{v:11.6f}" for v in vals))
    print(f"    max rel. dev. of the omega->0 column from Levy: {worst:.2e}")
    return worst


def run_figure():
    dos = FlatDOS(Wc=LAM_CUT)
    fr = np.linspace(0.05, 0.995, 40)
    lev = [levy_ratio(f * TC) for f in fr]
    w_small = [our_ratio(1e-4, f * TC)[0] for f in fr]
    w_fin = [our_ratio(0.05, f * TC)[0] for f in fr]
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    ax[0].plot(fr, lev, "k", lw=2.0, label="Levy 1963 Eq.(8):  $2f_0(\\varepsilon_0(T))$")
    ax[0].plot(fr, w_small, "o", ms=3.5, color="tab:blue", label="this work, $\\hbar\\omega=10^{-4}t$")
    ax[0].plot(fr, w_fin, "--", lw=1.2, color="tab:red", label="this work, $\\hbar\\omega=0.05t$")
    ax[0].set_xlabel("T / Tc")
    ax[0].set_ylabel("$\\alpha_s/\\alpha_n = W/W_N$")
    ax[0].set_title("s-wave reduction (flat DOS, $\\mu=0$, $\\Delta^s=0$, $V_{SO}=0$)")
    ax[0].legend(fontsize=8)
    ws = np.array([0.02, 0.005, 0.001, 2e-4, 5e-5, 1e-5])
    dev = [abs(our_ratio(w, 0.5 * TC)[0] - levy_ratio(0.5 * TC)) / levy_ratio(0.5 * TC) for w in ws]
    ax[1].loglog(ws, dev, "k-o", ms=4)
    ax[1].set_xlabel("$\\hbar\\omega$ [t]")
    ax[1].set_ylabel("rel. dev. from Levy")
    ax[1].set_title("convergence to Levy's limit as $\\omega\\to0$ (T/Tc=0.5)")
    ax[1].grid(True, which="both", alpha=0.3)
    fig.suptitle("Cross-check of the numerical procedure against Levy, Phys. Rev. 131, 1497 (1963)")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig("figures/levy1963_swave_check.png", dpi=150)
    plt.close(fig)
    print("figure: figures/levy1963_swave_check.png")


def regime_demo():
    """
    T5b  WHICH assumption produces 2f(Delta)?  Write the most general kernel-factorised
    ratio (the structure shared by Levy's Boltzmann solution and by any q->0 golden rule)

        L[K] = 2 INT_Delta^inf K(E) (-df/dE) dE  /  INT_-inf^inf K(|eps|) (-df_N/deps) deps

    and probe it:  K=1 (Levy's ql<1, omega*tau<<1 angular average; our C_I*Jacobian=1
    identity) gives 2f(Delta) exactly; any K with energy structure on the scale T or
    Delta destroys it.  Also compare, at FINITE omega, our E'=E+hw-resolved ratio against
    Levy's prescription (occupation factors at E only), to delimit the benchmark scope.
    """
    from scipy.integrate import quad
    print("=" * 96)
    print("T5b  regime analysis: kernel energy-independence is the operative assumption")
    print("=" * 96)
    T = 0.5 * TC
    D = gap_BCS(T)
    fp = lambda E: np.exp(E / T) / (T * (np.exp(E / T) + 1) ** 2)   # -df/dE > 0
    lev = 2.0 / (np.exp(D / T) + 1.0)
    for name, K in (("K=1 (Levy kernel / our identity)", lambda E: 1.0 + 0 * E),
                    ("K=E/Delta", lambda E: E / D),
                    ("K=1+0.1(E/Delta)^2", lambda E: 1.0 + 0.1 * (E / D) ** 2)):
        num = 2 * quad(lambda E: K(E) * fp(E), D, D + 60 * T, limit=400)[0]
        den = quad(lambda e: K(abs(e)) * fp(e), -60 * T, 60 * T, limit=400)[0]
        print(f"    {name:32s}  L[K] = {num/den:.6f}   vs 2f(Delta) = {lev:.6f}"
              f"   dev = {abs(num/den-lev)/lev:.2e}")
    print("    -> only an energy-INDEPENDENT kernel yields Levy's ratio; in Levy's regime")
    print("       that constancy comes from the angular average at ql<1, omega*tau<<1,")
    print("       in ours from the exact identity C_I(E,E)*(E/ell)*(ell/E) = 1 at q->0.")
    print()
    print("    finite-omega delimitation (T/Tc=0.5): our ratio resolves E'=E+hw, Levy's")
    print("    transport treatment keeps occupations at E; they differ at finite hw:")
    for w in (1e-4, 0.05, 0.7):
        r, _ = our_ratio(w, T)
        print(f"      hbar*omega = {w:<7.4f} t :  this work = {r:.5f}   Levy prescription = {lev:.5f}"
              f"   dev = {abs(r-lev)/lev:.2e}")
    print("    -> the benchmark is exact as omega->0 (physical ultrasound: hw/Delta ~ 1e-4);")
    print("       at the paper's theoretical hw = 0.7t ~ Delta the two prescriptions are")
    print("       different calculations, and Levy is NOT the applicable benchmark there.")


if __name__ == "__main__":
    run_table()
    regime_demo()
    run_figure()
