"""
energy_space.py -- the eta = +1 quasiparticle SCATTERING channel, computed in energy space with
both branches sigma, sigma' explicit, for three choices of the coherence factor.

---------------------------------------------------------------------------------------------
SET-UP
---------------------------------------------------------------------------------------------
Band energy eps, chemical potential mu, xi = eps - mu, gap Dbar(eps) = D0 - b*eps with
b = DS/(4t).   Quasiparticle energy

    E(eps) = sqrt( xi^2 + Dbar^2 ).

Squaring out:

    E^2 = (1 + b^2) eps^2 - 2 (mu + b D0) eps + (mu^2 + D0^2)

        = kappa^2 (eps - eps*)^2 + Emin^2,

    kappa = sqrt(1 + b^2),   eps* = (mu + b D0)/(1 + b^2),   Emin = |D0 - b mu| / kappa.

Inverting gives the TWO branches (this is where sigma enters):

    ell      = sqrt(E^2 - Emin^2)                        (>= 0)
    eps_s(E) = eps* + s * ell / kappa,                   s = +1, -1
    xi_s(E)  = eps_s(E) - mu
    Dbar_s(E)= D0 - b * eps_s(E)

and, since both branches must be counted, the quasiparticle DOS of one branch is

    nu_s(E) = N(eps_s(E)) * |d eps_s/dE| = N(eps_s(E)) * E / (kappa * ell),   E > Emin.

---------------------------------------------------------------------------------------------
THE THREE COHERENCE FACTORS  (eta = +1 throughout)
---------------------------------------------------------------------------------------------
    (a) signed xi xi'   :  C = 1/2 [ 1 + ( xi_s xi_s'   - Dbar_s Dbar_s' ) / (E E') ]   <- physically correct
    (b) |xi||xi'|       :  C = 1/2 [ 1 + ( |xi_s||xi_s'| - Dbar_s Dbar_s' ) / (E E') ]
    (c) paper's final   :  C = 1/2 [ 1 + ( |xi_s||xi_s'| + Dbar_s Dbar_s' ) / (E E') ]

---------------------------------------------------------------------------------------------
THE INTEGRAL
---------------------------------------------------------------------------------------------
Energy conservation on the QUASIPARTICLE energies (not on eps):

    E' = E + hbar*omega.

Scattering channel only (no pair creation):

    S_sc = sum_{s,s'} integral_{Emin}^{Emax} dE  nu_s(E) * nu_s'(E') * C(E,s;E',s')
                                        * [ f(E) - f(E') ],        f(E) = 1/(exp(E/T)+1)

    Emax is set by the band: eps_s(E) must satisfy |eps| <= Wband.

Normal state of the SAME model (Delta -> 0, so Emin -> 0, kappa -> 1, eps_s -> eps, nu_s -> N):

    S_N = integral dE N(E) N(E + hbar*omega) [ f_N(E) - f_N(E + hbar*omega) ],
    f_N(x) = 1/(exp((x - mu_N)/T)+1)   with the normal-state mu_N(T) at the same filling.

    W / W_N = S_sc / S_N        (the prefactors 2*pi*omega*B^2 and the k-sum normalisation cancel)

---------------------------------------------------------------------------------------------
QUADRATURE
---------------------------------------------------------------------------------------------
The substitution  ell = sqrt(E^2 - Emin^2)  (i.e. E = sqrt(Emax_eff^2 + ell^2)) removes the
1/sqrt(E^2 - Emin^2) behaviour of nu_s at the gap edge:  dE = ell dEll / E, hence
nu_s dE = N(eps_s(E)) dEll / kappa, which is smooth.  In the double-density integrand int dE nu_s(E) nu_s'(E+hw) C [f-f'] only ONE factor absorbs dE->dEll; the second stays per unit E':  nu_s'(E') = N(eps_s'(E')) E'/(kappa ell').  (Using N/kappa for both suppresses the integrand by ell'/E' -- the defect found by the Levy-1963 s-wave check.)
"""
from __future__ import annotations

import numpy as np

from parameters_ua import T_HOP, OMEGA, ETA
from dos_ua import NormalDOS


# ---------------------------------------------------------------------------------------------
# branch algebra
# ---------------------------------------------------------------------------------------------
def qp_constants(D0: float, DS: float, mu: float, t: float = T_HOP):
    b = DS / (4.0 * t)
    kappa = np.sqrt(1.0 + b * b)
    estar = (mu + b * D0) / (1.0 + b * b)
    Emin = abs(D0 - b * mu) / kappa
    return dict(b=b, kappa=kappa, estar=estar, Emin=Emin)


def branch(ell, sign, P, mu, D0, t=T_HOP):
    """eps_s, xi_s, Dbar_s at ell = sqrt(E^2 - Emin^2) on branch s = sign."""
    eps = P["estar"] + sign * ell / P["kappa"]
    return eps, eps - mu, D0 - P["b"] * eps


def E_of_ell(ell, Emin):
    return np.sqrt(ell ** 2 + Emin ** 2)


def fermi_qp(E, T):
    """Quasiparticle occupation 1/(exp(E/T)+1); no chemical potential (E measured from condensate)."""
    if T <= 1e-12:
        return np.zeros_like(np.asarray(E, float))
    x = np.asarray(E, float) / T
    return np.where(x > 500.0, 0.0, 1.0 / (np.exp(np.minimum(x, 500.0)) + 1.0))


def fermi_normal(eps, T, mu):
    if T <= 1e-12:
        return np.where(np.asarray(eps, float) < mu, 1.0, 0.0)
    x = (np.asarray(eps, float) - mu) / T
    return np.where(x > 500.0, 0.0, np.where(x < -500.0, 1.0,
                                             1.0 / (np.exp(np.clip(x, -500, 500)) + 1.0)))


def coherence(E, xi, Dbar, Ep, xip, Dbarp, variant: str) -> np.ndarray:
    """eta = +1 coherence factors.  variant: 'a_signed'=(a), 'b_abs'=(b), 'c_paper'=(c)."""
    EE = np.maximum(E * Ep, 1e-300)
    if variant in ("a_signed", "signed"):
        X = xi * xip                    # KEEP THE SIGN of xi  (requirement 1)
        D = -Dbar * Dbarp               # eta = +1
    elif variant in ("b_abs", "abs_minus"):
        X = np.abs(xi) * np.abs(xip)
        D = -Dbar * Dbarp
    elif variant in ("c_paper", "abs_plus"):
        X = np.abs(xi) * np.abs(xip)
        D = +Dbar * Dbarp               # the sign the repository/paper uses
    elif variant in ("d_std", "std"):
        # the textbook BCS Case-I coherence factor for quasiparticle SCATTERING:
        #   M = u_k u_k' - v_k v_k'  with u >= 0 and v = sign(Dbar)*sqrt((1-xi/E)/2)
        #   -> M^2 = 1/2 [1 - Dbar Dbar'/(E E')],  NO xi xi' term.
        # This is the only form with C -> 1 at the gap edge AND C -> 1/2 at Delta = 0.
        X = 0.0 * xi * xip
        D = -Dbar * Dbarp
    elif variant in ("case2_signed", "mb"):
        # Case II (electromagnetic / NMR): |u u' + v v'|^2 = 1/2[1 + (xi xi' + Dbar Dbar')/EE']
        X = xi * xip
        D = +Dbar * Dbarp
    else:
        raise ValueError(variant)
    return 0.5 * (1.0 + (X + D) / EE)


VARIANTS = (("a_signed", "(a) signed xi xi'"),
            ("b_abs", "(b) |xi||xi'|"),
            ("c_paper", "(c) paper's final"),
            ("d_std", "(d) textbook Case I, no xi xi'"))

VARIANTS3 = VARIANTS[:3]        # the three the question asks to compare


# ---------------------------------------------------------------------------------------------
# the scattering integral, branch-resolved, in ell
# ---------------------------------------------------------------------------------------------
def E_max_branch(P, Wb, s):
    """Largest quasiparticle energy on branch s whose band energy is still inside the band."""
    dl = max(Wb - s * P["estar"], 0.0)          # |eps_s - eps*| <= Wb - s*eps*
    return np.hypot(dl * P["kappa"], P["Emin"])


def pair_range(dos: NormalDOS, P, omega, s, sp):
    """
    E-range of the LOWER state for the branch pair (sigma = s, sigma' = sp).

        E >= Emin                       (both branches exist)
        E <= E_max(s)                   (lower state on branch s inside the band)
        E + omega <= E_max(sp)          (upper state on branch sp inside the band)

    E_max(+1) and E_max(-1) differ: eps* ~ mu sits near the TOP of the band here (the
    repository's filling n = 1.885 is almost full, 2 = full), so the particle-like branch
    s = +1 has a short range and the hole-like branch s = -1 a long one.  Using the symmetric
    bound kappa*(Wb - |eps*|) for both branches discards the hole-like branch entirely and
    gives S = 0 -- that is the trap.
    """
    Emin = P["Emin"]
    lo = Emin
    hi = min(E_max_branch(P, dos.Wband, s), E_max_branch(P, dos.Wband, sp) - omega)
    return lo, max(hi, lo)


def ell_limits(dos: NormalDOS, P, Emin, omega):
    """Generous outer range of ell = sqrt(E^2 - Emin^2) covering every branch pair."""
    hi = max(pair_range(dos, P, omega, s, sp)[1] for s in (+1, -1) for sp in (+1, -1))
    ell_hi = np.sqrt(max(hi ** 2 - Emin ** 2, 0.0))
    return 0.0, ell_hi


def S_scattering(dos: NormalDOS, T: float, omega: float, D0: float, DS: float, mu: float,
                 variant: str, nell: int = 6000, t: float = T_HOP, branches=(+1, -1)):
    """
    S_sc for one coherence-factor variant, with the sum over sigma, sigma' explicit.

    For each branch pair the integral over the lower-state quasiparticle energy E is done in
    ell = sqrt(E^2 - Emin^2), because

        nu_s(E) dE = N(eps_s(E)) * E/(kappa*ell) dE  and  dE = ell dEll / E
        =>  nu_s(E) dE = N(eps_s(E)) dEll / kappa,

    i.e. the 1/sqrt(E^2 - Emin^2) BCS singularity of the quasiparticle DOS is removed exactly.

    Returns (S_sc, {(sigma, sigma'): contribution}).
    """
    P = qp_constants(D0, DS, mu, t)
    Emin = P["Emin"]
    kappa = P["kappa"]
    parts = {}
    total = 0.0
    for sig in branches:                       # sigma  : branch of the LOWER state
        for sigp in branches:                  # sigma' : branch of the UPPER state
            E_lo, E_hi = pair_range(dos, P, omega, sig, sigp)
            if E_hi <= E_lo:
                parts[(sig, sigp)] = 0.0
                continue
            l_lo = np.sqrt(max(E_lo ** 2 - Emin ** 2, 0.0))
            l_hi = np.sqrt(max(E_hi ** 2 - Emin ** 2, 0.0))
            if l_hi <= l_lo:
                parts[(sig, sigp)] = 0.0
                continue
            dEll = (l_hi - l_lo) / nell
            ell = l_lo + (np.arange(nell) + 0.5) * dEll       # midpoint rule
            E = E_of_ell(ell, Emin)
            Ep = E + omega
            ellp = np.sqrt(np.maximum(Ep ** 2 - Emin ** 2, 0.0))
            eps_b, xi_b, Db = branch(ell, sig, P, mu, D0, t)
            eps_a, xi_a, Da = branch(ellp, sigp, P, mu, D0, t)
            nu_b = dos.N(eps_b) / kappa                       # LOWER state: nu dE = N dEll/kappa
            nu_a = dos.N(eps_a) * Ep / (kappa * np.maximum(ellp, 1e-12))  # UPPER state: per dE'
            C = coherence(E, xi_b, Db, Ep, xi_a, Da, variant)
            thermal = fermi_qp(E, T) - fermi_qp(Ep, T)        # >= 0
            val = float((nu_b * nu_a * C * thermal).sum() * dEll)
            parts[(sig, sigp)] = val
            total += val
    return total, parts


def S_normal(dos: NormalDOS, T: float, omega: float, mu_N: float, neps: int = 20000):
    """Normal-state absorption of the same model: S_N = int dE N(E) N(E+omega) [f_N(E) - f_N(E+omega)]."""
    E = np.linspace(-dos.Wband, dos.Wband, neps)
    dE = E[1] - E[0]
    Emid = E + 0.5 * dE
    val = dos.N(Emid) * dos.N(Emid + omega) * (fermi_normal(Emid, T, mu_N)
                                               - fermi_normal(Emid + omega, T, mu_N))
    return float(val.sum() * dE)


def ratio_all(dos: NormalDOS, st: dict, omega: float, variant_list=None, nell: int = 6000,
              mu_for_WN: str = "normal", t: float = T_HOP):
    """
    W/W_N for the three variants at one thermodynamic state st (from model.SCState.at).
    Returns a dict with S_sc per variant, the (sigma,sigma') decomposition, S_N and the ratios.
    """
    T, D0, DS, mu = st["T"], st["D0"], st["DS"], st["mu"]
    mu_N = st["mu_N"]
    out = dict(T=T, frac=st["frac"], D0=D0, DS=DS, Deff=st["Deff"], mu=mu, mu_N=mu_N)
    P = qp_constants(D0, DS, mu, t)
    out["Emin"] = P["Emin"]
    out["2Emin"] = 2.0 * P["Emin"]
    out["kappa"] = P["kappa"]
    out["estar"] = P["estar"]
    S_N = S_normal(dos, T, omega, mu_N if mu_for_WN == "normal" else mu)
    out["S_N"] = S_N
    vl = variant_list or [v for v, _ in VARIANTS]
    for v in vl:
        S, parts = S_scattering(dos, T, omega, D0, DS, mu, v, nell=nell, t=t)
        out[v] = S
        out[v + "_ratio"] = S / S_N if S_N != 0 else np.nan
        out[v + "_parts"] = parts
    return out


# ---------------------------------------------------------------------------------------------
# the pair-creation channel (needed for the Delta -> 0 sum rule; EXCLUDED from the answer above)
# ---------------------------------------------------------------------------------------------
def C_pair(E, xi, Dbar, Ep, xip, Dbarp, variant: str) -> np.ndarray:
    """
    Pair creation:  E + E' = hbar*omega,  C_p = 1 - C  for the same Bogoliubov matrix element
    (verified numerically: M_scat^2 + M_pair^2 = 1).

        (a) 1/2[1 - ( xi xi' - Dbar Dbar')/(E E')]
        (d) 1/2[1 + Dbar Dbar'/(E E')]
    """
    EE = np.maximum(E * Ep, 1e-300)
    if variant in ("a_signed", "signed"):
        num = -(xi * xip) + Dbar * Dbarp
    elif variant in ("b_abs", "abs_minus"):
        num = -np.abs(xi) * np.abs(xip) + Dbar * Dbarp
    elif variant in ("c_paper", "abs_plus"):
        num = -np.abs(xi) * np.abs(xip) - Dbar * Dbarp
    elif variant in ("d_std", "std"):
        num = +Dbar * Dbarp
    elif variant in ("case2_signed", "mb"):
        num = -(xi * xip) - Dbar * Dbarp
    else:
        raise ValueError(variant)
    return 0.5 * (1.0 + num / EE)


def S_paircreation(dos: NormalDOS, T: float, omega: float, D0: float, DS: float, mu: float,
                   variant: str, nell: int = 6000, t: float = T_HOP, branches=(+1, -1)):
    """
    S_pair = sum_{sig,sig'} int dE_b nu_s(E_b) nu_s'(E_a) C_p [1 - f(E_a) - f(E_b)],
    E_a = hbar*omega - E_b,  E_b in [E_min, hbar*omega/2].   Zero when hbar*omega < 2 E_min.
    """
    P = qp_constants(D0, DS, mu, t)
    Emin = P["Emin"]
    kappa = P["kappa"]
    if omega <= 2.0 * Emin:
        return 0.0, {}
    parts = {}
    total = 0.0
    for sig in branches:
        for sigp in branches:
            lo = Emin
            hi = min(0.5 * omega, E_max_branch(P, dos.Wband, sig), E_max_branch(P, dos.Wband, sigp))
            if hi <= lo:
                parts[(sig, sigp)] = 0.0
                continue
            l_lo = np.sqrt(max(lo ** 2 - Emin ** 2, 0.0))
            l_hi = np.sqrt(max(hi ** 2 - Emin ** 2, 0.0))
            if l_hi <= l_lo:
                parts[(sig, sigp)] = 0.0
                continue
            dEll = (l_hi - l_lo) / nell
            ell = l_lo + (np.arange(nell) + 0.5) * dEll
            E = E_of_ell(ell, Emin)
            Ep = omega - E
            ellp = np.sqrt(np.maximum(Ep ** 2 - Emin ** 2, 0.0))
            eps_b, xi_b, Db = branch(ell, sig, P, mu, D0, t)
            eps_a, xi_a, Da = branch(ellp, sigp, P, mu, D0, t)
            nu_b = dos.N(eps_b) / kappa                       # LOWER state: nu dE = N dEll/kappa
            nu_a = dos.N(eps_a) * Ep / (kappa * np.maximum(ellp, 1e-12))  # UPPER state: per dE'
            Cp = C_pair(E, xi_b, Db, Ep, xi_a, Da, variant)
            occ = 1.0 - fermi_qp(E, T) - fermi_qp(Ep, T)
            parts[(sig, sigp)] = float((nu_b * nu_a * Cp * occ).sum() * dEll)
            total += parts[(sig, sigp)]
    return total, parts
