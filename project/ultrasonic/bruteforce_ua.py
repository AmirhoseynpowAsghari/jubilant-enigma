"""
bruteforce.py -- an INDEPENDENT evaluation of the same scattering sum, done directly on the
lattice: every (k, s) state gets its own quasiparticle energy and its own branch label, the two
momentum sums are unrestricted (the paper's dirty-limit / constant-B approximation), and energy
conservation is imposed by a Gaussian representation of delta(E_a - E_b - hbar*omega).

Nothing here uses the analytic branch inversion eps_s(E) = eps* + s sqrt(E^2 - Emin^2)/kappa, the
Jacobian |d eps/dE|, or the change of variables to energy space.  It therefore validates
  * the branch algebra and the Jacobian in energy_space.py,
  * the DOS normalisation,
  * the Delta -> 0 sum rule.

    S_sc = sum_{k,s} sum_{k',s'} C(E_a, sigma_a; E_b, sigma_b) [f(E_b) - f(E_a)]
                                                    delta(E_a - E_b - hbar*omega)

Because xi and Dbar depend on k only through eps_s(k), and eps_s(k) is fixed by (E, s), the
coherence factor depends on (E_b, s_b, E_a, s_a) alone -- that is the whole content of the
energy-space reduction.  Here we verify it by keeping the k-labels explicitly.

Discretisation: E is binned with width de, delta is represented by exp(-(x/w)^2)/(w sqrt(pi)).
The exact answer is recovered by w -> 0 at fixed de; both errors are O(w^2) and O(de^2) and are
Richardson-extrapolated.
"""
from __future__ import annotations

import numpy as np

from kgrid_ua import kgrid
from parameters_ua import T_HOP
from absorption import fermi_qp


def lattice_states(VSO: float, Nk: int, D0: float, DS: float, mu: float, t: float = T_HOP):
    """
    All (k, helicity) states of the quarter zone: band energy eps, xi, Dbar, quasiparticle
    energy E, and the BRANCH LABEL sigma.

    IMPORTANT: sigma is NOT the Rashba helicity.  sigma is the label of the branch of the
    inversion  E -> eps_sigma(E) = eps* + sigma*sqrt(E^2 - Emin^2)/kappa,  i.e. it records
    which of the two band energies with the same quasiparticle energy a state belongs to.
    Since E(eps) is a single valley with its minimum at eps*, that label is

        sigma = sign(eps - eps*)      ( = sign(xi) when D0 = 0 )

    and every state, of either helicity, carries one of the two values.  Grouping the lattice
    states by the helicity instead mixes the two branches inside one group and makes the
    coherence factor wrong (it is exactly wrong in the Delta -> 0 limit, where it produces
    spurious cross-branch weight); that was the bug in the first version of this cross-check.
    """
    e0, R = kgrid(Nk)
    eps = np.concatenate([(e0 - 2.0 * (+1.0) * VSO * R).ravel(),
                          (e0 - 2.0 * (-1.0) * VSO * R).ravel()])
    hel = np.concatenate([np.full(e0.size, +1.0), np.full(e0.size, -1.0)])
    xi = eps - mu
    Dbar = D0 - DS * eps / (4.0 * t)
    E = np.hypot(xi, Dbar)
    b = DS / (4.0 * t)
    kappa = np.sqrt(1.0 + b * b)
    estar = (mu + b * D0) / (1.0 + b * b)
    sig = np.where(eps >= estar, +1.0, -1.0)       # branch of the E -> eps inversion
    # The midpoint grid covers ONE QUARTER of the zone, and eps_s(k) is invariant under
    # kx -> -kx and ky -> -ky, so the quarter zone represents the full zone exactly once.
    # The k-average <...>_k is therefore (1/Nk^2) * sum over this list -- NOT 1/(2 Nk^2).
    # With that weight, sum over both helicities of <1> = 2, which matches the convention
    # integral N(eps) deps = 2 used for the DOS.  (Using 1/(2 Nk^2) halves the DOS and is
    # exactly the factor-of-4 error this cross-check had before it was fixed: S ~ weight^2.)
    return dict(eps=eps, sig=sig, hel=hel, xi=xi, Dbar=Dbar, E=E, estar=estar,
                nstates=e0.size, nlist=eps.size)


def _coh(x, xp, D, Dp, E, Ep, variant):
    EE = np.maximum(E * Ep, 1e-300)
    if variant in ("a_signed", "signed"):
        num = x * xp - D * Dp
    elif variant in ("b_abs", "abs_minus"):
        num = np.abs(x) * np.abs(xp) - D * Dp
    elif variant in ("c_paper", "abs_plus"):
        num = np.abs(x) * np.abs(xp) + D * Dp
    elif variant in ("d_std", "std"):
        num = D * Dp
    elif variant in ("case2_signed", "mb"):
        num = x * xp + D * Dp          # Case-II (eta = -1) kernel, cf. case2_checks.py T6
    else:
        raise ValueError(variant)
    return 0.5 * (1.0 + num / EE)


def S_sc_bf(st, T, omega, variant, Egrid, nbin, w):
    """Cleaner driver: bins over E on Egrid = (Emin_grid, Emax_grid), sums all 4 (sigma,sigma')."""
    Emin_grid, Emax_grid = Egrid
    de = (Emax_grid - Emin_grid) / nbin
    Ebin = Emin_grid + (np.arange(nbin) + 0.5) * de
    E = st["E"]
    idx = ((E - Emin_grid) / de)
    inside = (idx >= 0) & (idx < nbin)
    idx = np.clip(idx.astype(np.int64), 0, nbin - 1)
    f = fermi_qp(Ebin, T)
    thermal = f[:, None] - f[None, :]
    delta = np.exp(-(((Ebin[None, :] - Ebin[:, None] - omega) / w) ** 2)) / (w * np.sqrt(np.pi))
    base = thermal * delta
    pref = 1.0 / st["nstates"]          # <...>_k ; the helicity sum is explicit below
    total = 0.0
    for sb in (+1, -1):
        mb = (st["sig"] == sb) & inside
        if not mb.any():
            continue
        cnt_b = np.bincount(idx[mb], minlength=nbin).astype(float)
        xi_b = np.bincount(idx[mb], weights=st["xi"][mb], minlength=nbin) / np.maximum(cnt_b, 1)
        D_b = np.bincount(idx[mb], weights=st["Dbar"][mb], minlength=nbin) / np.maximum(cnt_b, 1)
        for sa in (+1, -1):
            ma = (st["sig"] == sa) & inside
            if not ma.any():
                continue
            cnt_a = np.bincount(idx[ma], minlength=nbin).astype(float)
            xi_a = np.bincount(idx[ma], weights=st["xi"][ma], minlength=nbin) / np.maximum(cnt_a, 1)
            D_a = np.bincount(idx[ma], weights=st["Dbar"][ma], minlength=nbin) / np.maximum(cnt_a, 1)
            C = _coh(xi_b[:, None], xi_a[None, :], D_b[:, None], D_a[None, :],
                     Ebin[:, None], Ebin[None, :], variant)
            total += float((cnt_b[:, None] * cnt_a[None, :] * C * base).sum())
    return total * pref ** 2


def S_pair_bf(st, T, omega, variant, Egrid, nbin, w):
    """
    Lattice (k-space) pair-creation sum, mirror of S_sc_bf:
        S_pair = sum_{k,s} sum_{k',s'} [1 - C(E_a,s_a; E_b,s_b)] [1 - f(E_a) - f(E_b)]
                 delta(E_a + E_b - hbar*omega)
    i.e. the PP channel of Sardwal et al. arXiv:2608.19079 Eqs. (A6)-(A9) in the clean
    golden-rule limit.  States with E > hbar*omega cannot participate (E_a + E_b = omega),
    so Egrid must cover [E_min, omega].

    The double sum runs over ORDERED pairs; the physical channel counts unordered pairs
    once (the energy-space S_paircreation integrates E_b in [E_min, omega/2]), hence the
    factor 1/2.  NOTE: the integrand carries TWO BCS edge singularities, nu(E) nu(w-E),
    so the Gaussian smearing error is non-analytic in w and Richardson extrapolation is
    invalid here: use w <= 0.01 directly (w=0.005 agrees with energy space to ~0.3 %).
    """
    Emin_grid, Emax_grid = Egrid
    de = (Emax_grid - Emin_grid) / nbin
    Ebin = Emin_grid + (np.arange(nbin) + 0.5) * de
    E = st["E"]
    idx = ((E - Emin_grid) / de)
    inside = (idx >= 0) & (idx < nbin)
    idx = np.clip(idx.astype(np.int64), 0, nbin - 1)
    f = fermi_qp(Ebin, T)
    thermal = 1.0 - f[:, None] - f[None, :]
    delta = np.exp(-(((Ebin[None, :] + Ebin[:, None] - omega) / w) ** 2)) / (w * np.sqrt(np.pi))
    base = thermal * delta
    pref = 1.0 / st["nstates"]
    total = 0.0
    for sb in (+1, -1):
        mb = (st["sig"] == sb) & inside
        if not mb.any():
            continue
        cnt_b = np.bincount(idx[mb], minlength=nbin).astype(float)
        xi_b = np.bincount(idx[mb], weights=st["xi"][mb], minlength=nbin) / np.maximum(cnt_b, 1)
        D_b = np.bincount(idx[mb], weights=st["Dbar"][mb], minlength=nbin) / np.maximum(cnt_b, 1)
        for sa in (+1, -1):
            ma = (st["sig"] == sa) & inside
            if not ma.any():
                continue
            cnt_a = np.bincount(idx[ma], minlength=nbin).astype(float)
            xi_a = np.bincount(idx[ma], weights=st["xi"][ma], minlength=nbin) / np.maximum(cnt_a, 1)
            D_a = np.bincount(idx[ma], weights=st["Dbar"][ma], minlength=nbin) / np.maximum(cnt_a, 1)
            C = _coh(xi_b[:, None], xi_a[None, :], D_b[:, None], D_a[None, :],
                     Ebin[:, None], Ebin[None, :], variant)
            total += float((cnt_b[:, None] * cnt_a[None, :] * (1.0 - C) * base).sum())
    return 0.5 * total * pref ** 2


def S_sc_bf_from_state(st, omega, variant, VSO=None, Nk=320, nbin=3000):
    """
    Brute-force lattice value of S_sc for a state dict from physics_ua.SCState.at().
    Richardson extrapolation in the Gaussian delta width (as in validate_ua.T3).
    """
    import numpy as _np
    from absorption import qp_constants
    vso = st.get("VSO", VSO)
    lst = lattice_states(vso, Nk, st["D0"], st["DS"], st["mu"])
    P = qp_constants(st["D0"], st["DS"], st["mu"])
    Eg = (P["Emin"] * 0.995, float(lst["E"].max()) * 1.002)
    ws = _np.array([0.08, 0.05, 0.03])
    vals = _np.array([S_sc_bf(lst, st["T"], omega, variant, Eg, nbin, w) for w in ws])
    A = _np.vstack([_np.ones_like(ws), ws ** 2]).T
    coef, *_ = _np.linalg.lstsq(A, vals, rcond=None)
    return float(coef[0])
