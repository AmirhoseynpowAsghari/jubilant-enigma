"""
solvers.py -- T=0 grid-search solution, Tc bisection, and warm-started T sweeps.
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
Verbatim logic of the repository's stage-1 solver (branch master, commit cda0d56),
parameterised over (U, Delta_t, n, V_SO) so that every parameter set of arXiv:2607.02812
can be solved with the code that produced the paper's own stage-1 numbers.
"""
from __future__ import annotations

import json
import time
import numpy as np
from scipy.optimize import fsolve, brentq

from parameters_tc import (T, Nk_DEF, BISECTION_TOL, BISECTION_ITER,
                           GAP_THRESHOLD_FRAC, MIN_GAP_FOR_SC)
from gap_equations import make_residual, find_mu_normal, effective_gap
from kgrid import build_kgrid

def find_T0_solution(n_fill, epsilon_ks, term, residual):
    mu0 = find_mu_normal(n_fill, np.inf, epsilon_ks)
    best_sol, best_gap = None, 0.0
    for D0i in (0.1, 0.5, 1.0, 2.0, 5.0):
        for DSi in (0.1, 0.5, 1.0, 5.0, 10.0):
            for dmu in (-2.0, -0.5, 0.0, 0.5, 2.0):
                try:
                    sol, _, ier, _ = fsolve(residual, [D0i, DSi, mu0 + dmu],
                                            args=(np.inf, n_fill, epsilon_ks, term),
                                            full_output=True)
                    if ier == 1:
                        D0c, DSc, muc = abs(sol[0]), abs(sol[1]), sol[2]
                        Deff = effective_gap(D0c, DSc, epsilon_ks)
                        if Deff > best_gap:
                            best_gap, best_sol = Deff, (D0c, DSc, muc)
                except Exception:
                    pass
    if best_sol is None or best_gap < MIN_GAP_FOR_SC:
        return 0.0, 0.0, mu0, 0.0
    D0, DS, mu = best_sol
    return D0, DS, mu, effective_gap(D0, DS, epsilon_ks)


SEED_FRACS = (1.0, 0.8, 0.6, 0.4, 0.2, 0.05)
SEED_FRACS_FULL = tuple(round(0.05 * i, 3) for i in range(1, 21))


def _gap_at_T(TT, guess, n_fill, epsilon_ks, term, residual, t0=None, full=False):
    """Effective gap at temperature TT (max-gap over converged starts).

    The trivial solution (D0, DS) = 0 is an EXACT root of the residual at every T, and a
    single warm-started fsolve can slide into it even deep inside the SC state (the
    nontrivial root's basin is narrow at intermediate T for some (U, V_SO)).  Such a slide
    used to cap the Tc bisection at the slide temperature and produced notches in Tc(U)
    (e.g. U=17.5, V_SO=0.6: Tc 1.376 instead of ~2.55).  Multi-start over scaled copies of
    the T=0 root -- same max-gap philosophy as find_T0_solution -- removes the artifact;
    `full=True` widens the seed grid (used as escalation in the bisection).
    """
    seeds = [guess]
    if t0 is not None:
        for f in (SEED_FRACS_FULL if full else SEED_FRACS):
            seeds.append([t0[0] * f, t0[1] * f, t0[2]])
    best = (0.0, guess)
    for g in seeds:
        try:
            s, _, ier, _ = fsolve(residual, g, args=(1.0 / TT, n_fill, epsilon_ks, term),
                                  full_output=True)
        except Exception:
            continue
        if ier != 1:
            continue
        D0, DS = abs(s[0]), abs(s[1])
        Deff = effective_gap(D0, DS, epsilon_ks)
        if Deff > best[0]:
            best = (Deff, [D0, DS, s[2]])
    return best


def solve_set(U, Delta_t, n_fill, V_SO, Nk=Nk_DEF):
    """Tc, D0(0), DS(0), Deff(0), mu(0) for one parameter set (repo logic)."""
    eps, term = build_kgrid(V_SO, U, Delta_t, Nk)
    residual = make_residual(Delta_t)
    D0_T0, DS_T0, mu_T0, Deff_T0 = find_T0_solution(n_fill, eps, term, residual)
    if Deff_T0 < MIN_GAP_FOR_SC:
        return dict(U=U, Delta_t=Delta_t, n=n_fill, V_SO=V_SO,
                    Tc=0.0, D0_0=0.0, DS_0=0.0, Deff_0=0.0, mu_0=mu_T0)
    threshold = GAP_THRESHOLD_FRAC * Deff_T0
    t0 = [D0_T0, DS_T0, mu_T0]
    T_max = 0.57 * Deff_T0 * 3.0
    Deff, sol = _gap_at_T(T_max, t0, n_fill, eps, term, residual, t0=t0)
    for _ in range(10):
        if Deff <= threshold:
            break
        T_max *= 1.5
        Deff, sol = _gap_at_T(T_max, sol, n_fill, eps, term, residual, t0=t0)
    T_low, T_high = 0.0, T_max
    sol = list(t0)
    for _ in range(BISECTION_ITER):
        T_mid = 0.5 * (T_low + T_high)
        Deff, out = _gap_at_T(T_mid, sol, n_fill, eps, term, residual, t0=t0)
        if Deff <= threshold and (T_high - T_low) > 0.1:
            # narrow-basin guard: retry this midpoint with the full seed grid before
            # accepting "no SC here" (prevents trivial-root slides from capping Tc)
            Deff2, out2 = _gap_at_T(T_mid, sol, n_fill, eps, term, residual,
                                    t0=t0, full=True)
            if Deff2 > Deff:
                Deff, out = Deff2, out2
        if Deff > threshold:
            T_low, sol = T_mid, out
        else:
            T_high = T_mid
        if T_high - T_low < BISECTION_TOL:
            break
    Tc = 0.5 * (T_low + T_high)
    return dict(U=U, Delta_t=Delta_t, n=n_fill, V_SO=V_SO,
                Tc=Tc, D0_0=D0_T0, DS_0=DS_T0, Deff_0=Deff_T0, mu_0=mu_T0)


def sweep_T(set_dict, fracs, Nk=Nk_DEF):
    """
    Warm-started fsolve sweep of (D0, DS, mu) over T = frac * Tc.
    Returns lists aligned with fracs.  (Used for the auxiliary panels of the
    regenerated figures; the W/W_N evaluation itself uses the repository's
    rigid gap shape phi(T), see wwn_eta1/model.py.)
    """
    eps, term = build_kgrid(set_dict["V_SO"], set_dict["U"], set_dict["Delta_t"], Nk)
    residual = make_residual(set_dict["Delta_t"])
    guess = [set_dict["D0_0"], set_dict["DS_0"], set_dict["mu_0"]]
    out = []
    for f in fracs:
        TT = f * set_dict["Tc"]
        if TT <= 1e-9:
            D0, DS, mu = guess
        else:
            Deff, guess = _gap_at_T(TT, guess, set_dict["n"], eps, term, residual)
            D0, DS, mu = guess
        out.append(dict(frac=f, T=TT, D0=abs(D0), DS=abs(DS), mu=float(mu),
                        Deff=effective_gap(abs(D0), abs(DS), eps)))
    return out
