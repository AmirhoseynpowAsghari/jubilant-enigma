"""
gap_equations.py -- the three self-consistency equations (number equation included).
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

from parameters_tc import T, Nk_DEF
from kgrid import build_kgrid

def find_mu_normal(n_target, beta, epsilon_ks):
    """Return a normal-state chemical-potential seed for filling n in [0, 2].

    `epsilon_ks` contains both helicity bands, so the filling is twice the occupied
    fraction at T=0.  At finite beta, use the same number convention as make_residual:
    n = 1 - mean(tanh(beta * (epsilon - mu) / 2)).

    The T=0 equation is a step function on a finite k grid and generally has no exact
    root for arbitrary filling.  A quantile gives a robust seed without pretending
    that the discrete step function can be solved continuously.
    """
    epsilon_ks = np.asarray(epsilon_ks, dtype=float)
    if epsilon_ks.size == 0 or not np.all(np.isfinite(epsilon_ks)):
        raise ValueError("epsilon_ks must be a non-empty finite array")
    if not np.isfinite(n_target) or not 0.0 <= n_target <= 2.0:
        raise ValueError(f"n_target must lie in [0, 2], got {n_target!r}")

    if np.isinf(beta):
        if n_target == 0.0:
            return float(np.min(epsilon_ks) - 1e-12)
        if n_target == 2.0:
            return float(np.max(epsilon_ks) + 1e-12)
        return float(np.quantile(epsilon_ks, n_target / 2.0))
    if not np.isfinite(beta) or beta <= 0.0:
        raise ValueError(f"beta must be positive or +inf, got {beta!r}")

    def eq(mu):
        x = np.clip(beta * (epsilon_ks - mu) / 2.0, -350.0, 350.0)
        return 1.0 - float(np.tanh(x).mean()) - n_target

    margin = max(15.0, 20.0 / beta)
    lo, hi = float(np.min(epsilon_ks) - margin), float(np.max(epsilon_ks) + margin)
    return float(brentq(eq, lo, hi, xtol=1e-12))


def make_residual(Delta_t):
    def residual(x, beta, n_fill, epsilon_ks, term):
        D0, DS, mu = x
        Dbar = D0 - (DS / (4.0 * T)) * epsilon_ks
        xi = epsilon_ks - mu
        E = np.clip(np.sqrt(Dbar ** 2 + xi ** 2), 1e-12, None)
        tanh_f = np.ones_like(E) if np.isinf(beta) else np.tanh(beta * E / 2.0)
        g_E = tanh_f / (2.0 * E)
        Ntot = float(epsilon_ks.size // 2)   # repo: Nk**2, sum runs over 2*Nk**2 band elements
        r1 = D0 + (1.0 / (2.0 * Ntot)) * np.sum(term * Dbar * g_E)
        r2 = DS + (1.0 / (2.0 * Ntot)) * np.sum(8.0 * Delta_t * Dbar * g_E)
        r3 = (1.0 / (2.0 * Ntot)) * np.sum((xi / E) * tanh_f) - (1.0 - n_fill)
        return [r1, r2, r3]
    return residual


def effective_gap(D0, DS, epsilon_ks):
    Dbar = D0 - (DS / (4.0 * T)) * epsilon_ks
    return float(np.mean(np.abs(Dbar)))


