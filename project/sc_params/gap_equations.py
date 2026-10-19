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
    def eq(mu):
        if np.isinf(beta):
            return np.mean(epsilon_ks < mu) - n_target
        x = np.clip(beta * (epsilon_ks - mu), -500, 500)
        return 1.0 / (1.0 + np.exp(x)).mean() - n_target
    lo, hi = np.min(epsilon_ks) - 15.0, np.max(epsilon_ks) + 15.0
    try:
        return brentq(eq, lo, hi)
    except Exception:
        return float(np.mean(epsilon_ks))


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


