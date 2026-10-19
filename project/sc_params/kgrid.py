"""
kgrid.py -- k-space grid and pairing kernel for arbitrary (U, Delta_t, V_SO).
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

def build_kgrid(V_SO, U, Delta_t, Nk=Nk_DEF):
    kx = np.linspace(0, np.pi, Nk)
    ky = np.linspace(0, np.pi, Nk)
    KX, KY = np.meshgrid(kx, ky)
    sqrt_sin = np.sqrt(np.sin(KX) ** 2 + np.sin(KY) ** 2)
    cos_sum = np.cos(KX) + np.cos(KY)
    s_k_factor = 0.5 * (np.cos(KX) + np.cos(KY))
    sp = np.array([1.0, -1.0])
    eps = np.zeros((Nk, Nk, 2))
    eps[:, :, 0] = -2.0 * T * cos_sum - 2.0 * V_SO * sqrt_sin
    eps[:, :, 1] = -2.0 * T * cos_sum + 2.0 * V_SO * sqrt_sin
    term = (U + 8.0 * Delta_t * s_k_factor[..., None]
            + 4.0 * (Delta_t / T) * V_SO * sp[None, None :] * sqrt_sin[..., None])
    return eps, term


