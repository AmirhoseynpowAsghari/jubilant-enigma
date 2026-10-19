"""
parameters_tc.py -- solver tolerances and the registry of parameter sets.
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

T = 1.0
Nk_DEF = 80
BISECTION_TOL = 0.001
BISECTION_ITER = 35
GAP_THRESHOLD_FRAC = 0.01
MIN_GAP_FOR_SC = 1e-6

# ---------------------------------------------------------------------------
# registry of parameter sets (arXiv:2607.02812 groups A-D and the V_SO scans).
# each entry: name, U/t, Delta_t/t, filling n, V_SO  (V_SO=None -> scan)
# ---------------------------------------------------------------------------
SETS = [
    dict(group="A", label="U=-t",    U=-1.0,  Delta_t=0.0, n=0.8,   V_SO=0.0),
    dict(group="A", label="U=-2t",   U=-2.0,  Delta_t=0.0, n=0.8,   V_SO=0.0),
    dict(group="A", label="U=-3t",   U=-3.0,  Delta_t=0.0, n=0.8,   V_SO=0.0),
    dict(group="A", label="U=-4t",   U=-4.0,  Delta_t=0.0, n=0.8,   V_SO=0.0),
    dict(group="B", label="U=75t",   U=75.0,  Delta_t=4.5, n=1.885, V_SO=0.0),
    dict(group="B", label="U=115t",  U=115.0, Delta_t=4.5, n=1.885, V_SO=0.0),
    dict(group="B", label="U=125t",  U=125.0, Delta_t=4.5, n=1.885, V_SO=0.0),
    dict(group="C", label="Dt=3.5",  U=75.0,  Delta_t=3.5, n=1.885, V_SO=0.8),
    dict(group="C", label="Dt=4.0",  U=75.0,  Delta_t=4.0, n=1.885, V_SO=0.8),
    dict(group="C", label="Dt=4.5",  U=75.0,  Delta_t=4.5, n=1.885, V_SO=0.8),
    dict(group="D", label="VSO=0.0", U=95.0,  Delta_t=4.5, n=1.885, V_SO=0.0),
    dict(group="D", label="VSO=0.4", U=95.0,  Delta_t=4.5, n=1.885, V_SO=0.4),
    dict(group="D", label="VSO=0.8", U=95.0,  Delta_t=4.5, n=1.885, V_SO=0.8),
]
MASTER = dict(U=95.0, Delta_t=4.5, n=1.885)          # the paper's master set
SCAN_U = (95.0, 125.0, -75.0)                        # U values of the V_SO scans
VSO_GRID = list(np.linspace(0.0, 0.8, 25))
