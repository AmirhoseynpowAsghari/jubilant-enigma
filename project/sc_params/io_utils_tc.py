"""
io_utils_tc.py -- json I/O for stage-1 results.
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

import os

HERE = os.path.dirname(os.path.abspath(__file__))


def save(name, obj):
    path = os.path.join(HERE, "results", name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1)
    return path


def load(name):
    with open(os.path.join(HERE, "results", name)) as fh:
        return json.load(fh)
