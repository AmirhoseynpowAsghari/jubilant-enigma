"""
main_tc.py -- stage-1 entry point: python main_tc.py [groups|scan|uscans|all].
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

import sys
import scanning_tc as sc

what = sys.argv[1] if len(sys.argv) > 1 else "groups"
if what in ("groups", "all"):
    sc.run_groups()
if what in ("scan", "all"):
    sc.run_master_scan()
if what in ("uscans", "all"):
    sc.run_u_scans()
print("stage-1 done")
