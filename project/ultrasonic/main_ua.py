"""
main_ua.py -- stage-2 entry point: python main_ua.py [master|groups|scans|all].
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
"""
from __future__ import annotations

import json
import os
import time
import numpy as np

import sys
import scanning_ua as su
import sc_loader, io_utils_ua as io

what = sys.argv[1] if len(sys.argv) > 1 else "master"
if what in ("master", "all"):
    su.run_master()
if what in ("groups", "all"):
    import json, os
    gp = json.load(open(os.path.join(os.path.dirname(__file__), "..", "sc_params",
                                     "results", "groups_stage1.json")))
    N = {"A": 0.8, "B": 1.885, "C": 1.885, "D": 1.885}
    su.run_groups(gp, lambda s: N[s["group"]])
if what in ("scans", "all"):
    import json, os
    sc = json.load(open(os.path.join(os.path.dirname(__file__), "..", "sc_params",
                                     "results", "u_scans_stage1.json")))
    su.run_vso_scans(sc)
print("stage-2 done")
