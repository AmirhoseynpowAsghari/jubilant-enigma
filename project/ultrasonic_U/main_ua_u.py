"""
main_ua_u.py -- U-family entry point (groups A and B are U scans).
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
"""
from __future__ import annotations

import json
import os
import time
import numpy as np

import json, os, sys
import scanning_ua_u as su

gp = json.load(open(os.path.join(os.path.dirname(__file__), "..", "sc_params",
                                 "results", "groups_stage1.json")))
N = {"A": 0.8, "B": 1.885, "C": 1.885, "D": 1.885}
sel = [p for p in gp if p["set"]["group"] in ("A", "B")]
out = su.run_groups(sel, lambda s: N[s["group"]])
for sw in out:
    s = sw["set"]
    su_io = __import__("io_utils_ua_u")
    su_io.save(f"temp_sweep_U_{s['U']:.2f}.json", sw)
print("ultrasonic_U done")
