"""
sc_loader.py -- load stage-1 (sc_params) output: Tc, D0(0), DS(0), Deff(0).
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
This file keeps the repository's module name and role; the physics inside is the
validated one (see project/lecture_note_WWN_corrected.md and README.md).
"""
from __future__ import annotations

import json
import numpy as np

_STAGE1 = None


def stage1(root: str = None):
    """Load sc_params/results/tc_vs_vso.json (the master V_SO scan)."""
    global _STAGE1
    if _STAGE1 is None:
        import os
        root = root or os.path.join(os.path.dirname(__file__), "..", "sc_params")
        with open(f"{root}/results/tc_vs_vso.json") as fh:
            _STAGE1 = json.load(fh)
    return _STAGE1


def sc_params_at(VSO: float, root: str = None):
    """Tc, D0(0), DS(0), Delta_eff(0) at the nearest V_SO of the master scan."""
    r = stage1(root)["results"]
    V = np.asarray(r["V_SO"], float)
    i = int(np.argmin(np.abs(V - VSO)))
    return dict(VSO=float(V[i]), Tc=float(r["Tc"][i]), D0_0=float(r["D0_T0"][i]),
                DS_0=float(r["DS_T0"][i]), Deff_0=float(r["Deff_T0"][i]),
                meta=stage1(root)["meta"])
