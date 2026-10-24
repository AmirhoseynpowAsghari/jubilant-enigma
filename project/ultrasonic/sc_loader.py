"""
sc_loader.py -- load stage-1 (sc_params) output: Tc, D0(0), DS(0), Deff(0).
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation

AUDIT FIXES (2026-10-24)
  F-09  sc_params_at no longer SNAPS to the nearest grid value: the requested V_SO must
        match a finite grid entry to within VSO_TOL, and entries flagged unconverged
        (Tc=null, F-08) are refused with a descriptive error.
  F-10  the module-level cache is keyed by the RESOLVED stage-1 root, so two different
        roots can no longer silently share one cached scan.
"""
from __future__ import annotations

import json
import os
import numpy as np

VSO_TOL = 1e-6            # max |V_SO(requested) - V_SO(grid)| accepted (F-09)

_STAGE1 = {}              # resolved-root -> parsed json (F-10)


def _resolve_root(root):
    root = root or os.path.join(os.path.dirname(__file__), "..", "sc_params")
    return os.path.abspath(root)


def stage1(root: str = None):
    """Load sc_params/results/tc_vs_vso.json (the master V_SO scan), cached per root."""
    key = _resolve_root(root)
    if key not in _STAGE1:
        with open(os.path.join(key, "results", "tc_vs_vso.json")) as fh:
            _STAGE1[key] = json.load(fh)
    return _STAGE1[key]


def sc_params_at(VSO: float, root: str = None, tol: float = VSO_TOL):
    """Tc, D0(0), DS(0), Delta_eff(0) at EXACTLY the requested V_SO of the master scan.

    Raises ValueError if VSO is not finite, is not on the grid within `tol`, or if the
    matching stage-1 entry is flagged unconverged / non-superconducting (Tc null or <= 0).
    """
    VSO = float(VSO)
    if not np.isfinite(VSO):
        raise ValueError(f"V_SO must be finite (got {VSO!r})")
    data = stage1(root)
    r = data["results"]
    V = np.asarray(r["V_SO"], float)
    good = np.isfinite(V)
    if not good.any():
        raise ValueError("stage-1 scan contains no finite V_SO values")
    d = np.abs(V - VSO)
    d[~good] = np.inf
    i = int(np.argmin(d))
    if not (d[i] <= tol):
        raise ValueError(
            f"V_SO={VSO!r} is not on the stage-1 grid of {_resolve_root(root)} "
            f"(nearest grid value {V[i]!r} at distance {d[i]:.3e} > tol={tol:g}); "
            f"refusing to snap (audit F-09). Grid has {int(good.sum())} entries.")
    Tc = r["Tc"][i]
    if Tc is None or not np.isfinite(float(Tc)) or float(Tc) <= 0.0:
        raise ValueError(
            f"stage-1 entry V_SO={V[i]!r} has Tc={Tc!r} (flagged unconverged or "
            f"non-superconducting, audit F-08); refusing to build an SC state from it.")
    return dict(VSO=float(V[i]), Tc=float(Tc), D0_0=float(r["D0_T0"][i]),
                DS_0=float(r["DS_T0"][i]), Deff_0=float(r["Deff_T0"][i]),
                meta=data.get("meta", {}))
