"""
io_utils_ua.py -- json/csv I/O for stage-2 results.
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation

F-17 fix (2026-10-24): writers run with allow_nan=False and NaN/inf floats are sanitised
to JSON null BEFORE the dump (with a warning count), so a silent NaN can never enter a
shipped JSON; every run also stores results/PROVENANCE.json (git head, UTC time, versions,
numerics configuration).
"""
from __future__ import annotations

import json
import os
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

_NAN_COUNT = 0


def _sanitize(o):
    """Recursively replace non-finite floats by None; convert numpy scalars."""
    global _NAN_COUNT
    if isinstance(o, dict):
        return {k: _sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_sanitize(v) for v in o]
    if isinstance(o, (bool, np.bool_)):
        return bool(o)
    if isinstance(o, (int, np.integer)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        f = float(o)
        if np.isfinite(f):
            return f
        _NAN_COUNT += 1
        return None
    return o


def save(name, obj):
    global _NAN_COUNT
    _NAN_COUNT = 0
    path = os.path.join(HERE, "results", name)
    clean = _sanitize(obj)
    if _NAN_COUNT:
        print(f"  WARNING (F-17): {_NAN_COUNT} non-finite float(s) sanitised to null "
              f"while writing {name}", flush=True)
    with open(path, "w") as fh:
        json.dump(clean, fh, indent=1, allow_nan=False)
    return path


def load(name):
    with open(os.path.join(HERE, "results", name)) as fh:
        return json.load(fh)


def save_csv(name, header, rows):
    path = os.path.join(HERE, "results", name)
    with open(path, "w") as fh:
        fh.write(",".join(header) + "\n")
        for r in rows:
            fh.write(",".join(str(x) for x in r) + "\n")
    return path


def provenance(extra=None):
    """Run provenance block (F-17): git head, UTC time, versions, plus caller extras."""
    import subprocess
    import sys
    import datetime
    try:
        head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=HERE,
                                       text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        head = "unknown"
    d = dict(git_head=head,
             utc=datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
             python=sys.version.split()[0], numpy=np.__version__)
    try:
        import scipy
        d["scipy"] = scipy.__version__
    except Exception:
        pass
    if extra:
        d.update(extra)
    return d
