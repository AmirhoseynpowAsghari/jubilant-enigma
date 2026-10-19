"""
io_utils_ua.py -- json/csv I/O for stage-2 results.
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
"""
from __future__ import annotations

import json
import os
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def save(name, obj):
    path = os.path.join(HERE, "results", name)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1)
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
