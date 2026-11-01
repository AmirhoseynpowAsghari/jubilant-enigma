#!/usr/bin/env python3
"""Exploratory T=0 root diagnostics for intermediate fillings n=1.4 and 1.5.

This is an additive scan: it does not overwrite the existing n=1.885 scans or claim
thermodynamic stability.  For each point it records the largest-Deff root found by the
existing multi-start solver, its number/gap residuals, and the diagnostic bare-band
window for mu.  Free-energy ranking is intentionally not attempted here because the
repository does not yet implement the HHM mean-field thermodynamic potential.

Run from project/ultrasonic:
    PYTHONPATH=../sc_params python3 scan_filling_roots.py
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCP = os.path.join(HERE, "..", "sc_params")
sys.path.insert(0, SCP)

import solvers
from gap_equations import find_mu_normal, make_residual
from kgrid import build_kgrid

FILLINGS = (1.4, 1.5)
U_VALUES = tuple(float(x) for x in np.arange(15.0, 40.0001, 5.0))
VSO_VALUES = tuple(float(x) for x in np.arange(0.60, 0.9001, 0.05))
DELTA_T = 4.5
NK = 80
HIGH_NK = 160
HIGH_U_VALUES = (15.0, 25.0, 40.0)
HIGH_VSO_VALUES = (0.60, 0.80, 0.90)


def evaluate(U: float, n_fill: float, vso: float, nk: int = NK) -> dict:
    eps, kernel = build_kgrid(vso, U, DELTA_T, Nk=nk)
    residual = make_residual(DELTA_T)
    D0, DS, mu, Deff = solvers.find_T0_solution(n_fill, eps, kernel, residual)
    x = [D0, DS, mu]
    r = np.asarray(residual(x, np.inf, n_fill, eps, kernel), dtype=float)
    band = [float(np.min(eps)), float(np.max(eps))]
    normal_seed = find_mu_normal(n_fill, np.inf, eps)
    record = dict(
        U=float(U), Delta_t=float(DELTA_T), n=float(n_fill), VSO=float(vso), Nk=int(nk),
        root_method="existing T=0 multi-start; largest effective-gap root found",
        D0_0=float(D0), DS_0=float(DS), Deff_0=float(Deff), mu_0=float(mu),
        residual_max=float(np.max(np.abs(r))),
        number_residual=float(r[2]),
        band_window=band,
        mu_inside_band=bool(band[0] <= mu <= band[1]),
        mu_excess=float(max(mu - band[1], band[0] - mu, 0.0)),
        normal_mu_seed=float(normal_seed),
        has_nonzero_gap=bool(Deff > solvers.MIN_GAP_FOR_SC),
        free_energy=None,
        stability="not assessed: HHM mean-field free energy not implemented",
    )
    return record


def main() -> None:
    rows = [evaluate(U, n, v) for n in FILLINGS for v in VSO_VALUES for U in U_VALUES]
    high_rows = [evaluate(U, n, v, nk=HIGH_NK)
                 for n in FILLINGS for v in HIGH_VSO_VALUES for U in HIGH_U_VALUES]
    out_json = os.path.join(HERE, "results", "filling_roots_n140_n150.json")
    out_csv = os.path.join(HERE, "results", "filling_roots_n140_n150.csv")
    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    payload = dict(
        provenance=dict(
            scan="T=0 intermediate-filling root diagnostic",
            fillings=list(FILLINGS), U_values=list(U_VALUES), VSO_values=list(VSO_VALUES),
            Delta_t=DELTA_T, Nk=NK, policy=solvers.MU_WINDOW_POLICY,
            high_resolution=dict(Nk=HIGH_NK, U_values=list(HIGH_U_VALUES),
                                 VSO_values=list(HIGH_VSO_VALUES)),
            method="existing multi-start solver; only its largest-Deff root is recorded",
            stability="not assessed; no free-energy comparison is performed",
            note="mu outside the bare-band window is a diagnostic, not a rejection criterion",
            utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        ),
        records=rows,
        high_resolution_rechecks=high_rows,
    )
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, allow_nan=False)
    fields = [
        "U", "Delta_t", "n", "VSO", "Nk", "D0_0", "DS_0", "Deff_0", "mu_0",
        "residual_max", "number_residual", "band_min", "band_max", "mu_inside_band",
        "mu_excess", "normal_mu_seed", "has_nonzero_gap", "free_energy", "stability",
    ]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                **{k: row.get(k) for k in fields if k not in ("band_min", "band_max")},
                "band_min": row["band_window"][0],
                "band_max": row["band_window"][1],
            })
    print(f"Wrote {out_json}")
    print(f"Wrote {out_csv}")
    for n in FILLINGS:
        subset = [r for r in rows if r["n"] == n]
        high_subset = [r for r in high_rows if r["n"] == n]
        n_out = sum(not r["mu_inside_band"] for r in subset)
        n_out_high = sum(not r["mu_inside_band"] for r in high_subset)
        print(f"n={n:g}: {len(subset)} coarse points, {n_out} out-of-band; "
              f"{len(high_subset)} high-resolution points, {n_out_high} out-of-band")


if __name__ == "__main__":
    main()
