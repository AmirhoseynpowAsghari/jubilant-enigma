#!/usr/bin/env python3
"""Reproduce the HHM 2018 Fig. 8 parameter set at n=1.875.

Run from any working directory with the repository requirements installed:
    python project/sc_params/benchmarks/hhm_1801_08988_n1875/compare_fig8.py

This is an additive benchmark. It does not alter the production filling or overwrite the
project's existing result files. The paper's linearized Tc determinant is evaluated directly;
the repository's nonlinear 1%-gap threshold is included only as a comparison.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
SC_PARAMS = REPO / "project" / "sc_params"
sys.path.insert(0, str(SC_PARAMS))

from gap_equations import make_residual  # noqa: E402
from kgrid import build_kgrid  # noqa: E402
from parameters_tc import T as HOPPING  # noqa: E402
from solvers import _gap_at_T, solve_set  # noqa: E402

U = 90.0
DELTA_T = 4.5
FILLING = 1.875
VSO_VALUES = (0.0, 0.2, 0.4, 0.6, 0.8)
LINEARIZED_NK = (80, 160, 320)
NONLINEAR_NK = 160
LOW_T = 0.01


def paper_linearized_tc(vso: float, nk: int) -> float:
    """Solve HHM's linearized gap equations and number equation for Tc.

    The paper presents a 3x3 determinant in (Delta^0, Delta^s, Delta^(x-y)).
    Substituting its Eq. (35), Delta^(x-y)=(V_SO/2t) Delta^s, gives the equivalent
    2x2 system below, using the repository's HHM kernel and midpoint k grid.
    """
    eps, kernel = build_kgrid(vso, U, DELTA_T, Nk=nk)
    nstates_per_helicity = nk * nk

    def determinant(temp: float) -> float:
        # At Tc the linearized gap is zero, so Eq. (38) is the normal-state number equation.
        def number_residual(mu: float) -> float:
            z = np.clip((eps - mu) / (2.0 * temp), -350.0, 350.0)
            return 1.0 - float(np.tanh(z).mean()) - FILLING

        mu = brentq(
            number_residual,
            float(eps.min() - 20.0),
            float(eps.max() + 20.0),
            xtol=1e-12,
        )
        xi_abs = np.abs(eps - mu)
        g = np.tanh(xi_abs / (2.0 * temp)) / (2.0 * np.maximum(xi_abs, 1e-300))
        n = nstates_per_helicity
        t = HOPPING

        # Linearize the two independent HHM residuals at Delta=0:
        # [D0, DS]^T -> M [D0, DS]^T = 0.
        matrix = np.array(
            [
                [
                    1.0 + np.sum(kernel * g) / (2.0 * n),
                    -np.sum(kernel * eps * g) / (8.0 * t * n),
                ],
                [
                    4.0 * DELTA_T * np.sum(g) / n,
                    1.0 - DELTA_T * np.sum(eps * g) / (t * n),
                ],
            ],
            dtype=float,
        )
        return float(np.linalg.det(matrix))

    # Search for all determinant zero crossings and retain the (unique here) transition root.
    temps = np.geomspace(0.01, 2.0, 401)
    values = [determinant(float(temp)) for temp in temps]
    roots = []
    for left, right, f_left, f_right in zip(
        temps[:-1], temps[1:], values[:-1], values[1:]
    ):
        if f_left * f_right < 0.0:
            roots.append(brentq(determinant, float(left), float(right), xtol=1e-12))
    if len(roots) != 1:
        raise RuntimeError(
            f"Expected one linearized Tc root for VSO={vso}, Nk={nk}; found {roots}"
        )
    return float(roots[0])


def main() -> None:
    output = HERE / "fig8_n1875_results.csv"
    fields = [
        "VSO_over_t",
        "Tc_linearized_Nk80_over_t",
        "Tc_linearized_Nk160_over_t",
        "Tc_linearized_Nk320_over_t",
        "Tc_repo_1pct_threshold_Nk160_over_t",
        "Emin_at_T001_Nk160_over_t",
        "2Emin_over_Tc_linearized",
        "mu0_inside_band",
    ]
    rows = []

    for vso in VSO_VALUES:
        linearized = {
            nk: paper_linearized_tc(vso, nk) for nk in LINEARIZED_NK
        }

        # Existing implementation, at the same filling and Nk, retained as a comparator.
        stage1 = solve_set(
            U=U,
            Delta_t=DELTA_T,
            n_fill=FILLING,
            V_SO=vso,
            Nk=NONLINEAR_NK,
            check_conv=False,
        )
        if not stage1.get("Tc_bracketed", False) or stage1["Tc"] <= 0.0:
            raise RuntimeError(f"The repository Tc threshold did not bracket at VSO={vso}")

        eps, kernel = build_kgrid(vso, U, DELTA_T, Nk=NONLINEAR_NK)
        residual = make_residual(DELTA_T)
        t0 = [stage1["D0_0"], stage1["DS_0"], stage1["mu_0"]]
        _, low_t_solution = _gap_at_T(
            LOW_T,
            t0,
            FILLING,
            eps,
            kernel,
            residual,
            t0=t0,
            full=True,
        )
        d0, ds, mu = map(float, low_t_solution)
        b = ds / (4.0 * HOPPING)
        emin = abs(d0 - b * mu) / np.sqrt(1.0 + b * b)
        gap_ratio = 2.0 * emin / linearized[NONLINEAR_NK]
        window = stage1.get("mu_window") or {}

        rows.append(
            {
                "VSO_over_t": vso,
                "Tc_linearized_Nk80_over_t": linearized[80],
                "Tc_linearized_Nk160_over_t": linearized[160],
                "Tc_linearized_Nk320_over_t": linearized[320],
                "Tc_repo_1pct_threshold_Nk160_over_t": float(stage1["Tc"]),
                "Emin_at_T001_Nk160_over_t": float(emin),
                "2Emin_over_Tc_linearized": float(gap_ratio),
                "mu0_inside_band": bool(window.get("inside", False)),
            }
        )

    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {output}")
    for row in rows:
        print(
            "VSO/t={VSO_over_t:.1f}: Tc(linear)={Tc_linearized_Nk320_over_t:.9f}, "
            "Tc(1% proxy)={Tc_repo_1pct_threshold_Nk160_over_t:.9f}, "
            "Emin(T=.01)={Emin_at_T001_Nk160_over_t:.7f}, "
            "2Emin/Tc={2Emin_over_Tc_linearized:.6f}".format(**row)
        )


if __name__ == "__main__":
    main()
