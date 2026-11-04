#!/usr/bin/env python3
"""Regression tests for the optional fully self-consistent HHM finite-T solver."""
from __future__ import annotations

import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SCP = os.path.join(ROOT, "project", "sc_params")
BENCH = os.path.join(SCP, "benchmarks", "hhm_1801_08988_n1875")
ULTRA = os.path.join(ROOT, "project", "ultrasonic")
sys.path.insert(0, ULTRA)
sys.path.insert(0, SCP)
sys.path.insert(0, BENCH)

from hhm_finite_t import (HHMFiniteTemperatureSolver, NoAdmissibleRoot,
                          linearized_tc)
from solvers import solve_set_core
from compare_fig8 import paper_linearized_tc
from self_consistent_hhm import SelfConsistentHHMState


class TestHHMFiniteTemperature(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.params = dict(U=90.0, Delta_t=4.5, n_fill=1.875, V_SO=0.4, Nk=80)
        cls.solver = HHMFiniteTemperatureSolver(**cls.params)

    def test_linearized_tc_matches_fig8_benchmark(self):
        out = linearized_tc(**self.params)
        ref = paper_linearized_tc(0.4, 80)
        self.assertAlmostEqual(out["Tc"], ref, delta=2e-10)
        self.assertLess(abs(out["determinant"]), 1e-9)
        self.assertEqual(len(out["roots"]), 1)

    def test_zero_temperature_root_matches_existing_equations(self):
        old = solve_set_core(90.0, 4.5, 1.875, 0.4, Nk=80)
        new = self.solver.t0
        self.assertIsNotNone(new)
        for k, old_k in (("D0", "D0_0"), ("DS", "DS_0"),
                         ("mu", "mu_0"), ("Deff", "Deff_0")):
            self.assertAlmostEqual(new[k], old[old_k], delta=2e-7)
        self.assertLess(new["residual_inf"], 1e-10)

    def test_continuation_residuals_and_linearized_limit(self):
        path = self.solver.solve_fractions([0.0, 0.5, 0.9, 0.99, 0.999])
        requested = [r for r in path if r["requested"]]
        self.assertEqual(len(requested), 5)
        self.assertTrue(all(np.isfinite(r["D0"]) and np.isfinite(r["DS"])
                            and np.isfinite(r["mu"]) for r in path))
        self.assertLess(max(r["residual_inf"] for r in path), 1e-7)
        self.assertTrue(all(r["mu_window"]["inside"] for r in path))
        # At Tc-, the nonlinear gap direction approaches the null direction of M(Tc).
        near = requested[-1]
        mat = np.asarray(self.solver.tc_info["matrix"], float)
        _, _, vh = np.linalg.svd(mat)
        null = vh[-1]
        direction = np.array([near["D0"], near["DS"]])
        cross = abs(float(direction[0] * null[1] - direction[1] * null[0]))
        self.assertLess(cross / np.linalg.norm(direction), 0.03)

    def test_alternative_initial_root_selection_is_explicit(self):
        with self.assertRaises(IndexError):
            HHMFiniteTemperatureSolver(**self.params, initial_root_index=1)

    def test_ultrasonic_adapter_exposes_self_consistent_state_and_normal_mu(self):
        model = SelfConsistentHHMState(90.0, 4.5, 1.875, 0.4, Nk=80,
                                       mu_policy="relaxed")
        states = model.at_many([0.0, 0.5])
        self.assertEqual([s["frac"] for s in states], [0.0, 0.5])
        for state in states:
            for key in ("D0", "DS", "mu", "mu_N", "Tc", "residual_inf"):
                self.assertTrue(np.isfinite(state[key]))
            self.assertEqual(state["method"], "self_consistent_HHM")
            self.assertLess(state["residual_inf"], 1e-7)

    def test_strict_is_an_acceptance_filter_not_a_different_equation(self):
        relaxed = HHMFiniteTemperatureSolver(30.0, 4.5, 1.885, 0.0, Nk=80,
                                             mu_policy="relaxed")
        strict = HHMFiniteTemperatureSolver(30.0, 4.5, 1.885, 0.0, Nk=80,
                                            mu_policy="strict")
        self.assertIsNotNone(relaxed.t0)
        self.assertFalse(relaxed.t0["mu_window"]["inside"])
        self.assertIsNone(strict.t0)
        with self.assertRaises(NoAdmissibleRoot):
            strict.solve_fractions([0.0, 0.5])
        self.assertAlmostEqual(strict.Tc, relaxed.Tc, delta=1e-12)


if __name__ == "__main__":
    unittest.main(verbosity=2)
