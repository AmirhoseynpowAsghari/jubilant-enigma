#!/usr/bin/env python3
"""Regression checks for the postprocessed W/W_N temperature maxima."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ULTRA = ROOT / "project" / "ultrasonic"
sys.path.insert(0, str(ULTRA))

from extract_wwn_temperature_maxima import (  # noqa: E402
    FRACS, VSO_VALUES, RESULTS, _load_sources, _normalise_records, extract,
)


class TestWWNTemperatureMaxima(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw = _load_sources(RESULTS / "uvso_families.json",
                            RESULTS / "strict_policy" / "strict_family_data.json")
        cls.records = []
        for policy, data in raw.items():
            cls.records.extend(_normalise_records(data, policy))
        cls.peaks, cls.envelope = extract(cls.records)

    def test_one_peak_record_per_policy_soc_temperature(self):
        self.assertEqual(len(self.peaks), 2 * len(VSO_VALUES) * len(FRACS))
        for r in self.peaks:
            self.assertIsNotNone(r["U_peak"])
            self.assertIsNotNone(r["W_WN_peak"])
            self.assertGreater(r["n_valid_U"], 0)
            self.assertGreaterEqual(r["U_peak"], r["U_min"])
            self.assertLessEqual(r["U_peak"], r["U_max"])

    def test_strict_peaks_are_admissible_and_relaxed_peaks_are_diagnostics(self):
        strict = [r for r in self.peaks if r["policy"] == "strict"]
        relaxed = [r for r in self.peaks if r["policy"] == "relaxed"]
        self.assertTrue(all(r["mu0_inside_peak"] is True for r in strict))
        self.assertTrue(all(r["mu_T_inside_peak"] is True for r in strict))
        self.assertTrue(all(r["mu0_inside_peak"] is False for r in relaxed))

    def test_expected_strict_peak_tracks_from_source_grid(self):
        expected = {0.0: 52.5, 0.2: 54.5, 0.4: 61.0, 0.6: 72.0, 0.8: 86.5}
        for vso, expected_u in expected.items():
            q = [r for r in self.peaks if r["policy"] == "strict" and r["V_SO"] == vso]
            self.assertEqual(len(q), len(FRACS))
            for r in q:
                if not (vso == 0.8 and r["frac"] == 0.999):
                    self.assertEqual(r["U_peak"], expected_u)

    def test_upper_envelope_is_separate_reduction(self):
        self.assertGreater(len(self.envelope), len(self.peaks))
        self.assertTrue(all(r["n_valid_temperatures"] >= 1 for r in self.envelope))
        self.assertTrue(all(r["frac_at_max"] in FRACS for r in self.envelope))


if __name__ == "__main__":
    unittest.main(verbosity=2)
