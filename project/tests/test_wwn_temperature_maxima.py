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
from plot_temperature_peak_analysis import load_peak_rows


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

    def test_two_panel_export_uses_same_strict_peaks_and_flags_boundaries(self):
        source_csv = ULTRA / "analysis" / "wwn_temperature_maxima" / "temperature_specific_maxima.csv"
        rows = load_peak_rows(source_csv, policy="strict")
        self.assertEqual(len(rows), len(VSO_VALUES) * len(FRACS))
        self.assertTrue(all(r["policy"] == "strict" for r in rows))
        self.assertTrue(all(r["U_peak_over_t"] is not None for r in rows))
        self.assertTrue(all(r["W_WN_peak"] is not None for r in rows))
        self.assertTrue(all(r["at_admissible_domain_boundary"] is False for r in rows))
        extracted = {(r["V_SO"], r["frac"]): r for r in self.peaks if r["policy"] == "strict"}
        for r in rows:
            ref = extracted[(r["V_SO_over_t"], r["T_over_Tc"])]
            self.assertEqual(r["U_peak_over_t"], ref["U_peak"])
            self.assertEqual(r["W_WN_peak"], ref["W_WN_peak"])

    def test_relaxed_two_panel_data_are_separate_and_flag_scan_endpoints(self):
        source_csv = ULTRA / "analysis" / "wwn_temperature_maxima" / "temperature_specific_maxima.csv"
        rows = load_peak_rows(source_csv, policy="relaxed")
        self.assertEqual(len(rows), len(VSO_VALUES) * len(FRACS))
        self.assertTrue(all(r["policy"] == "relaxed" for r in rows))
        # A relaxed root has no band-admissibility cutoff: leave that field blank and flag
        # only whether the peak sits on the sampled U interval's edge.
        self.assertTrue(all(r["at_admissible_domain_boundary"] is None for r in rows))
        self.assertEqual(sum(r["at_sampled_U_domain_boundary"] is True for r in rows), 17)
        self.assertTrue(all(r["U_peak_over_t"] is not None and r["W_WN_peak"] is not None
                            for r in rows))


if __name__ == "__main__":
    unittest.main(verbosity=2)
