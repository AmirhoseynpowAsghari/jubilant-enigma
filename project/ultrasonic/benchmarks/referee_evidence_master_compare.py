"""Optional, non-destructive legacy-fit vs self-consistent HHM master comparison.

Writes only to results/referee_evidence/ (a new, separate evidence directory). It does not
regenerate or overwrite the production sweeps. Both gap methods are passed through the same
NormalDOS and absorption.ratio_all implementation, isolating the temperature-state model
change while keeping the response formula, frequency, and numerical integration settings fixed.

Run from any directory with repository requirements installed:
    python project/ultrasonic/benchmarks/referee_evidence_master_compare.py
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ULTRA = HERE.parent
ROOT = ULTRA.parents[1]
SC_PARAMS = ROOT / "project" / "sc_params"
sys.path.insert(0, str(ULTRA))
sys.path.insert(0, str(SC_PARAMS))

from absorption import ratio_all  # noqa: E402
from dos_ua import NormalDOS  # noqa: E402
from parameters_ua import N_TARGET, OMEGA, NELL  # noqa: E402
from physics_ua import SCState  # noqa: E402
from sc_loader import sc_params_at  # noqa: E402
from self_consistent_hhm import SelfConsistentHHMState  # noqa: E402

U = 95.0
DELTA_T = 4.5
FILLING = 1.885
SOC_VALUES = (0.0, 0.4082, 0.8)
# Match scanning_ua.FRACS exactly: 49 points from 0.02 through 0.98, then 0.99, 0.995, 0.999.
FRACTIONS = tuple(np.concatenate([np.linspace(0.02, 0.98, 49), [0.99, 0.995, 0.999]]))
OUTPUT = ULTRA / "results" / "referee_evidence"
CSV_PATH = OUTPUT / "master_fit_vs_hhm.csv"
JSON_PATH = OUTPUT / "master_fit_vs_hhm_summary.json"


def main() -> None:
    if CSV_PATH.exists() or JSON_PATH.exists():
        raise FileExistsError(
            f"Refusing to overwrite evidence output(s) in {OUTPUT}; move them aside explicitly first."
        )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    summary: dict[str, dict] = {}

    for vso in SOC_VALUES:
        # Use the existing production stage-1 record for the legacy-fit path.
        p = sc_params_at(vso)
        legacy = SCState.from_params(
            vso, p["Tc"], p["D0_0"], p["DS_0"], p["Deff_0"],
            n_tgt=FILLING, Nk_mu=300,
        )
        # Separate optional model: HHM linearized Tc + nonlinear finite-T gap/number equations.
        hhm = SelfConsistentHHMState(
            U, DELTA_T, FILLING, vso, Nk=80, mu_policy="relaxed",
            tc_scan_max=8.0, tc_scan_points=401, initial_root_index=0,
        )
        dos = NormalDOS(vso)
        soc_rows = []

        # Ascending fractions allow the HHM adapter to continue the selected root in T.
        hhm_states = hhm.at_many(FRACTIONS)
        for frac, hh_state in zip(FRACTIONS, hhm_states):
            fit_state = legacy.at(float(frac))
            fit_ratio = ratio_all(
                dos, fit_state, OMEGA, variant_list=("a_signed",), nell=NELL
            )["a_signed_ratio"]
            hhm_ratio = ratio_all(
                dos, hh_state, OMEGA, variant_list=("a_signed",), nell=NELL
            )["a_signed_ratio"]
            record = {
                "V_SO_over_t": float(vso),
                "T_over_Tc": float(frac),
                "Tc_legacy_fit_over_t": float(legacy.Tc),
                "Tc_HHM_linearized_over_t": float(hhm.Tc),
                "T_legacy_fit_over_t": float(fit_state["T"]),
                "T_HHM_over_t": float(hh_state["T"]),
                "D0_legacy_fit_over_t": float(fit_state["D0"]),
                "D0_HHM_over_t": float(hh_state["D0"]),
                "DS_legacy_fit_over_t": float(fit_state["DS"]),
                "DS_HHM_over_t": float(hh_state["DS"]),
                "Deff_legacy_fit_over_t": float(fit_state["Deff"]),
                "Deff_HHM_over_t": float(hh_state["Deff"]),
                "muSC_legacy_fit_over_t": float(fit_state["mu"]),
                "muSC_HHM_over_t": float(hh_state["mu"]),
                "muN_legacy_fit_over_t": float(fit_state["mu_N"]),
                "muN_HHM_over_t": float(hh_state["mu_N"]),
                "mu_HHM_inside_bare_band": bool(hh_state["mu_window"]["inside"]),
                "HHM_residual_inf": float(hh_state["residual_inf"]),
                "Wsc_over_WN_legacy_fit": float(fit_ratio),
                "Wsc_over_WN_HHM": float(hhm_ratio),
                "relative_response_change_percent": float(
                    100.0 * (hhm_ratio - fit_ratio) / fit_ratio
                ),
            }
            rows.append(record)
            soc_rows.append(record)

        summary[f"{vso:.4f}"] = {
            "n_points": len(soc_rows),
            "Tc_legacy_fit_over_t": float(legacy.Tc),
            "Tc_HHM_linearized_over_t": float(hhm.Tc),
            "Tc_relative_difference_percent": float(100.0 * (hhm.Tc - legacy.Tc) / legacy.Tc),
            "max_abs_response_difference": float(max(
                abs(r["Wsc_over_WN_HHM"] - r["Wsc_over_WN_legacy_fit"]) for r in soc_rows
            )),
            "max_abs_relative_response_change_percent": float(max(
                abs(r["relative_response_change_percent"]) for r in soc_rows
            )),
            "max_abs_D0_difference_over_t": float(max(
                abs(r["D0_HHM_over_t"] - r["D0_legacy_fit_over_t"]) for r in soc_rows
            )),
            "max_abs_DS_difference_over_t": float(max(
                abs(r["DS_HHM_over_t"] - r["DS_legacy_fit_over_t"]) for r in soc_rows
            )),
            "max_HHM_residual_inf": float(max(r["HHM_residual_inf"] for r in soc_rows)),
            "all_HHM_mu_in_band": bool(all(r["mu_HHM_inside_bare_band"] for r in soc_rows)),
        }
        print(
            f"V_SO/t={vso:.4f}: Tc fit={legacy.Tc:.9f}, Tc HHM={hhm.Tc:.9f}; "
            f"max |response relative change|={summary[f'{vso:.4f}']['max_abs_relative_response_change_percent']:.3f}%",
            flush=True,
        )

    with CSV_PATH.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "description": (
            "Optional comparison only; does not overwrite production outputs. Both paths use "
            "the same signed scattering ratio, DOS construction, omega, filling, and quadrature."
        ),
        "parameters": {"U_over_t": U, "Delta_t_over_t": DELTA_T, "n": FILLING,
                       "V_SO_over_t": list(SOC_VALUES), "hbar_omega_over_t": float(OMEGA)},
        "temperature_fractions": list(FRACTIONS),
        "n_points_per_SOC": len(FRACTIONS),
        "legacy_temperature_method": "rigid fitted phi(T); production stage-1 Tc",
        "self_consistent_temperature_method": "HHM linearized determinant Tc; nonlinear gap/number continuation",
        "response_method_both": "absorption.ratio_all, variant a_signed, S_sc/S_N",
        "quadrature_nell": int(NELL),
        "self_consistent_Nk": 80,
        "mu_policy": "relaxed (roots flagged; no stability claim)",
        "csv": CSV_PATH.name,
        "summary_by_SOC": summary,
        "caveat": (
            "This compares two state-temperature prescriptions within the existing scalar, "
            "unrestricted-transition response model. It does not validate the acoustic vertex, "
            "finite-q transport, pair-channel total, or thermodynamic stability."
        ),
    }
    JSON_PATH.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {CSV_PATH}")
    print(f"Wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
