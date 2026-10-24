#!/usr/bin/env bash
# Regeneration driver (audit F-05/F-16).  Errors are FATAL: `set -euo pipefail`, no `|| true`.
# Stage order is the dependency order of the corrected pipeline (see CORRECTION_LOG.md).
set -euo pipefail
R="$(cd "$(dirname "$0")" && pwd)/project"
step() { echo "== $* ==" ; }

step "stage 1: Tc / gap / DS / mu on the midpoint grid (F-02, F-08)"
( cd "$R/sc_params" && python3 main_tc.py all )

step "stage 2: production W/W_N"
( cd "$R/ultrasonic" && python3 main_ua.py all && python3 run_production.py && python3 export_results.py && python3 plotting_ua.py )

step "validation (exit codes are gating)"
make -C "$(dirname "$0")" test

step "companion families, SOC and peak studies"
( cd "$R/ultrasonic" && python3 make_uvso_families.py && python3 make_delta_families.py \
    && python3 make_soc_figures.py && python3 make_soc_explicit.py && python3 make_soc_enhancement_hist.py \
    && python3 param_scans_peaks.py && python3 param_scan_features.py && python3 make_param_figures.py )

echo "regen done: $(date -u +%FT%TZ)"
