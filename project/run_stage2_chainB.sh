#!/usr/bin/env bash
# Stage-2 chain B: companion families, SOC decomposition, peak study, U-tree.  Strict.
set -euo pipefail
cd "$(dirname "$0")/ultrasonic"
run() { echo "== [$(date -u +%H:%M:%S)] $* =="; "$@"; echo "   rc=$?"; }
run python3 make_uvso_families.py
run python3 make_delta_families.py
run python3 make_soc_figures.py
run python3 make_soc_explicit.py
run python3 make_soc_enhancement_hist.py
run python3 make_soc_teaching_figs.py
run python3 param_scans_peaks.py
run python3 param_scan_features.py
run python3 make_param_figures.py
cd ../ultrasonic_U
run python3 sync_u_tree.py
run python3 main_ua_u.py
echo "CHAIN_B_DONE"
