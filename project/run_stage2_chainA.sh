#!/usr/bin/env bash
# Stage-2 chain A: production W/W_N, exports, validation gates.  Strict: stops at first failure.
set -euo pipefail
cd "$(dirname "$0")/ultrasonic"
run() { echo "== [$(date -u +%H:%M:%S)] $* =="; "$@"; echo "   rc=$?"; }
run python3 main_ua.py all
run python3 run_production.py
run python3 export_results.py
run python3 plotting_ua.py
run python3 validate_ua.py
run python3 kernel_bdg_check.py
run python3 case2_checks.py
run python3 levy1963_check.py
run python3 dualspace_check.py
echo "CHAIN_A_DONE"
