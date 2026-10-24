#!/usr/bin/env bash
# FULL REGENERATION of every numerical product of the 2026-10-24 correction campaign, in
# dependency order.  Strict: the first failing step stops the run (set -euo pipefail).
# Usage:  bash project/run_full_regeneration.sh  [log file]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
LOG="${1:-$ROOT/regen_full.log}"
exec > >(tee "$LOG") 2>&1
step() { echo; echo "== [$(date -u +%H:%M:%S)] $*"; }
step "0. calibration of the DOS smoothing (F-03)"
( cd "$HERE/ultrasonic" && python3 calibrate_smooth.py )
step "1. stage 1: Tc / gap / DS / mu (midpoint grid, Nk->2Nk checks, window policy)"
( cd "$HERE/sc_params" && python3 main_tc.py all )
step "2. stage 2: production W/W_N, exports, validation gates"
( cd "$HERE/ultrasonic" && python3 main_ua.py all && python3 run_production.py && python3 export_results.py && python3 plotting_ua.py )
( cd "$HERE/ultrasonic" && python3 validate_ua.py && python3 kernel_bdg_check.py && python3 case2_checks.py > /dev/null \
    && python3 levy1963_check.py > /dev/null && python3 dualspace_check.py > /dev/null )
step "3. companion families, SOC decomposition, peak study"
( cd "$HERE/ultrasonic" && python3 make_uvso_families.py && python3 make_delta_families.py \
    && python3 make_soc_figures.py && python3 make_soc_explicit.py && python3 make_soc_enhancement_hist.py \
    && python3 make_soc_teaching_figs.py && python3 param_scans_peaks.py && python3 param_scan_features.py \
    && python3 make_param_figures.py )
step "4. U-family (groups A, B): sync physics modules, run, plot"
( cd "$HERE/ultrasonic_U" && python3 sync_u_tree.py && python3 main_ua_u.py && python3 plotting_ua_u.py )
step "5. before/after tables (from archive vs regenerated outputs)"
( cd "$ROOT" && python3 tools/before_after.py )
step "6. fix-verification suite (one executable test per finding)"
( cd "$HERE/tests" && python3 test_fixes.py )
step "7. paper inputs refreshed (sha256 verified) and paper figures regenerated"
( cd "$HERE/paper_correction" && python3 refresh_inputs.py && python3 make_replaced_figs.py )
step "8. manifest of every regenerated JSON"
( cd "$ROOT" && python3 tools/manifest.py )
echo; echo "FULL REGENERATION COMPLETE  $(date -u +%FT%TZ)"
