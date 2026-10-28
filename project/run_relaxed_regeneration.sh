#!/usr/bin/env bash
# D1-only regeneration for the adopted relaxed chemical-potential window.
# Deliberately omits the independent DOS-smoothing calibration (D2) from
# run_full_regeneration.sh. Regenerates all solver-dependent stage-1/stage-2,
# companion-family, U-family, comparison, and paper-figure outputs, then validates
# and writes the JSON manifest. Historical archives are not modified.
# Usage: bash project/run_relaxed_regeneration.sh [log file]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
LOG="${1:-${TMPDIR:-/tmp}/jubilant-enigma-relaxed-regeneration.log}"
exec > >(tee "$LOG") 2>&1
step() { echo; echo "== [$(date -u +%H:%M:%S)] $*"; }
step "1. stage 1: Tc, gaps, chemical potentials and policy metadata"
( cd "$HERE/sc_params" && python3 main_tc.py all && python3 plotting_tc.py )
cp "$HERE/sc_params/figures/Tc_vs_VSO.png" "$ROOT/SOCfigures/08_SOC_Tc_vs_VSO_stage1.png"
cp "$HERE/sc_params/figures/gap_vs_T_comparison.png" "$ROOT/SOCfigures/30_gap_vs_T_stage1.png"
step "2. stage 2: production, exports, figures and validation gates"
( cd "$HERE/ultrasonic" && python3 main_ua.py all && python3 run_production.py && python3 export_results.py && python3 plotting_ua.py )
( cd "$HERE/ultrasonic" && python3 validate_ua.py && python3 kernel_bdg_check.py && python3 case2_checks.py > /dev/null \
    && python3 levy1963_check.py > /dev/null && python3 dualspace_check.py > /dev/null )
step "3. companion families and derived SOC/parameter scans"
( cd "$HERE/ultrasonic" && python3 make_uvso_families.py && python3 make_delta_families.py \
    && python3 make_soc_figures.py && python3 make_soc_explicit.py && python3 make_soc_enhancement_hist.py \
    && python3 make_soc_teaching_figs.py && python3 param_scans_peaks.py && python3 param_scan_features.py \
    && python3 make_param_figures.py )
cp "$HERE/ultrasonic/figures/WWN_vs_TTc_and_Delta_Ufamilies.png" "$ROOT/SOCfigures/31_WWN_vs_TTc_and_Delta_Ufamilies.png"
cp "$HERE/ultrasonic/figures/WWN_vs_TTc_and_Delta_VSOfamilies.png" "$ROOT/SOCfigures/32_WWN_vs_TTc_and_Delta_VSOfamilies.png"
cp "$HERE/ultrasonic/figures/WWN_vs_U_families_VSO.png" "$ROOT/SOCfigures/33_WWN_vs_U_families_VSO.png"
cp "$HERE/ultrasonic/figures/WWN_vs_VSO_families_U.png" "$ROOT/SOCfigures/34_WWN_vs_VSO_families_U.png"
cp "$HERE/ultrasonic/figures/soc_explicit_vertex.png" "$ROOT/SOCfigures/35_soc_explicit_vertex.png"
cp "$HERE/ultrasonic/figures/soc_enhancement_histogram.png" "$ROOT/SOCfigures/36_SOC_enhancement_histogram.png"
step "4. U-family outputs"
( cd "$HERE/ultrasonic_U" && python3 sync_u_tree.py && python3 main_ua_u.py && python3 plotting_ua_u.py )
step "5. before/after tables"
( cd "$ROOT" && python3 tools/before_after.py )
step "6. regression and relaxed-policy tests"
( cd "$HERE/tests" && python3 test_fixes.py )
step "7. refresh paper inputs, replacement figures, addendum figures and PDFs"
( cd "$HERE/paper_correction" && python3 refresh_inputs.py && python3 make_replaced_figs.py \
    && python3 make_addendum_figs.py && python3 build_corrected_pdf.py && python3 build_fa_pdf.py )
step "8. regenerate manifest for refreshed JSON outputs"
( cd "$ROOT" && python3 tools/manifest.py )
echo; echo "RELAXED-POLICY REGENERATION COMPLETE  $(date -u +%FT%TZ)"
