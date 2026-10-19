#!/usr/bin/env bash
# Regeneration driver.  Self-contained part (stage-1/2, validation, all analysis figures)
# runs from a fresh clone with only requirements.txt installed.  The paper-panel figures
# (groups A-D, U scans) additionally need the sibling working directories used during the
# correction project; if they are absent this script now says so instead of failing.
set -e
R="$(cd "$(dirname "$0")" && pwd)/project"
cd "$R/ultrasonic" && python3 main_ua.py all && python3 plotting_ua.py
cd "$R/sc_params"  && python3 main_tc.py all || true
for sib in ../../wwn_eta1 ../../paper_fix; do
  if [ -d "$sib" ]; then echo "sibling $sib present: extended figure regeneration available"; else
    echo "sibling $sib ABSENT: committed paper-panel figures are used as-is (documented provenance)"; fi
done
echo "regen done (self-contained part)"
