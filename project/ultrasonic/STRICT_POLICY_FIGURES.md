# Strict-policy ultrasonic figures

The strict-policy outputs are separate from the archived relaxed-policy results. No relaxed JSON, CSV, or PNG is overwritten. The strict data use the analytic bare-band interval and reject an out-of-band chemical potential at stage 1 or stage 2; rejected ratios are stored as `null` and plotted as gaps, never as zero. This is the professor-recommended physical-acceptance criterion for this ultrasonic analysis; it is not a general theorem or a free-energy comparison between in-band branches.

## Figure set

Generated images are under `figures/strict_policy/`:

- `WWN_vs_U_families_VSO.png`
- `WWN_vs_VSO_families_U.png`
- `WWN_vs_TTc_and_Delta_Ufamilies.png`
- `WWN_vs_TTc_and_Delta_VSOfamilies.png`
- `extrema_by_U.png`
- `soc_scans_by_U.png`
- `extrema_tracks_by_soc.png`
- `refined_extrema_zoom_U.png`
- `refined_SOC_finite_differences.png`
- `strict_policy_zoom.png` (the earlier targeted strict-policy diagnostic)

The corresponding new family and extrema datasets/provenance are in `results/strict_policy/`. The targeted zoom dataset, `refined_extrema_zoom_strict.json`, predates the family-wide refresh and is read-only; the two strict refinement figures are rebuilt from it and the strict family data. The original files in `figures/` and `results/` retain their relaxed-policy content.

## Regeneration sources

- `regenerate_strict_family_figures.py`: re-solves unique stage-1 family points, reuses archived ratios only when the strict root matches the relaxed root, and otherwise re-integrates strict stage-2 states.
- `regenerate_strict_scan_figures.py`: strict extrema scan and screened fixed-temperature curves; companion fixed-U SOC profiles are based on the strict family data.
- `regenerate_strict_refinement_figures.py`: strict-only targeted refinement and SOC finite-difference visualizations. Missing endpoints are not bridged.

Each generator records grid/resolution and policy details in its JSON provenance. The solver's strict criterion enforces band admissibility but does not rank multiple in-band roots by free energy.
