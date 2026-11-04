# Optional self-consistent HHM ultrasonic comparison

`compare_ultrasonic.py` compares `legacy_fit` and `self_consistent_HHM` states by passing
both to the unchanged `absorption.ratio_all` response kernel. The fitted-gap state is sampled
at the same absolute temperature as the HHM state (its own legacy Tc is used to convert that
absolute temperature to a fraction). The normal-state reference chemical potential is solved
at the same filling for each state model. No response formula was modified.

Default export for `U=30,45,75`, `V_SO=0,0.4,0.8`, and fractions `0.5,0.9` is in `results/`.
The script's SOC points at `U=75` and its frequency/DOS/quadrature settings are finite
sampling choices, not a complete scan. A focused superconducting-grid check for `U=75` is
stored under `results/nk80`, `results/nk160`, and `results/nk320` (the DOS grid and response
quadrature are fixed across those runs).

## Results and acceptance interpretation

- At `Nk=80`, the optional solver found a relaxed residual-verified T=0 root for all nine
  `(U,V_SO)` pairs. The strict bare-band policy accepted only `U=75, V_SO=0` and `0.4`;
  it rejected the other seven pairs. In particular, the strict primary analysis has no
  accepted `U=75, V_SO=0.8` branch, so no strict-window high-SOC downturn can be claimed
  from this run. Relaxed roots remain diagnostics and are not a stability ranking.
- For strict accepted `U=75` points, at `T/Tc=0.5`, the self-consistent response ratios are
  `0.0425806` (`V_SO=0`) and `0.0489022` (`V_SO=0.4`); at `T/Tc=0.9` they are `0.429735`
  and `0.485625`. The corresponding matched-absolute-T legacy-fit values are `0.0442067`,
  `0.0508450`, `0.434543`, and `0.490550`. Thus the response changes by about `-3.7%` to
  `-3.8%` at half Tc and `-1.0%` to `-1.1%` near 0.9 Tc at these strict points.
- Relaxed `U=75,V_SO=0.8` values are `0.0599497` and `0.599451` at the two fractions;
  both exceed the relaxed `V_SO=0.4` values in this limited scan. This neither verifies nor
  rules out a downturn at larger SOC or under a different root; the strict policy excludes
  the `V_SO=0.8` root here.
- Across all 22 accepted model/fraction rows, self-consistent versus legacy-fit `W/W_N`
  differences range from roughly `-6.91%` to `-1.00%`. The largest absolute ratio
  difference is about `0.0136`. The CSV carries exact row-level values, policies, residuals,
  and both chemical potentials; the plot shows both state prescriptions.
- The `U=75` focused grid study has maximum nonlinear residual below `4.4e-15` for
  `Nk=80,160,320`. For matched policy/SOC/fraction rows available at all three resolutions,
  `Tc`, coupled gaps, `mu_SC`, and response ratios agree to the reported floating-point
  precision in this run. The `V_SO=0.8` strict root remains rejected at all three resolutions.
  This is a point-specific convergence result, not a universal grid-convergence guarantee.

## Stability and published comparisons

`strict` and `relaxed` are acceptance policies, not different equations. The strict window is
used as the primary criterion for this ultrasonic analysis per the project-specific guidance;
relaxed roots are retained for diagnosis only. A strict-policy rejection is not a general
thermodynamic theorem. The calculations verify residuals and grid consistency but do not
rank roots by free energy. No fixed-filling free-energy functional has been implemented.

No digitized published numerical table for these full `W/W_N` curves is included. These
results compare the optional self-consistent pathway against the repository's own legacy-fit
path, not against a claimed exact reproduction of a published plot. The `U=75` SOC downturn
remains unresolved under the strict criterion because the high-SOC candidate is not accepted.
