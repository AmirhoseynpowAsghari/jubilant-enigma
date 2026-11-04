# Optional self-consistent HHM benchmark workflow

This directory is deliberately separate from `hhm_1801_08988_n1875/`; it does not rewrite
its Fig. 8 outputs. The optional solver is `project/sc_params/hhm_finite_t.py`, and the
ultrasonic adapter is `project/ultrasonic/self_consistent_hhm.py`. Production defaults,
`physics_ua.SCState`, and response formulas remain unchanged.

## Reproduce the exports

From the repository root (with the normal Python dependencies installed):

```bash
python project/sc_params/benchmarks/self_consistent_hhm/run_hhm_grid.py \
  --output project/sc_params/benchmarks/self_consistent_hhm/results \
  --nk 80 160 320 --vso 0 0.2 0.4 0.6 0.8

python project/ultrasonic/benchmarks/self_consistent_hhm/compare_ultrasonic.py \
  --output project/ultrasonic/benchmarks/self_consistent_hhm/results \
  --U 30 45 75 --vso 0 0.4 0.8
```

The first runner exports determinant and legacy 1%-gap-proxy transition temperatures,
nonlinear coupled `(D0, DS, mu)` continuation, residuals, band-window diagnostics, normal
chemical potentials, and fitted-gap versus self-consistent chemical potentials. It defaults
to the HHM Fig. 8 filling `n=1.875` and parameters `U=90t, Delta_t=4.5t`. After running it,
`python project/sc_params/benchmarks/self_consistent_hhm/export_gap_shapes.py` exports
per-helicity complete-grid statistics and representative momentum samples of
`bar-Delta_{k,s}(T)`. The second runner uses the ultrasonic benchmark filling `n=1.885`, the
existing `absorption.ratio_all` kernel, and reports strict-window and relaxed roots separately. Strict-window response rows are the
primary acceptance data for the specific ultrasonic analysis described by the user.

Both tools can be narrowed with command-line arguments. `Nk` refers to the superconducting
HHM midpoint k-grid. The ultrasonic DOS resolution and response quadrature are controlled
separately by `--nk-dos` and `--nell`. Outputs are benchmark artifacts, not production inputs.

## Evidence and interpretation

**VERIFIED:** the optional path calculates `Tc` from the HHM normal-state number equation
and the linearized gap determinant; nonlinear continuation solves all three finite-T
unknowns and independently checks the coupled-equation residual. The tests compare the
linearized root to the existing Fig. 8 benchmark implementation, compare its T=0 root to the
existing HHM equations, verify a near-`Tc` gap-direction check, and distinguish strict from
relaxed band policies. The adapter passes its state into the unchanged response kernel and
solves `mu_N(T)` at the selected filling.

**VERIFIED, LIMITED GRID:** an exploratory point at `U=90t, n=1.875, Delta_t=4.5t,
V_SO=0, Nk=80` gave determinant `Tc=0.3273299275` and legacy 1%-gap proxy
`Tc=0.3271524221` (relative difference about `-0.0542%`). At `V_SO=0.4`, `Nk=80`, the
linearized result agrees with the independently implemented Fig. 8 determinant helper to
within `2e-10`; this point is also part of the automated tests. These are point checks, not
broad grid-convergence results.

**UNVERIFIED:** full quantitative reproduction of Fig. 8's temperature-dependent curves.
The paper's plotted curves are not accompanied here by digitized numerical data for
`D0(T)`, `DS(T)`, `mu(T)`, or `bar-Delta_{k,s}(T)`. The outputs enable an apples-to-apples
visual comparison, but a curve match must not be claimed until reference data are digitized
with uncertainty and compared quantitatively. `D0(0)` and `DS(0)` likewise need published
numerical values or a documented digitization before being presented as exact matches.

**LIKELY / REQUIRES RUN-SPECIFIC CHECK:** increasing `Nk` should reduce midpoint-grid error,
but no convergence conclusion should be generalized until the 80/160/320 exports are run
and examined for each parameter point. Continuation residuals establish equation solving,
not uniqueness; seed, path, and alternative-root sensitivity must be assessed from the
candidate diagnostics before relying on a selected branch.

**NOT ESTABLISHED:** thermodynamic stability. The solver's T=0 maximum-`Deff` convention
and finite-T nearest-predictor continuation are explicit numerical branch-selection rules,
not stability criteria. An out-of-band chemical potential is retained as a relaxed
diagnostic, not automatically labeled unphysical in general. The ultrasonic strict-window
criterion is a context-specific acceptance policy. No fixed-filling free-energy comparison
has been implemented, so neither a large effective gap nor a small residual is evidence of
thermodynamic preference.

The `legacy_fit` and `self_consistent_HHM` methods are intentionally distinct. Their
chemical potentials and ultrasonic ratios should be compared at the same requested
parameters, filling, response-kernel settings, and temperature fraction. The adapter does
not alter the legacy method or the response formulas.
