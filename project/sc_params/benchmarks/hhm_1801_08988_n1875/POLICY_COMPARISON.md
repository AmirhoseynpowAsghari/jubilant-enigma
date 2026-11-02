# HHM Fig. 8: paired relaxed/strict solver audit

This is an **additive benchmark** for the HHM parameter set in Fig. 8 of Hutchinson, Hirsch, and Marsiglio, [arXiv:1801.08988](https://arxiv.org/abs/1801.08988). It does not change `parameters_tc.py`, the production policy, or the existing Fig. 8 benchmark outputs.

Run from this directory (with project requirements installed):

```bash
python compare_policy_results.py
```

The script writes `policy_comparison.json`, `policy_tc_comparison.csv`, `policy_gap_functions.csv`, and `policy_comparison.png`. It tests the published set

\[
U=90t,\quad \Delta_t=4.5t,\quad n=1.875,\quad V_{SO}/t=0,0.2,0.4,0.6,0.8.
\]

The Fig. 8 benchmark explicitly overrides the production registry values; it does not change the repository's `n=1.885`, `U=95t` defaults.

## What is being compared

The policies are **chemical-potential admissibility rules**, not alternative gap equations:

- `relaxed`: accept residual-verified roots and record whether \(\mu\) is outside the bare band.
- `strict`: only accept roots with
  \[
  -E_{\rm edge}-10^{-8}t\le \mu\le E_{\rm edge}+10^{-8}t,
  \qquad E_{\rm edge}=4t\sqrt{1+\tfrac12(V_{SO}/t)^2}.
  \]

For each policy the audit compares the nonlinear stage-1 solution at \(N_k=160\), its \(N_k=320\) convergence check, and the gap components \(D_0(T),D_s(T),\Delta_{\rm eff}(T)\) at six common temperatures \(T/T_c^{\rm HHM}=0,0.25,0.5,0.75,0.9,0.99\). Gap points are solved and residual-checked at both \(N_k=80\) and 160. The fixed temperatures are fractions of the **paper's linearized** \(T_c\), so the two policies are compared at the same temperatures.

A separate \((U,n)=(46.5t,1.885)\), \(V_{SO}=0\) control is included solely to verify that the two modes actually differ when a root is just outside the band. It is not an HHM Fig. 8 benchmark point.

## Results

| \(V_{SO}/t\) | HHM linearized \(T_c/t\), \(N_k=320\) | Relaxed nonlinear 1% \(T_c/t\), \(N_k=160\) | Strict nonlinear 1% \(T_c/t\), \(N_k=160\) |
|---:|---:|---:|---:|
| 0.0 | 0.327329928 | 0.327154503 | 0.327154503 |
| 0.2 | 0.351847675 | 0.351647799 | 0.351647799 |
| 0.4 | 0.427122438 | 0.426720610 | 0.426720610 |
| 0.6 | 0.553334524 | 0.553185108 | 0.553185108 |
| 0.8 | 0.726188453 | 0.726150722 | 0.726150722 |

At all five Fig. 8 points, every tested \(T=0\) and finite-temperature chemical potential is in band. Consequently, **strict and relaxed are numerically identical** here: the maximum recorded strict–relaxed difference in \(T_c,D_0(0),D_s(0),\Delta_{\rm eff}(0),\mu(0)\), and the sampled gap-function values is zero at the stored precision. All 120 gap samples (two policies × five SOC values × six temperatures × two grids) pass the residual and band checks; the largest residual is \(1.18\times10^{-8}\), below the \(10^{-6}\) root tolerance. The largest absolute 80-to-160 gap-sample change among \(D_0,D_s,\Delta_{\rm eff},\mu\) is \(5.2\times10^{-6}t\). All stage-1 \(N_k=160\to320\) convergence checks pass. The direct linearized determinant's \(N_k=80,160,320\) spread is at most \(1.6\times10^{-15}t\).

The nonlinear solver's \(T_c\) is defined by the **1%-gap threshold**
\[
\Delta_{\rm eff}(T_c)=0.01\,\Delta_{\rm eff}(0),
\]
not by HHM's linearized instability condition. Its maximum difference from the direct linearized \(T_c\) in this table is \(4.02\times10^{-4}t\) (about \(0.0941\%\), at \(V_{SO}=0.4t\)). That small difference is between two Tc definitions; it is **not** a strict-versus-relaxed difference.

The policy control behaves as intended: at \(U=46.5t,n=1.885,V_{SO}=0\), relaxed returns \(T_c=0.99500094t\) with \(\mu=4.00058033t\), which exceeds the \(4t\) band edge by \(5.8033\times10^{-4}t\); strict returns no in-band \(T=0\) root. This demonstrates that the test harness distinguishes the policies when the constraint is active. It does not imply that this control point belongs to HHM Fig. 8.

## Conclusion and confidence

- **VERIFIED — for the tested Fig. 8 parameter set:** strict and relaxed produce the same self-consistent \(T_c\) proxy and sampled gap functions because the band constraint is inactive. This test cannot choose one policy over the other for Fig. 8.
- **VERIFIED — Tc method:** the direct implementation of HHM's linearized determinant is the method that computes the paper's stated \(T_c\) criterion. The solver's relaxed/strict switch only filters chemical potentials; neither switch turns its nonlinear 1%-gap threshold into the HHM determinant criterion.
- **LIKELY — qualitative agreement:** the independently calculated linearized \(T_c\) values and SOC trend are consistent with the Fig. 8 benchmark already documented in `README_FA.md`.
- **UNVERIFIED — exact curve reproduction:** the paper does not provide a digitized table of the plotted gap curves. No point-by-point equality to the author's plotted data is claimed. Neither policy compares free energies among multiple in-band solutions.

## Implementation references

- HHM linearized determinant and number equation: `compare_fig8.py::paper_linearized_tc`.
- Paired policies, gap residual checks, and outputs: `compare_policy_results.py`.
- Nonlinear gap and number residuals: `project/sc_params/gap_equations.py::make_residual`.
- Policy selection, root screening, Tc bisection, and analytic band edges: `project/sc_params/solvers.py::solve_set_core`, `_folded_root`, and `bare_band_window`.
- Existing automated tests for the relaxed out-of-band case and strict rejection: `project/tests/test_fixes.py::t_D1_relaxed_window` and `t_F20_strict_band`.
