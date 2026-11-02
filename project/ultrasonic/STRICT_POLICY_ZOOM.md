# Strict chemical-potential window: targeted low-U/high-SOC scan

## Run

From `project/ultrasonic`:

```bash
PYTHONPATH=../sc_params python3 strict_policy_zoom.py
```

The script applies `solvers.MU_WINDOW_POLICY = "strict"` locally. It deliberately does **not** change the historical global default or overwrite the relaxed-policy scan. It writes:

- `results/strict_policy/refined_extrema_zoom_strict.json`
- `results/strict_policy/refined_extrema_zoom_strict_summary.csv`
- `figures/strict_policy/strict_policy_zoom.png`

The left panel of the figure shows strict admissibility across the targeted grid; the right panel plots `W/W_N` only where a strict superconducting solution exists. Missing curves mean the strict solver found no admissible superconducting root at those parameter points; they are not zero-valued absorption data.

## Band condition

For the square-lattice Rashba dispersion used here, the continuum bare-band interval is

\[
\left[-4t\sqrt{1+\frac{V_{SO}^2}{2t^2}},\;
       +4t\sqrt{1+\frac{V_{SO}^2}{2t^2}}\right].
\]

Strict root acceptance uses this analytic interval, with a numerical tolerance of `1e-8 t`. This avoids using the narrower extrema of a finite midpoint grid as the physical band edges. Roots rejected by the strict condition are still recorded in `diagnostic_unrestricted` for audit, not used to compute physical `Tc` or `W/W_N` curves.

## Result for the current targeted window

At `n=1.885`, `Delta_t=4.5t`, `15 <= U/t <= 40` (`dU=0.5t`) and `0.60 <= V_SO/t <= 0.90` (`dV_SO=0.025t`), the strict solver found **0 in-band superconducting solutions out of 663 points**. Each point has an unrestricted diagnostic root, and all 663 such roots lie outside the analytic bare-band interval; their excess beyond the nearest edge ranges from approximately `0.958t` to `4.290t`. Thus the relaxed low-U/high-SOC peaks in this window are excluded from the strict physical curves. The strict W/W_N panel is intentionally empty rather than setting rejected points to zero.

This is a statement about the roots found by the current multistart solver and the advisor's strict admissibility rule. It is not an exhaustive theorem that no other in-band branch exists. Also, if strict roots were found at other parameters, this band criterion would not select among multiple in-band branches; that would still require the model's fixed-filling free-energy comparison.

## Files changed and safeguards

- `sc_params/solvers.py`: analytic band-window helper, strict filtering at zero and finite temperature, and unrestricted-root diagnostics when strict filtering rejects every root.
- `ultrasonic/refine_extrema_zoom.py`: preserve the rejected-root diagnostic in scan records.
- `ultrasonic/strict_policy_zoom.py`: separate reproducible strict scan and non-overwriting figure/data outputs.
- The default remains `relaxed` so existing historical outputs retain their provenance. The strict scan sets the policy explicitly and labels all outputs as strict.
