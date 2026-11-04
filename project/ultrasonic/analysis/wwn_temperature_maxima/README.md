# Temperature-specific maxima of the existing W/W_N versus U families

## Source and method

The requested source curves are the ones in `figures/WWN_vs_U_families_VSO.png`, generated
from `results/uvso_families.json`: `Delta_t=4.5t`, `n=1.885`, `hbar*omega=0.7t`,
`U=15..150t` in `0.5t` increments, `V_SO/t = 0, 0.2, 0.4, 0.6, 0.8`, and
`T/Tc = 0.3, 0.5, 0.7, 0.9, 0.999`. The source figure uses the **relaxed** chemical-potential
window. The strict-policy data from `results/strict_policy/strict_family_data.json` are
analyzed alongside it as the primary admissibility comparison; they are not silently mixed
with relaxed values.

For every `(policy, V_SO, T/Tc)` curve, the script takes the largest finite sampled `W/W_N`
over the existing U grid and records both the maximum and the U at which it occurs. No
interpolation or new response calculation is performed. Ties choose the smallest U. Strict
maxima only use finite strict-admissible samples. The script also exports a distinct
**pointwise temperature envelope**: at each U, take the largest `W/W_N` over the available
temperature samples. This second operation addresses the alternate reading “maximum over
T as a function of U”; it is not the same as the temperature-by-temperature peak extraction.

Reproduce from the repository root:

```bash
python project/ultrasonic/extract_wwn_temperature_maxima.py
```

The input family files and existing plots are read-only. Outputs are isolated in this
folder:

- `temperature_specific_maxima.csv`: peak value, peak U, temperature, SOC, policy, valid-U
  coverage, and root admissibility diagnostics (50 rows including both policies).
- `temperature_upper_envelope.csv`: per-U maximum over the five sampled temperatures.
- `temperature_specific_peak_values_vs_U.png`: one graph of each temperature's peak
  `(U_peak, (W/W_N)_peak)` coordinates; colors identify SOC, solid/open-marked lines are
  strict, dashed/filled-marked lines are relaxed, and marker shapes identify temperature.
- `temperature_envelope_max_vs_U.png`: the alternate pointwise maximum-over-temperature
  curves versus U, all SOC/policy curves on one graph.

## Results

Under the **strict** policy, the sampled peak U locations are essentially temperature
independent for each SOC: `52.5t` (`V_SO=0`), `54.5t` (`0.2t`), `61.0t` (`0.4t`),
`72.0t` (`0.6t`), and `86.5t` (`0.8t`; the `T/Tc=0.999` peak is at `86.0t`). The
corresponding peak values rise with temperature; for example, at `T/Tc=0.999` they range
from `0.97458` at zero SOC to `0.92616` at `V_SO=0.8t` (the intermediate-SOC peak values
are in the CSV).

In the original **relaxed** curves, most temperature-specific maxima sit at the lower edge
of the scan (`U=15t` or within a few grid steps). Every relaxed peak in this extraction is
marked `mu0_inside_peak=False`: these are out-of-band T=0 roots and are diagnostics, not the
primary ultrasonic acceptance results. The strict maxima therefore differ materially from
the largest values visible in the relaxed figure; compare policies rather than treating the
relaxed maxima as physical predictions for this ultrasonic analysis.

The maxima are grid-sampled over the pre-existing `0.5t` U spacing and finite temperature
list. They are not continuous-domain extrema, and the data do not resolve temperatures
between the five existing curves. The underlying family scan used the legacy fitted-gap
state model. It is not a recomputation with the optional self-consistent HHM pathway, nor a
thermodynamic stability ranking.

## Strict-only two-panel publication figure

Generate the professor-requested strict-admissible two-panel figure from the existing
maxima CSV without re-solving or mixing relaxed roots:

```bash
python project/ultrasonic/plot_temperature_peak_analysis.py
```

Outputs are new files; the earlier figures and source data remain unchanged:

- `temperature_peak_analysis_two_panel.png` (360 dpi)
- `temperature_peak_analysis_two_panel.pdf` (vector PDF)
- `temperature_peak_analysis_two_panel.csv` (strict-only, 25 SOC/temperature records)

Panel A is a scatter plot of each temperature curve's actual sampled maximum,
`(U_peak/t, (W/W_N)_peak)`, with exact, unjittered U positions. Color encodes SOC and marker
shape encodes temperature; no points are connected. Panel B connects those exact same 25
maxima against `T/Tc` to show their temperature dependence. The CSV is the shared source for
both panels and includes the strict admissible U interval and a boundary flag. The strict U
samples are contiguous at `0.5t` spacing. All 25 maxima are interior to their admissible U
intervals in the current data (`at_admissible_domain_boundary=False`); none is a boundary
maximum.

Scientifically, the strict peak location shifts upward with SOC: approximately `52.5t`,
`54.5t`, `61t`, `72t`, and `86.5t` for `V_SO/t=0,0.2,0.4,0.6,0.8`; only the last SOC at
`T/Tc=0.999` is `86t`. Peak magnitude rises with temperature at every SOC, approaching
`0.926–0.975` at `T/Tc=0.999`. Panel A directly answers the request by showing each
per-temperature maximum against its maximizing U. Panel B is complementary: it reveals how
those maxima change with temperature without implying that they share one U location.
