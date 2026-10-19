# Parameter scans of W/W_N(U, V_SO, Δ_t): peaks and non-monotonic features

Scripts: `ultrasonic/param_scans_peaks.py` (coarse fixed-T/Tc survey),
`ultrasonic/param_scan_features.py` (fixed-T ensembles, peak detection, refinement,
artifact checks), `ultrasonic/make_param_figures.py` (figures).
Data: `ultrasonic/results/param_scans_coarse.json`, `param_features.json`,
`param_features_summary.json`, `param_features_peaks.json`, logs `param_scans_log.txt`,
`param_features_log.txt`.  Reproduce: `python param_scans_peaks.py &&
python param_scan_features.py && python make_param_figures.py` (~15 min, 2 cores).

## 0. Honest summary

* **Primary channel (η=+1, signed ξξ′, scattering only — the standing convention of this
  project): NO peaks.**  Of 1320 curves (8 channels × 3 scans × families), 1292 are
  strictly monotone; 15 fail strict monotonicity only inside numerically flat tails
  (|values| < 10⁻⁶, e.g. the exponentially suppressed side); 13 prominent local maxima
  existed only in the now-withdrawn (S_sc+S_pair)/S_N combination.
* **WITHDRAWN (pair-dependent):** the previously reported sharp peaks of height 1.32-1.58
  in (S_sc+S_pair)/S_N arose from adding the pair channel to the numerator only; they are
  not quantitative results of this work (REPORT App. A2).
* **Case-II channel (η=−1 kernel): no resonance peaks** in the scanned ranges at fixed T
  (monotone); the Hebel–Slichter resonance of this model lives in the T-direction at
  fixed ω (already documented, T6), not in these parameter directions.
* Nothing was tuned to manufacture peaks: the peak mechanism below is a consequence of
  the model's own gap structure and of the pair channel that REPORT Appendix A derives.

## 1. Why the primary channel cannot peak (mechanism of the monotone curves)

The quasiparticle gap is E_min = |D0 − b·μ|/κ with b = Δˢ/4t.  The ridge E_min → 0
(D0 = b·μ, i.e. the effective gap Δ̄(ε) = D0 − b·ε vanishing at the Fermi energy) is
where one would naively expect a gapless peak.  The stage-1 solution shows that this
ridge **coincides with the collapse of the superconducting dome**: e.g. at (U=95,
V_SO=0) the ridge is reached at Δ_t ≈ 3.0–3.25 where T_c = 3·10⁻⁴ t, and at
(U≈180–200, Δ_t=4.5) where T_c = 4·10⁻⁴ t.  At fixed T/T_c both Δ(T) and T shrink
there, and the ratio dies smoothly (≈ 2f(Δ/T) with Δ/T finite); at fixed T the state
becomes normal and the ratio saturates at the normal golden-rule value.  Hence every
primary-channel curve is monotone or sigmoidal (dome crossing), see
`figures/param_scan_diagnostics.png`(a,b,d) and panels (a),(b) of the three scan figures.

## 2. Peak statement (corrected)

The primary (scattering-only) channel shows no peak in the temperature direction: 1292 of
1320 scanned curves are strictly monotone (the rest fail monotonicity only in numerically
flat tails), and its maxima over T/T_c are 0.400 / 0.505 / 0.761 at V_SO = 0 / 0.4082 / 0.8.
The 1.3-1.6 features previously reported arose from adding S_pair to the numerator only;
they are withdrawn.  One genuine scattering-only non-monotonicity survives and is retained:
a local maximum vs U at U* ~ 19-21 for V_SO <= 0.2 (heights 0.058/0.311/0.597/0.873/0.999 at
T/T_c = 0.3/0.5/0.7/0.9/0.999), from the competition of the falling Delta0/Tc with the
falling normal-state anchor R(hw, (T/Tc)Tc(U)) - see §7 addendum and
`figures/WWN_vs_U_families_VSO.png`.

## 3. Artifact checks

* Quadrature: doubling nell (20000→40000) changes the peak height by ≤ 2·10⁻⁵.
* Independent k-space evaluation (bruteforce lattice sums, N_k=320, Gaussian-δ widths
  0.08/0.05/0.03 with Richardson for the scattering leg; the new `S_pair_bf` at w=0.005
  for the pair leg — Richardson is invalid there because the pair integrand carries two
  BCS edge singularities ν(E)ν(ω−E), and the unordered-pair factor ½ is included):
  deviations 0.2–3.8 % for 12/13 peaks; the single 15 % outlier is the V_SO*=0.567 point
  where E_min is smallest and the lattice pair sum needs w < 0.005 (its *position* is the
  most stable of all, shift 0.003).
* The peak survives both refinements with sub-grid shifts; it is a property of the
  continuum integrals, not of the discretisation.

## 4. What was explored and where no peaks exist

* Ranges: U ∈ [50, 200] t (step 7.5–10), V_SO ∈ [0, 1.0] t (step 0.1, +0.4082),
  Δ_t ∈ [0.5, 6.0] t (step 0.25–0.5), at n = 1.885 (paper filling), ħω = 0.7 t,
  ensembles T/T_c = 0.5, 0.9 and fixed T = 0.1 t, 0.2 t; channels a_signed (primary),
  total (sc+pair), Case-II.  U = −75 excluded (unphysical strong-coupling root, errata E7).
* No peaks: primary channel in all three directions and both ensembles; Case-II channel
  in all three directions at fixed T; any channel at fixed T/T_c (the wedge is invisible
  there because T scales with T_c).
* Directions likely to show further structure if the study is extended: the T-direction
  at fixed parameters (Hebel–Slichter peak of the Case-II channel, T6); the ω-direction
  at fixed parameters (threshold at 2E_min); the filling n towards the van-Hove point
  n=1; finite τ/ql (Tsuneto-type kernel, outside the present model).

## 5. Figures

* `ultrasonic/figures/param_scan_U.png`, `param_scan_VSO.png`, `param_scan_Dt.png` —
  2×2 each: (a) primary T/T_c=0.5, (b) primary T=0.1t, (c) total T=0.1t with pair-open
  shading (SC region only) and peak stars, (d) total T=0.2t.
* `ultrasonic/figures/param_scan_diagnostics.png` — (a) T_c(U) dome with the T lines,
  (b) E_min(T=0.1t) vs ħω/2 threshold, (c) Case-II monotone, (d) primary T/T_c=0.9.

## 6. Standing convention

The in-range answer (2Δ > ħω, i.e. pair channel closed) remains the scattering-only
signed ratio; the total including S_pair is reported here *only* as the physical
attenuation of parameter regions where the pair channel is open, clearly labelled and
shaded in the figures — never merged into the in-range curves.

## 7. Addendum — extending U below 50: a primary-channel local maximum

The original U grid was [50, 200] t. Extending to U ∈ [15, 150] t
(`ultrasonic/make_uvso_families.py`, figure `WWN_vs_U_families_VSO.png`, data
`results/uvso_families.json`) reveals a **genuine local maximum of the primary
(scattering-only) ratio vs U** at U* ≈ 19–21 t for V_SO = 0–0.2 (and at T/T_c ≥ 0.9 also
for V_SO = 0.4), at every T/T_c:

| V_SO | T/T_c | 0.3 | 0.5 | 0.7 | 0.9 | 0.999 |
|---|---|---|---|---|---|---|
| 0 | W/W_N at U* | 0.058 | 0.311 | 0.597 | 0.873 | 0.999 |
| 0.2 | W/W_N at U* | 0.059 | 0.314 | 0.600 | 0.873 | 0.999 |

Mechanism (two opposing monotone factors at fixed T/T_c): (i) the suppression exponent
Δ⁰₀/T_c falls with U (1.10 at U=15 → 0.41 at U=25 → 0.06–0.11 at U=45–55), lifting the
ratio; (ii) the normal-state anchor R(ħω, T = (T/T_c)·T_c(U)) falls with U because T_c(U)
falls (ω/T grows), pulling it down. Their competition places an interior maximum at
U* ≈ 20. For V_SO ≥ 0.6 the SOC-driven T_c(U) is flat at small U, the two factors no
longer compete, and the curves are monotone (maximum at the U=15 endpoint).
The V_SO direction (U = 30…115, figure `WWN_vs_VSO_families_U.png`) and the Δ_t direction
remain monotone in the primary channel. The §0/§2 statement "no peaks in the primary
channel" therefore holds for the originally scanned range U ≥ 50 and for the V_SO/Δ_t
directions; the low-U maximum is an additional, now documented, feature of the U direction.
Flagged mean-field roots in these scans: U=55 at V_SO=0.4 and 0.6 (Δ⁰ ≪ Δ), U=150 at
V_SO=0 (dome edge, T_c=0.011 t).


## 8. Pair-channel decision (final): outlook only

The coherence factors satisfy C_I + C_pair = 1 per quasiparticle pair, so together the channels carry the total transition weight; because they occupy different sheets (E'-E = hw vs E+E' = hw) with different thermal weights, their sum is not bounded by and does not equal any single normal-state reference. This work therefore reports the scattering channel only; the pair channel is examined in the outlook without a quantitative claim.  All files that carried the withdrawn claims are listed in REFEREE_RESPONSE.md §6
(with per-file before -> after).  Figures SOCfigures/14, 15, 21-24 were moved to
SOCfigures/non_citable/; the peak panels survive there as historical record only.

