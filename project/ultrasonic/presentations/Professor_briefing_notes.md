# Presenter notes — U and spin–orbit coupling effects on \(W/W_N\)

**Audience:** Prof.
**Length:** about 7–9 minutes
**Language:** English
**Purpose:** show what the current computation does and does not support; do not present candidate relaxed roots as proven stable states.

## Slide 1 — Title

“I examined whether the calculated attenuation ratio has interior maxima or minima as we vary interaction strength \(U\) and spin–orbit coupling \(V_{SO}\). I have separated the fixed-reduced-temperature and fixed-absolute-temperature questions, then refined the low-\(U\), high-SOC region where candidate peaks appeared. The goal today is to distinguish reproducible numerical features from interpretations that still need validation.”

## Slide 2 — Definition and scan design

“The plotted quantity is \(R=W/W_N=S_{sc}/S_N\). A fixed \(T/T_c\) scan compares systems at the same fraction of their own transition temperature; a fixed absolute \(T\) scan instead changes whether a superconducting solution exists as \(U\) and \(V_{SO}\) change. These ensembles cannot be merged. The broad scan used \(U=15\) to \(150\) in \(2.5t\) steps and \(V_{SO}=0\) to \(0.9t\), plus \(0.4082t\). I then refined \(U=15\) to \(40\) at \(0.5t\) spacing and \(V_{SO}=0.60\) to \(0.90\) at \(0.025t\) spacing for \(T/T_c=0.3\) and \(0.5\).”

## Slide 3 — Broad U scan

“The full-range curves show why the extrema are hard to see on a single scale. At fixed \(T/T_c\), the broad trends dominate: most sampled minima sit at the high-\(U\) boundary, while many sampled maxima sit at the low-\(U\) boundary. At fixed absolute temperature, the superconducting profile ends where \(T_c\) falls to \(T\); the endpoint maximum is therefore a domain boundary, not an interior peak. The small numbers at low temperature should not be interpreted as exact zeros.”

## Slide 4 — Zoomed U profiles and U finite differences

“The zoom makes the low-\(U\) shoulder visible without smoothing or interpolation. At \(T/T_c=0.3\), the sampled maximum moves from around \(U=18\) at \(V_{SO}=0.60\)–\(0.75t\) to about \(U=31\)–\(34t\) from \(0.775t\) to \(0.90t\). That shift is abrupt between \(0.75t\) and \(0.775t\); I do not yet have evidence to assign a mechanism. At \(T/T_c=0.5\), the maximum stays at \(U=15t\) through \(V_{SO}=0.75t\), then becomes interior: \(U=15.5t\) at \(0.775t\), \(16t\) at \(0.80t\), and \(17.5t\) at \(0.90t\). The lower panels show raw adjacent-point differences, not a fitted derivative.”

## Slide 5 — SOC profiles at fixed U

“These panels show adjacent finite differences for five fixed-\(U\) profiles from the existing family, whose SOC grid is spaced by \(0.025t\). The slope changes are more visible than they are in the full-scale \(W/W_N\) plot, especially for \(T/T_c=0.3\). This figure is diagnostic: the family values are not newly recomputed at doubled SOC resolution, so the smaller sign changes remain candidates rather than established physical extrema.”

## Slide 6 — What does ‘cutoff’ mean?

“There are two distinct cutoffs. For the reduced-temperature curves, I can only quote a grid-based onset of an interior maximum. Under the current prominence screen, \(T/T_c=0.5\) first passes at \(V_{SO}=0.775t\), giving the sampled bracket \([0.75,0.775]t\). For \(T/T_c=0.3\), the maximum is interior in the zoom, but it does not pass that screen; therefore I do not claim a robust SOC cutoff there. Separately, at fixed absolute temperature I estimate the superconducting U endpoint from \(T_c(U,V_{SO})=T\). Those are linear crossing estimates between adjacent U-grid points, not exact continuum roots.”

## Slide 7 — Numerical checks

“The targeted computation completed 663 coarse points and 71 doubled-resolution points. Across the doubled Stage-1, chemical-potential, DOS, and absorption grids, the largest absolute change in \(W/W_N\) was about \(9.1\times10^{-7}\) at \(T/T_c=0.3\) and \(1.08\times10^{-5}\) at \(T/T_c=0.5\). The maximum \(T_c\) difference was \(3.33\times10^{-5}\). At 226 shared points with the prior family, the maximum ratio difference was zero. That last check verifies reproduction of overlapping data, not an independent convergence proof; and doubling everything together does not isolate each error source.”

## Slide 8 — Relaxed roots and the chemical-potential diagnostic

“I retained roots under the relaxed policy because rejecting them solely by a window check would hide candidate solutions. However, all 663 points in this targeted low-U region carry `mu_window.inside=False` in the Stage-1/zero-temperature diagnostic. That flag is not itself proof that a root is unphysical or unstable, but it is a serious qualification: the peak candidates are outside the checked band window. In the earlier broad scan, the sampled in-window \(R(U)\) profiles were monotonically decreasing for all tested reduced temperatures and SOC values. The HHM paper benchmark reproduced in the repository validates selected printed values and the band minimum; it does not establish the thermodynamic stability of every out-of-window root.”

## Slide 9 — Takeaway and discussion

“My defensible conclusion is that the relaxed calculation contains numerical interior maxima in a narrow low-U/high-SOC region, with a more persistent sampled feature at \(T/T_c=0.5\). The refined local maxima are stable to the combined resolution increase at the checked points, but all are outside the Stage-1 band-window diagnostic. The U-window prominence result at \(T/T_c=0.3\) is also sensitive to how the profile is windowed, so the earlier \([0.6,0.7]t\) ‘robust’ onset should not be treated as a confirmed cutoff. I propose to retain the relaxed roots as exploratory results and next test branch continuation, thermodynamic/free-energy stability, and finer SOC convergence before assigning physical significance. Does that validation route match your expectations?”

## Key equations

\[
R(U,V_{SO},T)=\frac{W}{W_N}=\frac{S_{sc}}{S_N},\qquad T=fT_c\ \text{(reduced-temperature ensemble)}.
\]

\[
U_c(T,V_{SO})\ \text{is estimated from}\ T_c(U_c,V_{SO})=T.
\]

## Questions to anticipate

**Why use the relaxed policy?**
It is a useful exploratory calculation that avoids imposing the solver’s band-window restriction as though it were already a physical stability criterion. Every root remains tagged by the diagnostic. Relaxed does not mean “all roots are stable.”

**Are the peaks artifacts?**
The checked maxima persisted in their three-point U neighborhoods when the listed grids were doubled, with the ratio changes reported above. That is encouraging numerical evidence for those sampled points, but the whole refined low-U region is outside the diagnostic window, and the check is not a free-energy or branch-stability test.

**Is \(0.775t\) an exact SOC threshold?**
No. It is the first grid point passing the stated prominence screen on a \(0.025t\) SOC grid; the bracket is \([0.75,0.775]t\). It is a sampled bracket only.

**Why did the \(T/T_c=0.3\) cutoff change from the first scan?**
The earlier screen used the full U-profile range, including the distant high-U baseline. In the low-U zoom the prominence drops below the same absolute threshold. Thus that screen is window-dependent here; the conservative conclusion is “interior sampled maximum, no robust cutoff established.”

## References / code provenance

- Hutchinson, Hirsch & Marsiglio, “Enhancement of superconducting \(T_c\) due to the spin-orbit interaction,” arXiv:1801.08988 (2018): https://arxiv.org/abs/1801.08988
- Repository benchmark statement: `project/DEFENSE_referee_QA.md:28`.
- Scan definitions and combined resolutions: `project/ultrasonic/refine_extrema_zoom.py:40–49`.
- Raw ratio evaluation: `project/ultrasonic/refine_extrema_zoom.py:61–65`.
- Local prominence screen: `project/ultrasonic/refine_extrema_zoom.py:124–149`; the underlying extrema helper is `project/ultrasonic/scan_extrema_vs_soc.py:119–140`.
- Doubled-resolution neighborhood comparisons: `project/ultrasonic/refine_extrema_zoom.py:159–209`.
- Raw curves and finite-difference figures: `project/ultrasonic/refine_extrema_zoom.py:220–275`.
- Diagnostic band window: `project/sc_params/solvers.py:96–100`.
- Full results: `project/ultrasonic/results/refined_extrema_zoom.json` and `project/ultrasonic/results/extrema_vs_soc_scan.json`.
