# Cover note — corrected version of arXiv:2607.02812

**File:** `arXiv_2607.02812_corrected.pdf` (24 pages = original 18 pages verbatim + 6 addendum pages)
**Companion files:** `figs/` (seven corrected figures, high resolution), `wwn_eta1/REPORT.md` (full
derivation and validation), `wwn_eta1/results.{json,csv}` (master-set numbers), the scripts
`stage1.py`, `run_groups.py`, `run_scans.py`, `make_paper_figs.py`, `build_pdf.py`.

## What was wrong (summary)

1. **Coherence-factor algebra (E1).** The final W/W_N ratio combined the coherence factor with
   |ξ|-type roots √(E²−Δ²) and the term +η Δ̄Δ̄′/EE′. The correct golden-rule combination for
   η = +1 (ultrasound, Case I), keeping the sign of ξ and imposing E′ = E + ħω on quasiparticle
   energies, is C_I = ½[1 + (ξξ′ − Δ̄Δ̄′)/(EE′)], summed explicitly over the two helicity branches.
2. **All W/W_N values (E2).** The reported flat 0.90–1.00 curves and tables are artefacts of (1) and
   of the normal-state reference (5). Corrected: W/W_N < 1 for every 0 < T < Tc and every parameter
   set, exponentially suppressed at low T, rising monotonically, saturating **below** 1 at Tc⁻
   (e.g. 0.402 / 0.509 / 0.762 at V_SO = 0 / 0.4082 / 0.8 for the master set U=95t, Δ_t=4.5t, n=1.885).
3. **The limit at Tc (E3).** The scattering channel alone returns R(hw,T_c) = 0.400 / 0.505 / 0.761 at T_c-; restoration of W/W_N -> 1 involves the pair-breaking channel, outside the scope of the quantitative result reported here. The coherence factors satisfy C_I + C_pair = 1 per quasiparticle pair, so together the channels carry the total transition weight; because they occupy different sheets (E'-E = hw vs E+E' = hw) with different thermal weights, their sum is not bounded by and does not equal any single normal-state reference. This work therefore reports the scattering channel only; the pair channel is examined in the outlook without a quantitative claim.

4. **The validity cut (E4).** 2Δ_eff > ħω is not the pair-breaking threshold; the true threshold is
   2E_min = ħω with E_min = |Δ⁰ − bμ|/√(1+b²). Inside 2Δ_eff > ħω the pair-breaking channel is
   identically zero, so the scattering-only result is complete in the window used.
5. **Normal-state reference and DOS (E5).** W_N must share the DOS with W, normalised to ∫N dε = 2
   with exact band support; the stage-2 DOS was mis-normalised (0.643) and leaked above the band
   edge; negative W, W_N were masked by abs().
6. **Parameter misstatements (E6, E7).** Group B's "Δ_t = 4" must be 4.5; Group D's auxiliary panels
   quoted Group A's U=−4t numbers; the U=−75t scan sits on an unphysical strong-coupling root
   (μ ≈ 35t, outside the band).

## What was changed in the PDF

* Original text: **untouched, verbatim** (pages 1–18).
* Pages 11–15: the seven wrong figures replaced by recomputations (tagged `[corrected]`), each page
  carrying a red errata note in the free band above the figure.
* Pages 19–24: Errata & Addendum — summary of corrections E1–E7, corrected formalism, corrected
  result tables (all 13 group sets + V_SO scans), validation suite (T0–T4, sum rules, brute-force
  lattice cross-check), reproduction guide, and three master-set detail figures.

## How the replacement numbers were obtained

Stage 1 (Tc, Δ⁰, Δ, μ for every set): verbatim replica of the repository solver
(`stage1.py`); it reproduces **every** Tc and Δ_eff printed in the paper to four decimals, which
also fixed the two parameter misstatements. Stage 2 (W/W_N): the validated energy-space pipeline
(`wwn_eta1/`, tests T0–T4, Δ=0 sum rule 100.0000 %, brute-force lattice agreement 0.7 % near Tc).
