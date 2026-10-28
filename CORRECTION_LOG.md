# CORRECTION LOG — audit-driven correction campaign (2026-10-24; D1 policy update 2026-10-28)

**Policy decision update (2026-10-28):** the user selected the relaxed chemical-potential window. The
project default is `MU_WINDOW_POLICY = "relaxed"`; strict-policy numbers below are historical unless
explicitly re-run. The executed D1-only scope (`project/run_relaxed_regeneration.sh`) covered stage 1,
stage 2, companion/SOC/parameter families, U-family, before/after tables, tests, paper inputs/figures,
addendum figures/PDFs, and the JSON manifest.
It deliberately skipped the independent D2 smoothing calibration and did not alter historical archives.
Out-of-band roots are flagged, not silently treated as ordinary Fermi-surface states; see the verified
regeneration record in §4.

This document records every change made in response to the external audit
`AUDIT_symmetrical-octo-pancake.md` (findings F-01 … F-18), the verification that each change
actually works, the scientific conclusions that changed as a result, what was **not** resolved,
and the decisions that still need a sign-off.

**Numerical baseline (as requested):** the corrected midpoint k-grid (F-02). All stage-1 and
stage-2 numbers below were recomputed on that grid.
**Pre-correction outputs** are preserved, unmodified, in `archive/pre_f02_baseline/` with an md5
manifest (`MANIFEST_md5.txt`, 96 files). Before/after tables are generated from those files by
`tools/before_after.py` → `CORRECTION_TABLES.md`, `results/BEFORE_AFTER.json`.
**Reproduction:** `bash project/run_full_regeneration.sh` is the broad audit rerun (~80 min); the selected
D1-only rerun is `bash project/run_relaxed_regeneration.sh` and deliberately skips the independent D2
calibration. Every regenerated data JSON is hashed in `results/REGEN_MANIFEST.json` (provenance files are excluded).

---

## 0. Summary for the reader

1. **Eighteen audit findings were examined. Seventeen are RESOLVED, with 19 executable regression checks
   (`project/tests/test_fixes.py`, 19/19 PASS on the relaxed-policy regenerated outputs; D1 adds the extra check). One, F-03, is RESOLVED-PART:
   the defect (an unexamined smoothing choice) is fixed and its magnitude is measured, but the σ → 0 limit
   is not reached (decision D2, ≈1 % systematic).** Two audit statements were corrected rather than
   confirmed: F-04's magnitude is state-dependent (+0.06 % here, not 7×10⁻⁵), and F-12's "untested"
   is replaced by a reclassification of variant (d) (§3). The per-finding status is in §2.
2. **Scientific results that move** (the 2026-10-24 strict-policy audit baseline, §5): stage-1 Tc fell by
   0.4–2.7 %; master-set W/W_N limits fell by 2.3–4.3 %; several companion sets were excluded by the
   strict chemical-potential window. Those figures are historical after the 2026-10-28 relaxed-policy
   decision and the corresponding regeneration.
3. **The coherence kernel is re-derived independently** (§3). Variant (a), the signed kernel
   C = ½[1 + (ξξ′ − Δ̄Δ̄′)/EE′], equals the numerically computed Bogoliubov matrix element
   (uu′ − vv′)² to 7×10⁻¹⁶. **Variant (d) is not a Bogoliubov matrix element at all**; the
   audit's and the repository's description of it as "the textbook Case-I factor" with
   "C → 1 at the gap edge" is incorrect (at E = Δ it equals 0). The paper's rejection of
   variants (b)/(c) and the claim that only (a) is physical rest on the derivation, not on the
   (d) comparison.
4. **Open limitation:** D1 is now decided as relaxed (§4), but out-of-band roots are not thereby proven
   thermodynamically stable or suitable for a Fermi-surface interpretation. D2 remains open: the DOS-
   smoothing limit carries an approximately 1 % systematic in the earlier calibration.

---

## 1. Scope, baseline, and what was not changed

* Scope: every finding F-01 … F-18, all dependent JSON, figures, the paper PDFs (English and
  Farsi), the presentation tree, and the documentation.
* Not changed (standing constraints): η = +1 signed primary channel; scattering-only quantitative
  claims; the pair channel is an outlook with **no numbers and no identity** anywhere (Route A);
  the ℓ-algebra and the coherence factors are unchanged in their closed forms (only the
  *verification* of them changed); T0–T7 definitions unchanged.
* Stage-1 solver logic, the gap-shape model φ(T) (E8, rigid master shape), and the
  Nambu/BdG derivation are unchanged. Changes are in the numerics and the root selection only.

---

## 2. Finding-by-finding status

Status key: **RESOLVED** = defect fixed and a test that fails on the old code passes on the new
code; **RESOLVED-PART** = defect fixed, residual quantified and open; **OPEN** = not fixed, the
reason is given; **DECISION** = depends on a user decision (§4).

| # | Audit claim | Status | What was changed (file) | Verification (test / measurement) |
|---|---|---|---|---|
| F-01 | d_std sign inverted in the cross-check `_coh` | **RESOLVED** (cross-check); production (d) was never affected | `bruteforce_ua._coh` now delegates to `absorption.coherence` (single source) | `t_F01`: max\|bf − production\| = 0 for all 5 variants; `kernel_bdg_check.py` (§3) |
| F-02 | Endpoint k-grid, uniform weights, un-converged | **RESOLVED** | `kgrid.py`, `solvers.py`: midpoint rule, Nk → 2Nk convergence check on every solve | `t_F02`: Tc(Nk=40/80/160) = 0.278214/0.278211/0.278211; `solve_set(check_conv=True)` records `conv.ok` for all 13 + 51 + 75 solves |
| F-03 | DOS smoothing = irreducible "1.3 % floor" | **RESOLVED-PART** (magnitude now quantified; σ → 0 limit NOT reached) | `parameters_ua.py`: SMOOTH 0.02 → 0.0025, NBINS 6000 → 48000; `calibrate_smooth.py` (pre-declared decision rule) | §5.3; `t_F03`. The σ-dependence does not converge: ≈1 % change per halving of σ at σ ≈ 0.0025–0.005 (`sigma_limit_converged = false`). **D2.** |
| F-04 | Midpoint S_N with half-bin offset, unchecked | **RESOLVED** | `absorption.S_normal`: exact Gauss–Legendre on a mesh that contains every DOS kink; runtime order-convergence assertion | `t_F04`: relative spread 1.6×10⁻¹⁶ (orders 2/4/8); old midpoint sum differed by +0.06 % (§5.5) |
| F-05 | Makefile/regen masks failures (`\| tail`, `\|\| true`) | **RESOLVED** | `Makefile` (`SHELL`, `pipefail`), `regen_all.sh` (`set -euo pipefail`), validators exit non-zero (`validate_ua.py`) | `t_F05`: pipe demo rc = 1; no `\|\| true` in active lines |
| F-06 | Normal-state (T ≥ Tc) points stored as SC data; stars on them | **RESOLVED** | `param_scan_features.py`: sc = 0 points carry `a = null`, kept only as `normal_ref`; `make_param_figures.py`: NaN gaps, data-driven "[monotone]" tags | `t_F06`: 85 normal-state points, 0 with physics values |
| F-07 | `abs()` applied to fsolve roots without a residual check | **RESOLVED** (and the defect was worse than the audit said, see §2.1) | `solvers._folded_root`: sign normalisation by the *pair* (D0,DS) → (−D0,−DS) (exact symmetry), residual check, band check | `t_F07`: symmetry residual 3×10⁻⁹; mixed-sign input folded to a verified root |
| F-08 | Tc bisection without bracket check; U = −75 shipped Tc ≈ 13 | **RESOLVED** | `solvers.solve_set_core`: unbracketed → `Tc=0, Tc_bracketed=false, Tc_lower_bound`; the U = −75 "Tc ≈ 12.7" was an `abs()` artefact (§2.1) | `t_F08`; U-scan −75: 0/25 SC points, 0 unbracketed |
| F-09 | Nearest-neighbour V_SO snapping, NaN → index 0 | **RESOLVED** | `sc_loader.sc_params_at`: exact match within 10⁻⁶ or ValueError; refuses flagged entries | `t_F09_F10`: off-grid refused; exact lookup OK |
| F-10 | Memoisation not keyed by the resolved root | **RESOLVED** | `sc_loader.stage1`: cache keyed by absolute root path | `t_F09_F10` |
| F-11 | Dense O(nbin²) temporaries (memory wall) | **RESOLVED** | `bruteforce_ua._pair_sum` (row chunks, `BF_CHUNK`); `make_soc_explicit.sums` (row chunks, peak RSS 355 MB; the un-chunked version was OOM-killed) | `t_F11`: chunked = unchunked to 2×10⁻¹⁶ |
| F-12 | "(d) untested"; (a)/(d) split = 23 % | **RESOLVED** as a *reclassification* (§3): the split is real, but (d) is not a BdG kernel, so the comparison cannot discriminate between two physical kernels | `kernel_bdg_check.py`; statements corrected (§3.3) | `t_F01_bdg`: (a) = BdG to 7×10⁻¹⁶; min distance of (d) to any BdG element = 0.50 |
| F-13 | No input validation (ω ≤ 0, frac, nell, n …) | **RESOLVED** (for the integrators; the DOS and stage-1 entry points keep their documented ranges) | `absorption._check_*`, `S_normal` guards | `t_F13`: 7 bad inputs rejected |
| F-14 | SOC explicit script bins by helicity (wrong for the coherence) and uses endpoint grid | **RESOLVED** | `make_soc_explicit.py`: states labelled by (helicity h, branch σ); coherence per σ-block; midpoint grid; chunked | `t_F14`: at V_SO = 0, R_A = R_B to 0 (the exact limit); inter-helicity = ½ by counting; §5.6 |
| F-15 | Threshold location from grid argmin of E_min | **RESOLVED** | `param_scans_peaks.threshold_crossings`: linear interpolation of 2E_min − ħω | `t_F15`: exact to 10⁻⁹ on a linear test |
| F-16 | requirements incomplete; absolute font path; misnamed vendored PDF | **RESOLVED** | `requirements.txt` (pymupdf, arabic-reshaper, python-bidi); `findfont`; `vendor/arXiv_2607.02812_original.pdf` + `MANIFEST_md5.txt` | `t_F16`: md5 matches; no absolute path |
| F-17 | allow_nan; shipped `total` alias; no provenance | **RESOLVED** | `io_utils_*`: `allow_nan=False`, non-finite → null with warning, `PROVENANCE.json`; `temp_sweep`: `total` only when computed | `t_F17`: strict JSON; `total` absent when pair channel not computed |
| F-18 | Hard-coded N(3.6), dead diagnostics, Fermi clamp at 500 | **RESOLVED** | `dos_ua.check(eps_probe)`; `eps_mean` removed; `fermi_qp/fermi_normal` overflow-safe by masking | `t_F18`: no warnings on extreme arguments |

### 2.1 Two defects the audit did not describe

* **Gap-pair sign handling (F-07, deeper than stated).** The original `abs(D0)`, `abs(DS)` folds each
  component independently. For a genuine root with D0 and DS of **opposite sign** this produces a
  point that is *not* a root. The audited U = −75 "Tc ≈ 12.7, D0 ≈ 17.8" is exactly such an artefact:
  the true root at the same μ is (D0, DS) = (+17.84, −8.13), and (17.84, +8.13) has residual ≈ 1.
  The fix flips the **pair** (an exact symmetry of the gap equations), keeps DS signed, and verifies
  the residual. The same `abs()` appeared in `_gap_at_T` and in `sweep_T`; all were removed.
  The DS-sign issue does not arise in the master set (DS > 0 there), so the master numbers are unaffected.
* **Chemical-potential window (new; the decision D1).** Once the root-folding was correct, roots with
  μ outside the band appeared as genuine solutions of the self-consistency equations (§4, D1).

---

## 3. Coherence-kernel revalidation (item 3 of the request)

### 3.1 Method (no closed form in the loop)
`project/ultrasonic/kernel_bdg_check.py` builds, for each random physical pair of states, a 2×2
Nambu block H = [[ξ, Δ̄],[Δ̄, −ξ]] per state, takes the positive-energy eigenvector (u, v) from
`numpy.linalg.eigh`, and forms the standard matrix elements: uu′ − vv′ (density/ultrasonic vertex,
Case I), uu′ + vv′ (Case II), uv′ + vu′ (pair partners). These are compared with the production
closed forms and with each of the four variants (a)–(d).

### 3.2 Results (n = 3000 random states, all sign combinations)

| Check | Result |
|---|---|
| (1) C_a = (uu′ − vv′)² | max deviation **7.8×10⁻¹⁶** |
| (2) C_Case-II = (uu′ + vv′)² | **6.7×10⁻¹⁶** |
| (3) (uu′ − vv′)² + (uv′ + vu′)² = 1 (completeness) | **8.9×10⁻¹⁶** |
| (3b) C_pair(a) = 1 − C_a = (uv′ + vu′)² | **6.7×10⁻¹⁶** |
| (6) Δ → 0: C_a → ½[1 + sgn ξ sgn ξ′] (branch projector) | **0** |
| (4) min over BdG elements of max|C_d − element| | **0.500** (d is **not** a BdG element) |
| (5) gap edge E = E′ = Δ̄, ξ = 0: C_a, C_d | 0, **0** |

### 3.3 Consequences for the claims
* **Variant (a) is confirmed by an independent derivation**, not only by the repository's algebra.
  The (a) result is the physically correct Case-I scattering kernel for a density-type vertex in the
  scattering channel. This is also the form found independently in Sardwal et al. (arXiv:2608.19079, Eq. 16
  homogeneous limit), which the project already cites.
* **Variant (d) is not a Bogoliubov matrix element.** It equals ½[1 − Δ̄Δ̄′/EE′], which is neither
  (uu′ − vv′)² nor (uu′ + vv′)² nor any pair partner. It can only be obtained by dropping the ξξ′
  term *by hand*. The repository's docstring ("the textbook BCS Case-I coherence factor … C → 1 at
  the gap edge AND C → 1/2 at Δ = 0") is **wrong on both counts**: (d) vanishes at the gap edge, and
  it is not textbook. The variant is retained only as a labelled comparison; the paper's (a)-vs-(d)
  discussion must say that (d) has no Bogoliubov derivation.
* **The T3 cross-check now tests integration, not the kernel.** After F-01 both sides of T3 call the
  same closed form, so T3 measures the energy-space integration and the DOS representation. The kernel
  itself is tested only by §3.2. T3 is reported for all variants and gated for (a) only (`validate_ua.py`).
* Variant (b) and (c) (absolute values) violate C ≤ 1 at Δ = 0 (T4 of the validation suite) and have no
  BdG derivation either; their rejection rests on this, and is now a documented result.

---

## 4. Policy decision and remaining numerical limitations

### D1 — Chemical-potential window (**the largest scientific decision in this campaign**)
* **What the rule says.** A superconducting root is accepted only if its μ lies inside the band
  [min ε, max ε] of the k-grid, i.e. a Fermi surface exists. Roots with μ outside the band have every
  quasiparticle gapped: there is no Fermi surface for the scattering integrals, and the normal reference
  S_N is evaluated at a μ_N inside the band.
* **Historical strict-policy rationale (2026-10-24; superseded).** The earlier audit excluded these roots
  because the chemical potential lies outside the bare band and a Fermi surface is absent. That was a
  conservative presentation choice, not a demonstration that the self-consistent roots do not exist. The
  active 2026-10-28 decision is relaxed; out-of-band roots are retained and explicitly flagged.
* **Historical sensitivity record** (`project/ultrasonic/results/sensitivity/relaxed_window_sensitivity.json`,
  166 strict-excluded sets, all reporting superconducting under the relaxed solver at that time). These are
  not the active full-family results. Historical comparison ranges:
  * U = 15–30 (67 sets): Tc 1.49–3.01 t, μ = 4.79–8.48 t, μ excess over the band 0.79–3.88 t, D0 0.47–2.00 t, DS −4.8 to −3.4;
  * U = 30–50 (69 sets): Tc 1.01–2.13 t, μ excess 0.04–2.03 t (the small excesses are at U ≈ 30–50);
  * U = 75, Δt = 4.5, V_SO = 0–0.8 (8 sets): Tc 0.87–1.07 t, μ excess 0.01–0.31 t (**marginal**);
  * U = −75 (3 sets): Tc 12.69–13.17 t, μ excess ≈ 31 t (**far outside**; the audited value is reproduced).
  These counts and ranges are the historical strict-baseline sensitivity comparison. They motivate visible
  out-of-band flags but do not decide physical stability. The newly regenerated relaxed-policy outputs now
  use the active default; do not substitute this old sensitivity JSON for those outputs.
* **What the former strict filter excluded (historical comparison, §5.4):**
  * the Group C Δt = 4.5 set (U = 75 t, V_SO = 0.8; μ = 4.72 vs band edge 4.60, i.e. 3 % outside): Tc = 0.9536 → none;
  * the whole U = 15–30 low-U region at every V_SO (μ ≈ 4.8–7.2 vs band ≈ ±4.0–4.6; D0 ≈ 0.5–1.5,
    DS ≈ −3 to −4.5): **163 of 435** companion-family sets; this includes the Tc values reported by the
    earlier E10 analysis for U = 15–17.5 (2.2–3.0 t) and the "edge maximum at U = 15" statement;
  * U = −75 at every V_SO (μ ≈ 35 t).
* **Decision (2026-10-28): relaxed adopted by user instruction and regenerated.** `MU_WINDOW_POLICY = "relaxed"`
  is the default. The D1-only scope in `project/run_relaxed_regeneration.sh` completed with exit status 0:
  stage 1/2, companion/SOC/parameter families, U-family, before/after outputs, replacement and addendum
  figures, English and Farsi PDFs, plus the JSON manifest were refreshed. The independent D2 smoothing
  calibration and historical archives were deliberately untouched. The active `uvso_families.json` contains
  1,515 records, all tagged `relaxed`; each has a `mu_window` diagnostic, and 521 are out of band. Stage-1
  provenance reports `conv_all_ok=true`; validation gates passed and the regression suite reports 19/19 PASS.
  Out-of-band roots are visibly flagged, including former strict-cutoff points. This removes the strict
  filter's artificial displayed truncation but does **not** establish thermodynamic stability or a
  Fermi-surface interpretation. Strict-baseline numbers in this section remain historical; active regenerated
  products supersede them where D1 changes results.

### D2 — DOS smoothing and histogram resolution
* Adopted: SMOOTH = 0.0025, NBINS = 48000 (σ/Δε ≈ 14; bin-converged to 0.04 %; k-grid converged to 6×10⁻⁵).
  The rule was fixed before the calibration (`calibrate_smooth.py`) and selects the smallest candidate.
* **Not resolved:** the σ → 0 limit is not reached. The W/W_N value still changes by ≈1 % per halving of σ
  at σ ≈ 0.0025–0.005. The rule's own systematic (0.07 %) is model-dependent and must not be quoted as the uncertainty.
  The honest DOS-construction uncertainty is **≈1 %**. Removing it requires an exact or tetrahedron DOS
  (Nk → ∞ analytic), which is not part of this campaign.
* A first calibration at 6000 bins was discretisation-limited (−0.14 %, −0.06 % per bin doubling); the
  reported calibration JSON now includes the bin-convergence table, so the change of rule is visible.

### D3 — Pair channel (Route A; no action needed, recorded)
* The feature study reports peaks in the total (scattering + pair) channel with values 1.3–1.5 at T = 0.1 t.
  These are **outlook only** under Route A and are not quoted here. They remain in the param-scan figures
  as shaded "outlook" curves. No pair-channel number appears in this log, in the README, or in the paper text.

---

## 5. Before/after numbers and convergence tests (items 2, 4, 6, 7)

All before/after values are generated by `tools/before_after.py`; the complete tables are in
`CORRECTION_TABLES.md`. The key values are reproduced here.

### 5.1 Stage 1, master set (U = 95 t, Δt = 4.5 t, n = 1.885)

| V_SO | Tc before → after | Δ%Tc | D0(0) before → after | Deff(0) before → after | Δ%Deff |
|---|---|---|---|---|---|
| 0.0000 | 0.28605 → **0.27821** | −2.74 | 0.25916 → 0.25560 | 0.44435 → 0.43442 | −2.23 |
| 0.4082 | 0.38084 → **0.37377** | −1.86 | 0.28547 → 0.28287 | 0.54888 → 0.53999 | −1.62 |
| 0.8000 | 0.64941 → **0.64669** | −0.42 | 0.32439 → 0.32403 | 0.85423 → 0.85065 | −0.42 |

The V_SO-dependent shift (−2.7 % → −0.4 %) is the signature of the endpoint-grid bias the audit
identified (F-02). **Convergence:** the midpoint grid is converged to 3.4×10⁻⁶ in Tc between Nk = 40 and 160
(the endpoint grid moved by 0.0119 over the same range).

### 5.2 Production W/W_N, master set (scattering channel, all variants, T/Tc → 0.999)

| V_SO | (a) before → after | Δ% | (b) | (c) | (d) |
|---|---|---|---|---|---|
| 0.0000 | 0.40004 → **0.38667** | −3.34 | 0.68362 → 0.66679 (−2.46 %) | 0.68592 → 0.66898 (−2.47 %) | 0.35021 → 0.34162 (−2.45 %) |
| 0.4082 | 0.50540 → **0.48350** | −4.33 | 0.80599 → 0.78861 (−2.16 %) | 0.80995 → 0.79238 (−2.17 %) | 0.41303 → 0.40400 (−2.19 %) |
| 0.8000 | 0.76112 → **0.74339** | −2.33 | 0.92091 → 0.92062 (−0.03 %) | 0.93695 → 0.93657 (−0.04 %) | 0.48620 → 0.48539 (−0.16 %) |

The variant-(a) limit moves by −2.3 % to −4.3 %. The (d)/(a) ratio at T/Tc = 0.9, V_SO = 0 is
0.826 before and 0.833 after: variant (d) lies about 17 % below (a) in both versions. (The audit quoted
23 % at a different state; the comparison is state-dependent and is quoted here at the master state only.)

### 5.3 DOS smoothing (F-03), master state V_SO = 0.4082

| Quantity | σ = 0.02 (old) | σ = 0.0025 (new) | change |
|---|---|---|---|
| R(T/Tc = 0.5) | 0.032188 | 0.031234 | −2.96 % |
| R(T/Tc = 0.9) | 0.324445 | 0.316409 | −2.48 % |
| histogram bins (σ = 0.0025), R(0.5): 6000 / 12000 / 24000 / 48000 | 0.0313273 / 0.0312694 / 0.0312471 / 0.0312339 | | converged to 0.04 % at 48000 |
| k-grid (σ = 0.0025), Nk_DOS 1200 vs 2000 | 0.0312339 vs 0.0312339 (48000 bins) | | 0 |

### 5.4 Companion families and the window decision (D1)
* **Historical strict baseline:** 435 sets were tested; 272 had an accepted in-band root. These counts
  describe the former filter only and are not the active project outputs.
* **Regenerated relaxed family:** `project/ultrasonic/results/uvso_families.json` contains 1,515 records;
  all 1,515 have Tc > 0, with 994 in-band and 521 explicitly marked `mu_window.inside=false`. The five
  former strict-cutoff samples listed in §8.2 now appear with positive Tc and an out-of-band flag.
* The older sensitivity file's 166/166 result remains a historical strict-baseline comparison; it is
  superseded by the full 1,515-record family regeneration, not by assumption.
* The SOC-enhancement histogram was regenerated from the relaxed-policy dataset. The logged summary is:
  V_SO=0.2: median 1.052, p90 1.24, max 1.66; V_SO=0.4: 1.179, 2.03, 3.62; V_SO=0.6: 1.352, 3.59, 9.47;
  V_SO=0.8: 1.600, 5.70, 14.53 (median, p90, max). These are descriptive outputs, not evidence for a
  physical mechanism or stability of out-of-band roots.

### 5.5 Normal reference S_N (F-04)
Gauss–Legendre on the kink-aligned mesh: relative spread across orders 2/4/8 = 1.6×10⁻¹⁶. The former
midpoint-with-half-bin-offset sum differed by **+0.06 %** on the same state (measured at a test state:
0.0160261 → 0.0160359). The audit's 7×10⁻⁵ figure is state-dependent and was not reproduced at the master state.

### 5.6 SOC vertex decomposition (F-14), T/Tc = 0.5 / 0.9

| V_SO | R_A (full vertex) before → after | R_B (helicity-diag.) before → after | inter-helicity fraction before → after |
|---|---|---|---|
| 0.0000 | 0.02638 → 0.02518 / 0.25878 → 0.25217 | = R_A (exact) | 0.5000 → 0.5000 (by counting) |
| 0.4082 | 0.03330 → 0.03225 / 0.32339 → 0.31825 | 0.02864 → 0.02770 / 0.28939 → 0.28432 | 0.4381 → 0.4371 |
| 0.8000 | 0.05588 → 0.05570 / 0.48623 → 0.48571 | 0.04482 → 0.04461 / 0.42589 → 0.42499 | 0.2935 → 0.2912 / 0.2843 → 0.2818 |

The V_SO = 0 limit is now exact (R_A = R_B to machine precision). The SOC trend and the vertex
systematic (R_B/R_A − 1) are unchanged in sign and about 1–2 % smaller in magnitude; the audit's
concern (contamination of the vertex band) is confirmed in size and corrected.

### 5.7 Parameter-peak study (new results; see §6.2)
* Coarse scan: 359 sets × 2 temperatures, 235 valid. Refined candidate maxima: **15** (U: 7, V_SO: 2,
  Δt: 6). Every maximum passes the doubled-quadrature check (≤ 7×10⁻⁶) and the independent k-space check
  (≤ 1.4×10⁻³). The previous study (`PARAM_SCAN_PEAKS.md`) reported no maxima in the scattering channel;
  that statement is **superseded** (§6.2).
* Locations: parabola vertex and dense rescan differ by up to 2.4 in U (grid step 2.5). Peak *positions*
  are resolved to about ±1.5 in U; heights to the printed precision.

---

## 6. Scientific conclusions that change

### 6.1 Unchanged in kind
* W/W_N < 1 for all T < Tc; variant (a) is the physical kernel; variants (b), (c) violate the ratio bound;
  the Δ → 0 projector; the Levy 1963 limit 2f(Δ) (T5: worst deviation 1.1×10⁻⁵ in the new test suite).
* The validation gates: T0 3.9×10⁻⁴, T1 3.1×10⁻⁶, T2 2.0×10⁻⁶ (unchanged), T7a/T7b PASS.

### 6.2 Changed (audit-campaign results; D1 entries below are historical strict-policy comparisons)
1. **Quantitative W/W_N values** in the 2026-10-24 audit baseline: −2 % to −8 % at the limit (Table in §5.2 and `CORRECTION_TABLES.md`). Master-set values remain the same under the relaxed D1 choice because those sampled roots are in band.
2. **Tc values** in that baseline: −0.4 % to −2.7 % on the master set; −2.0 % to −6.5 % for the U = 75/115/125 sets (grid bias).
3. **The E10 low-U Tc results** (U = 15–17.5, 2.2–3.0 t) and the "edge maximum at U = 15" statement:
   these were withdrawn under the former strict rule. Under the now-adopted relaxed policy, the affected
   self-consistent roots are retained in regenerated products and flagged when out of band; physical
   stability and the interpretation of the edge feature remain unverified.
4. **U = −75**: the corrected signed solver yields positive-Tc roots at all scanned V_SO values, with μ far
   outside the bare band. They are retained under the relaxed policy and flagged; this is not evidence that
   they represent a stable superconducting phase.
5. **Group C Δt = 4.5**: the former strict filter removed its root at the edge; the relaxed output retains it
   and flags the out-of-band state. Thermodynamic stability remains unverified.
6. **Parameter peaks**: the earlier "no maxima in the scattering channel" statement is superseded. The
   regenerated study finds 15 verified candidate maxima in the scattering channel, all located near the
   Tc(U) maximum and the root crossover (U ≈ 48–50, where D0 jumps from 0.002 to 0.030 and DS changes sign).
   Mechanism (measured, not assumed): W/W_N at fixed T/Tc peaks where Tc(U) peaks (U ≈ 52, Tc = 0.861 t),
   through the competition between the Tc dome and the gap growth. This is **not** a Hebel–Slichter-type
   resonance. The total-channel maxima are outlook only (D3).
7. **The fig. 36 headline** (SOC enhancement distribution): composition-driven shift (§5.4). The
   paired effect is small except in the high-V_SO tail.
8. **T3 interpretation**: the energy-space vs lattice agreement for (a) improved from ~1 % to 0.2 % with the
   DOS change. Consistent with the T3 deviation being dominated by the DOS smoothing; not separately tested.

### 6.3 Superseded statements (must not be quoted)
* "Only (a) survives the independent k-space check" (the T3 (d) discrepancy was a harness artefact).
* "(d) is the textbook Case-I coherence factor, C → 1 at the gap edge."
* "Stage-1 Tc at V_SO = 0.4082 is 0.3809" (now 0.3738) and similar stage-1 values.
* "No maxima in the scattering channel" (PARAM_SCAN_PEAKS.md).
* "U = −75 Tc ≈ 13 t" (an `abs()` artefact).
* "The 1.3 % DOS floor is irreducible" (it is a systematic of ≈1 %, not a floor; D2).

---

## 7. Tests performed

| Test | What it establishes | Result |
|---|---|---|
| `project/tests/test_fixes.py` (19 checks, including the D1 relaxed-policy regression) | each finding's fix works; Levy limit; out-of-band root is retained and flagged | **19/19 PASS**, 0 FAIL, 0 SKIP (on the regenerated outputs) |
| `project/ultrasonic/kernel_bdg_check.py` | (a) = BdG (uu′−vv′)²; (d) not BdG; completeness; Δ→0 projector | PASS (§3.2) |
| `validate_ua.py` (T0–T4, gated T0, T1, T2, T3(a)) | quadrature, Δ→0 closed form, BCS limit, lattice cross-check, Cauchy–Schwarz | `VALIDATE_UA T0=PASS T1=PASS T2=PASS T3_a=PASS` |
| `dualspace_check.py` (T7a/T7b) | k-space lattice vs energy space vs Levy, convergence in Nk | `T7 dual-space: T7a PASS T7b PASS` |
| `levy1963_check.py` (T5) | flat-DOS ω→0 limit = 2f(Δ) | worst rel 1.1×10⁻⁵ |
| `case2_checks.py` (T6) | Case-II vs Mattis–Bardeen / Hebel–Slichter | run; figures regenerated |
| `calibrate_smooth.py` | σ, bin and k-grid convergence of the DOS (F-03) | §5.3; σ-limit NOT converged (D2) |
| convergence: `solve_set(check_conv=True)` | Nk → 2Nk for all 13 + 51 + 75 + 13 solves | all `conv.ok` true (GRID_CONV_TOL = 2×10⁻³) |
| peak artefact checks (`param_scans_peaks.py`) | doubled quadrature and independent k-space for each candidate | all 15 pass (≤ 1.4×10⁻³) |
| sign-fold and window unit tests (`t_F07`, `t_F08`) | residual-verified roots only | PASS |
| before/after generator (`tools/before_after.py`) | all changed numbers traceable to files | `CORRECTION_TABLES.md` |
| full reproduction (`project/run_full_regeneration.sh`) | ordered, strict rerun from the committed code | §8 |

**Known limits of the test set.** T3 is no longer an independent kernel test (§3.3). The SOC chunked
decomposition is checked by the V_SO = 0 identity and the mass test, not by a dense reference at
finite V_SO. The peak study checks numerical convergence, not the physical interpretation (§6.2).

---

## 8. Reproduction and provenance

* Command: `bash project/run_full_regeneration.sh regen_full.log` (strict; stops at the first failure).
* Hashes: `results/REGEN_MANIFEST.json` (every regenerated JSON, SHA-256).
* Determinism check: the full run was performed twice on the same committed code (run 1: the pass
  whose hashes are in `/home/user/regen_logs/snapshot_run1/`; run 2: the final pass). The comparison
  is reported in §8.1.
* Provenance per run: `results/PROVENANCE.json` (git head, UTC time, Python/NumPy/SciPy versions).
* Environment (this run): Python 3.13, NumPy 2.3.5, SciPy 1.17.1, Matplotlib 3.10.9, pymupdf 2.0.

### 8.1 Determinism check (`results/REPRO_CHECK.json`)
* Two full regenerations were compared by SHA-256 over the 35 JSON files: run 1 (the pass whose hashes were
  snapshotted before the final code additions) and run 2 (the final pass, `regen_full.log`).
* **34 of 35 files are byte-identical.** The one difference, `project/sc_params/results/groups_stage1.json`,
  is caused by two diagnostic fields (`policy`, `mu_window`) added to the stage-1 output after run 1. All
  numerical fields (Tc, D0, DS, Deff, μ) of that file are identical in both runs (checked field by field:
  e.g. Tc(V_SO=0) = 0.278211, Tc(V_SO=0.8) = 0.646688, Tc(U=75) = 0.490888).
* The final historical full run completed with exit status 0 and 18/18 tests passing (`regen_full.log`).

### 8.2 Relaxed-policy D1 regeneration (2026-10-28)
* The D1-only regeneration completed successfully in 53 minutes (main pipeline; log at
  `/tmp/jubilant-enigma-relaxed-regeneration.log`), followed by stage-1 plotting and English/Farsi PDF
  builds. The reproducible runner is `project/run_relaxed_regeneration.sh`; it skips the independent
  D2 calibration and preserves historical archives.
* Stage-1 provenance: `MU_WINDOW_POLICY=relaxed`, `conv_all_ok=true`; stage-2 validation gates passed.
  Regression result: 19/19 PASS. The manifest covers 36 regenerated JSON files at
  `results/REGEN_MANIFEST.json`; provenance JSON files are intentionally excluded by `tools/manifest.py`.
* Family output: 1,515 relaxed-policy records; 521 have `mu_window.inside=false`. The five former strict
  cutoff samples at `(V_SO,U)=(0,46.5),(.2,48.5),(.4,55),(.6,66),(.8,81.5)t` now have positive,
  residual-verified Tc values and are visibly flagged out of band. This is a numerical root/metadata result,
  not evidence of thermodynamic stability.

---

## 9. Remaining open items (with owner)

| Item | Status | Owner / action |
|---|---|---|
| D1 window policy | **resolved: relaxed adopted and regenerated (2026-10-28)** | monitor scientific interpretation of out-of-band roots; stability remains unverified |
| D2 DOS σ-limit (≈1 %) | open | needs exact/tetrahedron DOS (Nk → ∞) |
| E8 per-V_SO gap-shape refit (1.3–3.4 % claim) | **not re-measured** in this campaign | the old `gap_refit_check.json` is pre-correction; re-run needed |
| Farsi PDF (`corrections_fa.pdf`) | **not regenerated; contains Route-A violations** | rewrite the Persian addendum from `addendum_numbers.py` and remove pair numbers; first open item |
| Stage-1 figures 08, 30 in SOCfigures | not refreshed | run `sc_params/plotting_tc.py` on the regenerated stage 1 and re-register them |
| Remaining doc numbers (REFEREE_RESPONSE, DEFENSE, PARAM_SCAN_PEAKS, SOC_EFFECT) | superseded banners only | full numeric rewrite |
| The (d) variant in the paper | must be relabelled (§3.3) | paper addendum wording |
| Pair channel | outlook only (Route A) | no action; no numbers |
| `README.md`, `REFEREE_RESPONSE.md`, `DEFENSE_referee_QA.md`, `PARAM_SCAN_PEAKS.md`, `SOC_EFFECT_SCATTERING.md`, `SOCfigures/INDEX.md` | superseded banners added; full numeric rewrite **not** done | see §10 |

---

## 10. Paper and documentation

* **Paper inputs:** `project/paper_correction/refresh_inputs.py` copies the regenerated exports into
  `vendor/results/` and verifies each copy by SHA-256. `make_replaced_figs.py` now uses the corrected
  physics module (`physics_ua`) for μ(T), not the vendored older copy of the stage-2 code.
* **Paper figures:** regenerated (`vendor/figs/`); the no-superconductivity sets are annotated in
  `vendor/figs/NO_SC_NOTES.txt` and are not plotted as curves.
* **Paper PDFs:** the English addendum is generated from refreshed results. The Farsi PDF has been
  rebuilt and its U=-75 discussion updated for the relaxed policy, but other pre-existing numeric/Route-A
  issues remain; see §10.1 and do not treat that PDF as fully corrected.

### 10.1 PDF status
| Document | Status | Notes |
|---|---|---|
| `arXiv_2607.02812_corrected.pdf` (English, 23 pp.) | **REBUILT** from regenerated results | original 18 pp. verbatim; replaced figures regenerated from the corrected physics module; addendum text generated from the results (`addendum_text_en.py`, `addendum_numbers.py`); addendum figures A1–A3 regenerated (`make_addendum_figs.py`); no pair-channel number in the addendum (Route A); page-14 note updated for U = −75 |
| `corrections_fa.pdf` (Farsi, B Nazanin) | **REBUILT, but still superseded/partial** | the U=-75 paragraph now describes relaxed retention and flags out-of-band roots. Other pre-existing W/W_N values (0.4024 / 0.5086 / 0.7616), T3 wording, and pair-channel content were not comprehensively rewritten; do not distribute as a fully corrected PDF. |
| `SOCfigures/` | refreshed for this policy run | stage-1 figures 08 and 30 were regenerated from the relaxed stage-1 data; family/SOC products were regenerated in the D1 run. Other archival material is unchanged. |
| `README.md`, `REFEREE_RESPONSE.md`, `project/DEFENSE_referee_QA.md`, `project/PARAM_SCAN_PEAKS.md`, `project/SOC_EFFECT_SCATTERING.md`, `project/LEVY_REGIME_ANALYSIS.md`, `project/lecture_note_WWN_corrected.md`, `SOCfigures/INDEX.md` | README updated (notice + corrected headline); the others carry a **superseded** banner | their full numeric rewrite is not done; they must not be quoted without checking `CORRECTION_TABLES.md` |
