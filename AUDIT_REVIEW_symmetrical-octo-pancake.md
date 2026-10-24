# AUDIT_REVIEW — evaluation of `AUDIT_symmetrical-octo-pancake.md` against HEAD `3b08785`

**Reviewer:** workspace agent (Arena). **Method:** every finding re-checked against the
*current* tree by (i) static read of the anchored code, (ii) independent dynamic probes
(kernel identities, Nk/σ/neps perturbation series, midpoint-grid stage-1 solve, JSON
censuses, Makefile exit-status demo). Probe outcomes are quoted per finding. Nothing in
this review modifies code; §6 lists the proposed remediation order.

**Calibration note.** The audit was run against a HEAD that already includes the Route-A
withdrawal, the vendored PDF stage and SOCfigures 36 (it cites them), i.e. essentially the
current tree. Its "verified sound" section (§15) matches our own referee audit
(`REFEREE_RESPONSE.md`) point for point, including the E10-class solver history. The audit
is competent and mostly right; the deviations found below are of three kinds: magnitudes
that are state-dependent, denominators/labels mis-stated, and one claim whose *narrative*
consequence is overstated.

---

## 1. Classification table

| # | Audit sev | Verdict here | Our severity | One-line basis |
|---|---|---|---|---|
| F-01 | S1 | **CONFIRMED** | S1 | `bruteforce_ua._coh` d_std = ½[1+ΔΔ′/EE′] vs `absorption.coherence` = ½[1−ΔΔ′/EE′]; probe: max\|C_abs−C_bf\| = 0 for a/b/c/case2, ≠0 **only** for d_std |
| F-02 | S1 | **CONFIRMED** | S1 | endpoint Nk series 0.293986/0.286054/0.282122 (40/80/160) reproduced exactly; midpoint grid converges by Nk=40 to 0.278211/0.255599/0.434424 (audit's numbers) |
| F-03 | S2 | **CONFIRMED (mechanism); optimum = needs-investigation** | S2 | σ sweep at master state: +1.05 % (0.010→0.020), +1.49 % (0.020→0.040) — audit's relative pattern exactly; absolute "σ*=0.010" rests on their lattice reference |
| F-04 | S2 | **CONFIRMED (mechanism); quoted magnitude not reproduced at master state** | S3 | neps 2e4/2e5/2e6 at master (0.4082, T/Tc=0.9): spread 8.4e-7 non-monotone (audit: 7.2e-5 at their state); half-bin offset present in code |
| F-05 | S2 | **CONFIRMED** | S2 | demo: recipe exiting 1 piped to `tail` → `make` exit 0; `regen_all.sh:9` `|| true` |
| F-06 | S2 | **CONFIRMED (both halves, with corrected denominators)** | S2 | 114 of 242 `temps` sub-rows (63/121 sets) are sc=0 (T>Tc) carrying a = 0.119–0.276; `make_param_figures.py` docstring l.7 still stars peaks in panels (c)/(d); figure is in the README gallery |
| F-07 | S3 | **CONFIRMED (latent)** | S3 | `abs()` at solvers.py:32,72 without residual re-check; masked by max-gap selection today |
| F-08 | S3 | **CONFIRMED** | S3 | escalation loop solvers.py:92-94 has no post-loop bracket check; `u_scans_stage1.json['-75']` ships Tc up to 13.17 t |
| F-09 | S3 | **CONFIRMED** | S3 | `sc_params_at` = bare argmin; NaN→index 0; no tolerance |
| F-10 | S3 | **CONFIRMED** | S4 | `stage1()` memoises without keying on `root` |
| F-11 | S3 | **CONFIRMED (limitation)** | S3 | dense O(nbin²) temporaries in `bruteforce_ua`; consistent with our own OOM history on this 2 GB host |
| F-12 | S3 | **PARTIAL — numbers confirmed, "untested" overstated** | S3 | (a)/(d) = 0.331862/0.255086 (−23.1 %) confirmed; but Levy T5, the Δ→0 projector structure and the Sardwal homogeneous limit still discriminate (a); what dies after F-01 is only the T3-based rhetoric |
| F-13 | S3 | **CONFIRMED** | S3 | no ω/frac/nell/n guards in `absorption.py`; audit's behaviour table plausible and spot-consistent with code |
| F-14 | S3 | **CONFIRMED** | S2 | `make_soc_explicit.py:90` bins by helicity `st["s"]`, bin-mean coherence, endpoint grid l.57 — the defect `bruteforce_ua`'s docstring warns about; contaminates SOCfigures 35/36 numbers at the 1–2.6 % level |
| F-15 | S3 | **CONFIRMED; impact reduced by Route A** | S4 | argmin threshold in `param_scans_peaks.py:203,236`; pair onset is a step (we measured S_pair≡0 in-window + jump at onset); affected figures are non-citable since Route A |
| F-16 | S3 | **CONFIRMED** | S3 | requirements.txt lacks pymupdf; TTF absolute path l.23; `vendor/src.tar.gz` is an 18-page PDF; xref/rect-based figure placement |
| F-17 | S4 | **PARTIAL** | S4 | `allow_nan` default and missing schema: true. `total`-alias code path exists (scanning_ua:60-62) but shipped `groups_WWN.json`, `ultrasonic_U/*` and `temp_sweep_*` at HEAD show total ≠ a_signed or no alias pair → latent, not shipped |
| F-18 | S4 | **CONFIRMED (cosmetic)** | S4 | `dos_ua.py:74` hard-codes N(3.6); `eps_mean` (l.55) unused; fermi clamp |x|=500 (l.101,108-109) |

**False positives: none wholesale.** The audit contains no invented defects; the four
"partial/not-reproduced" entries above are magnitude/denominator/narrative corrections,
not refutations.

---

## 2. Confirmed errors — detail, location, impact

### F-01 (S1) — d_std sign inversion in the independent cross-check
`bruteforce_ua.py:79-81` (`num = D * Dp`) vs `absorption.py:124` (`X=0, D=−Δ̄Δ̄′`).
Probe: on 5 000 random physical states the two implementations agree to 0.0 for
a_signed/b_abs/c_paper/case2_signed and disagree only for d_std (reflection about ½,
audit's `C_abs+C_bf≡1`). **Impact:** T3's d_std column (README row: 1.5e-2–4.5e-1) is a
harness artifact; the defence narrative "only (a) survives the independent k-space check"
loses that pillar (others survive, see F-12); `S_pair_bf` inherits the wrong C for d_std.
`dualspace_check.py` uses a_signed only → T7 unaffected.

### F-02 (S1) — stage-1 endpoint grid with uniform weights
`sc_params/kgrid.py:19-20` + `gap_equations.py:40-43`. Reproduced both columns: endpoint
Nk drift (Tc 0.293986→0.286054→0.282122 at 40/80/160, still moving at 320) and midpoint
convergence (0.278211/0.255599/0.434424 by Nk=40, stable to 1e-6). Stage 2
(`kgrid_ua.py:25`) already uses midpoints → the two stages integrate different measures.
**Impact:** every stage-1 number (Tc, Δ⁰₀, Δˢ₀, Δ_eff, μ) carries an O(Δk), V_SO-dependent
bias (−2.7 %→−0.4 %); headline triple 0.40004/0.50540/0.76112 → ≈0.3889/0.4966/0.7572;
SOC enhancement factors, all family figures, T0–T7 inputs and both corrected PDFs inherit
it. This is the single most consequential finding.

### F-03 (S2) — DOS smoothing width is a live systematic, not an irreducible floor
`parameters_ua.py:29` SMOOTH=0.02, `dos_ua.py:32,51`. Probe: R(σ) at the master state moves
+1.05 % / +1.49 % per σ-doubling around the shipped value — same relative pattern as the
audit. **Impact:** the T3/T7 "1.3 % floor" label is mis-classified; the bias is
parameter-dependent (audit: +0.3…+1.1 % across groups), so it leaks into U/V_SO trends.

### F-05 (S2) — `make test` cannot fail
Makefile pipes every validator to `tail`; no pipefail; `regen_all.sh:9` swallows status
with `|| true`. Demo reproduced (exit 0 with a recipe exiting 1). **Impact:** all
"PIPELINE OK" evidence, ours included, only proves the scripts ran, not that they passed.

### F-06 (S2) — withdrawn/normal-state content still shipped
Census: 114 sc=0 (T>Tc) sub-rows with a = 0.119–0.276 in `results/param_features.json`
(audit's "114 of 121 rows" is our 114 sub-rows / 63 sets — substance identical);
`make_param_figures.py` still draws peak stars in panels (c)/(d) (docstring l.7) and
`figures/param_scan_U.png` (README gallery) still shows them. **Impact:** a reader of the
figures/JSON alone sees retracted numbers and normal-state points labelled as SC data.

### F-07/F-08 (S3) — solver robustness
abs()-repair without residual check (solvers.py:32,72); Tc bisection without bracketing
guarantee (solvers.py:92-94) — the mechanism behind the shipped U=−75 Tc=13.17 t entries in
`u_scans_stage1.json`. Distinct from our E10 (trivial-root slide): E10 fixed *which root*
is found; F-08 is *no failure signal when no bracket exists*.

### F-09/F-10 (S3/S4) — silent plumbing
Nearest-neighbour V_SO snapping incl. NaN→0 (sc_loader.py:27-34); root-blind memoisation
(sc_loader.py:19-24).

### F-11 (S3) — brute-force memory wall
Dense nbin² temporaries; blocks any convergence study of the reference on this host.

### F-13 (S3) — no input validation
ω≤0, frac∉[0,1), nell=0, n_tgt out of range all produce silent wrong numbers or bare
scipy errors.

### F-14 (S2 here) — `make_soc_explicit` branch/helicity confusion
Bins by helicity with bin-mean coherence (l.74-121) + endpoint grid (l.57); its V_SO=0
self-check is degenerate. **Impact on our deliverables:** SOCfigures 35 and 36 and
`SOC_EFFECT_SCATTERING.md` §4 quote R_A/R_B and inter-helicity fractions from this script
at the 1–2.6 % contamination level — the same order as the vertex systematic it measures.
(Panel (b) of 36 reads `uvso_families.json`, not this script, and is unaffected.)

### F-16 (S3) — deliverable stage not reproducible as documented
requirements.txt (no pymupdf), absolute TTF path, `vendor/src.tar.gz` misnomer (18-page
PDF), xref-rect figure placement. README quick-start step 4 fails in a clean env.

### F-17/F-18 (S4) — hygiene
allow_nan writers, no schema/provenance; latent `total` alias; dead diagnostics; fermi
clamp.

---

## 3. False positives / not reproduced / overstated

1. **F-04 magnitude.** At the master state the neps spread is 8.4e-7 (relative 6e-5), not
   7.2e-5; the oscillation is state-dependent (Fermi-window position vs band edges). The
   *mechanism* (unchecked, half-bin-offset quadrature in the denominator of every number)
   stands; the severity is closer to S3/S4 than S2.
2. **F-06 denominators.** "114 of 121 rows … sc=0" conflates sub-rows and sets (114 sub-rows
   / 63 of 121 sets); the a-range 0.119–0.276 vs their 0.126–0.276. Substance unchanged.
3. **F-12 "untested modelling choice".** The 23 % (a)-(d) split is real, but (a) is not
   selected by T3 alone: Levy-1963 reduction (T5), the Δ→0 branch-projector structure, and
   the Sardwal et al. homogeneous limit (their Eq. 16 = our signed C_I) all discriminate in
   favour of (a). After F-01 the correct statement is "the k-space cross-check no longer
   *dis*favours (d)", not "nothing separates them".
4. **F-17 shipped alias.** At HEAD, `groups_WWN.json`, `ultrasonic_U/results/*` and
   `temp_sweep_*` do not carry total≡a_signed (with_pair=True runs); the alias is a latent
   code path, not shipped mislabelled data.

---

## 4. Unclear / needs investigation

1. **F-03 optimal σ and per-set bias table.** Requires calibrating against the
   Richardson-extrapolated lattice reference in *this* environment (and after any F-02
   change, since the reference itself moves). Relative σ-sensitivity confirmed; the
   "σ*=0.010, bias −0.006 %" claim is plausible but unverified here.
2. **F-01 post-fix T3 (d) number.** Audit says +0.47 %; must be re-measured here after the
   sign fix, and the README T3 row + DEFENSE wording restated accordingly.
3. **F-14 post-fix SOC numbers.** Whether ×1.70/×1.75 trends and the 0.500→0.294
   inter-helicity fractions survive the branch-correct binning at the same magnitudes.
4. **F-02 downstream blast radius.** A stage-1 grid fix moves the headline triple, T0–T7
   inputs, Levy/MB comparisons, both PDFs, the Farsi doc, SOC tables and every family
   figure. Scope = full regeneration + erratum (proposed E11) + headline restatement; a
   project-level decision, not a code patch.
5. **F-11 chunking threshold** on this 2 GB host (nbin at which the reference becomes
   runnable for a convergence study).

---

## 5. Interactions with this project's own errata history

* **F-08 vs E10:** same function, different defects. E10 (committed a7e0c1e) fixed the
  trivial-root *slide* via multi-start; F-08's missing *bracket* check remains and is
  independent (it fires when no root exists anywhere above, e.g. U=−75).
* **F-06 vs Route A:** the withdrawal moved SOCfigures 21–24 to `non_citable/`, but the
  starred `param_scan_U.png` remains in the README gallery and `param_features.json`
  still mixes sc=0 rows — the audit correctly catches the residue.
* **F-14 vs SOC deliverables:** SOCfigures 35/36 and `SOC_EFFECT_SCATTERING.md` §4 must be
  re-run from branch-correct code before any of their numbers are cited again.
* **F-02 vs everything:** our E8 (rigid φ(T)) and E10 notes quantify 1.3–3.4 % and local
  Tc effects; F-02's 2–3 % *global* stage-1 bias dominates both and was unknown to us.

---

## 6. Proposed remediation order (nothing applied yet)

1. **F-05** (30 min): `SHELL := /bin/bash`, `SHELLFLAGS := -o pipefail -c`, drop `|| true`;
   validators exit non-zero on FAIL. Restores meaning to every later "pipeline OK".
2. **F-01** (1 h): `num = -D * Dp` in `bruteforce_ua._coh`; delegate kernels to
   `absorption.coherence`; add `C_bf ≡ C_abs` unit test over a fixed state vector; re-run
   T3; restate README/DEFENSE (d)-variant wording (keep Levy/Δ→0/Sardwal grounds).
3. **F-03** (1 h + decision): set SMOOTH=0.010 *after* calibrating against the lattice
   reference here; publish residual as a ±0.8 % systematic instead of "floor".
4. **F-04** (1 h): Gauss–Legendre or analytic piecewise-linear integral for S_N;
   convergence assertion; move neps to parameters_ua.
5. **F-02** (**decision point**, then ~1 day): midpoint/weighted stage-1 grid + Nk
   convergence assertion inside `solve_set`; regenerate groups/tc_vs_vso/u_scans, all
   stage-2 JSONs, figures 00–36, T0–T7, both PDFs + Farsi; erratum E11 + headline
   restatement (0.3889/0.4966/0.7572 or 3-digit wording). Requires prof sign-off because
   every published digit moves by ≤3 %.
6. **Guards** (2 h): F-09 tolerance+finite check, F-13 input validation, F-17
   allow_nan=False + schema + omit `total` unless with_pair, F-10 cache key.
7. **Solver hygiene** (1 h): F-07 residual assert after abs(); F-08 raise/flag when the
   bracket is never established (retire the U=−75 rows or mark them invalid).
8. **F-14** (2 h): branch-σ binning + midpoint grid in make_soc_explicit; re-run 35/36 and
   SOC_EFFECT §4; keep the vertex caveat with refreshed numbers.
9. **F-06** (1 h): drop/NaN sc=0 rows; citable flag; regenerate param_scan_U/VSO/Dt
   panels (c)/(d) without stars; move them out of the README gallery or label withdrawn.
10. **F-16** (30 min): pymupdf in requirements; font via font_manager; rename
    src.tar.gz→source.pdf + hash manifest for REPL.
11. **F-11/F-15/F-18** (cleanup, 2 h): chunked brute force; threshold bisection +
    `threshold_kind:"step"` schema note; dead-diagnostics removal.

**Bottom line:** 14 of 18 findings confirmed as stated; 4 confirmed only in part
(F-04 magnitude, F-06 denominators, F-12 narrative, F-17 shipped-alias). The two S1 items
are real: F-01 invalidates one pillar of the variant-(d) rejection rhetoric (cheap fix),
and F-02 means every stage-1-derived digit in the project carries a 2–3 % grid bias
(expensive fix, needs a documented restatement). No finding contradicts the
scattering-only Route-A claim structure itself — but after F-02 the headline numbers must
be re-quoted.
