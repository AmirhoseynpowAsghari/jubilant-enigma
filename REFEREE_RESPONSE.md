# REFEREE_RESPONSE.md — replication audit of the corrected repository
(commit bbfaad1 mirror; manuscript arXiv:2607.02812 v1)

Every number below was re-measured in this workspace (scratch clone for destructive
commands; committed JSON read directly; independent implementations where required).
`PASS` = claim reproduces; `FAIL` = claim contradicted by measurement; `FIXED` = defect
found and corrected in this response.

## 1. Claim-by-claim table

| # | claim as written (file) | measured here | verdict |
|---|---|---|---|
| P0.1 | `main_ua.py master` reproduces committed `temp_sweep_VSO_*.json` | a_signed/b_abs/c_paper/d_std/pair bit-identical (maxdiff 0.0, 52 fracs × 3 V_SO) | PASS |
| P0.2 | Tc⁻ limits 0.40004/0.50540/0.76112 | 0.4000428/0.5054035/0.7611174 | PASS |
| P0.3 | stage-1 group B Tc=0.5009/0.1399/0.0861, Δ_eff=0.6272/0.2840/0.2070 | exact to 4 decimals (solve_set, Nk=80) | PASS |
| P0.4 | stage-1 and ℓ-substitution algebra sound | stage-1 exact; ℓ-algebra after T0 fix ≤3.9e-4 (was 0.30–0.35 with the committed single-Jacobian reference) | PASS / FIXED |
| P1 | "(W+W_p)/W_N = 1 at the 1e-4 level above the true threshold" (README:18, COVER_NOTE:20, PDF pp.19–24, plotting_ua.py:184) | max(total−1) above threshold: 0.3487 (V_SO=0, T/Tc=0.84), 0.4054 (0.4082, 0.98), 0.5139 (0.8, 0.98); total at first pair-open point 1.3487 / 1.2578 (0.92) / 1.5139; spike at 0.4082: total(0.96,0.98,0.99)=1.1956, 1.4054, 1.0885 — all reproduce the referee's table exactly | FAIL as written → FIXED (see §2) |
| P1b | quadrature excluded as cause | S_pair at nell=6k/20k/30k/120k: 0.9448278241/…8252/…8253/…8252 (stable ~1e-9) | CONFIRMED |
| P2 | peaks 1.32–1.58 in (S_sc+S_pair)/S_N vs established Case-I physics | peaks reproduce; decision and mechanism in §3 | FIXED (wording + reference) |
| P2b | REPORT §1.3 "0≤C_a≤1 ⇒ W/W_N≤1 guaranteed" | logically invalid (bound on factor ≠ bound on ratio with differently-measured reference); repo's own 1.514 is the counterexample | FAIL → FIXED (REPORT §1.3 + App. A2) |
| P3-T0 | README "2.1e-4"; committed test actual 2.99e-1/3.40e-1/3.50e-1 | reproduced; root cause confirmed (reference integral carried E/ℓ only); with BOTH Jacobians: 3.89e-4/1.81e-4/8.29e-5 (ω=0.05/0.01/0.002), converging as ω→0. Referee's rebuilt residuals (1.37e-2/6.38e-3/2.88e-3) do NOT reproduce with the correct double-Jacobian reference | FAIL (test) → FIXED; README now ≤3.9e-4 |
| P3-T4 | "0.0 % violations" for all four variants while README says (b),(c) break bound | reproduced; sampling is vacuous (|ξξ′|≤EE′ identically). New T4: structured-grid bound residual 0.00e+00 for (a); ratio bound at Δ=0, ħω=0.2, T=1: (a),(d)=0.950082≤1, (b),(c)=1.900165=2R>1 | VACUOUS → REPLACED |
| P3-T6(1) | 8.7e-16 kernel identity | tautology under the code's definitions (branch-summed C_II vs MB kernel are the same algebra); retained as consistency check, now labelled as such | RELABELLED |
| P3-T6(2) | "Hebel–Slichter peak" = σ₁/σ_N | NOT σ₁/σ_N: repo column = Case-II ultrasonic-style ratio at ω=1e-4 (below the δ-smearing scale). ω-scan at T/Tc=0.85: repo 1.7814/3.4740/5.2919/7.1148/8.9378 for ω=1e-2…1e-6; independent clean-limit Mattis–Bardeen with same gap: 2.0470/3.7906/7.1180/15.4708/39.9361 — BOTH diverge as ω→0 (clean limit has no relaxation), so the referee's "true σ₁/σ_N(ω→0) finite" is itself regularisation-dependent; at ω=0.7 my MB gives 0.701/0.817/0.923/1.006 vs repo 0.001/0.055/0.109/0.162 (different observable entirely) | MISLABELLED → RELABELLED qualitative |
| P3-T3 | README "0.4–1.6 %" | per variant: (a) 9.8e-3–1.3e-2, (b) 8.4e-4–7.3e-3, (c) 2.5e-4–6.9e-3, (d) 1.5e-2–4.5e-1 | PARTIAL → README per-variant |
| P4.1 | η/case convention undocumented (paper p.10 −ΔΔ′/EE′ vs integrand +ηΔΔ′/EE′) | table added (REPORT App. A3): printed factor = Case-I signed (correct, used); printed integrand with η=+1 = (uu′+vv′)² = Case II (wrong for ultrasound) | FIXED |
| P4.2 | GAP_FIT hard-coded 3.434/0.553 vs gap_vs_T_fit.json p=3.674,q=0.556 @V_SO=0 | confirmed; per-V_SO refit effect measured: 2.17/1.33 % (V_SO=0), 2.47/1.43 % (0.4082), 3.36/2.11 % (0.8) at T/Tc=0.5/0.9; erratum E8 + code comment added | FIXED (documented, quantified) |
| P4.3 | U=−75t unphysical root still in citable outputs | Tc≈12.7–13.2t, Δ⁰≈17.8t, μ≈35t confirmed; figure moved to `ultrasonic/figures/non_citable/`, removed from SOCfigures, erratum E9 | FIXED |
| P5 | absolute paths; broken quick start; schema drift; no requirements/CI; font licence; README drift | all reproduced (NameError sc_params_at AND band_edge; `main_ua.py all` dropped total/two_Emin/Deff/D0_0/DS_0/hbar_omega/n_target then plotting KeyError 'total'; regen_all.sh absolute paths; README said 31 PNGs/18 Q&A/2.1e-4/0.4–1.6 %) | FIXED (see §4) |

## 2. Sum rule: WITHDRAWN (Route A)

The previously printed W_N^QP sum-rule identity, with its 1e-4-level numerics, is
WITHDRAWN. No implementation of W_N^QP ever existed in the code, and values measured with
the formula as printed were 0.545 / 3.63 / 0.848 (V_SO = 0 / 0.4082 / 0.8 at T/Tc=0.5) and
0.594 / 1.81 / 0.626 (T/Tc=0.9) — never 1. What survives is the exact per-pair statement:

  C_I + C_pair = 1 per quasiparticle pair, so together the channels carry the total
  transition weight; because they occupy different sheets (E'-E = ħω vs E+E' = ħω) with
  different thermal weights ([f(E)-f(E+ħω)] vs [1-f(E)-f(ħω-E)]), their sum is not bounded
  by and does not equal any single normal-state reference.

Consequently this work reports the SCATTERING CHANNEL ONLY as its quantitative result; the
pair channel is examined in the outlook without numbers or identities. At Δ→0 the two
measures coincide and (S_sc+S_pair)/S_N = 100.0000 % (test T1) — a Δ→0 limit, not a
finite-Δ identity.

## 3. Pair-channel decision (Part 2): Route A — demote to outlook

Route B (implement the matched quasiparticle reference W_N^QP and add test T8 for the
identity) was REJECTED: the "identity" would compare three differently-measured objects,
and the referee is right that no such exact relation against a single normal-state
reference exists at finite Δ. Route A: withdraw all pair-channel-dependent claims.
Supporting measurements (all reproduce): S_pair ≡ 0 for 2E_min ≥ ħω (41/52, 45/52, 48/52
temperature points exactly zero at V_SO = 0 / 0.4082 / 0.8; first non-zero at T/Tc =
0.840 / 0.920 / 0.980), so the scattering-only result is COMPLETE inside the paper's own
window — not an approximation. The totals >1 (1.3487/1.4054/1.5139) and the 1.3–1.6 peaks
arose from adding S_pair to the numerator only; they are withdrawn as quantitative
results. Scattering-only maxima over the full sweep: (a) 0.4000/0.5054/0.7611,
(b) 0.6836/0.8060/0.9209, (c) 0.6859/0.8099/0.9370, (d) 0.3502/0.4130/0.4862 — the (a)
values are the headline. Files changed: full list with per-file before→after in §6.

## 4. Reproducibility patch (Part 5)

regen_all.sh rewritten relative-path with sibling-checkout guards; physics_ua.py quick
start fixed (sc_loader.sc_params_at, band_edge import); scanning_ua.temp_sweep now writes
the full committed schema (total, two_Emin, Deff, D0_0, DS_0, hbar_omega, n_target) —
`main_ua.py all` + `plotting_ua.py` verified end-to-end in a fresh clone; requirements.txt
(numpy>=2 for np.trapezoid, scipy, matplotlib); .gitignore; Makefile (test/master/figures);
FONT_NOTE.md for BNazanin.ttf (proprietary, not MIT, delete before redistribution);
README drift corrected (PNG count, Q&A count, T0/T3/T4/sum-rule lines). The PDF builder is
now fully relative: it vendors its inputs (vendor/src.tar.gz, vendor/figs/, vendor/eta1/)
and writes into project/paper_correction/results/ plus a copy next to itself; no absolute
paths remain in any .py/.sh.

## 5. Referee verdict

The corrected W/W_N (η=+1 signed, scattering-only, inside 2E_min>ħω) is defensible as it
stands: it reproduces Levy (3e-6), the Mattis–Bardeen kernel identity, the Sardwal-style
dual-formulation triangle (≤1.3 %, Nk-independent), and its stage-1 is exact; every defect
found by this audit was in "claims about" the computation (sum-rule wording, vacuous or
mislabelled tests, README drift, packaging), none in the primary channel's arithmetic, and
all are now corrected with measured numbers. The >1 totals and the 1.3–1.6 peaks are
WITHDRAWN together with the W_N^QP identity (Route A, §2–§3): the quantitative result is
the scattering channel alone, W/W_N < 1 everywhere inside 2Δ_eff > ħω, maxima 0.400 /
0.505 / 0.761 at V_SO = 0 / 0.4082 / 0.8; the pair channel is outlook only. The single experiment that would settle the remaining physics: a
threshold-crossing measurement of ultrasonic attenuation at ħω ≈ 2Δ (high-frequency
acoustic resonance or monochromatic phonons in a low-gap superconductor), where the
scattering-only and scattering+pair formulations differ in the onset shape; within theory,
the corresponding closure is a finite-q vertex with Rashba spinor overlaps (inter-helicity
weight ~q²v_F²/Δ_SO²), which would also close the 10–25 % vertex systematic of the SOC
magnitude at strong coupling.

## 6. Route-A execution: per-file before → after (commit e81a340 → fix)

Decision statement: **scattering-only quantitative claim; pair channel demoted to outlook;
W_N^QP identity withdrawn.** Stage-1 solver, the ℓ=√(E²−E_min²) algebra, the coherence
factors C_I/C_II and tests T0–T7 are untouched (they never depended on the pair channel).

| file | before | after |
|---|---|---|
| README.md:18 (headline) | "(W+W_p)/W_N = 1 exactly (1e-4 level)… total exceeds 1" | scattering-only headline: W/W_N < 1 inside window; maxima 0.4000/0.5054/0.7611; pair = outlook (no number, no identity) |
| README.md:28,94,99,116 | 36 PNGs; "sharp peaks 1.3–1.6"; sum-rule gallery line; T4 "0≤C_I≤1 pointwise" | 29 citable PNGs + non_citable/; peaks withdrawn (shallow U≈17.5–22.5 scattering-only feature kept, E10-corrected); outlook wording; T4 ratio-bound row; E8/E9 in errata |
| REFEREE_RESPONSE.md §2,§3,§5 | defended the identity and the >1 totals | WITHDRAWN; Route A vs rejected Route B recorded |
| project/paper_correction/COVER_NOTE.md item 3 | "sum rule holds exactly… total tends to 1 at Tc⁻" | verbatim Route-A T_c⁻ sentence + per-pair statement |
| project/paper_correction/build_corrected_pdf.py | W_N^QP identity strings (3 places); absolute /home/user paths | outlook wording (E3, §1 checks, Fig. A3 caption); E8/E9 appended; vendored relative inputs; output → results/ + local copy; PDF rebuilt (md5 differs from 47e6144b…) |
| project/ultrasonic/REPORT_WWN_corrected.md §1.3, App. A2 | "total 1.51 is the counterexample"; A2 = corrected sum rule | §1.3: ratio bound stated for scattering-only variants, >1 combination withdrawn; A2 rewritten as the withdrawal notice (measured 0.545/3.63/0.848, 0.594/1.81/0.626; S_pair≡0 counts; outlook) |
| project/ultrasonic/plotting_ua.py:184–186 | suptitle "above the true threshold it restores W/W_N → 1" | "above the true threshold the pair channel is an OUTLOOK (no quantitative claim)"; right-panel title "OUTLOOK only… not an identity" |
| project/ultrasonic/DERIVATION_HHM_scattering.md:250 | "(W+W_p)/W_N = 1 (1e-4 level, all T) holds exactly" | "C_I + C_pair = 1 per pair; no single-normal-state-reference identity; Δ→0 limit only (T1)" |
| project/ultrasonic/param_scan_features.py docstring | "(R2) … pair channel … TOTAL peaks" | marked WITHDRAWN as quantitative claim; shading = orientation only |
| project/ultrasonic/physics_ua.py __main__ | NameError: NormalDOS | `from dos_ua import NormalDOS` added; quick start runs |
| project/PARAM_SCAN_PEAKS.md §0,§2,§8 | "physical total peaks 1.32–1.58 … robust" | §0 bullet withdrawn; §2 = single honest paragraph (1292/1320 monotone; maxima 0.400/0.505/0.761; 1.3–1.6 withdrawn; shallow U≈17.5–22.5 feature kept, E10-corrected); §8 = Route-A decision |
| project/SOC_EFFECT_SCATTERING.md:53 | cross-ref to physical-total peaks | cross-ref replaced (scattering-only SOC figures unaffected) |
| project/DEFENSE_referee_QA.md:22,69–70,203,258 | B5 "sum rule (W+W_p)/W_N=1"; Q4 "identity holds exactly (1e-4 level)"; Q18c "physical total peaks"; oral Q6 "prove the sum rule" | B5 → C_I+C_pair=1 per pair, no identity; Q4 → outlook sentence; Q18c → withdrawal + surviving shallow U≈17.5–22.5 feature; Q6 → honest answer (per-pair conservation; no finite-Δ identity; T1 is the Δ→0 limit) |
| project/lecture_note_WWN_corrected.md:87,105 | "holds exactly (verified at the 1e-4 level)… W/W_N→1 at Tc⁻ once both channels kept"; "T0 2.1e-4; sum rule (5)" | Δ→0-limit-only wording + T_c⁻ sentence; T0 ≤3.9e-4 (double Jacobian); "not a finite-Δ identity" |
| SOCfigures/INDEX.md + files | 14,15,21–24 citable with identity/peak captions | moved to SOCfigures/non_citable/ (historical record); INDEX rows struck through; reading order updated; non_citable section added |


## 7. Post-audit finding (erratum E10): stage-1 T_c bisection artifact at low U

Found while densifying the family figures (33/34): isolated notches in W/W_N vs U at
U ≤ 20 (deepest at V_SO = 0.4–0.8). Diagnosis, all measured:

* The gaps (D0_0, DS_0, Deff_0) are SMOOTH across every notch; only T_c collapses
  (e.g. U=17.5, V_SO=0.6: Tc = 1.376 between neighbours 2.675 / 2.026).
* A multi-seed census of the gap equation at (17.5, 0.6) finds EXACTLY ONE nontrivial
  root at each T = 1.2…2.5 (Deff 2.457 → 0.504): no branch coexistence, no first-order
  transition — the branch is alive.
* Root cause: the trivial solution (D0,DS)=(0,0) is an exact root of the residual at
  every T. The old `_gap_at_T` ran ONE fsolve from the warm guess with no convergence /
  nontriviality check; at some bisection midpoints it slid into the trivial root
  (demonstrated: warm continuation dies at T=1.70 while the root exists at 2.5), and the
  bisection then capped T_high at the slide temperature forever. T = (T/T_c)·T_c too low
  ⇒ over-strong exponential suppression ⇒ the notch in the ratio.
* Fix (E10): `_gap_at_T` now takes the max-gap over converged multi-starts (scaled copies
  of the T=0 root, same philosophy as `find_T0_solution`), with a full 20-seed escalation
  at midpoints that report "no SC" while the bracket is wide (`sc_params/solvers.py`).
* Invariance proofs: master sets U=95: Tc 0.2861/0.3809/0.6494 unchanged; 50-point master
  V_SO scan max abs diff 5.6e-9; all sets with U ≥ 25 identical (tested across the
  parameter-scan and delta-family grids). Only U ≤ 20 moves (e.g. 1.185→2.183 at
  (15, 0), 1.376→2.533 at (17.5, 0.6)).
* Consequence: figures 31–34 regenerated (33/34 also on the dense grids); PARAM_SCAN §7,
  README gallery, INDEX and DEFENSE Q18c re-worded: with the corrected T_c(15) the ratio
  vs U is edge-maximal (U=15) for T/T_c ≤ 0.5 and only shallowly interior-maximal
  (U ≈ 17.5–22.5, +0.6 % at T/T_c=0.9) for T/T_c ≥ 0.7 — the formerly quoted pronounced
  U* ≈ 20 peak at every T/T_c was partly an E10 artifact.  Nothing else in the repository
  depends on T_c at U ≤ 20 (parameter scans start at U = 50; master/Levy/MB/T7 use U = 95).
