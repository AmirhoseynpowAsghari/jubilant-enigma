# Power absorbed by a Rashba superconductor — corrected reconstruction

> ## ⚠ Correction notice (2026-10-24) — read first
>
> This repository has been through a full audit-driven correction campaign
> (`AUDIT_symmetrical-octo-pancake.md`, findings F-01…F-18). **The numbers quoted anywhere in this
> README before this notice, and in any document listed as superseded below, are pre-correction.**
> * Full record, per-finding status, decisions and tests: **`CORRECTION_LOG.md`**.
> * Before/after numbers (generated from files): **`CORRECTION_TABLES.md`**, `results/BEFORE_AFTER.json`.
> * Pre-correction outputs, unmodified: `archive/pre_f02_baseline/` (md5 manifest included).
> * Reproduce everything: `bash project/run_full_regeneration.sh`; hashes in `results/REGEN_MANIFEST.json`.
>
> **Corrected headline (master set U=95t, Δt=4.5t, n=1.885, ħω=0.7t, variant (a), scattering channel):**
> lim W/W_N (T/Tc → 0.999) = **0.3867 / 0.4835 / 0.7434** at V_SO = 0 / 0.4082 / 0.8
> (previously 0.4000 / 0.5054 / 0.7611). Tc = **0.2782 / 0.3738 / 0.6467** (previously 0.2861 / 0.3808 / 0.6494).
>
> **Policy update (2026-10-28):** D1 now uses the **relaxed chemical-potential window**. Residual-verified
> roots are retained even when their chemical potential lies outside the bare band; outputs and figures
> must flag `mu_window.inside = false`. Such roots are algebraic mean-field solutions and are not, by this
> choice alone, established as stable Fermi-surface superconductors. Historical strict-policy results are
> documented in `CORRECTION_LOG.md` and must not be mixed with regenerated relaxed-policy outputs. **D2**
> remains open: the DOS smoothing σ→0 limit is not reached (≈1 % systematic). The coherence kernel (a) is
> confirmed by an independent Bogoliubov derivation; variant (d) is **not** a Bogoliubov matrix element.

> **Power absorbed in the presence of an external time-dependent perturbation**
> Ultrasonic attenuation ratio W/W_N for a spin–orbit-coupled Hubbard superconductor,
> recomputed from scratch, validated, and reconciled with arXiv:2607.02812.

This repository is a **drop-in reconstruction** of
[`AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation`](https://github.com/AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation)
(branch `master`): same stage layout, same module names and roles, same result-file
conventions — but every stage-2 quantity is computed with the **validated** pipeline and
every stage-1 number is reproducible with the repository's own solver logic, parameterised
over all parameter sets of the paper.

| | |
|---|---|
| Physics | η = +1 (ultrasound, Case I), scattering channel with signed coherence factor C_I = ½[1+(ξξ′−Δ̄Δ̄′)/(EE′)], energy conservation E′ = E+ħω on quasiparticle energies, explicit branch (σ) sums |
| Master set | U = 95t, Δ_t = 4.5t, n = 1.885, ħω = 0.7t, V_SO = 0 / 0.4082 / 0.8 |
| Headline result (corrected) | W/W_N < 1 for T < Tc (scattering-only quantitative claim); limits at T/Tc→0.999: 0.3867 / 0.4835 / 0.7434 at V_SO = 0 / 0.4082 / 0.8 (superseded: 0.4000 / 0.5054 / 0.7611); pair channel = outlook only (no number, no identity) |
| Defense package | `project/DEFENSE_referee_QA.md`: benchmarks B1-B6 vs accepted literature, 22 referee Q&A (Q1–Q18d), 12 oral-exam Q&A, honest limitations |
| Deliverables | `project/paper_correction/results/arXiv_2607.02812_corrected.pdf` (English, original text verbatim + replaced figures + addendum), `.../corrections_fa.pdf` (Farsi, B Nazanin, typeset equations) |

## Structure

```
.
├── README.md                          ← this file
├── LICENSE
├── SOCfigures/                        ← curated physics-figure collection (29 citable PNGs + INDEX.md; 6 pair-dependent/withdrawn panels in non_citable/):
│                                        SOC effect, branch/coherence structure, benchmarks,
│                                        peaks; reading order in SOCfigures/INDEX.md
└── project/
    ├── lecture_note_WWN_corrected.md  ← derivation, meaning of s vs σ, thresholds, errata
    ├── LEVY_REGIME_ANALYSIS.md        ← where ql<1, τ, ω=v_s q enter (and where they don't)
    ├── sc_params/                     ← STAGE 1 (Tc, Δ⁰, Δˢ, μ)  [repo solver, parameterised]
    │   ├── parameters_tc.py           ← tolerances + registry of the paper's parameter sets
    │   ├── kgrid.py  gap_equations.py  solvers.py  scanning_tc.py  main_tc.py
    │   ├── io_utils_tc.py  plotting_tc.py
    │   ├── results/   tc_vs_vso.json · groups_stage1.json · u_scans_stage1.json
    │   └── figures/   Tc_vs_VSO.png · gap_vs_T_comparison.png
    ├── ultrasonic/                    ← STAGE 2 (corrected W/W_N)
    │   ├── parameters_ua.py  kgrid_ua.py  dos_ua.py  physics_ua.py  sc_loader.py
    │   ├── absorption.py              ← coherence factors + branch-resolved integrals + W_p
    │   ├── bruteforce_ua.py           ← independent (k,helicity)×(k′,helicity) lattice sum
    │   ├── levy1963_check.py          ← s-wave reduction vs Levy, PR 131, 1497 (1963), Eq. (8)
    │   ├── case2_checks.py            ← Mattis-Bardeen 1958 kernels + Hebel-Slichter peak (T6)
    │   ├── dualspace_check.py         ← T7 k-space vs energy-space triangle (Sardwal et al., arXiv:2608.19079)
    │   ├── param_scans_peaks.py       ← coarse (U, V_SO, Delta_t) survey of W/W_N (peak hunt)
    │   ├── param_scan_features.py     ← fixed-T ensembles, peak detection/refinement/artifact checks
    │   ├── make_param_figures.py      ← publication figures of the parameter study
    │   ├── make_soc_figures.py        ← V_SO (spin-orbit) effect in the scattering channel
    │   ├── validate_ua.py             ← tests T0–T4
    │   ├── scanning_ua.py  main_ua.py  io_utils_ua.py  plotting_ua.py
    │   ├── REPORT_WWN_corrected.md    ← full report (derivation, tables, appendix A)
    │   ├── DERIVATION_HHM_scattering.md ← the scattering section derived equation-by-equation
    │   │                                  from the HHM Hamiltonian (arXiv:1801.08988)
    │   ├── results/   temp_sweep_VSO_*.json · vso_scan.json · groups_WWN.json · W_WN_table.csv
    │   └── figures/   (14 figures, gallery below)
    ├── ultrasonic_U/                  ← U-family stage (groups A, B are U scans)
    │   ├── *_u.py modules · main_ua_u.py · plotting_ua_u.py
    │   ├── results/   temp_sweep_U_*.json · groups_WWN_AB.json
    │   └── figures/   W_WN_U_temperature_sweep.png · W_WN_U_attractive_sweep.png · W_WN_vs_U.png
    └── paper_correction/              ← the arXiv:2607.02812 correction package
        ├── build_corrected_pdf.py  build_fa_pdf.py  eqrender.py  make_replaced_figs.py
        ├── COVER_NOTE.md
        ├── assets/BNazanin.ttf
        ├── vendor/    self-contained inputs: src.tar.gz (original source), figs/ (the seven
        │              replaced paper figures, regenerable via make_replaced_figs.py),
        │              eta1/ (addendum figures), eta1_code/ + results/ (stage-2 module & JSONs)
        ├── figures/   copy of vendor/figs (the seven replaced paper figures)
        └── results/   arXiv_2607.02812_corrected.pdf · corrections_fa.pdf
```

## Quick start

`pip install -r requirements.txt` (numpy>=2 required), then `make test` (validation suite),
`make master` + `make figures` (stage-2 + figures, all relative paths).

```bash
cd project/sc_params   && python main_tc.py all       # stage 1: groups + master scan + U scans
cd ../ultrasonic       && python main_ua.py all       # stage 2: T sweeps, group sweeps, V_SO scans
                       && python validate_ua.py       # tests T0–T4
                       && python plotting_ua.py       # figure gallery
cd ../ultrasonic_U     && python main_ua_u.py && python plotting_ua_u.py
cd ../paper_correction && python build_corrected_pdf.py && python build_fa_pdf.py
```

## Figure gallery (stage 2)

| figure | content |
|---|---|
| `ultrasonic/figures/W_WN_temperature_sweep.png` | corrected W/W_N(T/Tc), master set; scattering only vs + pair creation |
| `ultrasonic/figures/W_WN_coherence_variants.png` | the four coherence-factor choices (a)–(d); only (a) survives all limits |
| `ultrasonic/figures/branch_structure_sigma_vs_s.png` | **new** — what σ means: the two band-energy roots ε_σ(E) at fixed E, vs the helicity sum in N(ε) |
| `ultrasonic/figures/coherence_factors_per_branch.png` | **new** — C_I(E) per branch pair (σ,σ′); (b),(c) leave [0,1] as Δ→0 |
| `ultrasonic/figures/caseII_mattis_bardeen_hebel_slichter.png` | **new** — Case-II benchmarks: Hebel-Slichter coherence peak vs Levy monotone Case-I; Mattis-Bardeen threshold at 2Δ |
| `ultrasonic/figures/dual_space_k_vs_energy.png` | **new** — T7 dual-formulation benchmark: k-space lattice sum vs energy space vs Levy 2f(Δ), after Sardwal et al. arXiv:2608.19079 Fig. 3 |
| `ultrasonic/figures/param_scan_U.png`, `param_scan_VSO.png`, `param_scan_Dt.png` | **new** — parameter study: scattering-only ratio monotone in V_SO and Delta_t; shallow scattering-only interior max vs U at U~17.5-22.5 for T/T_c>=0.7 (edge max at U=15 below); previously reported 1.3-1.6 'total' peaks WITHDRAWN (pair channel demoted to outlook); report `project/PARAM_SCAN_PEAKS.md` |
| `ultrasonic/figures/param_scan_diagnostics.png` | **new** — dome T_c(x), pair threshold 2E_min=ħω, Case-II monotone |
| `ultrasonic/figures/soc_effect_scattering.png`, `soc_mechanism.png` | **new** — spin-orbit effect in the scattering channel: ×1.7–1.9 enhancement at fixed T/T_c, ×39 suppression at fixed T, sign change at T/T_c≈0.22; mechanism decomposition (DOS splitting, T_c & Δ⁰/T_c, branch weights) — note `project/SOC_EFFECT_SCATTERING.md` |
| `ultrasonic/figures/levy1963_swave_check.png` | **new** — s-wave reduction (Δˢ=0, V_SO=0, flat DOS, μ=0): ω→0 limit reproduces Levy 1963 Eq. (8), 2f(ε₀(T)), to 3e-6 |
| `ultrasonic/figures/pair_creation_channel.png` | **new** — the pair-creation channel: integrand on the shell E+E'=ħω per branch pair, and W_p/W_N opening exactly at 2E_min = ħω |
| `ultrasonic/figures/sum_rule_and_thresholds.png` | **new** — channels and the two thresholds (2Delta_eff vs 2E_min); pair channel shown as OUTLOOK only, no quantitative claim |
| `ultrasonic/figures/lattice_crosscheck.png` | **new** — energy-space pipeline vs brute-force lattice sum |
| `ultrasonic/figures/branch_decomposition.png` | (σ,σ′) decomposition of the scattering sum |
| `ultrasonic/figures/W_WN_limits.png` | Δ→0 and ω→0 limits vs analytic forms |
| `ultrasonic/figures/W_WN_vs_VSO*.png` | corrected V_SO scans (U = 95t, 125t, −75t) at fixed T/Tc |
| `ultrasonic/figures/W_WN_groups_*.png` | six-panel composites for groups C (Δ_t scan) and D (V_SO scan) |
| `ultrasonic_U/figures/W_WN_U_*.png`, `W_WN_vs_U.png` | U-family sweeps (groups A, B) and saturation value vs U |
| `sc_params/figures/Tc_vs_VSO.png`, `gap_vs_T_comparison.png` | stage-1 scans and gap functions |

## Validation summary

| test | statement | result |
|---|---|---|
| T0 | ℓ-substitution vs E-integral with BOTH Jacobians (E/κℓ)(E′/κℓ′) | ≤3.9e-4 rel. (3.9/1.8/0.8e-4 at ħω=0.05/0.01/0.002) |
| T1 | Δ→0 closed form R(ħω,T); sum rule at Δ=0 | 3.1e-6; 100.0000 % |
| T2 | ω→0 plateau = 2f(Δ) exactly = Levy 1963 Eq. (8) (s-wave reduction) | 2e-6 |
| T3 | brute-force lattice sum, σ = sign(ε−ε*), per variant | (a) 1.0–1.3 %; (b) ≤0.7 %; (c) ≤0.7 %; (d) 1.5–45 % (variant (d) not used for results) |
| T4 | structured-grid pointwise bound 0 ≤ C_I ≤ 1 for (a); ratio bound at Δ=0: (a),(d) → R ≤ 1, (b),(c) → 2R = 1.9002 > 1 (documented variant pathology) | ok |
| T5 | s-wave reduction vs Levy 1963: W/W_N(ω→0) = 2f(ε₀(T)) | 3e-6 |
| T6 | Case-II vertex: branch-summed kernels = Mattis-Bardeen 1958; Hebel-Slichter peak; 2Δ threshold | 8.7e-16 |
| T7 | dual formulation, Sardwal et al. protocol (arXiv:2608.19079): (i) k-space lattice sum ↔ (ii) energy space ↔ (iii) Levy closed form; Nk convergence | (i)-(ii) 0.2–1.3 % (Nk-independent floor); s-wave triangle ≤1.4 % |

## Errata to arXiv:2607.02812 (one line each; full text in the addenda)

E1 coherence-factor combination & energy conservation · E2 all W/W_N values (flat 0.90–1.00) ·
E3 limit at Tc: scattering channel saturates at R(ħω,Tc); W/W_N→1 would require the pair channel (outlook only) · E4 threshold 2Δ_eff vs 2E_min ·
E5 W_N/DOS normalisation, band-edge leak, abs() masking · E6 Δ_t typo (Group B), Group D
auxiliary panels copied from Group A · E7 U=−75t on an unphysical strong-coupling root · E8 rigid master gap shape phi(T) for all V_SO (per-V_SO refit changes W/W_N by 1.3-3.4 %) · E9 U=−75 outputs removed from citable figures. · E10 stage-1 Tc bisection: the single warm-started fsolve in `_gap_at_T` could slide into the trivial root (an exact solution at every T) and cap Tc at the slide temperature (notches at U ≤ 20 in the family figures); fixed by seeded multi-start in `sc_params/solvers.py`; master stage-1 and every set with U ≥ 25 unchanged (max diff 5.6e-9 over the 50-point V_SO scan); figures 31–34 regenerated.
