# SOCfigures — curated figure collection for the underlying physics

> **SUPERSEDED IN PART (2026-10-24; D1 update 2026-10-28).** The user selected the relaxed chemical-potential
> window. Figures 08, 30, and 31–36 were refreshed; 31–34 show relaxed-policy family data and mark out-of-band
> roots with open gray triangles. In those family plots, red squares mean `D0(0) < 0.05t`. Historical strict-cutoff
> counts and other pre-correction claims below are not current results. Read `CORRECTION_LOG.md` §4 and §8.2
> for verified D1 scope, outputs, and limitations; out-of-band roots are not thereby established as stable phases.


All figures are PNG, 150 dpi, regenerated from code in `project/` (nothing hand-drawn).
Theme: how spin–orbit coupling (V_SO), the correlated-hopping gap structure and the
scattering formalism shape the ultrasonic attenuation ratio W/W_N of the HHM model
(arXiv:1801.08988) as reconstructed in this repository (arXiv:2607.02812, corrected).

## A. Spin–orbit coupling — where it enters and what it does

| file | physics shown | regenerate with |
|---|---|---|
| `00_SOC_bands_rashba.png` | **new** — helicity bands ε_ks = −2t(cos kx+cos ky) ∓ 2V_SO√(sin²kx+sin²ky) along Γ–X–M–Γ for V_SO = 0, 0.4082, 0.8; μ(n=1.885) line and the two Fermi crossings per V_SO: the split Fermi surfaces behind every SOC effect | `ultrasonic/make_soc_teaching_figs.py` |
| `01_SOC_effect_scattering_channel.png` | the SOC effect **inside the scattering channel**: W/W_N vs T/T_c fans out with V_SO; vs V_SO at fixed T/T_c (×1.7–1.9); relative enhancement with its sign change at T/T_c≈0.22; fixed-T ensemble (×39 suppression) | `ultrasonic/make_soc_figures.py` |
| `02_SOC_mechanism_routes.png` | the four internal routes: van-Hove splitting of N(ε); T_c up ×2.27 & Δ⁰₀/T_c down ×0.55; (σ,σ′) branch-weight redistribution; E_min and μ moving | `ultrasonic/make_soc_figures.py` |
| `03_SOC_scan_master.png` | paper-regenerated V_SO scan of W/W_N (master set) at fixed T/T_c | `ultrasonic/plotting_ua.py` |
| `04,05_SOC_scan_U*.png` | the same scan at U = 95t and 125t. (The U=−75t panel was REMOVED from this collection as non-citable, errata E9; it lives in `ultrasonic/figures/non_citable/`.) | `ultrasonic/plotting_ua.py` |
| `07_SOC_groupD_VSO_scan.png` | Group-D composite of the paper (V_SO scan with auxiliary panels) | `ultrasonic/plotting_ua.py` |
| `08_SOC_Tc_vs_VSO_stage1.png` | stage-1: T_c(V_SO) — SOC strengthens pairing (0.286→0.649 t) | `sc_params/main_tc.py` |
| `09_SOC_branch_vs_helicity.png` | σ (branch index of E-space) vs s (helicity of k-space): why the two labelings differ once V_SO ≠ 0, with the (σ,σ′) weights | `ultrasonic/plotting_ua.py` |

## B. Quasiparticle branches, coherence factors, channels

| file | physics shown | regenerate with |
|---|---|---|
| `10_quasiparticle_branch_construction.png` | **new** — (a) E(ε)=√(ξ²+Δ̄²) with ε*, E_min and the two branch solutions ε_σ(E)=ε*+σ√(E²−E_min²)/κ at a sample E; (b) the signed coherence factor C_I(E,E+ħω) on each (σ,σ′) pair: (−,−) constructive, mixed pairs partially cancel | `ultrasonic/make_soc_teaching_figs.py` |
| `11_coherence_factors_per_branch.png` | C_I, C_pair, C_II per branch vs E — the sign structure η=+1 vs η=−1 | `ultrasonic/plotting_ua.py` |
| `12_coherence_variants_comparison.png` | W/W_N for the four coherence prescriptions (a) signed, (b) \|ξ\|\|ξ′\|, (c) paper's, (d) textbook — why (a) is the physical one | `ultrasonic/plotting_ua.py` |
| `13_branch_decomposition_of_W.png` | decomposition of S_sc into (σ,σ′) channel contributions vs T | `ultrasonic/plotting_ua.py` |
| ~~`14_pair_creation_channel_threshold.png`~~ → **non_citable/** | pair channel = OUTLOOK only (no quantitative claim; REPORT App. A2) | `ultrasonic/plotting_ua.py` |
| ~~`15_sum_rule_and_thresholds.png`~~ → **non_citable/** | the printed (W+W_p)/W_N = 1 identity is WITHDRAWN; panel survives as outlook with corrected titles | `ultrasonic/plotting_ua.py` |
| `16_analytic_limits.png` | analytic limits (ω→0 plateau 2f(Δ), Δ→0, T→T_c) against the numerics | `ultrasonic/plotting_ua.py` |

## C. Benchmarks (independent literature)

| file | physics shown | regenerate with |
|---|---|---|
| `17_benchmark_Levy1963.png` | s-wave reduction = Levy, PR 131, 1497 (1963), Eq. (8): α_s/α_n = 2f(ε₀(T)) to 3·10⁻⁶ | `ultrasonic/levy1963_check.py` |
| `18_benchmark_MattisBardeen_HebelSlichter.png` | Case-II kernel = Mattis–Bardeen 1958; Hebel–Slichter peak vs monotone Case I; 2Δ threshold | `ultrasonic/case2_checks.py` |
| `19_dual_space_k_vs_energy.png` | Sardwal-style triangle (arXiv:2608.19079): k-space lattice sum ↔ energy space ↔ Levy closed form | `ultrasonic/dualspace_check.py` |
| `20_lattice_crosscheck.png` | brute-force (k,helicity)² lattice sum vs the E-space pipeline (T3) | `ultrasonic/plotting_ua.py` |

## D. Parameter dependence, peaks, dome physics

| file | physics shown | regenerate with |
|---|---|---|
| ~~`21..23_peaks_vs_{U,VSO,Dt}.png`~~ → **non_citable/** | the 1.3–1.6 'physical total' peaks are WITHDRAWN (pair channel demoted; PARAM_SCAN_PEAKS §2, §8); scattering-only maxima 0.400/0.505/0.761 — use `33_WWN_vs_U_families_VSO.png` for the surviving shallow U≈17.5–22.5 feature (T/T_c≥0.7) | `ultrasonic/make_param_figures.py` |
| ~~`24_peaks_diagnostics_dome_threshold.png`~~ → **non_citable/** | peak diagnostics of the withdrawn combination (historical record only) | `ultrasonic/make_param_figures.py` |
| `25_Dt_scan_groupC.png` | Group-C Δ_t scan of the paper (regenerated) | `ultrasonic/plotting_ua.py` |
| `26..28_U_*.png` | U-family sweeps: production U scan, attractive-branch sweep, U-temperature composite | `ultrasonic_U/` plotting |
| `29_temperature_sweep_master.png` | master W/W_N vs T/T_c for the three V_SO of the paper (the baseline curves) | `ultrasonic/plotting_ua.py` |
| `30_gap_vs_T_stage1.png` | stage-1 gap functions Δ⁰(T), Δˢ(T) vs the fitted φ(T) used by stage-2 | `sc_params/main_tc.py` |

## Reading order (suggested for a first pass)

00 → 10 → 11 → 12 (formalism) → 29 → 01 → 02 (SOC) →
17 → 18 → 19 (benchmarks) → 31 → 32 → 33 → 34 (parameter families, scattering-only).

## non_citable/ (moved here under the Route-A withdrawal)

`14`, `15`, `21`–`24`: panels whose titles/claims depended on the withdrawn
(W+W_p)/W_N = 1 identity or on the withdrawn 1.3–1.6 'physical total' peaks.  Kept as a
historical record; NOT citable.  See `project/ultrasonic/REPORT_WWN_corrected.md` App. A2,
`project/PARAM_SCAN_PEAKS.md` §2/§8 and `REFEREE_RESPONSE.md` §6.

Companion texts: `project/SOC_EFFECT_SCATTERING.md`, `project/PARAM_SCAN_PEAKS.md`,
`project/LEVY_REGIME_ANALYSIS.md`, `project/ultrasonic/REPORT_WWN_corrected.md`,
`project/ultrasonic/DERIVATION_HHM_scattering.md`, `project/lecture_note_WWN_corrected.md`.

## E. Requested family grids (Δ_t = 4.5t, n = 1.885, ħω = 0.7t)

| file | physics shown | regenerate with |
|---|---|---|
| `31_WWN_vs_TTc_and_Delta_Ufamilies.png` | Refreshed relaxed-policy family curves; rows V_SO=0, 0.4, 0.8 and curves U=15…115. Open gray triangles flag out-of-band roots; red squares mean D0(0)<0.05t. These flags do not establish stability. | `ultrasonic/make_delta_families.py` |
| `32_WWN_vs_TTc_and_Delta_VSOfamilies.png` | Refreshed relaxed-policy family curves; rows U=35, 75, 115 and curves V_SO=0…0.9. Open gray triangles flag out-of-band roots; they are not stability markers. | `ultrasonic/make_delta_families.py` |

Data for both: `ultrasonic/results/delta_families.json` (per T/T_c: ratio, Δ⁰(T), E_min, T).

## F. Requested U- and V_SO-family scans (Δ_t = 4.5t, n = 1.885, ħω = 0.7t)

| file | physics shown | regenerate with |
|---|---|---|
| `33_WWN_vs_U_families_VSO.png` | W/W_N vs U=15…150 on the regenerated ΔU=0.5t grid; panels V_SO=0, 0.2, 0.4, 0.6, 0.8. Relaxed-policy out-of-band roots are marked with open gray triangles; red squares separately mean D0(0)<0.05t. Former strict-cutoff samples are retained, not shown as no-SC. | `ultrasonic/make_uvso_families.py` |
| `34_WWN_vs_VSO_families_U.png` | W/W_N vs V_SO=0…0.9 (ΔV_SO=0.025t); panels U=30, 45, 75, 95, 115. Open gray triangles flag relaxed-policy out-of-band roots; these are not indicators of phase stability. | `ultrasonic/make_uvso_families.py` |

Data for both: `ultrasonic/results/uvso_families.json`.

## G. Explicit vs indirect SOC content of the kernel

| file | physics shown | regenerate with |
|---|---|---|
| `35_soc_explicit_vertex.png` | (a) inter-helicity fraction of S_sc vs V_SO (0.50→0.29): the explicit kernel-level SOC signature; (b) isotropic vertex R_A vs strict q→0 helicity-diagonal R_B: SOC trend survives (×1.70 vs ×1.75); (c) signed effect R(V_SO)−R(0) vs T/T_c with sign change; (d) size of the vertex systematic | `ultrasonic/make_soc_explicit.py` |
| `36_SOC_enhancement_histogram.png` | regenerated on relaxed-policy data; the intentional D0(0)<0.05t root-switch subset and Tc=0 sets are excluded, while out-of-band roots are retained (not a stability claim). Current median enhancements for V_SO=0.2/0.4/0.6/0.8 are 1.052/1.179/1.352/1.600. | `ultrasonic/make_soc_enhancement_hist.py` |

Data: `ultrasonic/results/soc_explicit.json`; discussion `project/SOC_EFFECT_SCATTERING.md` §4.
