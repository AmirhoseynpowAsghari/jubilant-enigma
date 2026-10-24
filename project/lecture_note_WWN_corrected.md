# Lecture note — the corrected W/W_N for η = +1 (ultrasound, Case I)

> **SUPERSEDED IN PART (2026-10-24).** This document was written before the audit-driven correction
> campaign. Numbers, statements and figure references below are PRE-CORRECTION unless explicitly
> stated otherwise. The current status of every claim is in `CORRECTION_LOG.md` (§2 per-finding status,
> §5 before/after, §6.3 superseded statements). Do not quote figures from this file without checking
> them against `CORRECTION_TABLES.md` and `results/BEFORE_AFTER.json`.


Companion of `project/lecture_note_WWN_fluctuations.md` (original repository note).
Everything here is derived first and computed second; the code lives in
`project/sc_params` (stage 1) and `project/ultrasonic` (stage 2).
The model Hamiltonian is the one of Hutchinson–Hirsch–Marsiglio, arXiv:1801.08988
(hopping + on-site U + correlated hopping Δt + Rashba V_SO); the full
equation-by-equation derivation of the scattering section from that Hamiltonian is in
`project/ultrasonic/DERIVATION_HHM_scattering.md`.
Units: t = ħ = k_B = 1. Master set: U = 95t, Δ_t = 4.5t, n = 1.885, ħω = 0.7t.

## 1  Set-up and the two indices s and σ

**Bands (helicity index s = ±1).**  The Rashba coupling splits the normal state into two
helicity bands,

    ε_s(k) = −2t(cos kx + cos ky) − 2 s V_SO |sin k| ,      s = +1, −1 .

s is a property of the *normal* band structure (the spin texture). It enters the
computation only through the dispersions, the band edge
W_band = 4t√(1 + V_SO²/2t²), and the partial densities of states N_s(ε); the total normal
DOS is N(ε) = Σ_s N_s(ε), normalised to ∫N dε = 2 with exact band support.

**Gap (HHM form).**  Δ̄(ε) = Δ⁰ − (Δˢ/4t) ε, linear in the band energy; quasiparticle energy

    E(ε) = √( (ε−μ)² + Δ̄(ε)² ) .

**Branch index σ = ±1 (the one that causes confusion).**  For a *fixed* quasiparticle
energy E the relation above is a quadratic equation for the band energy ε, hence ε is
two-valued. With

    b = Δˢ/4t,   κ = √(1+b²),   ε* = (μ + bΔ⁰)/(1+b²),   E_min = |Δ⁰ − bμ|/κ,   ℓ = √(E² − E_min²),

the two solutions are

    ε_σ(E) = ε* + σ ℓ/κ ,      σ = +1, −1 ,

and the integration measure on each branch is ν_σ(E) dE = N(ε_σ(E)) dℓ/κ (the BCS
1/√(E²−E_min²) singularity cancels exactly). Physical meaning: σ = +1 is the
*particle-like* root (ξ = ε−μ > 0 away from the minimum, u² > v²), σ = −1 the *hole-like*
root; they merge at ε* where E = E_min. In the constant-gap limit (Δˢ = 0) this reduces to
the textbook ξ = ±√(E²−Δ²): **σ is exactly the sign of ξ that the paper replaced by |ξ|.**
See `ultrasonic/figures/branch_structure_sigma_vs_s.png`.

Why both indices appear: a quasiparticle state is labelled by (k, s); converting the k-sum
at fixed E into energy integrals gives Σ_s Σ_σ. The coherence factor depends on ξ, Δ̄, E —
functions of ε only — hence on (σ, σ′) but *not* on (s, s′); the helicity enters only
through the weight N_s(ε_σ). That is why the sums over channels in W and W_p run over
(σ, σ′) with Σ_s N_s(ε_σ) = N(ε_σ), and why the brute-force cross-check assigns
σ = sign(ε − ε*) per k-point (`ultrasonic/bruteforce_ua.py`).

## 2  Coherence factor for η = +1

The ultrasound vertex carries u u′ − η v v′ with η = +1 (Case I). Squaring and averaging
gives, with the sign of ξ kept and energy conservation imposed on the *quasiparticle*
energies (E′ = E + ħω),

    C_I(E, ξ; E′, ξ′) = ½ [ 1 + (ξ ξ′ − Δ̄ Δ̄′)/(E E′) ] .        (1)

Checks: 0 ≤ C_I ≤ 1 pointwise; at Δ = 0 it collapses to the branch projector
(1 if σ = σ′, 0 otherwise); the pair-creation element satisfies M_scat² + M_pair² = 1
(verified to 1e-16), i.e. C_pair = 1 − C_I.

## 3  The scattering channel (the in-window answer)

    W = Σ_{σ,σ′} ∫ dE ν_σ(E) ν_σ′(E+ħω) C_I [f(E) − f(E+ħω)] ,          (2)

with per-pair upper limits from the band support,
E ≤ min(E_max(σ), E_max(σ′) − ħω),  E_max(s) = √(κ²(W_band − s ε*)² + E_min²).
At ħω = 0.7t the (+,+) and (−,+) channels vanish identically.
Normal reference of the *same* model:

    W_N = ∫ dE N(E) N(E+ħω) [f_N(E) − f_N(E+ħω)] ,                      (3)

same DOS, μ_N from the normal-state number equation. W/W_N = S_sc/S_N (prefactors cancel).

## 4  Pair creation, thresholds, sum rule

    W_p = Σ_{σ,σ′} ∫ dE ν_σ(E) ν_σ′(ħω−E) (1−C_I) [1 − f(E) − f(ħω−E)] ,  (4)

non-zero only for ħω > 2 E_min(T) — the *true* kinematic threshold. The cut used in the
paper, 2Δ_eff(T) > ħω, is stricter and is **not** a threshold: inside it W_p ≡ 0, so the
scattering-only result (2) is complete there. Above the true threshold the sum rule

    (W + W_p)/W_N = 1                                                    (5)

is a Delta->0 limit only (T1): the two channels live on different sheets (E'-E=hw vs
E+E'=hw) with different thermal weights, so no finite-Delta identity with a single
normal-state reference exists (withdrawn, REPORT App. A2).  The scattering channel alone returns R(hw,T_c) = 0.400 / 0.505 / 0.761 at T_c-; restoration of W/W_N -> 1 involves the pair-breaking channel, outside the scope of the quantitative result reported here.
the scattering channel alone saturates at the band-structure value R(ħω,Tc) < 1,

    R(ħω,T) = (2/ħω) ∫₀^∞ [f(E) − f(E+ħω)] dE      (Δ→0 closed form).    (6)

## 5  Results (validated pipeline `ultrasonic/`)

* W/W_N < 1 for every 0 < T < Tc and every parameter set; exponential suppression at low
  T, monotone rise, saturation below 1 (tables in `ultrasonic/results/W_WN_table.csv`).
* Master set limits: 0.40004 / 0.50540 / 0.76112 at V_SO = 0 / 0.4082 / 0.8; flat-DOS
  R(0.7,Tc) = 0.49872 / 0.59351 / 0.74264 confirms the saturation values.
* Variants |ξ||ξ′| (b) and the paper's final combination (c) exceed 1 as Δ→0 (up to 1.97):
  unphysical for a scattering-only channel; the textbook Case-I form without ξξ′ (d) has a
  different limits: (d) gives f(Δ), half of Levy's 2f(Δ), which the signed factor reproduces exactly (T5).
* Thresholds (master): true 2E_min = ħω at T/Tc = 0.820 / 0.912 / 0.975; paper cut
  2Δ_eff = ħω at 0.740 / 0.840 / 0.936.

## 6  Validation suite (`ultrasonic/validate_ua.py`)

T0 quadrature vs adaptive (both Jacobians): <=3.9e-4. T1 Δ→0 closed form (6): 3.1e-6; two-channel sum = S_N at Δ=0:
100.0000 % (a Δ→0 limit, not a finite-Δ identity). T2/T5 ω→0 plateau = 2f(Δ) = Levy 1963 Eq. (8) to 2e-6. T3 independent brute-force lattice
sum over (k,helicity)×(k′,helicity) with σ = sign(ε−ε*): 0.4–1.6 % for 0.5 ≤ T/Tc ≤ 0.999.
T4 0 ≤ C_I ≤ 1 pointwise. T5 s-wave reduction reproduces Levy 1963 Eq. (8),
alpha_s/alpha_n = 2f(eps_0(T)), to 3e-6. T6 with the Case-II vertex (eta=-1) the same code
reproduces the Mattis-Bardeen 1958 kernels to 8.7e-16 and the Hebel-Slichter coherence peak.
Benchmarks B1-B6 and a referee Q&A are collected in `project/DEFENSE_referee_QA.md`;
the regime question (ql<1, tau, omega=v_s q) in the Levy benchmark is answered rigorously in
`project/LEVY_REGIME_ANALYSIS.md` (test T5b).
T7 (`dualspace_check.py`) solves the SAME model directly in k-space — explicit
(k,helicity)x(k',helicity') lattice double sum with a Gaussian delta and Richardson
extrapolation — following the three-way benchmark protocol of Sardwal-Kreisel-Mukherjee
(arXiv:2608.19079, Fig. 3: real space vs k space vs analytic): agreement with the
energy-space pipeline is 0.2-1.3 % at production hbar*omega=0.7t (independent of Nk), and
in the s-wave omega->0 limit the k-space/energy-space/Levy triangle closes at <=1.4 %
(figures/dual_space_k_vs_energy.png).

## 7  Errata to arXiv:2607.02812 (full text in `paper_correction/results/*.pdf`)

E1 coherence-factor combination and energy conservation (eq. (1) above);
E2 all printed/plotted W/W_N (flat 0.90–1.00) wrong;
E3 limit at Tc (scattering channel saturates at R(ħω,Tc); the paper's jump to 1 would need the pair channel — outlook only);
E4 threshold 2Δ_eff vs 2E_min;
E5 W_N/DOS normalisation, band-edge leak, abs() masking;
E6 Group B "Δ_t = 4" typo (4.5), Group D auxiliary panels copied from Group A;
E7 U=−75t scan sits on an unphysical strong-coupling root (μ ≈ 35t, Tc ≈ 12.7t).
