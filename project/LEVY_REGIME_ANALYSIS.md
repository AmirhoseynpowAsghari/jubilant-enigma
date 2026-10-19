# The regime question in the Levy (1963) benchmark — a rigorous answer

Addresses: *how is ql < 1 implemented? is τ or ℓ = v_F τ present? is ω = v_s q imposed? does the
calculation reproduce Levy's physical regime or only his limiting expression?*
Companion files: `levy1963_check.py` (tests T5, T5b), `DERIVATION_HHM_scattering.md`,
`DEFENSE_referee_QA.md` (B1, Part IV).

## 0. Honest summary (stated before any derivation)

1. **ql < 1 is NOT implemented in the code.** The code contains no wave vector q, no mean
   free path ℓ, no relaxation time τ, no sound velocity v_s. The ultrasound vertex is taken
   at q → 0 (the density operator ρ_q → ρ, helicity-diagonal), which is the common leading
   order of *all* regimes (ql ≪ 1 and clean q → 0 alike).
2. **No finite τ / ℓ.** The computation is the clean-limit Fermi golden rule (τ → ∞):
   W ∝ Σ |M|² (occupations) δ(E_f − E_i). Levy's τ enters his α_s and α_n only through an
   overall prefactor (ql)²/τ and an energy-independent angular kernel, both of which cancel
   in the *ratio* α_s/α_n — which is why a τ-free calculation can reproduce the ratio.
3. **ω = v_s q is NOT imposed.** ω is a free parameter and q → 0 is taken directly. Levy uses
   q = ω/v_s only to form the dimensionless combination ql = ωτ (v_F/v_s); note that
   ql < 1 ⇔ ωτ < v_s/v_F ≪ 1, i.e. his regime is a statement about τ at given ω — a
   statement a τ-free ratio never needs.
4. **The comparison is a valid benchmark of the universal (q→0, ω→0, kernel-constant) core
   of the Case-I machinery — and of nothing more.** Section 5 proves mathematically that
   Levy's ratio and ours are the same functional L[1] of the occupation factors, reached
   from opposite sides of the ql axis (his: collisional Boltzmann, ql<1; ours: clean golden
   rule, the BCS-type calculation Levy himself cites for ql>1). Section 6 delimits where
   the two prescriptions differ (finite ω, finite q, τ-dependent crossovers).

## 1. Audit: what the code actually contains

Present in `absorption.py` / `levy1963_check.py`: quasiparticle energies E_ks (HHM Eq. 26);
branch DOS ν_σ(E) = N(ε_σ)E/(κℓ) with the ℓ-measure; coherence factor C_I (η=+1) and
1−C_I; thermal factors [f(E)−f(E′)] with **E′ = E + ħω** and [1−f(E)−f(E′)]; normal
reference S_N with μ_N; flat-DOS/BCG s-wave reduction with BCS Δ(T).
Absent (grep-verifiable): τ, ℓ, q, v_s, any Boltzmann collision term, any angular kernel,
any impurity self-energy. The vertex is ρ_{q→0} (DERIVATION §5).

## 2. Levy's calculation in our notation

Boltzmann equation with collision drag (his Eqs. 1–4), shear or longitudinal wave,
q = ω/v_s, Meissner screening of electromagnetic fields, terms of order ωτ neglected:
the first-order distribution is f₁ ∝ (driving)/(1 + iωτ − i q v τ cosθ) · (dE/dε) df₀/dE, and
the cycle-averaged heat production after the angular integration is

    Q  ∝  INT_{ε₀(T)}^∞ dE  (−df₀/dE) · 𝒦(b, ql),        b = ωτ/(ql) = v_s/v_F ,

with 𝒦 a number depending on (ωτ, ql) but **not on E** (his Eqs. 5–6; the [1+ε₀²/ε²]^{1/2}
Jacobian is already absorbed into dE). The normal metal gives the same expression with
ε₀ = 0 and ∫_{−}^{∞} dε (−df₀/dε) = 1. Since α ∝ Q, and the prefactors (m, ρ, τ, (ql)², 𝒦)
are identical in the two states,

    α_s/α_n  =  [ 2 INT_{ε₀}^∞ (−f′) dE ] / [ INT_{−∞}^∞ (−f′) dε ]  =  2 f(ε₀(T)) ,   (Levy 8)

independent of ql, τ, b within the regime. Two structural features to note:
(i) the factor **2** is his observation that "the resulting integral is even in ε": the two
roots ε = ±√(E²−ε₀²) at fixed E — in our language, the branch sum Σ_σ in the s-wave case;
(ii) everything transport-related has cancelled; what survives is purely the
occupation-factor weight of the quasiparticle shell.

## 3. Our calculation in the same notation

Clean golden rule at q→0 (DERIVATION §§5–9), Case-I vertex:

    W = 2πB²ω Σ_{σ,σ′} INT dE ν_σ(E) ν_σ′(E+ω) C_I(σ,σ′) [f(E) − f(E+ω)] ,
    W_N = 2πB²ω INT dε N(ε) N(ε+ω) [f_N(ε) − f_N(ε+ω)] .

As ω→0: [f(E)−f(E+ω)] → ω(−f′); C_I(E,E) = ξ²/E² = ℓ²/E²; and with one measure per dE and
one per dℓ the integrand per branch collapses exactly, ν_σ(E)²C_I(E,E)(dE/dℓ)⁻¹… → N(ε_σ)²,
so that Σ_σ gives

    W/W_N  →  [ 2 INT_Δ^∞ (−f′) dE ] / [ INT_{−∞}^∞ (−f′) dε ]  =  2 f(Δ(T)) ,

verified numerically to 3·10⁻⁶ over T/Tc = 0.3–0.9 with convergence ∝ ω (test T5).

## 4. Correspondence table

| Levy (ql<1, ωτ≪1, Boltzmann + drag) | this work (clean, q→0, golden rule) |
|---|---|
| "integral even in ε" → factor 2 | branch sum Σ_σ over ε_σ(E) = ε* + σℓ/κ (s-wave: ε* = μ, ℓ = (E²−Δ²)) |
| angular kernel 𝒦(b, ql) constant in E | exact identity C_I(E,E)·(E/ℓ)·(ℓ/E) = 1 per branch |
| normal denominator ∫(−f′)dε = 1 | S_N/(N²ω) = 1 (flat DOS) |
| ħω neglected against Δ, T (never shifts E by ω) | ω→0 taken explicitly; E′ = E+ω retained at finite ω |
| prefactors (ql)²/τ, m, ρ cancel in ratio | prefactor 2πB²ω cancels in ratio |
| s-wave ε = ±√(E²−ε₀²) | HHM ε_σ = ε* + σℓ/κ, reduced to s-wave in T5 |

## 5. Mathematical demonstration that the same physical regime is reproduced (Q4)

Both computations reduce the ratio to the same functional

    L[K]  =  2 INT_Δ^∞ K(E)(−f′) dE  /  INT_{−∞}^∞ K(|ε|)(−f′) dε

evaluated at K ≡ 1: Levy because his angular kernel is E-independent at ql<1, ωτ≪1 and
cancels; we because the coherence×Jacobian identity makes the integrand E-weight-free at
q→0. Test T5b probes L[K] numerically:

    K = 1                 :  L[K] = 0.066128 = 2f(Δ)      dev 4.6·10⁻¹⁴
    K = E/Δ               :  L[K] = 0.209532              dev 2.17
    K = 1 + 0.1(E/Δ)²     :  L[K] = 0.075729              dev 0.145

so the *operative* assumption is kernel energy-independence on the thermal scale — and each
regime supplies it by its own mechanism (collisional angular averaging vs algebraic
identity). The agreement is therefore structural, not accidental. The bridge across the
ql axis is in Levy's own text: *"A similar result was obtained by BCS for longitudinal waves
for ql > 1 by computing the net rate of absorption of energy in the superconducting state
produced by direct absorption and induced emission of the imposed acoustic phonons"* —
i.e. the clean-limit golden rule (exactly our framework, with direct absorption minus
induced emission = [f(E)−f(E′)]) yields the same ratio from the ql>1 side. Our ω→0, q→0
limit sits in that universality class; Levy's ql<1 solution reaches it from the other side.

A third, modern bridge. Sardwal, Kreisel and Mukherjee (arXiv:2608.19079, 2026) compute the
attenuation from the phonon-number equation of motion in *real-space* BdG, in the clean
limit (mean free path larger than the phonon wavelength — the ql>1 side), with
chi_mn = (f(E_m)−f(E_n))/(hbar*omega + E_m − E_n + i*eta) and coherence-weighted matrix
element F_mn ∝ (u^m* u^n − v^n v^m*). Homogeneously this reduces to
(uu'−vv')^2 = ½[1+(xi xi' − Delta Delta')/(EE')] — our signed C_I — and their s-wave
benchmark reproduces the standard exponentially suppressed ratio, cross-checked against a
fully k-space formulation and the analytic Walker–Smith–Samokhin treatment (PRB 65, 014517
(2001)); their PP channel (pair creation/annihilation, weight 1−f_m−f_n) is "negligible for
homogeneous superconductors", the published counterpart of our S_pair ≡ 0 inside
2*Delta > hbar*omega. Collisional transport (Levy), clean golden rule (this work) and
real-space BdG (theirs) are therefore three formulations of one ratio, each validated
against the others; our T7 (`dualspace_check.py`) repeats their three-way protocol for this
model, closing the k-space/energy-space/Levy triangle at <=1.4 %.

## 6. Where the regimes differ — scope of the benchmark (T5b, second table)

At **finite ω** the prescriptions differ, because we resolve E′ = E + ħω in the occupations
and in the DOS/coherence structure while Levy's transport treatment keeps occupations at E:

    T/Tc = 0.5:   ħω = 1e-4 t : this work 0.06613 = Levy 0.06613   (dev 4·10⁻⁶)
                  ħω = 0.05 t : this work 0.05544 vs Levy 0.06613  (dev 0.16)
                  ħω = 0.70 t : this work 0.01030 vs Levy 0.06613  (dev 0.84)

Physical ultrasound has ħω/Δ ~ 10⁻⁴–10⁻³ (MHz–GHz vs a gap of order K), so for any real
experiment the ω→0 identity of the two regimes holds and Levy is the applicable benchmark.
The paper's ħω = 0.7t ~ Δ is a *theoretical* parameter outside the experimental ultrasound
window; there our curves are the clean-limit golden-rule (BCS-type) answer, and Levy's 2f
is not the correct comparison — a fact we state rather than hide. Finite q beyond the
vertex (Tsuneto's ql≫1 kernel, Pippard's transverse treatment) and τ-dependent crossovers
are likewise outside both this benchmark and the model's claims.

## 7. Direct answers

1. **Implementation of ql<1:** none; no q, ℓ in the code. The vertex is q→0, the shared
   leading order; Levy's ql<1 enters only through the E-independence (and eventual
   cancellation) of his angular kernel.
2. **τ / ℓ = v_F τ:** absent; clean golden rule. Justified for the *ratio* because in Levy's
   regime α_s and α_n carry identical (ql)²/τ·𝒦 prefactors, which cancel (his Eqs. 6–7).
3. **ω = v_s q:** not imposed; ω free, q→0 direct. Levy needs q=ω/v_s only to form ql, i.e.
   to bound τ — irrelevant to a τ-free ratio.
4. **Regime vs formula:** Sections 2–5 show both calculations reduce to L[1]; T5b shows
   L[K] ≠ 2f for any energy-dependent K, isolating the shared assumption; the BCS-ql>1
   citation in Levy's own paper places our clean golden rule in the same universality
   class. Section 6 delimits the difference at finite ω.

## 8. Validity statement (for referees)

The Levy comparison is a **valid and sharp benchmark of the universal core** of the
computation — the sign structure of the coherence factor, the DOS measure and its
Jacobian, the branch (σ) bookkeeping, the Fermi factors, and the q→0 vertex — because any
error in those ingredients breaks L[1] ≠ 2f(Δ). It is **not** a validation of finite-q,
finite-τ transport physics, which the model does not contain and does not claim
(DEFENSE_referee_QA.md, Part IV). Extending the pipeline with a Tsuneto-type finite-ql
kernel or a relaxation-time self-energy is well-defined future work, not a missing
ingredient of the present claims.

## 9. Reproduction

    python levy1963_check.py     # T5 (ratio vs Levy, omega-scan) and T5b (regime analysis)
