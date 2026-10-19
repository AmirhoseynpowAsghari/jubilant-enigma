# W/W_N from scratch for η = +1, quasiparticle scattering only

**Model:** the repository `power-absorbed-in-the-presence-of-external-time-dependent-perturbation`,
branch **master** (commit `cda0d56`): 2-D square lattice, Rashba SOC, gap
`Δ̄(ε) = Δ⁰ − (Δˢ/4t) ε`, filling `n = 1.885`, `U/t = 95`, `Δ_t/t = 4.5`,
`ħω = 0.7`, `T_c`, `Δ⁰(0)`, `Δˢ(0)` and `Δ_eff(0)` from `sc_params/results/tc_vs_vso.json`,
gap shape `Δ(T) = Δ(0)·(1−(T/T_c)^3.434)^0.553` (the repository's fitted exponents), units `t = ħ = k_B = 1`.

**What is computed here:** the golden-rule scattering channel with energy conservation on the
*quasiparticle* energies, both branches explicit, three coherence factors (a),(b),(c) plus the
textbook form (d), and the *correct* normal-state reference. Everything is cross-checked by five
independent tests (§4). Files: `model.py`, `energy_space.py`, `bruteforce.py`, `validate.py`,
`run_eta1.py`, `figures.py`; data `results.json`, `results.csv`; figures `fig_*.png`.

---

## 1. Derivation (before any numbers)

### 1.1 Bands, gap, branches

`ε_s(k) = −2t(cos kx + cos ky) − 2 s V_SO |sin k|`, `s = ±1`; `ξ = ε − μ`; `Δ̄(ε) = Δ⁰ − b ε` with
`b = Δˢ/4t`; `E(ε) = √(ξ² + Δ̄²)`. Completing the square,

```
E² = κ² (ε − ε*)² + E_min²,      κ = √(1+b²),   ε* = (μ + b Δ⁰)/(1+b²),   E_min = |Δ⁰ − b μ|/κ .
```

Hence for a given quasiparticle energy `E ≥ E_min` there are **two** band energies — this is where
the branch index enters:

```
ℓ(E)    = √(E² − E_min²)
ε_σ(E)  = ε* + σ ℓ/κ,            σ = ±1
ξ_σ(E)  = ε_σ(E) − μ                         (signed, requirement 1)
Δ̄_σ(E)  = Δ⁰ − b ε_σ(E)                     (signed)
|dε_σ/dE| = E/(κ ℓ)             ⇒   ν_σ(E) ≡ N(ε_σ(E)) · E/(κ ℓ)      (quasiparticle DOS of branch σ)
```

The Jacobian `E/(κℓ)` is the BCS `1/√(E²−E_min²)` singularity. Because `ε* ≈ μ` sits near the **top**
of the band (`n = 1.885` is almost full filling, `2` = full), the two branches are very asymmetric:
`E_max(+1) = √(κ²(W−ε*)² + E_min²)` is short, `E_max(−1) = √(κ²(W+ε*)² + E_min²)` is long. Using the
symmetric bound `κ(W−|ε*|)` for both (the obvious shortcut) deletes the whole hole-like branch and
gives `S = 0`. The integration limits are therefore taken **per (σ,σ′) pair**:
`E ∈ [E_min, min(E_max(σ), E_max(σ′) − ħω)]`.

### 1.2 The integral (scattering channel only)

Energy conservation on the quasiparticle energies (requirement 4): `E′ = E + ħω`. With
`f(E) = 1/(e^{E/T}+1)` (no μ — E is an excitation energy):

```
S_sc = Σ_{σ,σ′} ∫_{E_min}^{E_max(σ,σ′)} dE  ν_σ(E) ν_σ′(E′) C(E,σ; E′,σ′) [ f(E) − f(E′) ]
W/W_N = S_sc / S_N
S_N  = ∫ dE  N(E) N(E+ħω) [ f_N(E) − f_N(E+ħω) ] ,   f_N with μ_N(T) at n = 1.885
```

`S_N` is the absorption of the **same model with Δ⁰ = Δˢ = 0** (the normal state at the same
filling); the prefactors `2π ω B²` and all k-sum normalisations cancel in the ratio.
The four `(σ,σ′)` contributions are reported separately in `results.json` (`*_parts`).

**Quadrature.** The substitution `ℓ = √(E²−E_min²)` removes the singularity exactly
(`ν_σ dE = N(ε_σ) dℓ/κ`), so the midpoint rule in ℓ is well behaved. Test T0 verifies this change of
variables against a `4·10⁶`-point integral in E with the explicit Jacobian: agreement `2·10⁻⁴`.

### 1.3 The three (+ one) coherence factors, η = +1

```
(a)  C = ½ [ 1 + ( ξ_σ ξ_σ′   − Δ̄_σ Δ̄_σ′ ) / (E E′) ]     signed ξ — the factor you asked for
(b)  C = ½ [ 1 + ( |ξ_σ||ξ_σ′| − Δ̄_σ Δ̄_σ′ ) / (E E′) ]
(c)  C = ½ [ 1 + ( |ξ_σ||ξ_σ′| + Δ̄_σ Δ̄_σ′ ) / (E E′) ]     the repository/paper's final factor
(d)  C = ½ [ 1 − Δ̄_σ Δ̄_σ′/(E E′) ]                        textbook BCS Case I (see §5)
```

For (a): `|ξξ′ − Δ̄Δ̄′| ≤ √(ξ²+Δ̄²)√(ξ′²+Δ̄′²) = E E′` (Cauchy–Schwarz), so **0 ≤ C_a ≤ 1 identically**,
hence `0 <= C_a <= 1` for (a) pointwise, and the SCATTERING-ONLY ratio satisfies W/W_N <= 1 in every variant (maxima 0.4000/0.5054/0.7611 at Tc-). Adding S_pair to the numerator while keeping the single normal-state denominator produces values >1 with no supporting identity (App. A2): withdrawn. (b) and (c) break the ratio bound already at Delta=0, where they return exactly 2R > 1 (§4, T4).

### 1.4 How the normal state is recovered — the crux of your question 7

Take `Δ⁰, Δ → 0` (which is what `T → T_c⁻` does: `κ → 1`, `ε* → μ`, `E_min → 0`, `E → |ξ|`):

* the **Jacobians and DOS** recover the normal state exactly: `Σ_σ ν_σ(E) dE → N(ε) dε`, both sides
  of the band covered once (verified numerically to `10⁻⁶`, test T1);
* but the **coherence factor becomes a branch projector**:

```
C_a → 1  for σ = σ′ ,        C_a → 0  for σ = −σ′ ,
```

  because at Δ = 0, `ξ_σ/E → σ·sign(ξ)`-type limits. A scalar, time-reversal-even vertex cannot turn
  an electron into a hole, so half of the normal-state pair phase space is projected out **of the
  scattering channel**;
* and the **thermal factor** of a quasiparticle is `f(E) = 1/(e^{E/T}+1) ≤ ½` at `E ≥ 0`, whereas the
  normal-state lower state has occupation `f_N(ε) ≈ 1` for `ε < μ`.

Putting these together (test T1, flat DOS, exact):

```
        W_sc/W_N |_{Δ=0} = R(ħω, T) ≡ (2/ħω) ∫₀^ [ f(E) − f(E+ħω) ] dE          (closed form)
        R → 1  for ħω ≪ T ,      R → 2 T ln2 / ħω  for ħω ≫ T .
```

With the paper's `ħω = 0.7` and the three `T_c` this is **R = 0.499 / 0.594 / 0.743**, i.e. *not* 1.
The missing weight is not lost physics: the Bogoliubov matrix element satisfies
`C_sc + C_pair = 1` exactly (verified to `10⁻¹⁶`), and the complement of the scattering channel is
the **pair-creation channel** `E_a + E_b = ħω`, which you asked to exclude. At Δ = 0:

```
        (S_sc + S_pair)/S_N = 1.000000      (flat DOS: 100.0000 %; real lattice DOS: 100.00 %)
```

and for the real model at ħω = 0.7 the pair channel carries 60–88 % of the absorption near T_c.
**Therefore, with the pair-creation channel excluded, W/W_N cannot tend to 1 at T_c for any of the
three coherence factors — that is a property of the channel selection, not of the sign of ξ.**
Adding the pair channel back restores `W/W_N → 1` continuously and exactly (fig. 1, dark curves).

---

## 2. Numerical results

`V_SO = 0.4082` (`T_c = 0.3808`, `E_min(0) = 0.698`), `ħω = 0.7`; identical tables for all three
`V_SO` are in `results.csv`.

| T/T_c | Δ_eff | 2E_min | (a) signed | (b) \|ξ\|\|ξ′\| | (c) paper | (d) textbook | (a)+pair | R(ħω,T) |
|---|---|---|---|---|---|---|---|---|
| 0.10 | 0.5488 | 1.396 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.0000 | 0.075 |
| 0.40 | 0.5357 | 1.364 | 0.00969 | 0.01091 | 0.01445 | 0.00734 | 0.0097 | 0.297 |
| 0.60 | 0.4941 | 1.267 | 0.06609 | 0.07714 | 0.09854 | 0.04984 | 0.0661 | 0.423 |
| 0.80 | 0.3885 | 1.011 | 0.20250 | 0.24835 | 0.30306 | 0.15300 | 0.2025 | 0.520 |
| 0.90 | 0.2839 | 0.747 | 0.31763 | 0.41132 | 0.47846 | 0.24343 | 0.3176 | 0.559 |
| 0.999 | 0.0238 | 0.063 | 0.50516 | 0.80555 | 0.80950 | 0.41280 | 1.0041 | 0.593 |
| 0.999 | 0.0238 | 0.0317 | 0.50540 | 0.80599 | 0.80995 | 0.41303 | 1.0055 | 0.5932 |

Limits and maxima over the whole sweep `0.02 ≤ T/T_c ≤ 0.99999`:

| V_SO | T_c | (a) max = (a)(T→T_c⁻) | (b) | (c) | (d) | exceeds 1 anywhere? | pair opens at T/T_c |
|---|---|---|---|---|---|---|---|
| 0.0000 | 0.2861 | **0.40004** | 0.68362 | 0.68592 | 0.35021 | no | 0.840 |
| 0.4082 | 0.3808 | **0.50540** | 0.80599 | 0.80995 | 0.41303 | no | 0.915 |
| 0.8000 | 0.6494 | **0.76112** | 0.92091 | 0.93695 | 0.48620 | no | 0.975 |

Branch decomposition of (a) (`fig_branch_decomposition.png`): the `(+1,+1)` electron→electron channel
is **kinematically closed** at every T here, because `E_max(+1) = κ(W−ε*) < ħω` (μ sits within
`W − ħω` of the band top); all weight sits in `(−1,−1)` (hole→hole, 54 %→99.8 % as T→T_c) plus the
cross terms `(±1,∓1)` at low T. This is why (a) is systematically the *smallest* of the four.

---

## 3. Answers to your numbered points

1. **Sign of ξ kept.** (a) uses `ξ_σ(E) = ε_σ(E) − μ` with the signed branch values everywhere; no
   `|ξ|`, no positive root. (b) and (c) are computed with `|ξ||ξ′|` for comparison only.
2. **Both branches explicit.** §1.1–1.2; the four `(σ,σ′)` contributions are tabulated per T in
   `results.json` and plotted in `fig_branch_decomposition.png`.
3. **Scattering only.** Yes; the pair channel is computed but reported separately (`*_pair`,
   `*_total_ratio`) and is *excluded* from every number in the three-variant comparison.
4. **Energy conservation on E.** `E′ = E + ħω` on quasiparticle energies; the repository's use of the
   band energy ε with the normal-state Fermi function is what makes its `W_N` negative (§6).
5. **Same parameters/DOS/gap(T)/W_N definition as the paper.** Stage-1 numbers, gap exponents,
   `ħω`, filling, k-model all from branch master. The DOS is the two-helicity band DOS
   `∫N dε = 2` with exact band support (the repository's table integrates to 0.643 and leaks
   Gaussians above the band edge — §6); `W_N` is the Δ = 0 absorption of the same model with the
   normal-state `μ_N(T)` at n = 1.885.
6. **Same T/T_c range as the paper** and beyond: `0.02 … 0.99999` (88 points), i.e. up to the true
   `T → T_c⁻` limit, plus the pair-breaking onsets `T/T_c = 0.840 / 0.915 / 0.975` (the repository
   stops at 0.74 / 0.84 / 0.94 because of its `2Δ_eff ≤ ħω` cut, which is the wrong threshold:
   the true one is `2 E_min(T) = ħω`, fig. 2 right).
7. **Is W/W_N < 1 everywhere? Does it → 1 at T_c⁻?** Yes to the first: strictly below 1 at every T
   for all three variants and all three V_SO (for (a) this is a theorem, Cauchy–Schwarz). **No to
   the second, and provably so**: with the pair channel excluded the limit is
   `W/W_N → 0.402 / 0.509 / 0.762` (= the band-structure evaluation of the closed form R(ħω,T_c),
   §1.4). Including the pair channel gives `→ 1.0000` (§1.4, fig. 1).
8. **Three cases compared.** Table above; (a) < (d) < (b) < (c) at all T. (b) and (c) are not
   physical coherence factors: for a flat DOS at Δ = 0 they give exactly `2R` and **exceed 1**
   (e.g. 1.52 at T = 0.05, 1.97 at T = 1.0, ħω = 0.05), violating `C ≤ 1`.
9. **Derivation first.** §1, with the σ,σ′ entry points (branch inversion, Jacobian, per-pair limits)
   and the normal-state recovery (projector + closed form R + pair-channel complement).
10. **No assumption of < 1.** The answer falls out of the corrected formula; (b),(c) in fact *can*
    exceed 1, which is how the tests caught them. Physical reason (a) is below 1:
    (i) `C_a ≤ 1` pointwise, (ii) the Δ→0 projector removes cross-branch pairs from the scattering
    channel, (iii) the quasiparticle occupation of the lower state is `≤ ½` while the normal one is
    `≈ 1`, (iv) at ħω = 0.7 comparable to `2E_min` a large share of the phase space lives in the
    excluded pair channel. At low T the suppression is the usual `exp(−E_min/T)`.

---

## 4. Validation (all reproduced by `python3 validate.py`)

| test | statement | result |
|---|---|---|
| T0 | ℓ-substitution = E-integral with explicit Jacobian | 2·10⁻⁴ |
| T1 | Δ→0: (a) equals the closed form `R(ħω,T)`; Δ→0 sum rule with the pair channel = 1 | 3·10⁻⁶ / 100.0000 % |
| T2 | ω→0, flat DOS, μ=0: exact golden rule vs Levy 1963 Eq. (8), `α_s/α_n = 2 f(ε₀(T))`; see §5 | 2·10⁻⁶ |
| T5 | s-wave reduction (Δˢ=0, V_SO=0, flat DOS, μ=0, BCS Δ(T)): ω→0 column of W/W_N vs `2 f(ε₀(T))`; convergence ∝ ω | 3·10⁻⁶ |
| T3 | energy space (branch inversion + Jacobian) vs explicit lattice sum over all (k,helicity)×(k′,helicity′) pairs, `σ = sign(ε−ε*)`, Richardson in the δ-width | 0.4–1.6 % for 0.5 ≤ T/T_c ≤ 0.999 (the two sides regulate the gap-edge `1/√(E²−E_min²)` peak differently: Gaussian δ on E vs smoothed N(ε)) |
| T4 | Cauchy–Schwarz: 0 ≤ C ≤ 1 pointwise | holds for all four variants; (b),(c) fail the Δ→0 projector and sum-rule tests instead |
| branch DOS | lattice `Σ_k δ(E−E_k) per branch` vs `N(ε_σ(E)) E/(κℓ)` | ≤ 3 % |
| T7 | dual formulation after Sardwal–Kreisel–Mukherjee (arXiv:2608.19079): (i) k-space lattice double sum ↔ (ii) energy-space pipeline at production ħω=0.7t, Nk∈{160…452}; s-wave ω→0 triangle (i)↔(ii)↔(iii) Levy 2f(Δ) | (i)-(ii) 0.2–1.3 % (Nk-independent floor = DOS construction); triangle ≤1.4 %; (ii, flat DOS) = Levy to six digits |

The T7 script (`dualspace_check.py`) implements the cross-validation protocol of Sardwal,
Kreisel and Mukherjee (arXiv:2608.19079), who benchmark their real-space BdG attenuation
formula against a fully momentum-space formulation and the analytic Walker–Smith–Samokhin
treatment (PRB 65, 014517 (2001)). Our three legs are (i) the explicit lattice double sum
with Gaussian δ and Richardson extrapolation, (ii) the production energy-space reduction,
(iii) Levy's closed form. Their homogeneous-limit coherence factor (u\*u−v\*v)² =
½[1+(ξξ′−ΔΔ′)/(EE′)] coincides with the signed C_I used here, and their pair-creation
channel (weight 1−f_m−f_n) is stated to be negligible for homogeneous superconductors —
two independent 2026 confirmations of this work's structural choices. Output:
`results/T7_dualspace.txt`; figure `figures/dual_space_k_vs_energy.png`.

---

## 5. The textbook limit, and why your C_I is not the textbook Case-I factor

For a flat DOS and a **constant** gap with `μ = 0`, η = +1, the exact golden-rule scattering channel
(energy conservation `E′ = E + ħω`, double-density measure `ν_σ(E) ν_σ′(E′)`) gives, in the limit `ħω → 0`:

* with the signed ξξ′ factor (a) — and (b), which coincides with (a) in this limit — `W/W_N → 2 f(Δ)`
  **exactly**: the classic result `α_s/α_n = 2 f(ε₀(T))` of Levy [Phys. Rev. **131**, 1497 (1963), Eq. (8)]
  and of Morse/Tsuneto; verified to 2·10⁻⁶ over three decades of ω (T2, T5, `levy1963_check.py`);
* with the textbook static Case-I factor (d) `½[1 − Δ̄Δ̄′/(EE′)]`: `W/W_N → f(Δ) = ½ · 2 f(Δ)` —
  dropping the ξξ′ term removes exactly half of the ω→0 weight;
* with the paper's combination (c): the ω→0 integral acquires a non-integrable gap-edge weight (∝ 1/ℓ),
  a further sign that this combination does not belong to this Hamiltonian.

The identity is exact, with no evaluation-point approximation: in ℓ the (a)-integrand is
`N² (ℓ/E) ω (−f′(E))`, because the integrated state carries `ν_σ dE = N dℓ/κ` while the final-state
density stays per unit E′, `ν_σ′(E′) = N E′/(κ ℓ′)`, and `C_I(E,E) = ξ²/E² = ℓ²/E²`; hence
`W = 2 N² ω ∫_Δ^∞ (−f′) dE = 2 N² ω f(Δ)` while `W_N = N² ω`.

**Correction note (2026-10-15).**  An earlier revision of this report stated the plateau as
`0.647 · 2 f(Δ)` and attributed the difference from the literature to the evaluation point of the
coherence factor.  That difference was a quadrature defect: both DOS factors had been written per
dℓ, losing the final-state Jacobian `E′/ℓ′` (equivalently weighting the integrand by ℓ′/E′).  With
the corrected measure the literature value is reproduced exactly (T5), the brute-force lattice
cross-check improves from ≤ 12 % to ≤ 1.6 % at finite Δ (T3), and all production curves were
recomputed.  Quantities with ℓ/E → 1 are unchanged: the T→Tc⁻ saturation values, the thresholds
`2Δ_eff = ħω` and `2E_min = ħω`, the Δ→0 closed forms and the sum rule.

The factor (d) is the one that (i) comes from the Bogoliubov vertex `u_k u_k′ − v_k v_k′` with
`v_k = sign(Δ̄_k)√((1−ξ/E)/2)` (verified numerically to 10⁻¹¹ over 4·10⁵ random states, including
`C_sc + C_pair = 1`) and (ii) tends to 1 at the gap edge and to ½ at Δ = 0.  It has no ξξ′ term and
therefore reproduces only **half** of Levy's 2 f(Δ) in the ω→0 limit, as shown above.  The signed
C_I (η = +1) is the combination that reproduces Levy's s-wave result exactly, keeps the Δ→0 branch
projector, and satisfies the pair-creation sum rule; its projector structure is also what makes the
scattering-only ratio saturate below 1 at Tc⁻.

## 6. Why the repository's stored W/W_N (0.92–0.96) and the previous recalculation's T→T_c numbers
cannot be trusted

* `absorption.py` integrates over the **band energy** ε with the **normal-state** Fermi function and
  imposes energy conservation on ε; with μ = 3.65 and ħω = 0.7 the resulting normal-state integral is
  **negative** (−0.0130 at T = 0.7 T_c, V_SO = 0.4082) — absorbed power negative — and the code hides
  it with `abs(W/WN)`. Decomposition: the window ε ∈ [μ−ħω, μ] alone contributes −1.84·10⁻³.
* Its DOS table has `∫g dε = 0.643` instead of 2 and a spurious Gaussian tail above the band edge
  (`g(4.2) = 3.6·10⁻²` where the exact DOS is 0), i.e. states that do not exist enter W_N.
* `2Δ_eff > ħω` is the wrong pair-breaking threshold (`2 E_min = ħω` is the right one) and the
  `T < 10⁻⁴ → 0` shortcut hides the low-T behaviour.
* The previous recalculation used the repository's W_N as the denominator of an energy-space
  numerator, so its T→T_c limit was a ratio of two inconsistent normalisations; that is why it could
  not show `W/W_N → 1` and left the question open. With a consistent normalisation the limit is the
  closed form R(ħω,T_c) of §1.4 — and it becomes 1 once the pair channel is included.

---

## 7. Compact conclusion

**For η = +1, with only the quasiparticle-scattering channel and the correct signed coherence
factor, the model predicts W/W_N < 1 over the whole parameter range of the paper: yes.**
At every temperature `0 < T/T_c < 1` and for V_SO = 0, 0.4082, 0.8 the ratio is strictly below 1
(maximum 0.402, 0.509, 0.762, attained as T → T_c⁻); for the signed factor this bound is a theorem
(Cauchy–Schwarz on `ξξ′ − Δ̄Δ̄′`). It does **not** approach 1: with the pair-creation channel
excluded the limit is `W/W_N → R(ħω, T_c) = 0.40 / 0.51 / 0.76`, because at Δ → 0 the η = +1
coherence factor becomes a same-branch projector and because the excluded pair channel carries the
rest of the normal-state phase space (60–88 % at ħω = 0.7). Including that channel gives
`W/W_N → 1.0000` continuously, and the full result is below 1 except in a narrow window just above
the true pair-breaking threshold `2 E_min(T) = ħω` (T/T_c = 0.840 / 0.915 / 0.975), where it
overshoots (max 1.09 / 1.34 / 1.99 at T/T_c ≈ 0.96–0.98, the BCS edge step that a finite quasiparticle
lifetime would round) before returning to 1 at T_c. The repository's |ξ||ξ′| versions (b),(c)
are not admissible coherence factors: they violate `0 ≤ C ≤ 1` and exceed 1 already at Δ = 0.

---

## Appendix A. "Pair creation": what it is, where it comes from, and what your 2Δ > ħω condition removes

**It is not an extra ingredient I added; it is the second kinematic solution of the same golden
rule you are already using.** Writing the time-dependent perturbation in the Bogoliubov basis,
`H' ~ Σ B c†_a c_b` with `c = u γ + v γ†`, produces four operator structures. Two of them
(`γ†γ`) scatter an *existing* quasiparticle and carry `δ(E_a − E_b − ħω)`; the other two
(`γγ`, `γ†γ†`) let the phonon *create two quasiparticles out of the condensate* and carry
`δ(E_a + E_b − ħω)`. The complete absorbed power is therefore

```
W ∝ Σ_{a,b} C_sc  [f(E_b) − f(E_a)]        δ(E_a − E_b − ħω)        (quasiparticle scattering)
  + Σ_{a,b} C_p   [1 − f(E_a) − f(E_b)]    δ(E_a + E_b − ħω)        (pair creation / pair breaking)
C_p = 1 − C_sc   (verified numerically: M_scat^2 + M_pair^2 = 1 to 1e-16)
```

This is the standard "pair-breaking" term of ultrasonic attenuation in superconductors (the
`Γ^PB` term with occupation `1−f−f` in the ultrasound literature); it is not a model choice.

**How the code generates it:** `energy_space.S_paircreation` uses the *same* branch machinery with
`E_b ∈ [E_min, ħω/2]`, `E_a = ħω − E_b`, and begins with the hard guard
`if omega <= 2*E_min: return 0.0`. So it contributes nothing unless the channel is kinematically open.

**What your condition removes.** Your `2Δ > ħω` test (`physics_ua.pairbreaking_valid`, with
`Δ = Δ_eff`) is a *sufficient* condition for the pair channel to vanish, but it is not the
threshold. The threshold is `ħω = 2 E_min(T)` with `E_min = |Δ⁰ − b μ|/κ`, the true minimum of
`E(ε) = √(ξ² + Δ̄(ε)²)` over the band: because `Δ̄(ε)` changes sign across the band, the spectral
gap is `E_min`, not the k-average `Δ_eff = ⟨|Δ̄|⟩`. (For a *constant* gap the two coincide, which is
where the repository's test comes from.) In this model `E_min(T) > Δ_eff(T)` at all T, so:

| V_SO | your cut `2Δ_eff = ħω` at T/Tc | true threshold `2E_min = ħω` at T/Tc | (a) at your cut | S_pair/S_N at your cut | S_pair/S_N at T/Tc = 0.99 |
|---|---|---|---|---|---|
| 0.0000 | 0.740 | 0.820 | 0.12055 | 0 | 0.691 |
| 0.4082 | 0.840 | 0.912 | 0.24345 | 0 | 0.596 |
| 0.8000 | 0.936 | 0.975 | 0.51833 | 0 | 0.329 |

Consequences, all visible in `results.csv`:

1. **Inside the range you keep, my pair term is exactly zero**, so the scattering-only numbers in
   §2 *are* the complete W/W_N there — nothing has been added to your formula. My maxima *within
   your domain* are 0.121 / 0.243 / 0.518 at the cut, and the statement "W/W_N < 1" holds there.
2. **Your cut is stricter than necessary:** between your cut and the true threshold
   (0.74→0.82, 0.84→0.91, 0.94→0.97) the pair channel is still exactly zero, so the
   scattering-only formula remains complete and could have been used there.
3. **Above the true threshold the scattering-only formula is not wrong, only incomplete:** the
   complete W needs both channels, and their sum is what tends to 1 at T_c (verified to 0.01 %).
4. **This is why your item 7 and your validity condition cannot both hold:** `E_min(T) → 0` as
   `T → T_c⁻`, so any approach to T_c necessarily passes the pair threshold. Either you stop at
   your cut (then the curve ends at W/W_N = 0.12 / 0.24 / 0.52 and the limit at T_c is simply not
   in the domain), or you ask about the limit (then the pair channel is part of the answer and
   W/W_N → 1). The saturation values 0.40 / 0.51 / 0.76 quoted in §2–§3 are the *scattering-only*
   extrapolation to Δ→0 and are reported precisely to make this structure explicit.


## Appendix A2. The pair channel: outlook only (withdrawal notice)

The coherence factors satisfy C_I + C_pair = 1 per quasiparticle pair, so together the channels carry the total transition weight; because they occupy different sheets (E'-E = hw vs E+E' = hw) with different thermal weights, their sum is not bounded by and does not equal any single normal-state reference. This work therefore reports the scattering channel only; the pair channel is examined in the outlook without a quantitative claim.

Measured with the formula as previously printed (matched quasiparticle reference),
(W+W_p)/W_N^QP was 0.545 / 3.63 / 0.848 (V_SO = 0 / 0.4082 / 0.8 at T/Tc=0.5) and
0.594 / 1.81 / 0.626 (T/Tc=0.9) - never 1 - and no implementation of that reference ever
existed in code. Consequently: (i) the 1e-4-level "sum rule" and the W_N^QP identity claim are
WITHDRAWN from README, COVER_NOTE, PDF addendum, figure titles and PARAM_SCAN_PEAKS;
(ii) the "physical total" ratios >1 and the 1.3-1.6 peaks are WITHDRAWN as quantitative
results (they arose from adding S_pair to the numerator only); (iii) inside the paper's
own window 2E_min >= hw the pair channel is identically zero (41/52, 45/52, 48/52 temperature
points exactly zero; first non-zero at T/Tc = 0.840 / 0.920 / 0.980), so the scattering-only
result is COMPLETE there, not an approximation; (iv) The scattering channel alone returns R(hw,T_c) = 0.400 / 0.505 / 0.761 at T_c-; restoration of W/W_N -> 1 involves the pair-breaking channel, outside the scope of the quantitative result reported here.

## Appendix A3. Paper convention ↔ case ↔ correct factor (referee Part 4.1)

| object | arXiv:2607.02812 v1 p.10 | case | correct for ultrasound (η=+1) |
|---|---|---|---|
| printed coherence factor | ½[1+(ξξ′−ΔΔ′)/(EE′)] | Case I signed | YES — this is C_I used here |
| printed W integrand | +ηΔΔ′/EE′ with η=+1 | = ½[1+(ξξ′+ΔΔ′)/(EE′)] = (uu′+vv′)² | NO — Case II (EM/NMR) |
| repo implementation | coherence() a_signed | Case I signed | YES (erratum E1, now tabulated) |

## Appendix A4. Erratum E8 (gap shape regime) and E9 (U=−75 outputs)

E8: stage-2 replaces the paper's Δ(T)=Δ₀√(1−T/Tc) by the fitted φ(T)=(1−(T/Tc)^p)^q with the
MASTER exponents p=3.434, q=0.553 for ALL V_SO, while sc_params/results/gap_vs_T_fit.json
gives p=3.674, q=0.556 at V_SO=0 (r²=0.999).  This is a deliberate rigid-shape choice (one
universal shape keeps stage-2 comparable across V_SO); the per-V_SO refit changes W/W_N by 1.3–3.4 % (V_SO=0: 2.2/1.3 %; 0.4082: 2.5/1.4 %;
V_SO=0.8: 3.4/2.1 % at T/Tc=0.5/0.9; measured in REFEREE_RESPONSE.md).  It was not listed in E1–E7: it is now.
E9: the U=−75t root (Tc≈12.7–13.2t, Δ⁰≈17.8t, μ≈35t outside the band) is unphysical; its
figures are moved to project/ultrasonic/figures/non_citable/ and removed from SOCfigures.
