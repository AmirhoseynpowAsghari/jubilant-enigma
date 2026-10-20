# Referee & defense package — how to prove the computation is correct

Companion of `REPORT_WWN_corrected.md`, `DERIVATION_HHM_scattering.md`,
`levy1963_check.py`, `case2_checks.py`, `validate_ua.py`.
Everything below is reproducible by running the named scripts; no number in this document
is asserted without a command that regenerates it.

---

## Part I — The proof arsenal: benchmarks against accepted, famous results

The strongest possible answer to "why should we trust your numerics?" is: *the same code,
with no adjustable parameters, reproduces results that the community accepted decades
ago.* Five such benchmarks are implemented as permanent tests.

| # | benchmark (accepted literature) | what our code reproduces | precision | command |
|---|---|---|---|---|
| B1 | **Levy 1963**, Phys. Rev. 131, 1497, Eq. (8): ultrasound, ql<1, s-wave: α_s/α_n = 2f(ε₀(T)) | ω→0, q→0 limit of our golden-rule ratio, s-wave reduction (Δˢ=0, V_SO=0, flat DOS, μ=0, BCS Δ(T)); regime analysis in `LEVY_REGIME_ANALYSIS.md` | 3·10⁻⁶ over T/Tc=0.3–0.9, convergence ∝ ω | `python levy1963_check.py` |
| B2 | **Mattis & Bardeen 1958**, Phys. Rev. 109, 1553: the σ₁ kernels (EE′±Δ²)/(ℓℓ′) | branch-summed Case-II integrand of our pipeline equals the MB kernels (factor 2 = their particle-hole symmetrisation, cancels in ratios); MB threshold: absorption = 0 for ħω<2Δ, turns on above | 8.7·10⁻¹⁶ (machine precision), 2000 random states | `python case2_checks.py` |
| B3 | **Hebel & Slichter 1959**, Phys. Rev. 113, 1504: coherence peak of Case-II rates below Tc | same code, η=−1 vertex: ω→0 rate ratio peaks at ≈2.7 just below Tc and returns to 1 at Tc; η=+1 (our ultrasound Case I) stays monotone < 1 (Levy) | qualitative signature, exact shape | `python case2_checks.py` (figure `caseII_mattis_bardeen_hebel_slichter.png`) |
| B4 | **Hutchinson–Hirsch–Marsiglio 2018**, arXiv:1801.08988 (the model itself): kernel Eq. (20), gap linearity Eqs. (36–37), band minimum −4t√(1+V_SO²/2t²) | our stage-1 solver (same kernel) reproduces every printed Tc, Δ⁰, Δˢ of the paper's parameter sets to 4 decimals; band edge of `kgrid_ua` equals their exact minimum | 10⁻⁴ | `python main_tc.py groups` |
| B5 | BCS textbook limits: Δ→0 branch projectors, Cauchy–Schwarz 0≤C_I≤1; scattering-only ratio ≤1 in every variant; C_I+C_pair=1 per pair (different sheets: no normal-state identity) | tests T1 (Δ→0 limit), T4 (ratio bound) | 3·10⁻⁶ / 100.0000 % at Δ=0 | `python validate_ua.py` |
| B6 | internal independence: brute-force (k,helicity)×(k′,helicity) lattice sum vs energy-space pipeline | T3, T7a | 0.4–1.6 % for 0.5≤T/Tc≤0.999; T7a: 0.2–1.3 %, Nk-independent floor | `python validate_ua.py`, `python dualspace_check.py` |
| B7 | independent 2026 dual-formulation precedent + confirmation of the signed coherence factor: Sardwal–Kreisel–Mukherjee, arXiv:2608.19079 — real-space BdG ↔ k-space ↔ Walker–Smith–Samokhin analytic (PRB 65, 014517); their coherence factor (u\*u−v\*v)² = ½[1+(ξξ′−ΔΔ′)/(EE′)] (Eqs. 16/A3) = our signed C_I; their PP-channel weight (1−f_m−f_n) = our S_pair, "negligible for homogeneous superconductors" | T7 triangle | k↔E ≤1.3 %; triangle ≤1.4 %; analytic leg = Levy 3e-6 (T5) | `python dualspace_check.py` |

**How to present B1–B3 to a referee (one slide):** "Our machinery has three limits in which
it must reduce to celebrated results: (i) s-wave + ql<1 + ω→0 → Levy's 2f(ε₀); (ii) Case-II
vertex → the Mattis–Bardeen kernels identically; (iii) Case-II rate → the Hebel–Slichter
peak. All three are reproduced by the *unchanged* production code. A defect in the
coherence factor, the DOS measure, the energy conservation or the branch bookkeeping would
break at least one of them."

---

## Part II — Likely referee questions, with answers

**Q1. Why do you introduce a branch index σ in addition to the helicity s? Is it physical?**
A. The conversion of the k-sums into energy integrals requires counting states at fixed
quasiparticle energy E. The map ε→E(ε)=√((ε−μ)²+Δ̄(ε)²) is two-to-one (a valley with minimum
E_min at ε*), so states at energy E come from two disjoint band-energy pockets,
ε_σ(E)=ε*+σ√(E²−E_min²)/κ. σ labels them; in the BCS limit it is exactly sign(ξ), the ± the
textbooks hide in √(E²−Δ²). It is bookkeeping, but *necessary* bookkeeping: C_I contains
ξξ′, whose sign differs between channels (values 0.007…0.897 for the four (σ,σ′) pairs at
one (E,E′) in the master set), the Δ→0 limits are projectors δ_{σσ′} (scattering) and
δ_{σ,−σ′} (pair creation), and the band cutoffs E_max(σ) are per-branch. Figure
`branch_structure_sigma_vs_s.png`; derivation §8 of `DERIVATION_HHM_scattering.md`.

**Q2. Do you sum over helicities s, s′ anywhere? Your loops show only two indices.**
A. Yes, but that sum is done analytically once: ξ, Δ̄, E depend on a state only through its
band energy ε, so C_I is blind to s and the s-sum merely adds DOS weights,
Σ_s N_s(ε)=N(ε), which `dos_ua.NormalDOS` builds from both Rashba bands. The explicit
double loop is over (σ,σ′) (variables now named `sig, sigp`). The brute-force cross-check
keeps s explicit in the state list and bins by σ for the coherence factor; the two
formulations agree to ~1 % (T3).

**Q3. Energy conservation is imposed on quasiparticle energies E′=E+ħω. Justify.**
A. The δ-function comes from the time integral of the phonon field in the quasiparticle
(α) representation of the HHM mean-field Hamiltonian (their Eq. 23–27); the eigenenergies
there are E_ks. Imposing it on band energies ε (as the original manuscript did) mixes
states of different true excitation energy and, with μ≈3.6t, ħω=0.7t, makes the normal-state
integral negative — absorbed power negative — which the old code hid with abs(). Erratum E1.

**Q4. Why is W/W_N < 1 everywhere and why does it NOT reach 1 at Tc in your scattering-only
curves, while BCS intuition says attenuation ratios tend to 1?**
A. Two statements. (i) W/W_N<1 follows from 0≤C_I≤1 pointwise (Cauchy–Schwarz, T4) plus the
DOS restructuring. (ii) The limit at Tc⁻ of the *scattering channel alone* is the Δ→0 golden
rule with the branch projectors still active: it equals R(ħω,Tc)=(2/ħω)∫₀^∞[f−f′]dE weighted
by the band DOS — 0.4987/0.5935/0.7426 (flat-DOS values) — not 1. The full absorption
includes pair creation, which is examined as an OUTLOOK only: The scattering channel alone returns R(hw,T_c) = 0.400 / 0.505 / 0.761 at T_c-; restoration of W/W_N -> 1 involves the pair-breaking channel, outside the scope of the quantitative result reported here. The coherence factors satisfy C_I + C_pair = 1 per quasiparticle pair, so together the channels carry the total transition weight; because they occupy different sheets (E'-E = hw vs E+E' = hw) with different thermal weights, their sum is not bounded by and does not equal any single normal-state reference. This work therefore reports the scattering channel only; the pair channel is examined in the outlook without a quantitative claim.
Figures `sum_rule_and_thresholds.png` (outlook panel), `W_N_limits.png` (W_WN_limits.png).

**Q5. Your validity window is 2Δ_eff > ħω, yet you say pair creation matters. Contradiction?**
A. No: inside that window W_p ≡ 0 identically (its phase space opens only at
ħω > 2E_min(T) ≤ 2Δ_eff(T)), so the scattering-only result is *complete* there — the window
is safe. But the window is not a threshold: it closes at T/Tc = 0.740/0.840/0.936 while the
true threshold 2E_min=ħω sits at 0.820/0.912/0.975; any neighbourhood of Tc lies beyond it,
which is why the Tc⁻ limit needs the pair channel. Erratum E4; appendix A of the report.

**Q6. Why should the coherence factor keep the sign of ξξ′? The published formula uses |ξ|.**
A. Because C_I derives from squaring u u′ − η v v′ with the HHM coefficients
u v* = −Δ̄/(2E) (their Eq. 26): the ξξ′ term appears with its sign. Three independent
diagnostics fail when the sign is dropped: (i) the Δ→0 limit stops being the branch
projector and W/W_N overshoots 1 (up to 1.97); (ii) the sum rule (W+W_p)/W_N=1 is violated
(overshoots 1.0924/1.3386/1.9874); (iii) the s-wave ω→0 limit no longer equals Levy's
2f(ε₀) (the |ξ| variants give the wrong edge weight). Erratum E1; tests T1, T5.

**Q7. Is the deformation-potential vertex with η=+1 the right coupling for ultrasound?**
A. For ql<1 sound couples to the density/kinetic energy (strain modulates the hopping);
in the helicity basis that vertex is helicity-diagonal and, after the Bogoliubov
transformation, carries the Case-I combination u u′ − v v′ (η=+1). The Case-II combination
u u′ + v v′ belongs to electromagnetic/NMR probes; our code contains both and reproduces
the Mattis–Bardeen/Hebel–Slichter signatures with η=−1 (B2, B3) and Levy with η=+1 (B1) —
the vertex classification is therefore validated, not assumed.

**Q8. Your normal reference W_N: which DOS, which chemical potential?**
A. The same Hamiltonian at Δ⁰=Δˢ=0: the exact Rashba DOS N(ε)=Σ_sN_s(ε) of HHM Eq. (18) with
exact band support ±4t√(1+V_SO²/2t²) (normalised ∫N dε=2; the repository's stored table had
∫g dε=0.643 and a Gaussian leak above the band edge — erratum E5), and f_N at μ_N(T) from
the normal-state number equation at the same filling. Using μ_SC instead changes ratios by
O(Δ²) (column `*_ratio_muSC` in the results); the normal-state μ is the consistent choice.

**Q9. How accurate are the integrals? Why midpoint in ℓ?**
A. The substitution ℓ=√(E²−E_min²) removes the 1/√(E²−E_min²) BCS singularity exactly
(ν_σ dE = N(ε_σ)dℓ/κ), leaving a smooth integrand; midpoint on ℓ converges as h². T0 checks
it against adaptive quadrature (2·10⁻⁴); T1 against a closed form (3·10⁻⁶); T3 against an
entirely independent lattice sum (≤1.6 %). Residual T3 deviation at low T is the different
regularisation of the gap-edge peak (Gaussian δ on E vs smoothed N(ε)), not a convergence
failure.

**Q10. Why is the brute-force cross-check only ~1 % and not 10⁻⁶?**
A. It regulates the δ-function as a Gaussian of width w on E and extrapolates w→0
(Richardson), and bins states in E; near the gap edge the 1/√ singularity makes the
bin-average and the exact measure differ at O(1 %). The energy-space value is the
convergent one (T0/T1 prove its quadrature exact). The cross-check's role is to catch
*structural* errors (wrong index bookkeeping, wrong Jacobian) — which it did, historically.
T7a sharpens this: the ~1 % offset is *independent of Nk* (160→452), i.e. it is the
NormalDOS construction (Gaussian smoothing 0.02 t) seen from the two sides, not a k-grid
artefact; and the ω→0 s-wave leg (T7b) closes the k-space/energy-space/Levy triangle at
≤1.4 %. For scale: Sardwal–Kreisel–Mukherjee (arXiv:2608.19079, Fig. 3) report
real-space-vs-k-space agreement of the same order (a few %) for their d-wave panels.

**Q11. The published curves are flat at 0.90–1.00. Yours rise from ~0. Why?**
A. The flat curves cannot come from this Hamiltonian: any golden-rule ratio with
0≤C≤1 and a gapped DOS must vanish exponentially as T→0 (f(Δ) suppression), rise
monotonically, and saturate below 1 at Tc⁻ (Q4). Flat 0.9–1.0 behaviour is the signature of
the defects listed in errata E2/E5 (band-energy energy-conservation with abs(), misnormalised
DOS). Our curves reproduce Levy/Morse/Tsuneto physics in every limit (B1–B3).

**Q12. Group B caption says Δ_t = 4 but you use 4.5. Group D auxiliary panels look like
Group A's. Explain.**
A. Both are manuscript defects, documented as errata E6: the stage-1 numbers printed for
Group B are only reproducible with Δ_t=4.5t (Tc=0.5009/0.1399/0.0861), and the Group D
auxiliary μ(T), Δ(T) panels are verbatim copies of Group A's (U=−4t) solution. Our
reproduction uses the parameters that regenerate the printed data, and states the
discrepancy explicitly.

**Q13. The U=−75t scan looks anomalous. Is it physical?**
A. No — erratum E7. Solving the HHM gap+number equations at U=−75t, n=1.885 lands on an
unphysical strong-coupling root (μ≈35t, Tc≈12.7t, far outside the band); the published scan
sits on that root. Our reproduction flags it and the figure marks the branch.

**Q14. Mean-field only? No vertex corrections, no impurities, ql<1?**
A. Yes: the scope is the HHM mean-field quasiparticles + golden rule, i.e. exactly the
framework of Levy (ql<1, collision drag) and of the manuscript. Impurity self-energies and
ql>1 (Pippard/Tsuneto regime) are outside both; stating this is part of the error budget,
not a hidden assumption. For ql>1 the Pippard kernel replaces the constant one and the
ratio changes — our benchmarks deliberately sit in the ql<1 regime where Levy's result is
valid.

**Q15. Why does the T→Tc⁻ value depend on V_SO so strongly (0.40→0.76)?**
A. R(ħω,Tc) is a DOS-weighted average over the band: V_SO reshapes N(ε) (van Hove
structure, band edge) and moves μ; the saturation value is a band-structure fingerprint,
computable in closed form in the flat-DOS limit (0.4987/0.5935/0.7426 vs the lattice values
0.4000/0.5054/0.7611 at T/Tc=0.999). This sensitivity is a prediction, testable against
material-specific DOS.

**Q16. Numerical parameters: Nk, bins, nell — converged?**
A. DOS: Nk=1400, 6000 bins, smoothing 0.02t (∫N dε=2 to 10⁻⁸); μ: Nk=300 number equation;
integrals nell=6000 (production) with T0 showing 10⁻⁴-level quadrature error; brute force
Nk=320, nbin=3000 + Richardson. Doubling nell changes ratios <10⁻⁴.

**Q17. What breaks if Δ̄(ε) changes sign inside the band (Δ_eff vs nodes)?**
A. Nothing structurally: E_min=|Δ⁰−bμ|/κ can vanish while Δ⁰,Δˢ≠0 (gap nodes in ε); the
branch algebra, cutoffs and code handle E_min→0 continuously (the ℓ-substitution remains
valid; tests at small E_min included in T1/T5 sweeps). Physical interpretation changes
(node-line superconductivity), which the report flags where parameters approach it.

**Q18a. Where is ql<1 in your code? Do you have τ or ℓ? Do you impose ω=v_s q?**
A. Nowhere — and we say so explicitly (`LEVY_REGIME_ANALYSIS.md`). The code is the clean-limit
golden rule with a q→0 density vertex: no q, τ, ℓ, v_s. Levy's ratio is reproducible without
them because in his regime α_s and α_n carry identical prefactors (ql)²/τ and an
energy-independent angular kernel, all of which cancel in α_s/α_n; and because both his
collisional solution (ql<1, ωτ≪1) and our clean q→0 identity reduce the ratio to the same
universal functional L[1] = 2∫(−f′)dE/∫(−f′)dε = 2f(Δ). Levy's own paper notes that BCS
obtained the same ratio for longitudinal waves at ql>1 by exactly our method (direct
absorption minus induced emission), so the agreement is structural across the ql axis.
Test T5b shows the operative assumption (kernel energy-independence) and delimits the
scope: at finite ħω ~ Δ the two prescriptions differ (dev 0.84 at ħω=0.7t), where our
curves are the clean-limit answer and Levy is not the applicable comparison.

**Q18b. Your pipeline lives in energy space. Is the same model ever solved in k-space, and does it agree?**
A. Yes — twice, and a 2026 paper validates exactly this protocol. T3 (validate_ua.py) sums
all (k,helicity)×(k′,helicity′) lattice pairs with a Gaussian δ and Richardson
extrapolation: 0.4–1.6 % against the energy-space pipeline. T7 (dualspace_check.py)
upgrades this to the benchmark protocol of Sardwal–Kreisel–Mukherjee, arXiv:2608.19079,
Sec. III — (i) k-space lattice sum, (ii) energy-space reduction, (iii) analytic closed
form: at production ħω=0.7t the (i)–(ii) deviation is 0.2–1.3 % with an Nk-independent
floor; in the s-wave ω→0 limit the triangle (i)↔(ii)↔(iii)=Levy closes at ≤1.4 %, the
flat-DOS leg reproducing Levy to six digits (T5). The same paper independently confirms
the *structure*: its coherence-weighted electron–phonon matrix element (Eqs. 16, A3)
reduces homogeneously to (uu′−vv′)² = ½[1+(ξξ′−ΔΔ′)/(EE′)] — the signed Case-I factor —
and its pair-creation channel (Eqs. A6–A9, weight 1−f_m−f_n) is "negligible for
homogeneous superconductors", matching our S_pair ≡ 0 inside 2Δ > ħω.
Figure: `ultrasonic/figures/dual_space_k_vs_energy.png`.

**Q18c. Your professor expects peaks or non-monotonic features in W/W_N vs U, V_SO, Δ_t. Do they exist?**
A. In the primary channel (η=+1 signed, scattering only): **no** — 1292 of 1320 scanned
curves are strictly monotone (the rest fail monotonicity only in numerically flat tails).
The reason is structural: the naive peak candidate, the gapless ridge E_min = |D0−bμ|/κ → 0,
coincides with the collapse of the SC dome (T_c → 0), so the ratio dies smoothly instead of
peaking (`PARAM_SCAN_PEAKS.md` §1). Features of height 1.32–1.58 previously reported for (S_sc+S_pair)/S_N are WITHDRAWN
(pair channel demoted to outlook, REPORT App. A2); the only surviving non-monotonicity is
the shallow scattering-only interior maximum vs U at U ≈ 17.5–22.5 for V_SO ≤ 0.2 and T/T_c ≥ 0.7 (edge max at U=15 for lower T/T_c; E10-corrected),
from the competition of falling Δ⁰₀/T_c with the falling anchor R(ħω,(T/T_c)T_c(U)); at
U=125, Δ_t* ≈ 4.0–5.1 (grows with U and T, shrinks with V_SO). They survive 13/25-point
refinement (shifts ≤0.06 in Δ_t), doubled quadrature (≤2·10⁻⁵) and an independent k-space
lattice evaluation including the new pair sum S_pair_bf (0.2–3.8 %, one 15 % outlier at the
smallest E_min). Case-II channel: monotone in these directions (its Hebel–Slichter
resonance lives in the T-direction, T6). Figures `param_scan_{U,VSO,Dt}.png`; nothing was
tuned to create peaks, and the null result for the primary channel is reported as such.
Addendum: extending the U-grid below 50 (U ∈ [15,150], `WWN_vs_U_families_VSO.png`) adds ONE
primary-channel feature — a shallow interior maximum at U ≈ 17.5–22.5 for V_SO ≤ 0.2 and T/T_c ≥ 0.7 (heights 0.60/0.87
at T/T_c = 0.5/0.7/0.9) from the competition of the falling Δ⁰₀/T_c with the falling
normal-state anchor R(ħω, (T/T_c)T_c(U)); documented in `PARAM_SCAN_PEAKS.md` §7 addendum.

**Q18d. Where is the spin–orbit effect in the scattering channel? My professor insists it must be visible.**
A. It is, strongly, in four internal routes (`SOC_EFFECT_SCATTERING.md`,
`soc_effect_scattering.png`): at fixed T/T_c the signed scattering ratio rises monotonically
with V_SO — 0.400→0.761 at T/T_c=0.999, ×1.7 at 0.5 — because Rashba splitting (i) splits the
van-Hove structure of N(ε), (ii) raises T_c ×2.27 while lowering the coupling ratio Δ⁰₀/T_c
from 0.906 to 0.500 (weaker f(Δ/T) suppression at fixed T/T_c), (iii) moves the branch
structure ε_σ(E) and the (σ,σ′) weights ((−,−) 70→66 %), and (iv) changes the (ξ, D̄) pairs
entering C_I. The relative enhancement even changes sign at T/T_c≈0.22 (−70 % below, +90 %
above): deep in the dome the SOC-grown E_min (0.47→1.22 t) cuts the scattering phase space.
At fixed physical T the sign flips again (×39 suppression at T=0.2 t) because T_c doubles —
so any statement of "the SOC effect" must name its ensemble; both are shown.

**Q18. Reproducibility?**
A. Single command per stage (`main_tc.py all`, `main_ua.py all`, `validate_ua.py`,
`levy1963_check.py`, `case2_checks.py`, `dualspace_check.py`, `plotting_ua.py`); all results
stored as JSON/CSV in
`results/`; figures regenerated from stored results only. PDFs rebuilt by
`build_corrected_pdf.py` / `build_fa_pdf.py`.

---

## Part III — Oral-exam style questions your professor may ask (with model answers)

1. *"Derive C_I from the vertex in three lines."* — u u′ − η v v′ squared with
   u²=½(1+ξ/E), v²=½(1−ξ/E), uv=−Δ̄/2E ⇒ ½[1+(ξξ′−ηΔ̄Δ̄′)/(EE′)]; η=+1 for the density
   (Case-I) vertex. Pair channel: (uv′+ηvu′)² = 1−C_I. (DERIVATION §6.)
2. *"Show me by hand that ω→0 gives 2f(Δ)."* — integrand in ℓ: N²(ℓ/E)ω(−f′) because
   C_I(E,E)=ℓ²/E² cancels the Jacobians; ∫d → ∫dE N²ω(−f′) = 2N²ωf(Δ); W_N=N²ω. (Levy.)
3. *"Why does the pair channel vanish in your window but restore 1 at Tc?"* — phase space
   [E_min, ħω−E_min] empty while ħω<2E_min; sum rule is spectral-weight conservation
   (M_scat²+M_pair²=1), so weight merely migrates between channels as Δ opens.
4. *"What is η physically, and how would you measure the Case-I/II distinction?"* — relative
   sign of the anomalous (pair-fluctuation) coupling of the probe: density/strain probes
   (sound) vs vector probes (EM, NMR). Measured: ultrasound shows monotone suppression
   (Levy), EM/NMR show the Hebel–Slichter peak — our figure reproduces both.
5. *"Why intraband pairing only?"* — HHM Sec. II: interband pairs carry finite centre-of-mass
   momentum at zero field (FFLO-like); time-reversed intraband pairing is the Anderson
   prescription. The Rashba gauge phase s e^{−iθ_k} cancels in all rates.
6. *"Prove the sum rule analytically."* — At fixed (E,E′) shell, C + (1−C) = 1 per state
   pair, so the two channels carry the total transition weight; but they live on different
   sheets (E′−E=ħω vs E+E′=ħω) with different thermal weights, so no identity with a single
   normal-state reference follows — the previously printed W_N^QP "identity" is withdrawn
   (measured 0.545/3.63/0.848, never 1). What survives: at Δ→0 the measures coincide and
   (S_sc+S_pair)/S_N = 100.0000 % (T1); the quantitative result is scattering-only.
7. *"Why μ_N ≠ μ_SC, and which enters where?"* — particle number is fixed; the gap shifts the
   density of occupied states, so the chemical potentials differ by O(Δ²). W (SC states) uses
   μ_SC; W_N (normal states at same filling) uses μ_N. Consistency checked via the
   `mu_for_WN` option (differences <10⁻³ in ratios).
8. *"Where does your Tc come from?"* — the repository's own stage-1 solver (HHM Eqs. 32–38)
   re-implemented parameterised; reproduces printed values to 10⁻⁴ (B4). φ(T) fit shape is
   the repository's convention; limits at Tc⁻ are insensitive to its details (Δ→0).
9. *"What would change with a d-wave or nodal gap?"* — C_I keeps its form (functions of ε),
   but E_min→0 on node lines and the Δ→0 projectors are replaced by nodal averages; the
   branch machinery is unchanged. Out of scope here (HHM gap is extended s-wave).
10. *"Explain the 0.647 factor you once reported."* — historical quadrature defect (missing
    final-state Jacobian E′/ℓ′), found and fixed via the Levy benchmark; correction note in
    REPORT §5. Good example of why external benchmarks matter.
11. *"Is q→0 legitimate for sound (ql<1 but q≠0)?"* — yes: Levy's regime; the q-dependence
    enters only through the Pippard kernel for ql≳1; our ω→0 benchmark matches Levy exactly
    in this regime.
12. *"If I gave you experimental α_s/α_n data for Sr₂RuO₄ or a Rashba interface, what would
    you need?"* — material DOS and μ(T), gap amplitudes from fit, ħω, and the regime (ql);
    then the pipeline is parameter-free. Prediction: saturation value <1 at Tc⁻ for the
    scattering channel, sum-rule restoration above 2E_min.

---

## Part IV — Honest limitations (state them before the referee does)

1. Mean-field HHM quasiparticles + golden rule; no impurity self-energy, no vertex
   corrections, no ql>1 Pippard kernel (same scope as Levy 1963 and the manuscript).
2. Gap shape φ(T) is the repository's fitted form; Tc⁻ limits are Δ→0 quantities and do not
   depend on it, but mid-T curve shapes do (weakly).
3. Lattice vs energy-space agreement is 0.4–1.6 % at finite Δ (regularisation of the gap-edge
   peak), not machine precision; the energy-space side is the quadrature-exact one (T0/T1).
4. Group-B Δ_t typo and Group-D copied panels, U=−75t unphysical root: manuscript defects we
   document, not physics we defend.
5. The 2026-10-15 measure correction (final-state Jacobian) is documented openly in REPORT
   §5; all numbers after that date include it, and the Levy/Mattis–Bardeen/Hebel–Slichter
   benchmarks pass only with it.
6. The Levy benchmark tests the universal q→0, ω→0 core only: ql<1, finite τ and ω=v_s q are
   not implemented (no q, τ, ℓ in the code). At finite ħω ~ Δ the clean golden rule differs
   from Levy's transport prescription (T5b); there our results are the clean-limit
   (BCS-type) answer. Full regime analysis: `LEVY_REGIME_ANALYSIS.md`.
7. The q->0 vertex is taken isotropic (inter-helicity matrix elements = intra-helicity),
   as in the paper, Walker et al. and Sardwal et al.'s homogeneous limit. A strict q->0
   density/strain vertex is helicity-diagonal; recomputing with it (R_B,
   `soc_explicit_vertex.png`) shifts W/W_N by -12…-20 % at V_SO=0.8 but leaves the SOC
   trend (x1.70 vs x1.75 at T/T_c=0.5) intact: the SOC effect is physical, its magnitude
   carries a vertex systematic quantified in `SOC_EFFECT_SCATTERING.md` §4. A finite-q
   spinor-overlap vertex is future work.
8. The dual-formulation (k-space vs energy-space) agreement is 0.2–1.3 %; the floor is set
   by the Gaussian smoothing (0.02 t) of the NormalDOS construction, seen from the two
   sides (T7a: deviation independent of Nk, 160→452). A raw lattice-histogram DOS in both
   legs would push it lower — a numerical refinement, not a structural issue. T7.
