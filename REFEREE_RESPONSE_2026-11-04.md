# Evidence-backed response to the referee

**Subject:** Power absorption and ultrasonic attenuation in a Rashba spin–orbit-coupled superconductor
**Response date:** 5 November 2026
**Status:** Draft author response with new repository checks. This report distinguishes measured/code evidence from physical claims that remain unproved.

## Executive response

We thank the referee for identifying the main gap between a numerically evaluated absorption model and a quantitatively established physical ultrasonic response. We have now added executable evidence for the repository’s internal calculation and have completed a separate 52-temperature comparison of the legacy fitted-gap calculation against the optional self-consistent HHM solver for all three master SOC values.

That evidence supports a specific, narrower conclusion: **the repository is on a credible numerical path for calculating the stated finite-frequency, scalar-vertex quasiparticle-scattering model.** It does **not** prove that this model’s unrestricted transition vertex is the physical finite-q acoustic vertex, nor that its finite-frequency results describe experimental ultrasonic attenuation. We therefore accept **major revision** and will not describe the existing ratios as quantitatively validated material predictions until the finite-q response and transport regime are established.

The headline observable remains the scattering contribution \(W_{\rm sc}/W_N\), not total \(W/W_N\). Residual-verified gap roots are not thereby proven thermodynamically stable.

### Evidence labels

- **VERIFIED HERE:** command or calculation run in this workspace, with outputs identified below.
- **SOURCE EVIDENCE:** directly visible in the cited implementation/documentation; not by itself a validation of the physics assumption.
- **OPEN:** evidence does not yet establish the claim.

## 1. Evidence for the numerical road already implemented

### 1.1 Reproducible test run

The default environment initially lacked NumPy, so the first `make test` attempt stopped at `ModuleNotFoundError: No module named 'numpy'`. I installed dependencies outside the repository and then ran the complete suite:

```bash
PYTHONPATH=/tmp/arena-referee-deps make test
```

Environment used: Python 3.11.2; NumPy 2.4.6; SciPy 1.17.1; Matplotlib 3.11.2. **VERIFIED HERE:** the command exited 0 and printed `MAKE TEST: all checks passed`. In detail, the run included the independent BdG-kernel check; `test_fixes.py` (21 PASS, 0 FAIL); the HHM finite-temperature tests (6 PASS); the temperature-maxima tests (6 PASS); `validate_ua.py` (T0–T3 PASS); and the Case-II, Levy, and dual-space checks.

This is evidence for numerical consistency and regressions covered by those tests. It is not independent verification of every stored result or of the physical vertex.

### 1.2 Analytic reduction actually tested by the code

For \(b=\Delta^s/(4t)\), the model uses

\[
E(\epsilon)^2=(\epsilon-\mu)^2+(\Delta^0-b\epsilon)^2
=\kappa^2(\epsilon-\epsilon_*)^2+E_{\min}^2,
\]
\[
\kappa=\sqrt{1+b^2},\qquad
\epsilon_*={\mu+b\Delta^0\over1+b^2},\qquad
E_{\min}={|\Delta^0-b\mu|\over\kappa}.
\]

Therefore the two band-energy branches at quasiparticle energy \(E\) are

\[
\epsilon_\sigma(E)=\epsilon_*+\sigma {\sqrt{E^2-E_{\min}^2}\over\kappa},\quad
\nu_\sigma(E)=N(\epsilon_\sigma(E)){E\over\kappa\sqrt{E^2-E_{\min}^2}}.
\]

With \(\ell=\sqrt{E^2-E_{\min}^2}\), \(dE=\ell\,d\ell/E\), so
\[
\nu_\sigma(E)\,dE=N(\epsilon_\sigma(E))\,{d\ell\over\kappa}.
\]
This is the Jacobian cancellation implemented in [`absorption.py`](project/ultrasonic/absorption.py#L232-L279). The production integrand then evaluates

\[
S_{\rm sc}=\sum_{\sigma,\sigma'}\int dE\,
\nu_\sigma(E)\nu_{\sigma'}(E+\hbar\omega)C_I
[f(E)-f(E+\hbar\omega)],\qquad
{W_{\rm sc}\over W_N}={S_{\rm sc}\over S_N}.
\]

In code the upper energy is set to `Ep = E + omega`; all four energy-inversion branch pairs are summed, and the lower and upper DOS/Jacobian factors are treated separately. The normal denominator has its own \(\mu_N(T)\) and a kink-aligned Gauss–Legendre rule whose order is checked at runtime ([`absorption.py`](project/ultrasonic/absorption.py#L282-L353)).

For the scalar density-scattering assumption, the signed Case-I coherence is

\[
C_I=|uu'-vv'|^2
={1\over2}\left[1+{\xi\xi'-\bar\Delta\bar\Delta'\over EE'}\right].
\]

`coherence()` retains the sign of both \(\xi\) factors in variant `a_signed` ([`absorption.py`](project/ultrasonic/absorption.py#L160-L185)). The test independently diagonalizes BdG matrices and compares eigenvector matrix elements with these closed forms: maximum errors are \(7.772\times10^{-16}\) for Case I, \(6.661\times10^{-16}\) for Case II, and \(8.882\times10^{-16}\) for scattering-plus-pair matrix-element completeness. This is strong algebraic evidence **conditional on the scalar vertex**.

### 1.3 Independent numerical checks and their limits

The full run returned these quantitative checks:

| Check | Measured result in the executed suite | What it supports |
|---|---:|---|
| \(\ell\)-integral vs explicit energy integral (T0) | worst relative deviation \(3.89\times10^{-4}\) | branch/Jacobian implementation |
| Flat-DOS, \(\Delta\to0\) closed form (T1) | worst deviation \(3.05\times10^{-6}\) | signed kernel and thermal/normalization path in this limit |
| BCS, \(\omega\to0\): \(W/W_N\to2f(\Delta)\) (T2) | worst relative deviation \(2.01\times10^{-6}\) | low-frequency limiting algebra |
| Energy-space vs lattice sum for signed kernel (T3) | relative deviation \(2.03\times10^{-3}\) at \(T/T_c=0.9\), maximum \(2.59\times10^{-3}\) on the reported grid | implementation of the same unrestricted-transition model |
| Normal-reference quadrature | order 2/4/8 relative spread \(1.6\times10^{-16}\) | denominator quadrature at tested state |
| Levy reduced s-wave benchmark | worst relative deviation \(1.1\times10^{-5}\) in the PHYS check | low-frequency benchmark only |

The lattice check is useful but not an independent physical-vertex derivation: the source explicitly sums unrestricted \((\mathbf k,s)\times(\mathbf k',s')\) transitions and calls this the dirty-limit/constant-matrix-element model ([`bruteforce_ua.py`](project/ultrasonic/bruteforce_ua.py#L1-L18)). Consequently, its agreement with energy space validates the transformation and numerics within that assumption, not the assumption itself.

## 2. Point-by-point response to the referee

### Major concern 1 — Is the physical absorption vertex established?

**Answer: not yet.** We agree with the referee’s central objection. The derivation writes a density operator with phonon momentum, but the production energy-space integral has no phonon momentum and the lattice implementation sums over all initial/final momenta. The coherence-factor derivation establishes the result for the assumed scalar transition matrix element; it does not derive the acoustic vertex from the Rashba Hamiltonian.

The physically targeted calculation should start from a finite-q perturbation, for example
\[
H_{ep}=\sum_{\mathbf k,\mathbf q}g_{\mathbf q}(\mathbf k)
(b_{\mathbf q}+b^\dagger_{-\mathbf q})
 c^\dagger_{\mathbf{k+q}}c_{\mathbf k},
\]
and evaluate the response from the corresponding BdG matrix elements,
\[
W(\mathbf q,\omega)\propto \omega[-\operatorname{Im}\chi^R_{\Gamma\Gamma}
(\mathbf q,\omega)],\quad
\chi^R_{\Gamma\Gamma}\sim\sum_{m,n,\mathbf k}
{f(E_{m\mathbf k})-f(E_{n,\mathbf{k+q}})\over
\hbar\omega+E_{m\mathbf k}-E_{n,\mathbf{k+q}}+i0^+}
|\langle n,\mathbf{k+q}|\Gamma_{\mathbf q}|m,\mathbf k\rangle|^2.
\]
Here \(\Gamma_{\mathbf q}\) must include the strain/deformation-potential vertex, Rashba helicity spinors, and the Nambu coherence factors. For normalized Rashba spinors with unit texture \(\hat{\mathbf g}_{\mathbf k}\), the spin-only overlap obeys
\[
|\langle u_{s',\mathbf{k+q}}|u_{s,\mathbf k}\rangle|^2
={1+s s'\hat{\mathbf g}_{\mathbf{k+q}}\cdot\hat{\mathbf g}_{\mathbf k}\over2}.
\]
At \(q=0\), this becomes \(\delta_{s's}\) away from helicity degeneracies. It is not generally the unrestricted constant weight used by the current production integral.

There is also an exact conservation check. If the Hamiltonian conserves total particle number \(N\), then
\[
(E_n-E_m)\langle n|N|m\rangle
=\langle n|[H,N]|m\rangle=0.
\]
Thus at exactly \(q=0\), for a finite-frequency transition \(E_n-E_m=\hbar\omega\ne0\), the total-density matrix element vanishes. A broken-symmetry BdG calculation must implement a gauge-conserving response (including phase/collective and screening corrections as appropriate) to recover the corresponding conservation law; a bare bubble alone is not enough to settle it.

**Code evidence:** [`S_scattering`](project/ultrasonic/absorption.py#L232-L279) accepts `T, omega, D0, DS, mu`, not \(q\), \(v_s\), or \(\tau\). It sets `Ep = E + omega`, and sums energy branches. [`bruteforce_ua.py`](project/ultrasonic/bruteforce_ua.py#L1-L18) explicitly states that its \(\mathbf k,\mathbf k'\) sums are unrestricted. This confirms the referee’s identified gap.

**Evidence for the right numerical road:** conditional on this effective vertex, the BdG coherence test, Jacobian test, and independent lattice/energy-space agreement above show the response integrator is not arbitrary: it has a well-defined matrix-element origin in the scalar model and reproduces its limiting cases.

**Required next code, not claimed as implemented:** construct helicity eigenspinors on the same k-grid; evaluate the finite-q overlap and Nambu matrix element for each \((\mathbf k,s)\to(\mathbf{k+q},s')\); enforce \(E'-E=\hbar\omega\) in scattering and \(E+E'=\hbar\omega\) in pair creation; include the selected phonon polarization and screened/disorder vertex; then test the Ward/conservation limit. Until that implementation is independently checked, the physical vertex remains **OPEN**.

### Major concern 2 — Is the ultrasound transport regime implemented?

**Answer: no; the code is a finite-frequency quasiparticle-transition calculation, not a transport-regime solver.** The active parameter file sets `OMEGA = 0.7*T_HOP` ([`parameters_ua.py`](project/ultrasonic/parameters_ua.py#L13-L25)). There are no corresponding `q`, sound velocity \(v_s\), mean free path \(\ell\), or relaxation time \(\tau\) in the absorption API. In a relaxation-time description,
\[
q\ell=qv_F\tau={\omega\tau v_F\over v_s},
\]
so the condition \(q\ell\ll1\) cannot be evaluated without material parameters and a transport model. Setting a free \(\omega\) in a transition integral does not impose \(\omega=v_s q\).

**Positive but limited evidence:** the Levy check reduces the model to flat DOS, s-wave, and small \(\omega\), where the ratio approaches \(2f(\Delta)\); the executed test agrees to \(1.1\times10^{-5}\). This verifies a valuable low-frequency limit. It does not validate the finite \(\hbar\omega=0.7t\) response or any clean/dirty crossover.

**Required next step:** add explicit acoustic dispersion, \(q\), phonon polarization, and either a clean finite-q response or a justified disorder/relaxation treatment; report \(q\ell\), \(\omega\tau\), and the energy/frequency scale relative to the gap. Until then, we will qualify the outputs as model results and not call them material-realistic ultrasonic predictions.

### Major concern 3 — Does the reported ratio include all absorption channels?

**Answer: no. The quantitative headline is scattering-only.** In the code,
\[
S_{\rm sc}\propto C_I[f(E)-f(E')],\quad E'=E+\hbar\omega,
\]
whereas a candidate pair contribution uses
\[
S_{\rm pair}\propto C_{\rm pair}[1-f(E)-f(E')],\quad E+E'=\hbar\omega.
\]
The pair integral is kinematically closed for \(\hbar\omega\le2E_{\min}\). [`S_paircreation`](project/ultrasonic/absorption.py#L356-L423) implements this threshold check and returns zero below threshold. `ratio_all()` itself returns the scattering ratio and does not add the pair term ([`absorption.py`](project/ultrasonic/absorption.py#L329-L353)).

This code is useful as a candidate-channel calculation within the same assumed scalar vertex; it is not a derivation of the physical total acoustic response. Threshold accessibility near \(T_c\) does not prove a large rate. The pair matrix elements, phase space, physical vertex, screening, and a consistent normal-state normalization must all be checked before claiming total \(W/W_N\).

**Revision:** label tables and plotted outputs \(W_{\rm sc}/W_N\). Do not infer total absorption from the threshold or the per-transition identity \(C_{\rm sc}+C_{\rm pair}=1\); that identity is not a finite-gap sum rule for the integrated ratio.

### Major concern 4 — How large is the SOC vertex sensitivity?

**Answer: it is measurable and must be reported as a systematic, not dismissed.** The referee’s A/B comparison is also stored in [`soc_explicit.json`](project/ultrasonic/results/soc_explicit.json). Values extracted from that tracked result at \(T/T_c=0.9\) are:

| \(V_{\rm SO}/t\) | Isotropic A, \(R_A\) | Helicity-diagonal B, \(R_B\) | \(100(R_B/R_A-1)\) |
|---:|---:|---:|---:|
| 0 | 0.252165 | 0.252165 | 0.00% |
| 0.4082 | 0.318251 | 0.284325 | −10.66% |
| 0.8 | 0.485714 | 0.424987 | −12.50% |

These values are directly present in the repository artifact; the percentages are calculated from its values. At \(V_{\rm SO}=0\), A and B coincide, which is a useful control. At finite SOC the difference is not negligible. The positive SOC trend remains in both approximations, but that does not identify the physical vertex. The B calculation suppresses inter-helicity transitions according to a strict small-q prescription; it is still not a full finite-q response with phonon kinematics and gauge-conserving corrections.

**Revision:** show the A/B spread and call it a vertex-model sensitivity. The matrix elements must be derived from the actual acoustic coupling before one curve is selected as the physical prediction.

### Major concern 5 — Are the self-consistent roots stable phases?

**Answer: stability is not established.** A self-consistency solver finds stationary solutions, schematically \(\partial\Omega_{MF}/\partial\Delta_i=0\) together with the number equation. A small residual proves only that those equations are solved. It does not prove a local minimum or the globally preferred phase.

At fixed filling, competing solutions should be compared using the Helmholtz free energy
\[
F_{MF}(N,T)=\Omega_{MF}(\mu,T)+\mu N,
\]
with the same mean-field functional, including its interaction/double-counting terms, and with \(\mu\) chosen for the same \(N\). The Hessian and branch continuity should also be checked. We do not substitute \(|\Delta_{\rm eff}|\) for this comparison.

**Code evidence:** [`hhm_finite_t.py`](project/sc_params/hhm_finite_t.py#L38-L48) folds an exact global gap-sign symmetry and checks the residual; its module and solver documentation explicitly say root status is not a free-energy/stability judgment. The six executed HHM tests verify the linearized \(T_c\) benchmark, equation residuals, continuation, and adapter; they do not calculate \(F_{MF}\). Thus the repository demonstrates a careful route to finding roots but not their thermodynamic ranking.

**Revision:** call these “residual-verified self-consistent solutions” and retain the stability question as open until the fixed-filling free-energy comparison is implemented.

### Major concern 6 — Does a self-consistent finite-T solver change the production results?

**New evidence generated for this response:** yes, a full sampled comparison has now been run without overwriting production outputs.

The legacy path uses
\[
\Delta^0(T)=\Delta^0(0)\phi(T),\quad \Delta^s(T)=\Delta^s(0)\phi(T),\quad
\phi(T)=[1-(T/T_c)^{3.434}]^{0.553}.
\]
The optional HHM path determines \(T_c\) from \(\det M(T_c)=0\) and solves the nonlinear gap and number equations along a continuation path ([`hhm_finite_t.py`](project/sc_params/hhm_finite_t.py#L63-L103), [`self_consistent_hhm.py`](project/ultrasonic/self_consistent_hhm.py#L16-L63)). Both states were evaluated with the **same** DOS construction, frequency \(\hbar\omega=0.7t\), `absorption.ratio_all`, signed coherence variant, normal-state reference convention, and `nell=6000`. This isolates the temperature-state prescription within the existing response model.

I added the non-destructive runner [`referee_evidence_master_compare.py`](project/ultrasonic/benchmarks/referee_evidence_master_compare.py). It uses \(U=95t\), \(\Delta_t=4.5t\), \(n=1.885\), \(V_{\rm SO}/t=0,0.4082,0.8\), and the same 52 fractions as the production sweep. It writes a separate 156-row CSV (52 per SOC) and summary JSON under `project/ultrasonic/results/referee_evidence/`; it refuses to overwrite those files. The CSV contains \(T_c\), physical \(T\), \(\Delta^0\), \(\Delta^s\), effective gap, \(\mu_{SC}\), \(\mu_N\), root residual, band flag, and both \(W_{\rm sc}/W_N\) values at every sampled fraction.

**Measured results from the 52-point comparison:**

| \(V_{\rm SO}/t\) | \(T_c\), fitted legacy | \(T_c\), HHM linearized | \(\Delta T_c/T_c\) | max \(|\Delta R|\) across 52 points | max relative \(|\Delta R/R|\) | max HHM residual |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.278211481 | 0.278146613 | −0.02332% | 0.003610 | 3.602% | \(8.9\times10^{-16}\) |
| 0.4082 | 0.373772706 | 0.373640483 | −0.03538% | 0.004794 | 3.932% | \(7.8\times10^{-16}\) |
| 0.8 | 0.646687589 | 0.646676682 | −0.00169% | 0.010087 | 4.650% | \(9.7\times10^{-15}\) |

Here \(R=W_{\rm sc}/W_N\), and the relative change is \(100(R_{HHM}-R_{fit})/R_{fit}\). Maximum relative differences occur around \(T/T_c=0.36,0.38,0.40\) for SOC 0, 0.4082, 0.8 respectively. Near \(T/T_c=0.999\), the fitted/HHM ratios are 0.386668/0.385633, 0.483500/0.482067, and 0.743388/0.742425. All 156 HHM points in this comparison remained in-band; the maximum gap/number-equation residual was below \(1.0\times10^{-14}\). The two \(T_c\) definitions are very close for these sets, but the finite-T gaps and response are not identical; the response discrepancy reaches about 4.7% at intermediate temperature.

The full machine-readable results are [`master_fit_vs_hhm.csv`](project/ultrasonic/results/referee_evidence/master_fit_vs_hhm.csv) and [`master_fit_vs_hhm_summary.json`](project/ultrasonic/results/referee_evidence/master_fit_vs_hhm_summary.json). The script can be rerun with:

```bash
PYTHONPATH=<installed-dependencies> python3 project/ultrasonic/benchmarks/referee_evidence_master_compare.py
```

**Interpretation:** this is direct evidence that the repository’s optional solver can be coupled to the absorption pipeline, and it quantifies the legacy-fit sensitivity for the actual three master SOC sets. It answers the requested comparison within the current scalar response model. It does not resolve the physical vertex, finite-q transport, pair-channel completeness, or branch stability.

## 3. Reproducibility and technical corrections

The referee correctly noted that lower bounds are not exact dependency pins. `requirements.txt` has now been corrected to say “minimum supported versions,” and the F-16 regression-test label/output was changed from “pinned” to “listed in requirements”; the test checks package presence, not version locking. This is a documentation/test-description correction, not a numerical result.

The new comparison artifacts have these SHA-256 checksums:

```text
5e0825162ec9a5684345dc13604a4de9d5f981d14f80ffda831072752e61fb73  project/ultrasonic/benchmarks/referee_evidence_master_compare.py
ecd810bb56c8c58bf703b81736e7064e5122ff4327a9c017fc766441a9c47703  project/ultrasonic/results/referee_evidence/master_fit_vs_hhm.csv
4d38587cc123f4dcddf5a365307fd1d1ed2a752de0062a49fd33648bab89c09a  project/ultrasonic/results/referee_evidence/master_fit_vs_hhm_summary.json
```

These new files are a separate evidence snapshot and have not yet been added to the repository’s central regeneration manifest. The existing approximately 1% DOS-smoothing systematic remains an uncertainty; the small quadrature and kernel-test errors do not remove it.

## 4. Proposed revised scientific claim

> We calculate the normalized quasiparticle-scattering contribution \(W_{\rm sc}/W_N\) for a specified finite-frequency, unrestricted-transition model with a signed scalar coherence factor. The two-branch energy-space calculation, BdG matrix-element check, limiting-case tests, and independent lattice evaluation support the internal numerical implementation of this model. A separate 52-point-per-SOC comparison shows that the optional nonlinear HHM finite-temperature solver changes the legacy fitted response by at most about 4.7% on the sampled master sweeps. These checks do not derive or validate the finite-q electron–phonon vertex of a material-specific acoustic perturbation. The present implementation does not impose \(\omega=v_s q\), model \(q\ell\) or a relaxation-time crossover, or establish thermodynamic stability of all roots. Pair creation is not included in the quantitative headline. The results are therefore model-dependent quasiparticle-scattering calculations, not yet quantitatively validated predictions of experimental ultrasonic attenuation.

## 5. Referee decision and remaining work

The evidence establishes that the project is **on the right computational road for a clearly specified effective absorption model**: the equations are explicit, limiting cases are executable, key integrals agree across formulations, and the temperature-state prescription has now been compared point-by-point against an optional self-consistent solver. That is meaningful progress—not proof of physical completeness.

The decisive next work is:

1. Derive and implement the finite-q electron–phonon/Rashba vertex with the correct spinor and BdG matrix elements; verify the number-conservation/Ward limit.
2. Add the acoustic dispersion and transport/disorder parameters, or restrict the scientific claim to the phenomenological finite-frequency model.
3. Compare competing roots using the fixed-filling mean-field free energy.
4. Derive and validate the pair channel before claiming total absorbed power.
5. Quantify the DOS-smoothing uncertainty and add the new evidence outputs to a commit-specific provenance manifest.

**Final assessment:** internal numerical evidence is strong for the stated model; physical validation as real ultrasonic absorption remains open primarily because the finite-q vertex and transport regime are not implemented. No honest report can claim that the current repository proves the latter.

## Verification boundary

The test suite and three-SOC comparison were executed in this workspace on 5 November 2026. The comparison script and its two outputs are newly generated and separate from the pre-existing production data. The finite-q calculation, free-energy ranking, experiment comparison, and total-response validation were not run or claimed here.
