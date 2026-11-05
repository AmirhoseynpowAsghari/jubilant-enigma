# Finite-q acoustic coupling: read-only scientific assessment

**Branch scope examined:** `main`, commit `39621eead0cb4060dbafb335e8696d6065023941` (2026-11-05).  
**Code changes:** none.  
**Numerical finite-q predictions:** none; the microscopic coupling is not sufficiently specified to calculate them without inventing parameters.

## Executive conclusion

The referee’s objection is **valid for the numerical production calculation**. The derivation introduces a momentum-resolved density operator

$$\rho_{\mathbf q}=\sum_{\mathbf k,s}c^\dagger_{\mathbf{k+q},s}c_{\mathbf k,s},$$

which connects $\mathbf k$ to $\mathbf k+\mathbf q$. But the golden-rule equations and production implementation subsequently use independently summed initial and final states with only energy conservation. The production kernel has no $\mathbf q$, phonon branch, $v_s$, or $\mathbf k' = \mathbf k+\mathbf q$ restriction. Its lattice cross-check explicitly describes its momentum sums as unrestricted.

There is a second, related gap. The HHM Hamiltonian determines the electronic bands and pairing solution; it does not specify how those parameters change under an ionic displacement. A deformation-potential density coupling, strain-modulated hopping/correlated hopping, and strain-modulated Rashba coupling are distinct microscopic contributions. Their relative sizes and phases require derivatives of the electronic Hamiltonian with respect to strain or an underlying electron-ion potential. They cannot be inferred from the HHM electronic Hamiltonian alone.

Thus the current calculation is reproducible as a **constant-matrix-element, unrestricted-transition model**, but the repository does not establish that this is the finite-q acoustic response of the stated Rashba system. No corrected numerical calculation or revised SOC peak search is justified from the available physical inputs.

## 1. What the repository sources establish

### 1.1 HHM specifies electronic structure, not the phonon coupling

The repository’s transcription of HHM gives nearest-neighbor hopping, on-site interaction, correlated hopping, and Rashba SOC; see [`DERIVATION_HHM_scattering.md`, §1–2](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/DERIVATION_HHM_scattering.md#L10-L45). It explicitly says HHM stops at thermodynamics, $T_c$, and the gap; the ultrasound perturbation is an external addition ([§1 opening](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/DERIVATION_HHM_scattering.md#L1-L8)). The HHM parameter set does not supply an ionic displacement field, phonon eigenvectors, sound velocities, mass density, or strain derivatives of $t$, $\Delta_t$, or $V_{SO}$.

The original paper PDF vendored in the repository is the manuscript itself: `project/paper_correction/vendor/arXiv_2607.02812_original.pdf`. In its printed discussion on pp. 4–5, the transition matrix element is denoted $B(\mathbf k\alpha|\mathbf k'\alpha')$; the text says it averages over the directions of $\mathbf k,\mathbf k'$ and spin and then treats $B$ as constant. It replaces the momentum sums by a double energy/DOS integral. Those pages do not derive that constant from a phonon displacement or impose $\mathbf k'=\mathbf k+\mathbf q$.

### 1.2 The later derivation writes momentum selection, but does not retain it in the rate

The derivation adds a sound perturbation and defines $\rho_{\mathbf q}$ with the $\mathbf{k+q}$ electron operator ([§5, lines 102–114](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/DERIVATION_HHM_scattering.md#L102-L114)). It also writes $\mathbf k'=\mathbf k+\mathbf q$ in a pair-channel expansion ([lines 177–186](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/DERIVATION_HHM_scattering.md#L177-L186)). Those are operator-level statements.

However, the scattering golden-rule equation is written as $\sum_{\mathbf k,s}\sum_{\mathbf k',s'}$ with $\delta(E'-E-\hbar\omega)$ and no momentum delta or fixed-q condition ([§7, lines 155–169](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/DERIVATION_HHM_scattering.md#L155-L169)). The text later declares the two quasiparticle energies independently integrated over an energy shell and calls that a $q\to0$ phase-space convention ([lines 242–257](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/DERIVATION_HHM_scattering.md#L242-L257)); this is an assumption/convention, not an implementation of the operator’s momentum constraint.

### 1.3 Production code enforces energy conservation, not momentum conservation

In [`absorption.py::S_scattering`](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/absorption.py#L232-L279), the inputs are DOS, $T$, $\omega$, gaps, chemical potential, and energy-branch choices. It sets `Ep = E + omega`, loops over the two energy-inversion branches `sig, sigp`, and multiplies two DOS factors, the coherence factor, and $f(E)-f(E')$. There is no $\mathbf k$, $\mathbf q$, or phonon branch in this routine. `ratio_all` uses this integral and a normal-state DOS convolution ([lines 329–353](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/absorption.py#L329-L353)).

The lattice cross-check is independent of the energy inversion but not momentum-conserving. Its docstring says the two momentum sums are unrestricted and calls this the paper’s “dirty-limit / constant-B approximation” ([`bruteforce_ua.py`, lines 1–18](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/bruteforce_ua.py#L1-L18)). Its binned kernel uses a Gaussian for $E_a-E_b-\hbar\omega$ and multiplies all pairs of state-bin counts ([lines 101–125](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/bruteforce_ua.py#L101-L125)); the driver groups by energy branch, not by $\mathbf k$ and $\mathbf{k+q}$ ([lines 128–148](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/bruteforce_ua.py#L128-L148)). The Gaussian delta is an energy delta only.

The active parameter is a free `OMEGA = 0.7*T_HOP` ([`parameters_ua.py`, line 16](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/parameters_ua.py#L13-L19)); it is not related to $q$ by $\omega_{\mathbf q}=v_s q$.

### 1.4 What the repository quotes from the cited finite-q theory

[`dualspace_check.py`, lines 13–29](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/dualspace_check.py#L13-L29) transcribes Sardwal–Kreisel–Mukherjee’s attenuation formula and vertex as

$$F_{mn}(\mathbf q)=\sum_{ij\sigma}\delta t_{ij}\,\phi^*_{ij}(\mathbf q)
\left(u^{m*}_{i\sigma}u^n_{j\sigma}-v^n_{j\sigma}v^{m*}_{i\sigma}\right).$$

As quoted there, this contains the phonon-mode factor $\phi_{ij}(\mathbf q)$ and bond modulation $\delta t_{ij}$ explicitly; the homogeneous limit gives the familiar $u u'-v v'$ coherence factor. This supports the distinction between a finite-q bond vertex and a constant all-to-all $B$. It also makes clear which inputs a real-space implementation needs. The local repository contains this transcription, not the cited paper itself; I have not independently verified the external paper’s full derivation in this read-only pass. The same file says its analytic Walker–Smith–Samokhin comparison is a homogeneous-superconductor benchmark; that benchmark checks a limiting coherence/attenuation result, not the HHM-specific electron-ion derivatives. The cited expression therefore supports the structural point, but neither citation as represented locally supplies HHM-specific $\delta t_{ij}$ or phonon data.

## 2. Which acoustic coupling follows from the stated model?

**There is no unique answer from the supplied Hamiltonian.** A general linear strain coupling is obtained by differentiating the electronic Hamiltonian with respect to strain:

$$\delta H_{ep}=\sum_{ab}\frac{\partial H_{HHM}}{\partial\epsilon_{ab}}\epsilon_{ab}(\mathbf r,t),\qquad
\epsilon_{ab}=\frac12(\partial_a u_b+\partial_b u_a).$$

The derivative can contain distinct operators, schematically,

$$\delta H_{ep}=\delta H_t+\delta H_{\Delta_t}+\delta H_R+\delta H_{def}+\cdots,$$

where $\delta H_t$ modulates ordinary hopping, $\delta H_{\Delta_t}$ modulates correlated hopping, $\delta H_R$ modulates Rashba SOC, and $\delta H_{def}=D_{ab}\epsilon_{ab} n$ is a possible local deformation-potential term. These terms are not interchangeable. Their coefficients require microscopic input such as $\partial t/\partial\epsilon_{ab}$, $\partial\Delta_t/\partial\epsilon_{ab}$, $\partial V_{SO}/\partial\epsilon_{ab}$, and the phonon polarization/eigenvector. None is specified by the HHM electronic Hamiltonian in this repository.

For example, if one **chooses** a bond-length-dependent nearest-neighbor hopping model, then for a bond $\boldsymbol\delta$ and acoustic displacement $\mathbf u_{\mathbf q\lambda}$,

$$\delta t_{\boldsymbol\delta}(\mathbf q,\lambda)
\simeq 2i\frac{\partial t_{\boldsymbol\delta}}{\partial r_{\boldsymbol\delta}}
(\hat{\boldsymbol\delta}\!\cdot\mathbf e_{\mathbf q\lambda})
\sin\!\left(\frac{\mathbf q\!\cdot\boldsymbol\delta}{2}\right)
\sqrt{\frac{\hbar}{2\rho V\omega_{\mathbf q\lambda}}},$$

and the resulting tight-binding vertex contains a form factor proportional to

$$g_t(\mathbf k,\mathbf q,\lambda)\propto
-2\sum_{\boldsymbol\delta}\delta t_{\boldsymbol\delta}(\mathbf q,\lambda)
\cos[(\mathbf k+\mathbf q/2)\!\cdot\boldsymbol\delta].$$

This is momentum dependent and, for an acoustic mode, scales with the strain (at small $q$, $\sin(\mathbf q\cdot\boldsymbol\delta/2)\sim q$). It is not a constant density vertex. Here $\rho V$ denotes the mass in the chosen phonon normalization volume; a strictly two-dimensional phonon normalization would instead use areal mass density times area. This example is a derivation **conditional on a bond-length hopping model**; the derivatives and parameters needed to evaluate it are not available, so it is not adopted as a prediction here.

A scalar deformation potential $D_{ab}\epsilon_{ab}\,n$ is another possible model. It is a density vertex, but its value and tensor structure must come from a deformation-potential calculation/measurement. It does not follow merely from saying “strain modulates hopping.” A strain-dependent Rashba term is also possible because the Rashba field depends on inversion-breaking electric fields and bond geometry; its size and spin structure require $\partial\mathbf g_R/\partial\epsilon_{ab}$. The present inputs do not choose among these mechanisms or establish that one dominates.

### 2.1 Momentum and Rashba-helicity matrix elements

For the electronic Rashba Hamiltonian written in the repository, use spectral projectors

$$P_s(\mathbf k)=\frac12[1+s\hat{\mathbf g}_{\mathbf k}\cdot\boldsymbol\sigma],\qquad s=\pm1.$$

For a **scalar spin-independent** vertex $g(\mathbf k,\mathbf q)\sigma_0$, the helicity overlap is

$$|\langle u_{s',\mathbf{k+q}}|u_{s,\mathbf k}\rangle|^2
=\operatorname{Tr}[P_{s'}(\mathbf{k+q})P_s(\mathbf k)]
=\frac{1+ss'\hat{\mathbf g}_{\mathbf{k+q}}\cdot\hat{\mathbf g}_{\mathbf k}}2.$$

At $q=0$ this is $\delta_{s's}$ away from Rashba degeneracies. At small nonzero $q$, inter-helicity amplitude is generally $O(q)$ and its squared weight is $O(q^2)$, subject to degeneracies and the direction of $\mathbf q$. A modulation of SOC itself has a vertex proportional to a Pauli matrix and need not obey this scalar-vertex selection rule.

The derivation asserts helicity diagonality for its $q\to0$ density vertex ([lines 111–114](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/DERIVATION_HHM_scattering.md#L111-L114)). By contrast, the standard production energy kernel has no helicity label. The separate `make_soc_explicit.py` study calls its A vertex isotropic/all-helicity and its B vertex helicity diagonal ([lines 4–21, 84–90](https://github.com/AmirhoseynpowAsghari/jubilant-enigma/blob/main/project/ultrasonic/make_soc_explicit.py#L4-L21)); even B only filters helicity and still does not enforce $\mathbf k'=\mathbf k+\mathbf q$. Its reported A/B sensitivity is evidence that helicity weights matter, not a finite-q calculation.

## 3. What a finite-q BdG response would require

For a specified one-body acoustic perturbation $V(\mathbf k,\mathbf q,\lambda)$, define the BdG vertex $\Gamma_{\mathbf q\lambda}(\mathbf k)$ in the Nambu basis, including its particle and hole blocks. If $|\Phi_{m\mathbf k}\rangle$ are BdG eigenvectors, the transition matrix elements are

$$M_{nm}(\mathbf k,\mathbf q,\lambda)=
\langle\Phi_{n,\mathbf{k+q}}|\Gamma_{\mathbf q\lambda}(\mathbf k)|\Phi_{m,\mathbf k}\rangle.$$

This expression retains the helicity spinors and superconducting coherence factors together; one must not multiply by a scalar coherence factor and an independent spin factor unless the chosen vertex allows that factorization. Crystal momentum conservation is then explicit: in an extended-zone convention, $\mathbf k'=\mathbf k+\mathbf q+\mathbf G$ with reciprocal-lattice vector $\mathbf G$; equivalently, fold $\mathbf k+\mathbf q$ back into the first Brillouin zone and retain the band/folding label. For pair creation the two outgoing quasiparticles carry total crystal momentum $\mathbf q$ modulo $\mathbf G$.

For clarity, in the Nambu convention $\Psi_{\mathbf k}=(c_{\mathbf k\uparrow},c_{\mathbf k\downarrow},c^\dagger_{-\mathbf k\uparrow},c^\dagger_{-\mathbf k\downarrow})^T$, a one-body electron vertex $V(\mathbf k,\mathbf q)$ is lifted to a BdG vertex $\Gamma_{\mathbf q}(\mathbf k)$ by including the corresponding particle and hole blocks. For a real, spin-independent scalar potential, the Nambu structure is $\Gamma=g\,\tau_z\otimes\sigma_0$ (with the finite-q momentum arguments on the two blocks). Thus $M_{nm}=\Phi^\dagger_{n,\mathbf{k+q}}\Gamma_{\mathbf q}(\mathbf k)\Phi_{m,\mathbf k}$ contains both the spin/helicity overlaps and the $u/v$ coherence amplitudes. For a more general strain-modulated hopping or Rashba term, one must use its actual matrix-valued $V$; it need not reduce to $g\tau_z$.

For positive-energy quasiparticle branches, the scattering contribution has the structure

$$P_{sc}(\mathbf q,\lambda)=\frac{2\pi}{\hbar}(\hbar\omega_{\mathbf q\lambda})
\sum_{\mathbf k,m,n}|M_{nm}|^2
[f(E_{m\mathbf k})-f(E_{n,\mathbf{k+q}})]
\delta(E_{n,\mathbf{k+q}}-E_{m\mathbf k}-\hbar\omega_{\mathbf q\lambda}),$$

and pair creation has terms of the form

$$P_{pair}(\mathbf q,\lambda)\propto\sum_{\mathbf k,m,n}|M^{pair}_{nm}|^2
[1-f(E_{m\mathbf k})-f(E_{n,\mathbf{q-k}})]
\delta(E_{m\mathbf k}+E_{n,\mathbf{q-k}}-\hbar\omega_{\mathbf q\lambda}).$$

The overall normalization depends on phonon-mode normalization and the microscopic vertex; these equations specify the transition structure, not a numerical rate. For a long-wavelength acoustic branch, the dispersion must be imposed, $\omega_{\mathbf q\lambda}=v_{s,\lambda}|\mathbf q|$ (or the appropriate anisotropic dispersion), rather than treating $\omega$ and $q$ independently. The calculated $W/W_N$ must use the same $q$, branch, polarization, normalization, and screening prescription in the superconducting and normal states.

For a scalar, intraband density vertex in the appropriate small-q limit, the BdG matrix element reduces to the familiar $u u'-v v'$ structure and hence the Case-I coherence factor quoted in the repository. That reduction checks the algebra **conditional on the vertex**; it does not remove the finite-q kinematics or supply the missing deformation-potential magnitude.

## 4. Conservation, screening, disorder, and collective response

At exactly $q=0$, a scalar density perturbation is proportional to total particle number $N$. For a number-conserving Hamiltonian, $[H,N]=0$, so

$$ (E_n-E_m)\langle n|N|m\rangle=\langle n|[H,N]|m\rangle=0.$$

It cannot cause a finite-frequency transition between different-energy exact states. In a superconducting mean-field BdG calculation, gauge invariance and charge conservation require the appropriate phase/collective and vertex corrections; a bare density bubble need not satisfy the Ward identity by itself. In a charged system, long-range Coulomb screening also modifies the scalar density response. Whether these effects are numerically large for a particular acoustic mode depends on the material, $q$, frequency, screening, and regime; those inputs are absent here.

Disorder is **not automatically implied** by using a constant matrix element. An impurity-assisted/dirty-limit calculation can relax electronic momentum selection after disorder averaging, but it requires a disorder correlator, scattering rate/mean free path, the regime (e.g. $q\ell$), and consistent self-energy and vertex corrections. The `bruteforce_ua.py` phrase “dirty-limit / constant-B approximation” is a label for the implemented unrestricted sum, not a derivation of those conditions. No $\tau$, $\ell$, impurity correlator, or impurity vertex correction appears in the production calculation.

## 5. Consequences for current results and SOC-dependent maxima

- The current $W/W_N$ values are **numerical results for the unrestricted, constant-matrix-element transition model**. Internal energy-space/lattice agreement validates numerical reduction under that same assumption; it does not validate the microscopic acoustic vertex.
- The original paper’s constant-$B$ angular/spin average is itself a phenomenological approximation. It may be used as such, but the later density-operator derivation does not establish it as equivalent to a clean finite-q acoustic vertex.
- The repository’s A/B helicity comparison shows that a change of helicity selection can alter the response; at $V_{SO}/t=0.8$, $T/T_c=0.9$, the stored isotropic and helicity-diagonal values are about 0.4857 and 0.4250. This is a useful sensitivity indication, but B is not a full momentum-conserving calculation.
- Therefore the existing SOC-dependent peak positions and shifts are not validated predictions for the correct acoustic coupling. They remain properties of the stated phenomenological model until recalculated with specified $V(\mathbf k,\mathbf q,\lambda)$, $\omega_{\mathbf q\lambda}$, and any required screening/disorder response.

No corrected finite-q code or numerical comparison is included: the repository does not specify the microscopic deformation potentials, strain derivatives of hopping/correlated hopping/Rashba SOC, phonon velocities and polarizations, mass density, screening, or disorder regime. Choosing values or a “dominant” vertex would introduce assumptions the task explicitly forbids.

## 6. Resolved and unresolved

**Resolved by source inspection**

1. The derivation defines $\rho_{\mathbf q}$ with $\mathbf{k+q}$ and mentions $\mathbf k'=\mathbf k+\mathbf q$ for pair creation.
2. The production absorption integral does not retain that momentum relation; its delta is energy-only after reduction, and its lattice cross-check is explicitly unrestricted.
3. The production coherence factor is the scalar Case-I BdG factor; it has no explicit Rashba helicity overlap. A separate helicity-only comparison exists but is not momentum resolved.
4. The original manuscript itself describes averaging its transition matrix element over angle and spin and treating it as constant before replacing sums by a DOS convolution.

**Unresolved / required before a physical finite-q prediction**

1. Which microscopic electron-ion coupling applies: local deformation potential, ordinary hopping modulation, correlated-hopping modulation, Rashba-SOC modulation, or a measured/derived combination.
2. The strain derivatives and phonon mode data that set the vertex and acoustic dispersion.
3. Whether screening, collective modes, disorder, and conserving vertex corrections are required and their parameter values.
4. The finite-q BdG response, its Ward/conservation limit, numerical convergence, and the normal-state normalization with the same vertex.
5. Whether the previously reported SOC-dependent maxima survive that response.

**Bottom line:** do not interpret the present production curves as a verified finite-q ultrasonic absorption calculation. Keep them identified as model-dependent unrestricted-transition results. A corrected numerical calculation should begin only after the missing microscopic and transport inputs are supplied or independently derived.
