# Editable Physical Review-style manuscript draft

## Files

- `draft_prb.tex` — APS REVTeX 4.2 manuscript source, structured as a Physical Review B working draft.
- `draft_references.bib` — bibliography entries used in the source.
- `draft_prb.pdf` — generated reading proof; the REVTeX source remains authoritative and editable.

The author names and affiliation follow the original manuscript in `../vendor/arXiv_2607.02812_original.pdf`; please confirm the author list, corresponding-author details, acknowledgments, and funding before submission.

## Scientific scope of this draft

This is a corrected **model-paper draft**, not a claim that the project has completed a microscopic finite-wave-vector ultrasonic calculation. It deliberately calls the computed observable a constant-vertex quasiparticle-scattering ratio and states that the acoustic vertex is unspecified by the HHM electronic Hamiltonian. Its SOC trends are conditional on the assumed vertex and supplied gap inputs.

The manuscript corrects or qualifies these original claims/formulas:

- Uses the signed Case-I factor derived from the stated Bogoliubov matrix element, with energy conservation imposed on quasiparticle energies.
- Keeps the two energy-inversion roots and the Jacobian in the branch-resolved quasiparticle DOS.
- Distinguishes the helicity label from the energy-branch label.
- Identifies the actual pair-channel threshold in this model as `hbar*omega = 2 E_min`, not `2 Delta_eff` in general.
- Does not call `W/W_N` a superfluid weight or a material-specific attenuation coefficient.
- Does not equate strain-modulated hopping with a scalar density vertex without a microscopic derivation.
- Does not claim the constant-vertex, unrestricted momentum sum is a clean finite-q acoustic response.
- Identifies the stage-one nonlinear gap-threshold `Tc` proxy separately from the HHM linearized-determinant criterion, and treats the rigid fitted gap shape as a modeling choice.

The full earlier finite-q audit is in `../ultrasonic/FINITE_Q_ACOUSTIC_COUPLING_ASSESSMENT.md`.

## Figures included by reference

No image was copied or changed. The `.tex` file points directly to existing repository figures:

1. `ultrasonic/figures/strict_policy/WWN_vs_TTc_and_Delta_VSOfamilies.png` — strict-window temperature/gap families.
2. `ultrasonic/figures/strict_policy/WWN_vs_VSO_families_U.png` — strict-window SOC families.
3. `ultrasonic/figures/branch_decomposition.png` — energy-root branches and the model-specific spectral threshold.
4. `ultrasonic/figures/dual_space_k_vs_energy.png` — numerical cross-check within the unrestricted model.
5. `ultrasonic/figures/soc_explicit_vertex.png` — helicity-vertex sensitivity diagnostic, explicitly labeled as not finite-q.

Missing or rejected strict-policy samples are described as gaps, not zeros. The legacy figures in `paper_correction/figures/` were not substituted for these current strict-policy plots.

## Build

For a publication-layout build, compile the editable source from this directory with REVTeX 4.2 and a working LaTeX/BibTeX installation:

```sh
latexmk -pdf -interaction=nonstopmode draft_prb.tex
```

The checked-in PDF is a reading proof generated in this sandbox through Pandoc-to-Typst because no TeX engine is installed here. It was visually checked for the shared affiliation, figure captions, and bibliography. This conversion is not a substitute for compiling the REVTeX source: equation/layout details can differ, so use the `.tex` and `.bib` files as the authoritative manuscript and rebuild them with REVTeX before submission. Figure paths are relative and assume compilation from this directory.

## Open items for author review

1. Confirm the author list, correspondence, acknowledgments, and funding.
2. Decide whether the paper should remain a clearly labeled phenomenological model study or wait for a specified electron--phonon vertex and finite-q calculation.
3. Verify references and any claims against the cited source papers before submission; the current draft cites only the HHM and Levy works represented in the repository.
4. Review the selected strict-policy figures and captions. They are already-generated outputs and have not been regenerated or modified for this draft.
