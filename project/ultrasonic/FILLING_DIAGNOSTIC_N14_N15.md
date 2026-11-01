# Intermediate-filling diagnostic: n = 1.4 and 1.5

## Purpose

This additive diagnostic asks whether the T=0 self-consistency solver finds nontrivial roots with a chemical potential outside the sampled bare Rashba-band interval at two intermediate-high fillings. It holds the other parameters fixed at `Delta_t=4.5`, samples `U/t = 15, 20, 25, 30, 35, 40` and `V_SO/t = 0.60, 0.65, ..., 0.90`, and uses the midpoint `Nk=80` grid. These points lie in the low-U/high-SOC region where the earlier relaxed-policy scan reported out-of-band roots.

Run from `project/ultrasonic`:

```bash
PYTHONPATH=../sc_params python3 scan_filling_roots.py
```

The script writes new files only:

- `results/filling_roots_n140_n150.json`
- `results/filling_roots_n140_n150.csv`

## Current result

On the listed 84-point grid (42 points per filling), the largest-Deff root found by the existing multistart solver was inside the sampled bare-band interval at all 84 points. The largest reported residual was `1.33e-8` for `n=1.4` and `1.74e-8` for `n=1.5`. The nonzero effective gaps ranged from `2.50` to `5.51` at `n=1.4` and `2.45` to `5.34` at `n=1.5` (in units of `t`). The JSON also contains 18 selected rechecks at `Nk=160` (nine per filling). At those shared points, the largest coarse-to-fine `Deff_0` change was `4.08e-5` for `n=1.4` and `3.43e-5` for `n=1.5`; the largest `mu_0` changes were `1.17e-8` and `3.88e-9`, respectively. All 18 high-resolution selected roots were also in-band.

This is evidence that the particular selected roots in this intermediate-filling pilot do not reproduce the out-of-band flag seen in the earlier `n=1.885` targeted scan. It does **not** establish that every branch is in-band or that the selected branch has the lowest free energy.

It does not overwrite the existing n=1.885 scans.

## What is fixed in this change

`sc_params/gap_equations.py::find_mu_normal` now uses the project filling convention, in which `epsilon_ks` contains both helicity bands and `0 <= n <= 2`. At finite beta it solves `n = 1 - mean(tanh(beta*(epsilon-mu)/2))`; at beta=infinity it uses the `n/2` quantile of the finite-grid band energies as a seed. The zero-temperature number function is a step function on a finite grid and may not have an exact root for an arbitrary target filling, so the quantile is used only as an initial guess.

A regression check covers both finite-beta and zero-temperature seeds at `n=1.4` and `n=1.5`.

## Important interpretation limits

- Each row records the largest-effective-gap T=0 root found by the existing multi-start solver. It is not a complete enumeration of all self-consistent roots.
- `mu_inside_band` is diagnostic only. An out-of-band value is neither automatically rejected nor declared physical.
- This scan does **not** calculate the HHM Helmholtz/free energy, compare competing roots, or identify the thermodynamically stable branch. The repository does not currently implement the exact mean-field thermodynamic potential needed for that ranking. Consequently, `stability` is explicitly marked `not assessed` and the scan cannot answer which point is the true equilibrium state.
- The scan does not compute `W/W_N`; doing so before choosing the thermodynamic branch would only propagate an unranked candidate state.

The next step for stability is to derive and implement the grand potential from the paper's decoupled mean-field Hamiltonian (including all constant/double-counting terms), enumerate all roots at each fixed `(U, V_SO, Delta_t, n, T)`, convert each to the canonical free energy `F = Omega + mu*N`, and compare at the same filling and temperature. Do not infer stability from the residual norm or band-window flag alone.
