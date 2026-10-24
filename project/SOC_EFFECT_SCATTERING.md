# The spin–orbit (V_SO) effect in the scattering channel

> **SUPERSEDED IN PART (2026-10-24).** This document was written before the audit-driven correction
> campaign. Numbers, statements and figure references below are PRE-CORRECTION unless explicitly
> stated otherwise. The current status of every claim is in `CORRECTION_LOG.md` (§2 per-finding status,
> §5 before/after, §6.3 superseded statements). Do not quote figures from this file without checking
> them against `CORRECTION_TABLES.md` and `results/BEFORE_AFTER.json`.


Companion to `ultrasonic/make_soc_figures.py` (data: `ultrasonic/results/soc_scan.json`);
figures `ultrasonic/figures/soc_effect_scattering.png` and `soc_mechanism.png`.
Channel: η=+1 signed coherence factor, **scattering only**, ħω = 0.7 t, n = 1.885,
U = 95 t (and 125 t in panels b/d), Δ_t = 4.5 t.

## 1. The effect is large, monotone in V_SO at fixed T/T_c, and changes sign in T

| V_SO/t | T_c/t | Δ⁰₀/T_c | R(T/T_c=0.5) | R(0.9) | R(0.999) | R(T=0.1t) | R(T=0.2t) |
|---|---|---|---|---|---|---|---|
| 0    | 0.2861 | 0.906 | 0.02628 | 0.26098 | 0.40004 | 0.00414 | 0.10348 |
| 0.2  | 0.3082 | 0.863 | 0.02820 | 0.28023 | 0.42920 | 0.00279 | 0.08388 |
| 0.4082 | 0.3809 | 0.750 | 0.03342 | 0.33187 | 0.50542 | 0.00071 | 0.04170 |
| 0.6  | 0.4921 | 0.626 | 0.03912 | 0.39476 | 0.58964 | 0.00008 | 0.01362 |
| 0.8  | 0.6494 | 0.500 | 0.04601 | 0.47571 | 0.76111 | 0.00000 | 0.00266 |

* **Fixed T/T_c:** W/W_N grows monotonically with V_SO by ×1.75 (T/T_c=0.5), ×1.82 (0.9),
  ×1.90 (0.999) from V_SO=0 to 0.8 — the curves fan out visibly at every temperature
  (figure 1a,b,c).
* **Sign change:** the relative enhancement (R(V_SO)−R(0))/R(0) is *negative* (down to
  −70 %) for T/T_c ≲ 0.22 and positive (up to +90 %) above (figure 1c). Deep in the SC
  state the SOC-grown gap minimum E_min (0.47→1.22 t at T/T_c=0.5) cuts the thermal
  scattering phase space; near T_c the fallen coupling ratio Δ⁰₀/T_c (0.906→0.500)
  weakens the f(Δ/T) suppression and wins. The crossing at T/T_c ≈ 0.22 is itself a
  spin–orbit signature.
* **Histogram view (figure 36, `make_soc_enhancement_hist.py`):** over all (U, T/T_c)
  points per V_SO level (root-switch sets excluded), the per-point enhancement
  R(V_SO)/R(0) has median 1.06 / 1.21 / 1.40 / 1.67 and p90 1.26 / 2.17 / 3.72 / 5.73 for
  V_SO = 0.2 / 0.4 / 0.6 / 0.8; the fraction of points enhanced by >5 % grows 0.58 → 0.90.
  Per U-bin: median enhancement ≈1 for U ≤ 30 at T/T_c ≥ 0.7 (saturated plateau) but
  ×1.5–×14 in the decay tail (U ≥ 75); at T/T_c = 0.3 the suppression is active already at
  U ≈ 30, so the SOC leverage appears there too.  Each +0.2 step of V_SO raises the ratio
  at ~96–100 % of the points (median ×1.05–×1.25 per step).
* **Fixed physical T:** the effect is even stronger and opposite in sign — at T=0.2 t the
  ratio falls ×39 (0.103→0.0027) because T_c doubles (0.286→0.649 t) and the system sits
  deeper inside the dome (figure 1d). Any comparison statement must name its ensemble.

## 2. Four internal routes (figure 2)

1. **Band splitting** ε_ks = ε_k − 2 s V_SO|sin k|: the van-Hove peak of N(ε) splits into
   two peaks moving apart with a dip at ε=0 (fig. 2a), redistributing DOS weight around μ.
2. **Self-consistent scales:** T_c rises ×2.27 while Δ⁰₀/T_c falls ×0.55 (fig. 2b); at fixed
   T/T_c this directly weakens the exponential suppression, at fixed T it strengthens it.
3. **Branch structure:** ε_σ(E) = ε* + σ√(E²−E_min²)/κ moves (ε*, κ, E_min), redistributing
   the (σ,σ′) weights of the scattering sum — the (−,−) share falls 70→66 % while (−,+)
   grows 30→34 % between V_SO=0 and 0.8 (fig. 2c); E_min and μ both rise monotonically
   (fig. 2d).
4. **Gap anisotropy sampling:** D̄(ε) = Δ⁰ − (Δ_t/4t)ε is evaluated on split helicity
   sheets, so the coherence factor C_I = ½[1+(ξξ′−D̄D̄′)/(EE′)] sees different (ξ, D̄) pairs.

All four act *inside* the scattering channel; none is an artifact of the normal-state
reference (the same N(ε) and μ_N conventions are used for W and W_N at each V_SO).

## 3. Relation to earlier figures

The regenerated paper figures `W_WN_vs_VSO*.png` (U = 95, 125, −75 at fixed T/T_c) and the
Group-D panels show the same monotone enhancement; this note adds the fixed-T ensemble,
the relative-enhancement sign change, and the mechanism decomposition. The parameter
study (`PARAM_SCAN_PEAKS.md`) confirms the V_SO dependence of the peak positions of the
withdrawn pair-dependent peak study; the scattering-only SOC figures (01, 02, 32, 34) are unaffected.


## 4. What is EXPLICIT and what is INDIRECT in the SOC effect (vertex caveat)

The kernel C_I(xi,xi',Dbar,Dbar'; E,E') contains no V_SO: SOC reaches W/W_N only through
band splitting (DOS, Fermi volumes, phase space), the self-consistent scales (T_c, D0, DS,
mu), and E_min/thresholds. Consequently, in representations normalised by those scales
(T/T_c, Delta_eff(T)/Delta_eff(0)) the SOC content partly collapses out -- this is a
property of the model, not of the plotting.

Explicit kernel-level signatures that DO exist (figure `soc_explicit_vertex.png`,
SOCfigures/35, data `results/soc_explicit.json`, helicity-resolved k-space lattice sums):

* **helicity-channel redistribution:** the inter-helicity fraction of S_sc falls from 0.500
  (degenerate bookkeeping at V_SO=0) to 0.438 (V_SO=0.4082) and 0.29 (V_SO=0.8) -- SOC moves
  scattering weight between helicity channels inside the kernel;
* **vertex systematic:** with a strict q->0 vertex (helicity-diagonal, Rashba spinor overlap
  -> delta_ss'), R_B vs our isotropic-vertex R_A: identical at V_SO=0 (sanity), -12 % to
  -20 % at V_SO=0.8; yet the SOC trend SURVIVES: R_B(0.8)/R_B(0) = 1.70 at T/T_c=0.5
  (vs 1.75 for R_A). The effect is therefore physical and not manufactured by the vertex
  approximation, but ~a quarter of its magnitude at strong SOC rides on inter-helicity
  channels that a strict q->0 vertex suppresses.
* **Caveat / future work:** a finite-q vertex with the full spinor overlap (weight
  ~ q^2 v_F^2 / Delta_SO^2 for inter-helicity processes) would interpolate between R_B and
  R_A; ql-dependence is outside the present model (see LEVY_REGIME_ANALYSIS.md). The
  isotropic-vertex practice is the same as the paper's, Walker-Smith-Samokhin's and the
  homogeneous limit of Sardwal et al. (arXiv:2608.19079).
