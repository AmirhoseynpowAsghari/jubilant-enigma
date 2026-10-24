"""
parameters_ua.py -- physical and numerical constants (stage 2).
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
This file keeps the repository's module name and role; the physics inside is the
validated one (see project/lecture_note_WWN_corrected.md and README.md).
"""
from __future__ import annotations

import json
import numpy as np

T_HOP  = 1.0
HBAR   = 1.0
KB     = 1.0
OMEGA  = 0.7 * T_HOP       # hbar*omega of the ultrasound
B_COUP = 1.0               # coupling prefactor B (cancels in W/W_N)
N_TARGET = 1.8850          # filling of the master parameter set
ETA    = +1.0              # case I (ultrasound): coherence factor u u' - eta v v'

GAP_FIT_P = 3.434          # MASTER rigid gap shape phi(T) = (1-(T/Tc)^p)^q for ALL V_SO.
# Deliberate model choice (REPORT App. A4 / erratum E8): sc_params/results/gap_vs_T_fit.json
# gives per-V_SO fits (p,q)=(3.674,0.556) at V_SO=0; using per-V_SO refits changes W/W_N
# by 1.3-3.4 % at T/Tc=0.5/0.9 (quantified in REFEREE_RESPONSE.md).
GAP_FIT_Q = 0.553

Nk_DOS   = 1400            # k-grid of the normal DOS histogram
NBINS    = 48000           # DOS histogram bins: sigma/de >= 14; bin-converged to 0.04 % (F-03)
SMOOTH   = 0.0025          # DOS Gaussian smoothing (units of t), adopted by the pre-declared rule of calibrate_smooth.py (F-03)
                           # (calibrate_smooth.py, results/calibrate_smooth*.json): the old
                           # 0.02 carried a 2.2 % bias relative to 0.005 and the dependence
                           # on sigma has NOT converged (see CORRECTION_LOG.md).
Nk_MU    = 300             # k-grid of the number-equation solve for mu(T)
NELL     = 6000            # midpoint points per (sigma, sigma') branch integral
NEPS_NORMAL = 4          # Gauss-Legendre nodes per kink-aligned cell of S_normal (F-04)
SN_RTOL     = 1e-8         # accepted agreement of the order-n and order-(n+2) S_normal rules (F-04); achieved values are typically <1e-12, the 1e-10 setting was too tight for the tiny-omega Levy integrand (rel 1.9e-10)
BF_CHUNK    = 512          # row-chunk of the brute-force dense pair sums (F-11 memory bound)
