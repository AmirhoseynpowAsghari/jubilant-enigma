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
NBINS    = 6000            # DOS bins
SMOOTH   = 0.02            # Gaussian smoothing of the DOS (units of t)
Nk_MU    = 300             # k-grid of the number-equation solve for mu(T)
NELL     = 6000            # midpoint points per (sigma, sigma') branch integral
