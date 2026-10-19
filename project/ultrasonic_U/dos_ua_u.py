"""
dos_ua.py -- normal-state DOS of the two helicity bands (normalised, exact band support).
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
This file keeps the repository's module name and role; the physics inside is the
validated one (see project/lecture_note_WWN_corrected.md and README.md).
"""
from __future__ import annotations

import json
import numpy as np

from parameters_ua_u import T_HOP
from kgrid_ua_u import kgrid, band_edge

# ---------------------------------------------------------------------------------------------
# normal-state DOS  N(eps) = sum_s < delta(eps - eps_s(k)) >_k ,  normalised to  integral = 2
# ---------------------------------------------------------------------------------------------
class NormalDOS:
    """
    Histogram + Gaussian-smoothed normal DOS of the two Rashba helicity bands.

    Differences from the repository's project/ultrasonic/dos_ua.py (both are defects there):
      * normalisation: the repository's prefactor 1/(pi^1.5 sigma t) with norm = dk^2/(2pi)^2
        gives integral N deps = 0.643 instead of 2.  Here it is renormalised to exactly 2.
      * band support: the repository interpolates its table on [-6, 6] and returns a nonzero
        (spurious) Gaussian tail above the band edge.  Here N(eps) = 0 exactly for |eps| > W_band.
    Both matter for W_N; the ratio W/W_N is what we report, and it is only meaningful if the two
    integrals use the *same* N.
    """

    def __init__(self, VSO: float, Nk: int = 1400, nbins: int = 6000, smooth: float = 0.02,
                 t: float = T_HOP):
        self.VSO, self.t, self.smooth = VSO, t, smooth
        self.Wband = band_edge(VSO, t)
        edges = np.linspace(-self.Wband, self.Wband, nbins + 1)
        de = edges[1] - edges[0]
        self.de = de
        e0, R = kgrid(Nk)
        cnt = np.zeros(nbins)
        wsum = np.zeros(nbins)          # sum of eps inside each bin (for the exact mean eps)
        for s in (+1.0, -1.0):
            e = (e0 - 2.0 * s * VSO * R).ravel()
            idx = np.searchsorted(edges, e) - 1
            inside = (idx >= 0) & (idx < nbins)
            idx = idx[inside]
            cnt += np.bincount(idx, minlength=nbins)
            wsum += np.bincount(idx, weights=e[inside], minlength=nbins)
        from scipy.ndimage import gaussian_filter1d
        raw = cnt / (Nk ** 2 * de)                       # quarter-zone midpoint grid -> full zone
        dens = gaussian_filter1d(raw, smooth / de, mode="constant")
        dens *= 2.0 / (dens.sum() * de)                  # exact normalisation: integral N deps = 2
        self.grid = 0.5 * (edges[1:] + edges[:-1])
        self.dens = dens
        self.eps_mean = np.where(cnt > 0, wsum / np.maximum(cnt, 1), self.grid)
        # quadrature grid for the "state b" sum (midpoint rule, weight = N deps)
        self.eps_i = self.grid.copy()
        self.w_i = self.dens * de
        # for the brute-force k-space cross-check: the exact bin index of every k point
        self._kdata = (e0, R, edges)

    def N(self, eps):
        eps = np.asarray(eps, float)
        out = np.interp(eps, self.grid, self.dens, left=0.0, right=0.0)
        return np.where(np.abs(eps) <= self.Wband, out, 0.0)     # exact band support

    def bin_of(self, eps):
        """bin index of a band energy (clipped); used by the brute-force check."""
        edges = self._kdata[2]
        return np.clip(np.searchsorted(edges, eps) - 1, 0, len(self.grid) - 1)

    def check(self):
        tot = self.w_i.sum()
        return dict(integral=tot, Wband=self.Wband, N_at_mu_region=float(self.N(3.6)),
                    nbins=len(self.grid), de=self.de)


