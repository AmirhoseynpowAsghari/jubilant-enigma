"""
physics_ua.py -- gap shape phi(T), chemical potential from the number equation, and the
superconducting state at a given T/Tc (stage 2).
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
"""
from __future__ import annotations

import json
import numpy as np

from parameters_ua_u import (T_HOP, KB, N_TARGET, GAP_FIT_P, GAP_FIT_Q, Nk_MU)
from kgrid_ua_u import kgrid
import sc_loader_u as sc_loader


def phi_gap(T: float, Tc: float, p: float = GAP_FIT_P, q: float = GAP_FIT_Q) -> float:
    """The repository's fitted gap shape (master physics_ua.Delta_T with Delta0 = 1)."""
    if T >= Tc:
        return 0.0
    return (1.0 - (T / Tc) ** p) ** q


def electron_density(mu: float, D0: float, DS: float, T: float, VSO: float,
                     Nk: int = 300, t: float = T_HOP) -> float:
    """
    n(mu) = 1 - (1/2) < xi/E tanh(E / 2T) >_k summed over both helicities
    (identical to project/ultrasonic/physics_ua.compute_electron_density, except that the gap is
    the energy-dependent Dbar(eps) = D0 - DS eps/4t rather than one constant).
    """
    e0, R = kgrid(Nk)
    tot = 0.0
    for s in (+1.0, -1.0):
        eps = e0 - 2.0 * s * VSO * R
        xi = eps - mu
        Dbar = D0 - DS * eps / (4.0 * t)
        E = np.hypot(xi, Dbar)
        th = np.tanh(E / (2.0 * KB * T)) if T > 1e-12 else np.ones_like(E)
        tot += np.mean(xi / np.maximum(E, 1e-300) * th)
    return 1.0 - 0.5 * tot


def solve_mu(T: float, D0: float, DS: float, VSO: float, n_tgt: float = N_TARGET,
             Nk: int = 300, bracket=(-100.0, 100.0)) -> float:
    """mu(T) such that n(mu) = n_tgt.  D0 = DS = 0 gives the normal-state mu_N(T)."""
    from scipy.optimize import brentq
    f = lambda mu: electron_density(mu, D0, DS, T, VSO, Nk=Nk) - n_tgt
    return brentq(f, *bracket, xtol=1e-12)


class SCState:
    """
    D0(T), DS(T), Delta_eff(T), mu_SC(T), mu_N(T) at a given T (fraction of the repository's Tc).

    mu_SC : number equation with the energy-dependent gap  (self-consistent with the gap used below)
    mu_N  : number equation with D0 = DS = 0              (the normal state of the same model)
    """

    def __init__(self, VSO: float = None, Nk_mu: int = 300,
                 repo_root: str = None,
                 params: dict = None, n_tgt: float = None):
        if params is None:
            p = sc_loader.sc_params_at(VSO, repo_root)
            self.meta = p["meta"]
        else:                            # explicit stage-1 numbers (paper parameter sets)
            p = params
            self.meta = {"source": "explicit params"}
        self.VSO = p["VSO"]
        self.Tc, self.D0_0, self.DS_0, self.Deff_0 = p["Tc"], p["D0_0"], p["DS_0"], p["Deff_0"]
        self.n_tgt = N_TARGET if n_tgt is None else n_tgt
        self.Nk_mu = Nk_mu

    @classmethod
    def from_params(cls, VSO, Tc, D0_0, DS_0, Deff_0, n_tgt=N_TARGET, Nk_mu=300):
        """SC state built from explicit stage-1 numbers instead of the repository json."""
        return cls(params=dict(VSO=VSO, Tc=Tc, D0_0=D0_0, DS_0=DS_0, Deff_0=Deff_0),
                   n_tgt=n_tgt, Nk_mu=Nk_mu)

    def at(self, frac: float, mu_mode: str = "consistent"):
        T = frac * self.Tc
        phi = phi_gap(T, self.Tc)
        D0, DS = phi * self.D0_0, phi * self.DS_0
        Deff = phi * self.Deff_0
        mu_N = solve_mu(T, 0.0, 0.0, self.VSO, n_tgt=self.n_tgt, Nk=self.Nk_mu)
        if mu_mode == "normal":          # W_N convention: use mu_N everywhere
            mu = mu_N
        else:                            # repository convention: mu from the SC number equation
            mu = solve_mu(T, D0, DS, self.VSO, n_tgt=self.n_tgt, Nk=self.Nk_mu)
        return dict(frac=frac, T=T, phi=phi, D0=D0, DS=DS, Deff=Deff, mu=mu, mu_N=mu_N)


if __name__ == "__main__":
    for v in (0.0, 0.4082, 0.8):
        p = sc_params_at(v)
        print(f"V_SO={p['VSO']:.4f}  Tc={p['Tc']:.4f}  D0(0)={p['D0_0']:.4f}  "
              f"DS(0)={p['DS_0']:.4f}  Deff(0)={p['Deff_0']:.4f}  band edge={band_edge(v):.4f}")
    d = NormalDOS(0.4082, Nk=400)
    print("DOS check:", d.check())
