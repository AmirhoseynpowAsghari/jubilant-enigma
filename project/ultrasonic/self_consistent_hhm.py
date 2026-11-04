"""Optional ultrasonic-state adapter for a fully self-consistent HHM gap.

``physics_ua.SCState`` remains the production-compatible ``legacy_fit`` implementation.
This adapter is opt-in: it solves D0(T), DS(T), and mu_SC(T) from the coupled HHM gap and
number equations, and computes mu_N(T) from the same fixed-filling number equation.  It
returns state dictionaries accepted directly by ``absorption.ratio_all``; no ultrasonic
response formula is changed here.
"""
from __future__ import annotations

from typing import Iterable

from hhm_finite_t import HHMFiniteTemperatureSolver, NoAdmissibleRoot


class SelfConsistentHHMState:
    """Opt-in SC state with the HHM linearized Tc and nonlinear finite-T self-consistency.

    Args are explicit so the legacy registry/default filling cannot be changed implicitly.
    The chemical-potential policy is applied as an acceptance filter only; the relaxed
    root is retained in solver diagnostics if strict mode has no admissible branch.
    """

    method = "self_consistent_HHM"

    def __init__(self, U: float, Delta_t: float, n_fill: float, V_SO: float,
                 Nk: int = 80, mu_policy: str = "relaxed", tc_scan_max: float = 8.0,
                 tc_scan_points: int = 401, initial_root_index: int = 0):
        self.U = float(U)
        self.Delta_t = float(Delta_t)
        self.n_fill = float(n_fill)
        self.V_SO = float(V_SO)
        self.Nk = int(Nk)
        self.mu_policy = mu_policy
        self.solver = HHMFiniteTemperatureSolver(
            self.U, self.Delta_t, self.n_fill, self.V_SO, Nk=self.Nk,
            mu_policy=mu_policy, tc_scan_max=tc_scan_max,
            tc_scan_points=tc_scan_points, initial_root_index=initial_root_index)
        self.Tc = self.solver.Tc
        self._cache: dict[float, dict] = {}

    def at_many(self, fractions: Iterable[float]) -> list[dict]:
        """Return states for ascending fractions; continuation checkpoints are cached too."""
        fracs = [float(f) for f in fractions]
        if not fracs or any(f < 0.0 or f >= 1.0 for f in fracs):
            raise ValueError("fractions must be a nonempty sequence in [0,1)")
        if any(b <= a for a, b in zip(fracs, fracs[1:])):
            raise ValueError("fractions must be strictly increasing")
        path = self.solver.solve_fractions(fracs)
        for r in path:
            rec = dict(frac=float(r["T"] / self.Tc), T=float(r["T"]),
                       D0=float(r["D0"]), DS=float(r["DS"]),
                       Deff=float(r["Deff"]), mu=float(r["mu"]),
                       mu_N=float(r["mu_N"]), phi=None,
                       Tc=float(self.Tc), method=self.method,
                       Tc_method="HHM_linearized",
                       mu_window=dict(r["mu_window"]),
                       residual_inf=float(r["residual_inf"]),
                       branch_selection=r["branch_selection"],
                       requested=bool(r.get("requested", False)),
                       mu_policy=self.mu_policy)
            self._cache[round(rec["frac"], 12)] = rec
        return [self._nearest_cached(f) for f in fracs]

    def _nearest_cached(self, frac: float) -> dict:
        key = round(float(frac), 12)
        if key in self._cache:
            return dict(self._cache[key])
        if not self._cache:
            raise KeyError("no temperatures have been solved yet; call at_many first")
        actual = min(self._cache, key=lambda k: abs(k - key))
        if abs(actual - key) > 1e-10:
            raise KeyError(f"requested fraction {frac} was not solved")
        return dict(self._cache[actual])

    def at(self, frac: float) -> dict:
        """Convenience method for one fraction (builds a low-to-high continuation path)."""
        frac = float(frac)
        if frac < 0.0 or frac >= 1.0:
            raise ValueError("frac must lie in [0,1)")
        if round(frac, 12) in self._cache:
            return dict(self._cache[round(frac, 12)])
        return self.at_many([0.0, frac] if frac > 0.0 else [0.0])[-1]

    def attenuation(self, dos, frac: float, omega: float, nell: int = 6000) -> dict:
        """Evaluate the existing signed-coherence scattering formula for this state.

        The unchanged response kernel is ``absorption.ratio_all``; its normal reference
        uses this adapter's same-filling ``mu_N(T)``.
        """
        from absorption import ratio_all
        return ratio_all(dos, self.at(frac), omega, variant_list=("a_signed",), nell=nell)


def make_sc_state(gap_method: str, *, legacy_params: dict | None = None,
                  V_SO: float | None = None, U: float | None = None,
                  Delta_t: float | None = None, n_fill: float | None = None,
                  Nk: int = 80, Nk_mu: int = 300, mu_policy: str = "relaxed",
                  initial_root_index: int = 0):
    """Factory keeping the two temperature models explicitly selectable.

    ``legacy_fit`` delegates unchanged to ``physics_ua.SCState`` and requires its usual
    explicit stage-1 params (or a resolvable V_SO).  ``self_consistent_HHM`` requires all
    four physical parameters and uses the linearized HHM Tc.
    """
    if gap_method == "legacy_fit":
        from physics_ua import SCState
        if legacy_params is not None:
            return SCState(params=legacy_params, Nk_mu=Nk_mu, n_tgt=n_fill)
        return SCState(VSO=V_SO, Nk_mu=Nk_mu, n_tgt=n_fill)
    if gap_method == "self_consistent_HHM":
        if None in (U, Delta_t, n_fill, V_SO):
            raise ValueError("self_consistent_HHM requires U, Delta_t, n_fill, and V_SO")
        return SelfConsistentHHMState(U, Delta_t, n_fill, V_SO, Nk=Nk,
                                      mu_policy=mu_policy,
                                      initial_root_index=initial_root_index)
    raise ValueError("gap_method must be 'legacy_fit' or 'self_consistent_HHM'")
