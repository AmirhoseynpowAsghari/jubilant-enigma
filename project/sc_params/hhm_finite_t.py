"""Optional HHM-compatible, fully self-consistent finite-temperature solver.

The existing production solver and its 1%-of-zero-temperature-gap Tc proxy are left
unchanged.  This module provides a separate path that (i) locates Tc from HHM's
linearized gap determinant, (ii) solves the nonlinear gap and number equations at each
T below that Tc, and (iii) follows a selected nontrivial root continuously from T=0.

Root status means residual-verified only.  It is not a free-energy/stability judgment.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from scipy.optimize import brentq, fsolve

from gap_equations import effective_gap, find_mu_normal, make_residual
from kgrid import build_kgrid
from solvers import BAND_EDGE_TOL, bare_band_window

RESIDUAL_TOL = 1e-7
NONTRIVIAL_GAP_TOL = 1e-10


class NoAdmissibleRoot(RuntimeError):
    """No residual-verified nontrivial root met the selected mu-window policy."""


def _max_residual(residual, x, beta, n_fill, eps, term) -> float:
    try:
        r = np.asarray(residual(list(x), beta, n_fill, eps, term), dtype=float)
    except Exception:
        return float("inf")
    return float(np.max(np.abs(r))) if np.all(np.isfinite(r)) else float("inf")


def _fold_and_verify(residual, raw, beta, n_fill, eps, term, tol=RESIDUAL_TOL):
    """Use HHM's exact global gap-sign symmetry, then independently check all equations."""
    x = np.asarray(raw, dtype=float).copy()
    if x.shape != (3,) or not np.all(np.isfinite(x)):
        return None
    if x[0] < 0.0:
        x[:2] *= -1.0
    rn = _max_residual(residual, x, beta, n_fill, eps, term)
    if rn > tol:
        return None
    return x, rn


def _deduplicate(roots, candidate, atol=2e-6):
    x, _ = candidate
    for old in roots:
        y = old["x"]
        scale = np.maximum(1.0, np.maximum(np.abs(x), np.abs(y)))
        if np.max(np.abs(x - y) / scale) <= atol:
            old["starts"] += 1
            old["residual_inf"] = min(old["residual_inf"], candidate[1])
            return
    roots.append(dict(x=x.copy(), residual_inf=float(candidate[1]), starts=1))


def linearized_tc(U: float, Delta_t: float, n_fill: float, V_SO: float,
                  Nk: int = 80, T_min: float = 0.01, T_max: float = 8.0,
                  n_scan: int = 401) -> dict:
    """Find HHM's linearized instability temperature(s) from det M(T)=0.

    The HHM relation for the third gap component is already incorporated in
    ``build_kgrid``.  The remaining two linearized residuals in (D0, DS) form M.
    The normal-state number equation is solved independently at each temperature.
    All sign-changing determinant roots are returned; the highest positive root is
    the normal-state instability reported as Tc.
    """
    if Nk < 4 or not (0.0 < T_min < T_max) or n_scan < 20:
        raise ValueError("require Nk>=4, 0<T_min<T_max, and n_scan>=20")
    eps, kernel = build_kgrid(V_SO, U, Delta_t, Nk)
    nstates = float(Nk * Nk)

    def matrix_and_mu(temp):
        beta = 1.0 / temp

        def number_eq(mu):
            z = np.clip(beta * (eps - mu) / 2.0, -350.0, 350.0)
            return 1.0 - float(np.tanh(z).mean()) - n_fill

        mu = brentq(number_eq, float(eps.min() - max(20.0, 20.0 * temp)),
                    float(eps.max() + max(20.0, 20.0 * temp)), xtol=1e-12)
        ab_xi = np.abs(eps - mu)
        g = np.tanh(ab_xi / (2.0 * temp)) / (2.0 * np.maximum(ab_xi, 1e-300))
        mat = np.array([
            [1.0 + np.sum(kernel * g) / (2.0 * nstates),
             -np.sum(kernel * eps * g) / (8.0 * nstates)],
            [4.0 * Delta_t * np.sum(g) / nstates,
             1.0 - Delta_t * np.sum(eps * g) / nstates],
        ], dtype=float)
        return mat, float(mu)

    def determinant(temp):
        return float(np.linalg.det(matrix_and_mu(float(temp))[0]))

    ts = np.geomspace(T_min, T_max, int(n_scan))
    ds = np.asarray([determinant(float(t)) for t in ts])
    roots = []
    for a, b, fa, fb in zip(ts[:-1], ts[1:], ds[:-1], ds[1:]):
        if fa == 0.0:
            root = float(a)
        elif fa * fb < 0.0:
            root = float(brentq(determinant, float(a), float(b), xtol=1e-12))
        else:
            continue
        if not roots or abs(root - roots[-1]) > 1e-8 * max(1.0, root):
            roots.append(root)
    if ds[-1] == 0.0:
        roots.append(float(ts[-1]))
    root_records = []
    for temp in roots:
        mat, mu = matrix_and_mu(temp)
        root_records.append(dict(Tc=float(temp), mu_normal=float(mu),
                                 determinant=float(np.linalg.det(mat)),
                                 matrix=mat.tolist()))
    if not root_records:
        return dict(Tc=None, roots=[], scan=dict(T_min=T_min, T_max=T_max,
                                                 n_scan=int(n_scan), determinant_min=float(ds.min()),
                                                 determinant_max=float(ds.max())))
    selected = max(root_records, key=lambda r: r["Tc"])
    return dict(Tc=selected["Tc"], mu_normal=selected["mu_normal"],
                determinant=selected["determinant"], matrix=selected["matrix"],
                roots=root_records,
                scan=dict(T_min=T_min, T_max=T_max, n_scan=int(n_scan),
                          determinant_min=float(ds.min()), determinant_max=float(ds.max())))


@dataclass
class HHMFiniteTemperatureSolver:
    """Finite-T nonlinear continuation at fixed (U, Delta_t, n, V_SO).

    ``mu_policy`` is explicit and independent of the module-level production switch:
    ``relaxed`` retains every residual-verified root; ``strict`` selects only roots whose
    mu is inside the analytic bare-band interval.  T=0 branch selection follows the
    existing solver's maximum-Deff convention among the candidates found.  Thereafter,
    the root nearest the continuation predictor is followed.  The full candidate list is
    retained in each returned record.
    """
    U: float
    Delta_t: float
    n_fill: float
    V_SO: float
    Nk: int = 80
    mu_policy: str = "relaxed"
    residual_tol: float = RESIDUAL_TOL
    tc_method: str = "HHM_linearized"
    tc_scan_max: float = 8.0
    tc_scan_points: int = 401
    initial_root_index: int = 0
    seed_d0: tuple = (0.1, 0.5, 1.0, 2.0, 5.0)
    seed_ds: tuple = (-10.0, -5.0, -2.0, -0.5, 0.5, 2.0, 5.0, 10.0)
    seed_dmu: tuple = (-4.0, -2.0, 0.0, 2.0, 4.0)

    def __post_init__(self):
        if self.Nk < 4:
            raise ValueError("Nk must be >=4")
        if self.mu_policy not in ("relaxed", "strict"):
            raise ValueError("mu_policy must be 'relaxed' or 'strict'")
        if self.tc_method != "HHM_linearized":
            raise ValueError("this optional solver defines Tc by HHM_linearized")
        self.eps, self.term = build_kgrid(self.V_SO, self.U, self.Delta_t, self.Nk)
        self.residual = make_residual(self.Delta_t)
        self.band = bare_band_window(self.eps, self.V_SO)
        self.tc_info = linearized_tc(
            self.U, self.Delta_t, self.n_fill, self.V_SO, Nk=self.Nk,
            T_max=self.tc_scan_max, n_scan=self.tc_scan_points)
        self.Tc = self.tc_info["Tc"]
        self.t0_roots = self._find_t0_roots()
        valid = [r for r in self.t0_roots if self._policy_accepts(r["x"][2])
                 and r["Deff"] > NONTRIVIAL_GAP_TOL]
        if not valid:
            self.t0 = None
        else:
            valid.sort(key=lambda r: r["Deff"], reverse=True)
            if not 0 <= int(self.initial_root_index) < len(valid):
                raise IndexError(f"initial_root_index={self.initial_root_index} is outside "
                                 f"the {len(valid)} accepted nontrivial T=0 roots")
            self.t0 = valid[int(self.initial_root_index)]

    def _mu_window(self, mu):
        lo, hi = self.band
        excess = float(max(float(mu) - hi, lo - float(mu), 0.0))
        return dict(mu=float(mu), band=[float(lo), float(hi)], inside=excess <= BAND_EDGE_TOL,
                    excess=excess, edge_tolerance=BAND_EDGE_TOL)

    def _policy_accepts(self, mu):
        return self.mu_policy == "relaxed" or self._mu_window(mu)["inside"]

    def _record(self, x, residual_inf, T, starts=1):
        D0, DS, mu = map(float, x)
        return dict(T=float(T), D0=D0, DS=DS, mu=mu,
                    Deff=effective_gap(D0, DS, self.eps),
                    mu_window=self._mu_window(mu), residual_inf=float(residual_inf),
                    accepted=bool(self._policy_accepts(mu)), starts=int(starts))

    def _find_t0_roots(self):
        mu0 = find_mu_normal(self.n_fill, np.inf, self.eps)
        roots = []
        for d0 in self.seed_d0:
            for ds in self.seed_ds:
                for dmu in self.seed_dmu:
                    guess = [float(d0), float(ds), float(mu0 + dmu)]
                    try:
                        raw, _, ier, _ = fsolve(
                            self.residual, guess,
                            args=(np.inf, self.n_fill, self.eps, self.term),
                            full_output=True, xtol=1e-10, maxfev=4000)
                    except Exception:
                        continue
                    if ier != 1:
                        continue
                    checked = _fold_and_verify(self.residual, raw, np.inf, self.n_fill,
                                               self.eps, self.term, self.residual_tol)
                    if checked is None:
                        continue
                    # At T=0 the exactly normal solution has a finite-grid plateau in mu;
                    # discard its near-zero gap copies rather than misreport them as branches.
                    if effective_gap(checked[0][0], checked[0][1], self.eps) <= NONTRIVIAL_GAP_TOL:
                        continue
                    _deduplicate(roots, checked)
        out = []
        for r in roots:
            x = r["x"]
            q = self._record(x, r["residual_inf"], 0.0, r["starts"])
            q["nontrivial"] = True
            out.append(q | dict(x=x))
        out.sort(key=lambda r: r["Deff"], reverse=True)
        return out

    def _roots_at(self, temp: float, seeds: Iterable[Iterable[float]]):
        roots = []
        beta = 1.0 / float(temp)
        for guess in seeds:
            try:
                raw, _, ier, _ = fsolve(
                    self.residual, list(guess),
                    args=(beta, self.n_fill, self.eps, self.term),
                    full_output=True, xtol=1e-10, maxfev=4000)
            except Exception:
                continue
            if ier != 1:
                continue
            checked = _fold_and_verify(self.residual, raw, beta, self.n_fill,
                                       self.eps, self.term, self.residual_tol)
            if checked is None:
                continue
            _deduplicate(roots, checked)
        out = []
        for r in roots:
            x = r["x"]
            q = self._record(x, r["residual_inf"], temp, r["starts"])
            q["nontrivial"] = q["Deff"] > NONTRIVIAL_GAP_TOL
            out.append(q | dict(x=x))
        return out

    def normal_mu(self, temp: float) -> float:
        if temp <= 0.0:
            return find_mu_normal(self.n_fill, np.inf, self.eps)
        return find_mu_normal(self.n_fill, 1.0 / float(temp), self.eps)

    def solve_path(self, temperatures: Iterable[float], include_diagnostics: bool = True) -> list[dict]:
        """Solve an ascending temperature path, starting from the selected T=0 branch.

        Intervals are internally subdivided to at most 0.05 Tc.  The continuation uses
        the previous solution and a linear predictor from the previous two points, plus
        scaled T=0 seeds and the linearized null direction near Tc.  The root nearest the
        predictor is followed; it is not reselected by maximum gap at every temperature.
        Returned records include intermediate continuation checkpoints and flag whether
        each was explicitly requested.
        """
        ts = [float(t) for t in temperatures]
        if not ts or any(not np.isfinite(t) or t < 0.0 for t in ts):
            raise ValueError("temperatures must be a nonempty finite sequence with T>=0")
        if any(b <= a for a, b in zip(ts, ts[1:])):
            raise ValueError("temperatures must be strictly increasing")
        if self.Tc is None:
            raise RuntimeError("HHM linearized determinant has no sign-changing Tc root")
        if any(t >= self.Tc for t in ts):
            raise ValueError("solve_path accepts only T<Tc; Tc itself is the normal-state limit")
        if self.t0 is None:
            raise NoAdmissibleRoot(
                f"no nontrivial T=0 root accepted under mu_policy={self.mu_policy}; "
                f"unrestricted candidates={self.t0_roots}")

        t0x = np.asarray([self.t0["D0"], self.t0["DS"], self.t0["mu"]], dtype=float)
        t0rec = {k: v for k, v in self.t0.items() if k != "x"}
        t0rec.update(T=0.0, mu_N=self.normal_mu(0.0), Tc=self.Tc,
                     Tc_method=self.tc_method,
                     branch_selection="max_Deff_at_T0_then_continuation",
                     requested=(0.0 in ts),
                     candidates=[{k: v for k, v in q.items() if k != "x"}
                                 for q in self.t0_roots] if include_diagnostics else [])
        path = [t0rec]
        tmax_step = 0.05 * self.Tc
        for target in ts:
            if target == 0.0:
                continue
            start_t = float(path[-1]["T"])
            if target <= start_t:
                continue
            nsub = max(1, int(np.ceil((target - start_t) / tmax_step)))
            subtemps = np.linspace(start_t, target, nsub + 1)[1:]
            for j, temp in enumerate(subtemps):
                prev = path[-1]
                prevx = np.asarray([prev["D0"], prev["DS"], prev["mu"]], dtype=float)
                predictor = prevx.copy()
                if len(path) >= 2:
                    prev2 = path[-2]
                    prev2x = np.asarray([prev2["D0"], prev2["DS"], prev2["mu"]], dtype=float)
                    dt_prev = float(prev["T"] - prev2["T"])
                    if dt_prev > 0.0:
                        ratio = (float(temp) - float(prev["T"])) / dt_prev
                        predictor = prevx + ratio * (prevx - prev2x)
                else:
                    amp = np.sqrt(max(1.0 - float(temp) / self.Tc, 1e-4))
                    predictor[:2] *= amp
                seeds = [predictor, prevx]
                for frac in (0.9, 0.7, 0.5, 0.3, 0.15, 0.05):
                    seeds.append([t0x[0] * frac, t0x[1] * frac, float(predictor[2])])
                # A small-gap seed in the linearized null direction helps close to Tc.
                mat = np.asarray(self.tc_info["matrix"], dtype=float)
                _, _, vh = np.linalg.svd(mat)
                null = vh[-1]
                seeds.append([float(null[0]) * 0.02, float(null[1]) * 0.02,
                              float(predictor[2])])
                candidates = self._roots_at(float(temp), seeds)
                nontrivial = [q for q in candidates if q["nontrivial"] and q["accepted"]]
                if not nontrivial:
                    raise NoAdmissibleRoot(
                        f"no nontrivial residual-verified root at T={temp:.9g} under "
                        f"mu_policy={self.mu_policy}; candidates={candidates}")
                scale = np.asarray([max(abs(prevx[0]), 0.1), max(abs(prevx[1]), 0.1), 1.0])
                chosen = min(nontrivial, key=lambda q: float(np.linalg.norm(
                    (q["x"] - predictor) / scale)))
                is_requested = any(abs(float(temp) - q) <= 1e-12 * max(1.0, self.Tc) for q in ts)
                rec = {k: v for k, v in chosen.items() if k != "x"}
                rec.update(mu_N=self.normal_mu(float(temp)), Tc=self.Tc,
                           Tc_method=self.tc_method,
                           branch_selection="nearest_continuation_predictor",
                           requested=bool(is_requested),
                           candidates=[{k: v for k, v in q.items() if k != "x"}
                                       for q in candidates] if include_diagnostics else [])
                path.append(rec)
        return path

    def solve_fractions(self, fractions: Iterable[float], include_diagnostics: bool = True) -> list[dict]:
        fracs = [float(f) for f in fractions]
        if any(f < 0.0 or f >= 1.0 for f in fracs):
            raise ValueError("fractions must lie in [0, 1)")
        if any(b <= a for a, b in zip(fracs, fracs[1:])):
            raise ValueError("fractions must be strictly increasing")
        return self.solve_path([f * self.Tc for f in fracs], include_diagnostics)


def compare_nonlinear_tc_proxy(U: float, Delta_t: float, n_fill: float, V_SO: float,
                               Nk: int = 80) -> dict:
    """Compare HHM determinant Tc with the repository's separate 1%-gap-threshold proxy."""
    from solvers import solve_set_core

    hh = linearized_tc(U, Delta_t, n_fill, V_SO, Nk=Nk)
    legacy = solve_set_core(U, Delta_t, n_fill, V_SO, Nk=Nk)
    tc_hhm = hh["Tc"]
    tc_proxy = legacy.get("Tc")
    return dict(Tc_HHM_linearized=tc_hhm,
                Tc_legacy_1pct=(None if tc_proxy is None else float(tc_proxy)),
                delta_Tc=(None if tc_hhm is None or tc_proxy is None else float(tc_proxy - tc_hhm)),
                relative_delta=(None if tc_hhm is None or tc_proxy is None or tc_hhm == 0.0
                                else float((tc_proxy - tc_hhm) / tc_hhm)),
                legacy_stage1=legacy, linearized=hh)
