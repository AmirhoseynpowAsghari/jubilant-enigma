"""
solvers.py -- T=0 grid-search solution, Tc bisection, and warm-started T sweeps.
Part of the corrected reconstruction of
  AmirhoseynpowAsghari/power-absorbed-in-the-presence-of-external-time-dependent-perturbation
Stage-1 logic of the repository (branch master, commit cda0d56), parameterised over
(U, Delta_t, n, V_SO), with the corrections of the 2026-10-24 audit campaign:

  F-02  the k-grid is the MIDPOINT rule (kgrid.py); solve_set(check_conv=True) re-solves
        at 2*Nk and records the difference, so grid convergence is measured, not assumed.
  F-07  roots are sign-normalised by flipping the PAIR (D0,DS)->(-D0,-DS) (an exact symmetry
        of the gap equations) so that D0 >= 0, with DS kept SIGNED.  The original abs() of
        each component independently turns a genuine root into a non-root whenever D0 and DS
        have opposite signs (the audited U=-75 "Tc=12.7" entries are of this kind).  Every
        candidate is residual-checked after normalisation and must have mu inside the band.
  F-08  the Tc bisection only runs when the bracket is VERIFIED: the escalation loop must
        find T_max with Deff(T_max) <= threshold.  If it never does, the set is returned
        with Tc=None, Tc_bracketed=False and Tc_lower_bound=T_max instead of a fabricated
        "Tc" at the top of the (arbitrary) escalation ladder -- the defect that shipped
        Tc(U=-75) up to 13.17 t without any convergence signal.

The E10 multi-start escalation of _gap_at_T (trivial-root-slide guard) is unchanged.
"""
from __future__ import annotations

import json
import time
import numpy as np
from scipy.optimize import fsolve, brentq

from parameters_tc import (T, Nk_DEF, BISECTION_TOL, BISECTION_ITER,
                           GAP_THRESHOLD_FRAC, MIN_GAP_FOR_SC, GRID_CONV_TOL,
                           MU_WINDOW_POLICY)
from gap_equations import make_residual, find_mu_normal, effective_gap
from kgrid import build_kgrid

RESID_TOL = 1e-6          # max-|residual| a folded/re-polished root must satisfy (F-07)


def _residual_norm(residual, x, beta, n_fill, epsilon_ks, term):
    """max |r_i| of the three self-consistency residuals at x = [D0, DS, mu]."""
    try:
        r = residual(list(x), beta, n_fill, epsilon_ks, term)
    except Exception:
        return float("inf")
    r = np.asarray(r, float)
    if not np.all(np.isfinite(r)):
        return float("inf")
    return float(np.max(np.abs(r)))


def _folded_root(residual, sol, beta, n_fill, epsilon_ks, term, mu_band=None):
    """Sign-normalise a converged fsolve output and verify it is a root (F-07 / F-08).

    The residual is invariant under the EXACT symmetry (D0, DS, mu) -> (-D0, -DS, mu)
    (r1, r2 odd in the gap pair, r3 even).  The physical normalisation is therefore
    D0 >= 0 with DS kept SIGNED: the pair is flipped as a whole.  Taking abs() of D0 and DS
    independently (the original repository logic) produces a point that is generally NOT a
    root when D0 and DS have opposite signs -- the U=-75 "Tc=12.7" entries of the audited
    tree are exactly such artefacts.  The candidate must pass the residual check after
    normalisation and, if `mu_band` is given, its chemical potential must lie inside the
    band [mu_band[0], mu_band[1]] (a Fermi surface must exist for the filling).
    """
    x = [sol[0], sol[1], sol[2]]
    if x[0] < 0.0:
        x = [-x[0], -x[1], x[2]]
    cand = None
    if _residual_norm(residual, x, beta, n_fill, epsilon_ks, term) <= RESID_TOL:
        cand = x
    else:
        try:
            s2, _, ier2, _ = fsolve(residual, x, args=(beta, n_fill, epsilon_ks, term),
                                    full_output=True)
        except Exception:
            return None
        if ier2 != 1:
            return None
        x2 = [s2[0], s2[1], s2[2]]
        if x2[0] < 0.0:
            x2 = [-x2[0], -x2[1], x2[2]]
        if _residual_norm(residual, x2, beta, n_fill, epsilon_ks, term) <= RESID_TOL:
            cand = x2
    if cand is None:
        return None
    if mu_band is not None and not (mu_band[0] <= cand[2] <= mu_band[1]):
        return None
    return cand


def _band(epsilon_ks):
    """Band window [min eps, max eps] of the k-grid, or None under the "relaxed" policy."""
    if MU_WINDOW_POLICY != "strict":
        return None
    return (float(np.min(epsilon_ks)), float(np.max(epsilon_ks)))


def mu_window(epsilon_ks, mu):
    """Diagnostic: is mu inside the band window (independent of the policy)?"""
    lo, hi = float(np.min(epsilon_ks)), float(np.max(epsilon_ks))
    return dict(mu=float(mu), band=[lo, hi], inside=bool(lo <= mu <= hi),
                excess=float(max(mu - hi, lo - mu, 0.0)))


def find_T0_solution(n_fill, epsilon_ks, term, residual):
    band = _band(epsilon_ks)
    mu0 = find_mu_normal(n_fill, np.inf, epsilon_ks)
    best_sol, best_gap = None, 0.0
    for D0i in (0.1, 0.5, 1.0, 2.0, 5.0):
        for DSi in (0.1, 0.5, 1.0, 5.0, 10.0):
            for dmu in (-2.0, -0.5, 0.0, 0.5, 2.0):
                try:
                    sol, _, ier, _ = fsolve(residual, [D0i, DSi, mu0 + dmu],
                                            args=(np.inf, n_fill, epsilon_ks, term),
                                            full_output=True)
                    if ier != 1:
                        continue
                    folded = _folded_root(residual, sol, np.inf, n_fill, epsilon_ks, term,
                                          mu_band=band)
                    if folded is None:
                        continue
                    D0c, DSc, muc = folded
                    Deff = effective_gap(D0c, DSc, epsilon_ks)
                    if Deff > best_gap:
                        best_gap, best_sol = Deff, (D0c, DSc, muc)
                except Exception:
                    pass
    if best_sol is None or best_gap < MIN_GAP_FOR_SC:
        return 0.0, 0.0, mu0, 0.0
    D0, DS, mu = best_sol
    return D0, DS, mu, effective_gap(D0, DS, epsilon_ks)


def diagnose_outside_band(n_fill, epsilon_ks, term, residual):
    """Diagnostic only (F-08 transparency): the best T=0 root WITHOUT the in-band requirement.

    Returned in the stage-1 JSON when the physical-window rule rejects every root, so that a
    reader can see what the unrestricted solver finds (it is NOT used for any Tc).
    """
    mu0 = find_mu_normal(n_fill, np.inf, epsilon_ks)
    best = None
    for D0i in (0.1, 0.5, 1.0, 2.0, 5.0):
        for DSi in (-10.0, -5.0, -2.0, -0.5, 0.5, 2.0, 5.0, 10.0):
            for dmu in (-2.0, -0.5, 0.0, 0.5, 2.0):
                try:
                    sol, _, ier, _ = fsolve(residual, [D0i, DSi, mu0 + dmu],
                                            args=(np.inf, n_fill, epsilon_ks, term),
                                            full_output=True)
                except Exception:
                    continue
                if ier != 1:
                    continue
                folded = _folded_root(residual, sol, np.inf, n_fill, epsilon_ks, term)
                if folded is None:
                    continue
                Deff = effective_gap(folded[0], folded[1], epsilon_ks)
                if best is None or Deff > best["Deff_0"]:
                    band = _band(epsilon_ks)
                    best = dict(D0_0=folded[0], DS_0=folded[1], mu_0=folded[2], Deff_0=Deff,
                                band=list(band), mu_inside_band=bool(band[0] <= folded[2] <= band[1]))
    return best


SEED_FRACS = (1.0, 0.8, 0.6, 0.4, 0.2, 0.05)
SEED_FRACS_FULL = tuple(round(0.05 * i, 3) for i in range(1, 21))


def _gap_at_T(TT, guess, n_fill, epsilon_ks, term, residual, t0=None, full=False):
    """Effective gap at temperature TT (max-gap over converged, residual-checked starts).

    The trivial solution (D0, DS) = 0 is an EXACT root of the residual at every T, and a
    single warm-started fsolve can slide into it even deep inside the SC state (the
    nontrivial root's basin is narrow at intermediate T for some (U, V_SO)).  Such a slide
    used to cap the Tc bisection at the slide temperature and produced notches in Tc(U)
    (e.g. U=17.5, V_SO=0.6: Tc 1.376 instead of ~2.55).  Multi-start over scaled copies of
    the T=0 root -- same max-gap philosophy as find_T0_solution -- removes the artifact;
    `full=True` widens the seed grid (used as escalation in the bisection).  Every accepted
    start passes the F-07 residual check after sign normalisation.
    """
    seeds = [guess]
    if t0 is not None:
        for f in (SEED_FRACS_FULL if full else SEED_FRACS):
            seeds.append([t0[0] * f, t0[1] * f, t0[2]])
    best = (0.0, guess)
    for g in seeds:
        try:
            s, _, ier, _ = fsolve(residual, g, args=(1.0 / TT, n_fill, epsilon_ks, term),
                                  full_output=True)
        except Exception:
            continue
        if ier != 1:
            continue
        folded = _folded_root(residual, s, 1.0 / TT, n_fill, epsilon_ks, term,
                              mu_band=_band(epsilon_ks))
        if folded is None:
            continue
        D0, DS, mu = folded
        Deff = effective_gap(D0, DS, epsilon_ks)
        if Deff > best[0]:
            best = (Deff, [D0, DS, mu])
    return best


def solve_set_core(U, Delta_t, n_fill, V_SO, Nk=Nk_DEF):
    """Tc, D0(0), DS(0), Deff(0), mu(0) for one parameter set on the midpoint grid.

    Returns Tc=0.0 + Tc_bracketed=False + Tc_lower_bound when the gap never falls below the
    threshold inside the escalation ladder (F-08): no Tc is invented; consumers treat Tc<=0
    as "no superconductivity" and the flag keeps the distinction visible in every JSON.
    """
    eps, term = build_kgrid(V_SO, U, Delta_t, Nk)
    residual = make_residual(Delta_t)
    base = dict(U=U, Delta_t=Delta_t, n=n_fill, V_SO=V_SO, Nk=int(Nk), grid="midpoint")
    D0_T0, DS_T0, mu_T0, Deff_T0 = find_T0_solution(n_fill, eps, term, residual)
    if Deff_T0 < MIN_GAP_FOR_SC:
        diag = diagnose_outside_band(n_fill, eps, term, residual)
        return dict(base, Tc=0.0, D0_0=0.0, DS_0=0.0, Deff_0=0.0, mu_0=mu_T0,
                    Tc_bracketed=True, no_inband_root=True, diagnostic_unrestricted=diag,
                    policy=MU_WINDOW_POLICY)
    threshold = GAP_THRESHOLD_FRAC * Deff_T0
    t0 = [D0_T0, DS_T0, mu_T0]
    T_max = 0.57 * Deff_T0 * 3.0
    Deff, sol = _gap_at_T(T_max, t0, n_fill, eps, term, residual, t0=t0)
    bracketed = Deff <= threshold
    for _ in range(10):
        if bracketed:
            break
        T_max *= 1.5
        Deff, sol = _gap_at_T(T_max, sol, n_fill, eps, term, residual, t0=t0)
        bracketed = Deff <= threshold
    if not bracketed:
        print(f"  WARNING (F-08): U={U}, Dt={Delta_t}, V_SO={V_SO}, Nk={Nk}: Deff(T) "
              f"survives to T_max={T_max:.4f}; Tc NOT bracketed -> flagged (Tc=None)",
              flush=True)
        # Tc is NOT established: report 0 (= "no SC" for every consumer, which already
        # excludes Tc<=0) together with the explicit flag and the lower bound (F-08).
        return dict(base, Tc=0.0, Tc_lower_bound=float(T_max), Tc_bracketed=False,
                    D0_0=D0_T0, DS_0=DS_T0, Deff_0=0.0, mu_0=mu_T0)
    T_low, T_high = 0.0, T_max
    sol = list(t0)
    for _ in range(BISECTION_ITER):
        T_mid = 0.5 * (T_low + T_high)
        Deff, out = _gap_at_T(T_mid, sol, n_fill, eps, term, residual, t0=t0)
        if Deff <= threshold and (T_high - T_low) > 0.1:
            # narrow-basin guard: retry this midpoint with the full seed grid before
            # accepting "no SC here" (prevents trivial-root slides from capping Tc)
            Deff2, out2 = _gap_at_T(T_mid, sol, n_fill, eps, term, residual,
                                    t0=t0, full=True)
            if Deff2 > Deff:
                Deff, out = Deff2, out2
        if Deff > threshold:
            T_low, sol = T_mid, out
        else:
            T_high = T_mid
        if T_high - T_low < BISECTION_TOL:
            break
    Tc = 0.5 * (T_low + T_high)
    return dict(base, Tc=Tc, Tc_bracketed=True,
                D0_0=D0_T0, DS_0=DS_T0, Deff_0=Deff_T0, mu_0=mu_T0,
                mu_window=mu_window(eps, mu_T0), policy=MU_WINDOW_POLICY)


def solve_set(U, Delta_t, n_fill, V_SO, Nk=Nk_DEF, check_conv=False):
    """solve_set_core plus the optional Nk -> 2Nk convergence check (F-02).

    check_conv=True re-solves on the doubled midpoint grid and records Tc/Deff at both
    resolutions plus an `ok` flag (|dTc|, |dDeff| <= GRID_CONV_TOL).  Production stage-1
    runs (scanning_tc) and the test suite enable it; the bulk parameter scans rely on the
    converged-grid evidence recorded there plus their own spot checks (documented policy,
    CORRECTION_LOG.md).
    """
    r = solve_set_core(U, Delta_t, n_fill, V_SO, Nk=Nk)
    if check_conv:
        r2 = solve_set_core(U, Delta_t, n_fill, V_SO, Nk=2 * int(Nk))
        conv = dict(Nk=int(Nk), Nk2=2 * int(Nk), Tc2=r2["Tc"], D0_02=r2["D0_0"],
                    Deff_02=r2["Deff_0"])
        if r["Tc_bracketed"] and r2["Tc_bracketed"]:
            conv["dTc"] = float(abs(r2["Tc"] - r["Tc"]))
            conv["dDeff_0"] = float(abs(r2["Deff_0"] - r["Deff_0"]))
            conv["ok"] = bool(conv["dTc"] <= GRID_CONV_TOL
                              and conv["dDeff_0"] <= GRID_CONV_TOL)
        else:
            conv["dTc"] = None
            conv["dDeff_0"] = None
            conv["ok"] = bool(r["Tc_bracketed"] == r2["Tc_bracketed"])
        r["conv"] = conv
        if not conv["ok"]:
            print(f"  WARNING (F-02): grid convergence check FAILED for U={U}, "
                  f"Dt={Delta_t}, V_SO={V_SO}: {conv}", flush=True)
    return r


def sweep_T(set_dict, fracs, Nk=Nk_DEF):
    """
    Warm-started fsolve sweep of (D0, DS, mu) over T = frac * Tc.
    Returns lists aligned with fracs.  (Used for the auxiliary panels of the
    regenerated figures; the W/W_N evaluation itself uses the repository's
    rigid gap shape phi(T), see wwn_eta1/model.py.)
    """
    eps, term = build_kgrid(set_dict["V_SO"], set_dict["U"], set_dict["Delta_t"], Nk)
    residual = make_residual(set_dict["Delta_t"])
    guess = [set_dict["D0_0"], set_dict["DS_0"], set_dict["mu_0"]]
    out = []
    for f in fracs:
        TT = f * set_dict["Tc"]
        if TT <= 1e-9:
            D0, DS, mu = guess
        else:
            Deff, guess = _gap_at_T(TT, guess, set_dict["n"], eps, term, residual)
            D0, DS, mu = guess
        out.append(dict(frac=f, T=TT, D0=float(D0), DS=float(DS), mu=float(mu),
                        Deff=effective_gap(D0, DS, eps)))
    return out
