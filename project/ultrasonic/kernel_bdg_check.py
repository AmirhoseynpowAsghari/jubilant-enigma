"""
kernel_bdg_check.py -- independent revalidation of the coherence kernels (audit F-01/F-12).

The production kernels (absorption.coherence) are CLOSED FORMS.  This script derives the
matrix elements NUMERICALLY, with no closed form in the loop:

  * each quasiparticle state is a 2x2 Nambu block  H = [[xi, Dbar], [Dbar, -xi]]
    (real Dbar, the on-site gap of that state's band energy);
  * the positive-energy eigenvector (u, v) is taken from numpy.linalg.eigh;
  * the scattering matrix elements of the two standard vertices are
        M_tau3 = u u' - v v'     (density / ultrasonic-type vertex,  Case I, eta=+1)
        M_1    = u u' + v v'     (Case II, eta=-1 / electromagnetic-NMR type)
    and the pair-creation partners are u v' + v u' (and its Case-II analogue).

Checked, on random physical states (|xi|, |xi'| arbitrary sign, Dbar of either sign):
  (1) M_tau3^2 == C_a  (a_signed, eta=+1)           -- the signed kernel IS the BdG one
  (2) M_1^2    == C_case2                              (Case II)
  (3) M_tau3^2 + (u v' + v u')^2 == 1                   (completeness: C_pair = 1 - C_a)
  (4) variant (d) [1/2(1 - Dbar Dbar'/EE')] equals NONE of the standard BdG elements
      (uu'-vv')^2, (uu'+vv')^2, (uv'+vu')^2, (uv'-vu')^2  -- it is not a BdG matrix element.
  (5) the gap-edge behaviour claimed for (d) is checked directly (E = E' = Dbar).
  (6) the Delta -> 0 limit of (a): C_a -> 1/2 [1 + sgn-projector] -- the branch projector.
Also the BCS limit: flat DOS, const gap, omega->0 gives W/W_N = 2 f(Delta)  (Levy 1963) is
re-derived for variant (a) and for variant (d) in levy1963_check / validate_ua (see T2, T5).
"""
from __future__ import annotations

import numpy as np

from absorption import coherence


def bdg_elements(xi, D, xip, Dp):
    """Numerical BdG eigenvectors -> standard matrix elements (vectorised over states)."""
    xi, D, xip, Dp = map(lambda a: np.asarray(a, float), (xi, D, xip, Dp))
    out = {}
    def eigvec(x, d):
        H = np.zeros(x.shape + (2, 2))
        H[..., 0, 0] = x; H[..., 1, 1] = -x
        H[..., 0, 1] = d; H[..., 1, 0] = d
        w, V = np.linalg.eigh(H)
        # eigh returns ascending eigenvalues: the positive-energy state is the LAST column
        vec = V[..., :, -1]
        E = w[..., -1]
        return vec[..., 0], vec[..., 1], E
    u, v, E = eigvec(xi, D)
    up, vp, Ep = eigvec(xip, Dp)
    out["E"], out["Ep"] = E, Ep
    out["M_tau3"] = u * up - v * vp          # (uu' - vv')
    out["M_1"] = u * up + v * vp             # (uu' + vv')
    out["P_tau3"] = u * vp + v * up          # (uv' + vu')   pair partner
    out["P_1"] = u * vp - v * up             # (uv' - vu')   pair partner, second vertex
    return out


def random_states(n, rng):
    xi = rng.normal(0, 1.5, n); xip = rng.normal(0, 1.5, n)
    D = rng.normal(0, 0.8, n); Dp = rng.normal(0, 0.8, n)
    return xi, D, xip, Dp


def run(n=4000, seed=7, verbose=True):
    rng = np.random.default_rng(seed)
    xi, D, xip, Dp = random_states(n, rng)
    b = bdg_elements(xi, D, xip, Dp)
    E, Ep = b["E"], b["Ep"]
    res = {}
    # (1) signed Case I
    Ca = coherence(E, xi, D, Ep, xip, Dp, "a_signed")
    res["(1) C_a == (uu'-vv')^2"] = float(np.max(np.abs(Ca - b["M_tau3"] ** 2)))
    # (2) Case II
    C2 = coherence(E, xi, D, Ep, xip, Dp, "case2_signed")
    res["(2) C_case2 == (uu'+vv')^2"] = float(np.max(np.abs(C2 - b["M_1"] ** 2)))
    # (3) completeness
    res["(3) (uu'-vv')^2 + (uv'+vu')^2 == 1"] = float(np.max(np.abs(
        b["M_tau3"] ** 2 + b["P_tau3"] ** 2 - 1.0)))
    Cp = 1.0 - Ca
    res["(3b) C_pair(a) == 1 - C_a == (uv'+vu')^2"] = float(np.max(np.abs(
        Cp - b["P_tau3"] ** 2)))
    # (4) variant (d) against every standard BdG element
    Cd = coherence(E, xi, D, Ep, xip, Dp, "d_std")
    cands = {"(uu'-vv')^2": b["M_tau3"] ** 2, "(uu'+vv')^2": b["M_1"] ** 2,
             "(uv'+vu')^2": b["P_tau3"] ** 2, "(uv'-vu')^2": b["P_1"] ** 2}
    dist_d = {k: float(np.max(np.abs(Cd - v))) for k, v in cands.items()}
    res["(4) min_k max|C_d - BdG_k|  (d is NOT a BdG element if >>0)"] = min(dist_d.values())
    res["(4a) per-candidate distances for (d)"] = dist_d
    # (5) gap-edge values for (d) vs (a)
    Dg = 0.37
    Cd_edge = coherence(np.array([Dg]), np.array([0.0]), np.array([Dg]),
                        np.array([Dg]), np.array([0.0]), np.array([Dg]), "d_std")
    Ca_edge = coherence(np.array([Dg]), np.array([0.0]), np.array([Dg]),
                        np.array([Dg]), np.array([0.0]), np.array([Dg]), "a_signed")
    res["(5) at gap edge E=E'=Dbar (xi=0): C_a, C_d"] = (float(Ca_edge[0]), float(Cd_edge[0]))
    # (6) Delta -> 0 : (a) becomes the branch projector  1/2 [1 + sgn(xi xi')] (for xi xi' != 0)
    xs = rng.normal(0, 1.0, 2000); xps = rng.normal(0, 1.0, 2000)
    z = np.zeros_like(xs)
    Ca0 = coherence(np.abs(xs), xs, z, np.abs(xps), xps, z, "a_signed")
    proj = 0.5 * (1.0 + np.sign(xs) * np.sign(xps))
    res["(6) Delta->0: max|C_a - 1/2(1+sgn xi sgn xi')|"] = float(np.max(np.abs(Ca0 - proj)))
    # (7) the production (d) sign now equals the BdG-independent convention used by
    #     bruteforce_ua after the F-01 fix (checked in test_fixes.py)
    if verbose:
        print("=" * 84)
        print("kernel_bdg_check: numerical BdG eigenvectors vs production closed forms")
        print("=" * 84)
        for k, v in res.items():
            if isinstance(v, dict):
                for kk, vv in v.items():
                    print(f"   {k}: {kk:>16s}  {vv:.3e}")
            elif isinstance(v, tuple):
                print(f"   {k}: {v[0]:.6f}  {v[1]:.6f}")
            else:
                print(f"   {k}: {v:.3e}")
    return res


if __name__ == "__main__":
    r = run()
    ok = (r["(1) C_a == (uu'-vv')^2"] < 1e-12 and r["(2) C_case2 == (uu'+vv')^2"] < 1e-12
          and r["(3) (uu'-vv')^2 + (uv'+vu')^2 == 1"] < 1e-12
          and r["(3b) C_pair(a) == 1 - C_a == (uv'+vu')^2"] < 1e-12
          and r["(6) Delta->0: max|C_a - 1/2(1+sgn xi sgn xi')|"] < 1e-12)
    print("RESULT:", "PASS (signed kernel (a) and Case II are BdG matrix elements)" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)
