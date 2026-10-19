"""
case2_checks.py -- test T6: the SAME pipeline reproduces the Mattis-Bardeen / Hebel-Slichter
Case-II classics (Mattis & Bardeen, Phys. Rev. 109, 1553 (1958); Hebel & Slichter, Phys. Rev.
113, 1504 (1959)), i.e. the two most famous accepted benchmarks of any coherence-factor
machinery, alongside the Case-I (ultrasound) Levy 1963 result of levy1963_check.py.

(1) KERNEL IDENTITY.  Summing our branch-resolved integrand over (sigma, sigma') with the
    Case-II vertex |uu' + vv'|^2 = 1/2[1 + (xi xi' + Dbar Dbar')/(EE')] gives EXACTLY the
    Mattis-Bardeen kernels:

        sum_{sig,sig'} (E/ell)(E'/ell') C_II        = 2 (E E' + Delta^2)/(ell ell')  (scattering)
        sum_{sig,sig'} (E/ell)(E'/ell') (1 - C_II)  = 2 (E E' - Delta^2)/(ell ell')  (pair breaking)
    (the factor 2 is MB's particle-hole symmetrisation, absorbed in their prefactor 2/hw)

    verified below to machine precision.  This is an algebraic proof that our numerical
    procedure IS the Mattis-Bardeen integral, written in the branch (sigma) basis.

(2) HEBEL-SLICHTER PEAK.  s-wave reduction (flat DOS, mu=0, BCS Delta(T)): the omega->0
    Case-II rate ratio diverges like sqrt(Delta/T) below Tc (the coherence peak of NMR T1
    and EM absorption), while the Case-I (eta=+1, ultrasound) ratio is Levy's monotone
    2 f(Delta(T)) < 1.  One code, both classic signatures.

(3) THRESHOLD.  At low T the Case-II absorption is zero for hbar*omega < 2 Delta and turns
    on above (pair breaking), the Mattis-Bardeen threshold.

Run:  python case2_checks.py
"""
from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from absorption import S_scattering, S_normal, coherence
from validate_ua import FlatDOS
import levy1963_check as L


def kernel_identity():
    print("=" * 96)
    print("T6a  branch-summed Case-II kernel vs Mattis-Bardeen 1958 kernels")
    print("=" * 96)
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(2000):
        D = rng.uniform(0.05, 1.5)
        E = D + rng.uniform(1e-6, 3.0)
        Ep = D + rng.uniform(1e-6, 3.0)
        l, lp = np.sqrt(E ** 2 - D ** 2), np.sqrt(Ep ** 2 - D ** 2)
        k_sc = k_pb = 0.0
        for sg in (+1, -1):
            xi, Db = sg * l, D            # xi_sigma = sigma*ell, Dbar = Delta (s-wave, mu=0)
            for sp in (+1, -1):
                xip, Dbp = sp * lp, D
                w = (E / l) * (Ep / lp)
                C = coherence(np.array([E]), np.array([xi]), np.array([Db]),
                              np.array([Ep]), np.array([xip]), np.array([Dbp]),
                              "case2_signed")[0]
                k_sc += w * C
                k_pb += w * (1.0 - C)
        mb_sc = 2.0 * (E * Ep + D ** 2) / (l * lp)
        mb_pb = 2.0 * (E * Ep - D ** 2) / (l * lp)
        worst = max(worst, abs(k_sc - mb_sc) / mb_sc, abs(k_pb - mb_pb) / max(abs(mb_pb), 1e-12))
    print(f"    max relative deviation from 2x MB kernels over 2000 random states: {worst:.2e}")
    print("    (the factor 2 is Mattis-Bardeen's particle-hole symmetrisation, folded into their")
    print("     prefactor 2/hbar*omega; it cancels identically in sigma_1/sigma_N)")
    return worst


def hebel_slichter():
    print("=" * 96)
    print("T6b  Hebel-Slichter coherence peak (Case II) vs Levy monotone suppression (Case I)")
    print("     QUALITATIVE kernel demonstration, NOT sigma_1/sigma_N: the clean-limit ratio")
    print("     diverges as omega->0 (no relaxation); fixed omega=1e-4 is below the delta-smearing")
    print("     scale, so heights are regularisation-dependent.  See REFEREE_RESPONSE.md T6.")
    print("=" * 96)
    dos = FlatDOS(Wc=L.LAM_CUT)
    fr = np.linspace(0.30, 0.995, 26)
    rI, rII = [], []
    for f in fr:
        T = f * L.TC
        D = L.gap_BCS(T)
        w = 1e-4
        SI, _ = S_scattering(dos, T, w, D, 0.0, 0.0, "a_signed", nell=20000)
        SII, _ = S_scattering(dos, T, w, D, 0.0, 0.0, "case2_signed", nell=20000)
        SN = S_normal(dos, T, w, 0.0)
        rI.append(SI / SN)
        rII.append(SII / SN)
        print(f"    T/Tc={f:.3f}   Case I (ultrasound) = {SI/SN:.5f}   Case II (EM/NMR) = {SII/SN:.4f}")
    return fr, np.array(rI), np.array(rII)


def threshold_curve():
    dos = FlatDOS(Wc=L.LAM_CUT)
    T = 0.2 * L.TC
    D = L.gap_BCS(T)
    ws = np.linspace(0.05, 1.2, 60)
    ys = []
    for w in ws:
        SII, _ = S_scattering(dos, T, w, D, 0.0, 0.0, "case2_signed", nell=8000)
        from absorption import S_paircreation
        P, _ = S_paircreation(dos, T, w, D, 0.0, 0.0, "case2_signed", nell=8000)
        SN = S_normal(dos, T, w, 0.0)
        ys.append((SII + P) / SN)
    return ws, np.array(ys), 2 * D


def run_figure():
    fr, rI, rII = hebel_slichter()
    ws, ys, twoD = threshold_curve()
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax[0].plot(fr, rII, "o-", color="tab:red", ms=3.5, lw=1.4,
               label="Case II ($\\eta=-1$, EM/NMR): Hebel-Slichter peak")
    ax[0].plot(fr, rI, "s-", color="tab:blue", ms=3.5, lw=1.4,
               label="Case I ($\\eta=+1$, ultrasound): Levy $2f(\\varepsilon_0)$")
    ax[0].axhline(1.0, color="0.5", lw=0.7)
    ax[0].axvline(1.0, color="0.4", ls=":", lw=0.8)
    ax[0].set_yscale("log")
    ax[0].set_xlabel("T / Tc")
    ax[0].set_ylabel("$\\omega\\to0$ rate ratio")
    ax[0].set_title("one code, both classic signatures (s-wave, flat DOS)")
    ax[0].legend(fontsize=8, loc="lower left")
    ax[1].plot(ws, ys, color="tab:red", lw=1.6)
    ax[1].axvline(twoD, color="k", ls=":", lw=1.2, label=f"$\\hbar\\omega = 2\\Delta = {twoD:.3f}$")
    ax[1].set_xlabel("$\\hbar\\omega$ [t]")
    ax[1].set_ylabel("Case-II absorption / normal")
    ax[1].set_title("Mattis-Bardeen threshold at $2\\Delta$ (T/Tc = 0.2)")
    ax[1].legend(fontsize=8)
    fig.suptitle("Benchmark T6: the branch-resolved pipeline reproduces Mattis-Bardeen 1958 "
                 "and Hebel-Slichter 1959")
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    fig.savefig("figures/caseII_mattis_bardeen_hebel_slichter.png", dpi=150)
    plt.close(fig)
    print("figure: figures/caseII_mattis_bardeen_hebel_slichter.png")


if __name__ == "__main__":
    kernel_identity()
    run_figure()
