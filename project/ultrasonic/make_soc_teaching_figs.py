"""
make_soc_teaching_figs.py -- two explanatory figures for the SOCfigures collection:
  00_SOC_bands_rashba.png                 Rashba band splitting along Gamma-X-M-Gamma,
                                          V_SO = 0 / 0.4082 / 0.8, with mu(T/Tc=0.5) and
                                          the two Fermi momenta (route 1 of the SOC effect)
  10_quasiparticle_branch_construction.png  how E = sqrt(xi^2 + Dbar^2) builds the two
                                          branches eps_sigma(E) = eps* + sigma*sqrt(E^2-Emin^2)/kappa,
                                          and what the signed coherence factor C_I(E,E+hw)
                                          does on each (sigma,sigma') channel
Writes directly into ../../SOCfigures/.
"""
from __future__ import annotations

import os
import sys
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))
OUT = os.path.join(HERE, "..", "..", "SOCfigures")
os.makedirs(OUT, exist_ok=True)

from physics_ua import SCState
from absorption import qp_constants, coherence
from parameters_ua import OMEGA

VSO_MU = {0.0: 3.32, 0.4082: 3.63, 0.8: 4.41}     # mu(0.5 Tc) from results/soc_scan.json


def bands_figure():
    n = 400
    seg = []
    lab = []
    # Gamma-(pi,0): (k,0); X-M: (pi,k); M-Gamma: (k,k)
    k = np.linspace(0, np.pi, n)
    paths = [np.stack([k, 0 * k], -1), np.stack([np.full(n, np.pi), k], -1),
             np.stack([np.pi - k, np.pi - k], -1)]
    names = [r"$\Gamma$", "X", "M", r"$\Gamma$"]
    K = np.concatenate(paths, 0)
    x = np.arange(len(K))
    cs = np.cos(K[:, 0]) + np.cos(K[:, 1])
    g = np.sqrt(np.sin(K[:, 0]) ** 2 + np.sin(K[:, 1]) ** 2)
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0), sharey=True)
    cols = ("tab:blue", "tab:red")
    for ax, v in zip(axes, (0.0, 0.4082, 0.8)):
        for s, c in zip((+1, -1), cols):
            ax.plot(x, -2 * cs - 2 * s * v * g, lw=1.6, color=c,
                    label=f"helicity s={s:+d}" if v == 0.0 else None)
        mu = VSO_MU[v]
        ax.axhline(mu, color="k", ls="--", lw=1.2, label=f"$\\mu$ = {mu:.2f}t")
        # Fermi points
        for s, c in zip((+1, -1), cols):
            e = -2 * cs - 2 * s * v * g
            cross = np.where(np.diff(np.sign(e - mu)) != 0)[0]
            ax.plot(x[cross], mu * np.ones_like(cross), "o", ms=6, mfc="none", mec=c, mew=1.6)
        ax.set_title(f"$V_{{SO}}$ = {v:g}t", fontsize=10)
        ax.set_xticks([0, n - 1, 2 * n - 1, 3 * n - 1])
        ax.set_xticklabels(names)
        ax.grid(alpha=0.3)
        ax.set_ylabel("$\\varepsilon_{k,s}/t$" if v == 0.0 else "")
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Rashba splitting of the two helicity bands and the two Fermi momenta "
                 "(n = 1.885, $\\mu$ at $T/T_c=0.5$); open circles: Fermi crossings",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(OUT, "00_SOC_bands_rashba.png"), dpi=150)
    plt.close(fig)
    print("wrote 00_SOC_bands_rashba.png")


def branch_figure():
    sc = SCState(params=dict(VSO=0.4082, Tc=0.3808, D0_0=0.2855, DS_0=1.1124,
                             Deff_0=0.5489), Nk_mu=200)
    st = sc.at(0.5)
    D0, DS, mu = st["D0"], st["DS"], st["mu"]
    P = qp_constants(D0, DS, mu)
    b, kap, eps_star, Emin = P["b"], P["kappa"], P["estar"], P["Emin"]
    eps = np.linspace(mu - 3.5, mu + 3.0, 800)
    xi = eps - mu
    Db = D0 - b * eps
    E = np.hypot(xi, Db)
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.6))
    ax = axes[0]
    ax.plot(eps, xi, ls=":", lw=1.4, label=r"$\xi(\varepsilon)=\varepsilon-\mu$")
    ax.plot(eps, Db, ls="--", lw=1.4, label=r"$\bar\Delta(\varepsilon)=\Delta^0-b\varepsilon$")
    ax.plot(eps, E, lw=2.0, label=r"$E(\varepsilon)=\sqrt{\xi^2+\bar\Delta^2}$")
    ax.plot(eps, -E, lw=1.0, color="0.7", label=r"$-E(\varepsilon)$ (hole branch)")
    ax.axhline(0, color="k", lw=0.8)
    ax.axvline(eps_star, color="0.5", ls=":", lw=1.0)
    ax.annotate(r"$\varepsilon^*$", (eps_star, 0), xytext=(4, 8), textcoords="offset points")
    ax.axhline(Emin, color="tab:red", ls="--", lw=1.2)
    ax.annotate(f"$E_{{min}}$ = {Emin:.3f}t", (mu - 3.3, Emin), xytext=(4, 4),
                textcoords="offset points", color="tab:red")
    Es = Emin + 0.3
    ell = np.sqrt(Es ** 2 - Emin ** 2)
    for sg, mk in ((+1, "o"), (-1, "s")):
        e_s = eps_star + sg * ell / kap
        ax.plot([e_s], [Es], mk, ms=9, mfc="none", mew=1.8,
                color="tab:green" if sg > 0 else "tab:orange")
        ax.annotate(f"$\\varepsilon_{{{sg:+d}}}(E)$", (e_s, Es), xytext=(6, -14),
                    textcoords="offset points")
    ax.axhline(Es, color="0.4", ls=":", lw=1.0)
    ax.set_xlabel(r"$\varepsilon/t$")
    ax.set_ylabel(r"$t$")
    ax.set_title("(a) branch construction at $T/T_c=0.5$ ($V_{{SO}}=0.4082$): every $E>E_{{min}}$\n"
                 "has two solutions $\\varepsilon_\\sigma(E)=\\varepsilon^*+\\sigma\\sqrt{E^2-E_{{min}}^2}/\\kappa$",
                 fontsize=9)
    ax.legend(fontsize=8, loc="upper center")
    ax.grid(alpha=0.3)
    ax = axes[1]
    Egrid = np.linspace(Emin + 1e-3, Emin + 1.2, 300)
    Ep = Egrid + OMEGA
    for (sa, sb), c in zip([(-1, -1), (-1, 1), (1, -1), (1, 1)],
                           ("tab:blue", "tab:orange", "tab:green", "tab:red")):
        l = np.sqrt(np.maximum(Egrid ** 2 - Emin ** 2, 0))
        lp = np.sqrt(np.maximum(Ep ** 2 - Emin ** 2, 0))
        xi_a = sa * l / kap
        xi_b = sb * lp / kap
        D_a = D0 - b * (eps_star + sa * l / kap)
        D_b = D0 - b * (eps_star + sb * lp / kap)
        C = coherence(Egrid, xi_a, D_a, Ep, xi_b, D_b, "a_signed")
        ax.plot(Egrid, C, lw=1.6, color=c, label=f"$(\\sigma,\\sigma')=({sa:+d},{sb:+d})$")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("$E/t$")
    ax.set_ylabel(r"$C_I(E,\,E+\hbar\omega)$")
    ax.set_title(f"(b) signed coherence factor per branch pair, $\\hbar\\omega={OMEGA}t$:\n"
                 "(-,-) constructive, (+,-)/( -,+) partially cancel -- the weights of S_sc",
                 fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "10_quasiparticle_branch_construction.png"), dpi=150)
    plt.close(fig)
    print("wrote 10_quasiparticle_branch_construction.png")


if __name__ == "__main__":
    bands_figure()
    branch_figure()
