"""
plotting_ua.py -- figure gallery of the corrected stage-2 (run: python plotting_ua.py).
Figures are written to figures/ with the repository's naming style.
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import io_utils_ua as io
from parameters_ua import OMEGA, T_HOP
from dos_ua import NormalDOS
from physics_ua import SCState
from absorption import qp_constants, branch, coherence, ratio_all, S_paircreation

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3, "grid.linewidth": 0.5,
                     "font.size": 9, "legend.fontsize": 8, "figure.facecolor": "white"})
COLV = {"a_signed": "tab:blue", "b_abs": "tab:orange", "c_paper": "tab:red", "d_std": "tab:green"}


def load_sweep(v):
    return io.load(f"temp_sweep_VSO_{v:.4f}.json")


def fig_temperature_sweep():
    """W/W_N vs T/Tc for the master set: corrected (a), scattering-only quantitative result.
    The pair channel is an OUTLOOK only (REPORT App. A2) and is not plotted here."""
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharey=True)
    for i, v in enumerate((0.0, 0.4082, 0.8)):
        r = load_sweep(v)
        f = np.array(r["fracs"])
        a = np.array(r["ratios"]["a_signed"])
        ax[i].plot(f, a, color="tab:blue", lw=1.8, label="scattering only (corrected, quantitative)")
        ax[i].axhline(1.0, color="0.5", lw=0.7)
        ax[i].axvline(1.0, color="0.4", ls=":", lw=0.8)
        ax[i].set_title(f"V_SO = {v:.4f},  Tc = {r['Tc']:.4f}")
        ax[i].set_xlabel("T / Tc")
        if i == 0:
            ax[i].set_ylabel("W / W_N")
            ax[i].set_ylim(0, 1.15)
    ax[0].legend(loc="upper left")
    fig.suptitle("Corrected ultrasonic attenuation ratio (master set U=95t, $\\Delta_t$=4.5t, n=1.885, $\\hbar\\omega$=0.7t)")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{FIG}/W_WN_temperature_sweep.png", dpi=150)
    plt.close(fig)


def fig_coherence_variants():
    """the four coherence-factor choices compared (why only (a) survives)."""
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0), sharey=True)
    for i, v in enumerate((0.0, 0.4082, 0.8)):
        r = load_sweep(v)
        f = np.array(r["fracs"])
        for key, lab in (("a_signed", "(a) signed $\\xi\\xi'$ (correct)"),
                         ("b_abs", "(b) $|\\xi||\\xi'|$"),
                         ("c_paper", "(c) paper's final"),
                         ("d_std", "(d) textbook Case I")):
            ax[i].plot(f, r["ratios"][key], color=COLV[key], lw=1.4, label=lab)
        ax[i].axhline(1.0, color="0.5", lw=0.7)
        ax[i].axvline(1.0, color="0.4", ls=":", lw=0.8)
        ax[i].set_yscale("log")
        ax[i].set_title(f"V_SO = {v:.4f}")
        ax[i].set_xlabel("T / Tc")
    ax[0].set_ylabel("W / W_N")
    ax[0].legend(loc="lower right")
    fig.suptitle("Coherence-factor variants: only the signed form (a) stays below 1 with a consistent $\\Delta\\to0$ limit")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{FIG}/W_WN_coherence_variants.png", dpi=150)
    plt.close(fig)


def fig_branch_structure():
    """NEW: what sigma means -- the two band-energy roots at fixed quasiparticle energy."""
    import sc_loader
    p = sc_loader.sc_params_at(0.4082)
    sc = SCState.from_params(p["VSO"], p["Tc"], p["D0_0"], p["DS_0"], p["Deff_0"])
    st = sc.at(0.5)
    dos = NormalDOS(sc.VSO)
    P = qp_constants(st["D0"], st["DS"], st["mu"])
    emin, kappa, estar = P["Emin"], P["kappa"], P["estar"]
    eps = np.linspace(-dos.Wband, dos.Wband, 1200)
    Dbar = st["D0"] - st["DS"] * eps / 4.0
    Ecurve = np.hypot(eps - st["mu"], Dbar)
    Esample = emin + 0.45 * (dos.Wband - emin)
    roots = [estar + s * np.sqrt(Esample ** 2 - emin ** 2) / kappa for s in (+1, -1)]
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
    ax[0].plot(eps, Ecurve, "k", lw=1.6, label="$E(\\varepsilon)=\\sqrt{(\\varepsilon-\\mu)^2+\\bar\\Delta(\\varepsilon)^2}$")
    ax[0].axhline(Esample, color="tab:red", ls="--", lw=1.0)
    ax[0].axvline(estar, color="0.4", ls=":", lw=0.9)
    ax[0].plot(roots, [Esample] * 2, "o", color="tab:red", ms=7)
    ax[0].annotate("$\\sigma=+1$ (particle-like root)", (roots[0], Esample),
                   textcoords="offset points", xytext=(-10, 14), ha="right", color="tab:red")
    ax[0].annotate("$\\sigma=-1$ (hole-like root)", (roots[1], Esample),
                   textcoords="offset points", xytext=(10, -20), color="tab:red")
    ax[0].annotate("$\\varepsilon^*$", (estar, 0.02), textcoords="offset points", xytext=(3, 30), color="0.3")
    ax[0].annotate("$E_{min}$", (estar, emin), textcoords="offset points", xytext=(6, 0), color="k")
    ax[0].plot([estar], [emin], "k*", ms=10)
    ax[0].set_xlabel("band energy $\\varepsilon$  [t]")
    ax[0].set_ylabel("quasiparticle energy E  [t]")
    ax[0].set_title("fixed E $\\Rightarrow$ two roots $\\varepsilon_\\sigma(E)=\\varepsilon^*+\\sigma\\ell/\\kappa$")
    ax[0].legend(loc="upper center", fontsize=7.5)
    ax[1].plot(dos.grid, dos.dens, "k", lw=1.4, label="$N(\\varepsilon)=\\sum_s N_s(\\varepsilon)$ (helicity sum)")
    for s, c in ((+1, "tab:blue"), (-1, "tab:orange")):
        e = estar + s * np.sqrt(Esample ** 2 - emin ** 2) / kappa
        ax[1].axvline(e, color=c, ls="--", lw=1.1)
        ax[1].annotate(f"$\\sigma={s:+d}$", (e, dos.dens.max() * 0.9), color=c,
                       textcoords="offset points", xytext=(4, 0))
    ax[1].axvline(-dos.Wband, color="0.5", ls=":", lw=0.8)
    ax[1].axvline(dos.Wband, color="0.5", ls=":", lw=0.8)
    ax[1].set_xlabel("band energy $\\varepsilon$  [t]")
    ax[1].set_ylabel("normal DOS  $N(\\varepsilon)$")
    ax[1].set_title("helicity s enters only through $N_s$; $\\sigma$ through the roots")
    ax[1].legend(loc="upper left", fontsize=7.5)
    fig.suptitle("The branch index $\\sigma$ (particle/hole root) vs the helicity index s (Rashba band)"
                 f"   [master set, V_SO=0.4082, T/Tc=0.5]")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(f"{FIG}/branch_structure_sigma_vs_s.png", dpi=150)
    plt.close(fig)


def fig_coherence_factors():
    """NEW: C_I(E) on each branch pair, for the four variants, at two temperatures."""
    import sc_loader
    p = sc_loader.sc_params_at(0.4082)
    sc = SCState.from_params(p["VSO"], p["Tc"], p["D0_0"], p["DS_0"], p["Deff_0"])
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3), sharey=True)
    for j, frac in enumerate((0.6, 0.999)):
        st = sc.at(frac)
        P = qp_constants(st["D0"], st["DS"], st["mu"])
        emin = P["Emin"]
        E = np.linspace(emin + 1e-6, emin + 2.5, 400)
        Ep = E + OMEGA
        for (s, sp), ls in (((-1, -1), "-"), ((-1, +1), "--"), ((+1, +1), ":")):
            eb, xib, Db = branch(np.sqrt(np.maximum(E ** 2 - emin ** 2, 0)), s, P, st["mu"], st["D0"])
            ea, xia, Da = branch(np.sqrt(np.maximum(Ep ** 2 - emin ** 2, 0)), sp, P, st["mu"], st["D0"])
            for key in ("a_signed", "b_abs", "c_paper"):
                C = coherence(E, xib, Db, Ep, xia, Da, key)
                ax[j].plot(E, C, ls=ls, color=COLV[key], lw=1.2,
                           label=f"{key}  ({s:+d},{sp:+d})" if j == 0 else None)
        ax[j].axhline(1.0, color="0.5", lw=0.7)
        ax[j].axhline(0.0, color="0.5", lw=0.7)
        ax[j].set_title(f"T/Tc = {frac}")
        ax[j].set_xlabel("E  [t]")
    ax[0].set_ylabel("coherence factor C(E, E+$\\hbar\\omega$)")
    ax[0].legend(loc="lower left", fontsize=7, ncol=1)
    fig.suptitle("Coherence factors per branch pair: (a) stays in [0,1]; (b),(c) leave the interval as $\\Delta\\to0$")
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    fig.savefig(f"{FIG}/coherence_factors_per_branch.png", dpi=150)
    plt.close(fig)


def fig_sum_rule_thresholds():
    """NEW: scattering + pair = W_N exactly, and the two thresholds."""
    r = load_sweep(0.4082)
    f = np.array(r["fracs"])
    a = np.array(r["ratios"]["a_signed"])
    pr = np.array(r["ratios"]["pair"])
    tot = np.array(r["ratios"]["total"])
    deff = np.array(r["Deff"])
    tem = np.array(r["two_Emin"])
    cut_user = f[np.argmin(np.abs(2 * deff - OMEGA))]
    cut_true = f[np.argmin(np.abs(tem - OMEGA))]
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    ax[0].plot(f, a, color="tab:blue", lw=1.7, label="scattering $W/W_N$")
    ax[0].plot(f, pr, color="tab:red", lw=1.7, label="pair creation $W_p/W_N$")
    ax[0].plot(f, tot, color="k", ls="--", lw=1.2, label="sum")
    ax[0].axhline(1.0, color="0.5", lw=0.7)
    ax[0].axvline(cut_user, color="tab:purple", ls=":", lw=1.2,
                  label=f"paper cut $2\\Delta_{{eff}}=\\hbar\\omega$ (T/Tc={cut_user:.3f})")
    ax[0].axvline(cut_true, color="tab:green", ls=":", lw=1.2,
                  label=f"true threshold $2E_{{min}}=\\hbar\\omega$ (T/Tc={cut_true:.3f})")
    ax[0].set_xlabel("T / Tc")
    ax[0].set_ylabel("ratio")
    ax[0].set_title("channels and thresholds (V_SO = 0.4082)")
    ax[0].legend(loc="center left", fontsize=7.5)
    ax[1].semilogy(f, np.abs(tot - 1.0) + 1e-12, color="k", lw=1.4)
    ax[1].set_xlabel("T / Tc")
    ax[1].set_ylabel("$|(W+W_p)/W_N - 1|$")
    ax[1].set_title("OUTLOOK only: (W+W_p)/W_N against the band reference;\nnot an identity (different sheets, REPORT App. A2)")
    fig.suptitle("Pair creation is identically zero inside the paper's window $2\\Delta_{eff}>\\hbar\\omega$;"
                 " above the true threshold the pair channel is an OUTLOOK (no quantitative claim)")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(f"{FIG}/sum_rule_and_thresholds.png", dpi=150)
    plt.close(fig)


def fig_lattice_crosscheck():
    """NEW: independent brute-force lattice sum vs the energy-space pipeline."""
    import bruteforce_ua as bf
    from absorption import S_scattering, S_normal
    import sc_loader
    p = sc_loader.sc_params_at(0.4082)
    sc = SCState.from_params(p["VSO"], p["Tc"], p["D0_0"], p["DS_0"], p["Deff_0"])
    dos = NormalDOS(sc.VSO)
    fracs = (0.3, 0.5, 0.7, 0.9, 0.99, 0.999)
    es, bfv = [], []
    for fr in fracs:
        st = sc.at(fr)
        r = ratio_all(dos, st, OMEGA, variant_list=["a_signed"])
        Sbf = bf.S_sc_bf_from_state(st, OMEGA, "a_signed", VSO=sc.VSO)
        es.append(r["a_signed_ratio"])
        bfv.append(Sbf / r["S_N"])
    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax[0].plot(fracs, es, "o-", color="tab:blue", label="energy-space pipeline")
    ax[0].plot(fracs, bfv, "s--", color="tab:red", ms=5, label="brute-force lattice sum")
    ax[0].set_xlabel("T / Tc")
    ax[0].set_ylabel("W / W_N")
    ax[0].set_title("V_SO = 0.4082")
    ax[0].legend()
    rel = np.abs(np.array(bfv) - np.array(es)) / np.maximum(np.array(es), 1e-30)
    ax[1].semilogy(fracs, rel, "k-o", ms=4)
    ax[1].set_xlabel("T / Tc")
    ax[1].set_ylabel("relative difference")
    ax[1].set_title("agreement 0.7% near Tc (gap-edge regularisation at low T)")
    fig.tight_layout()
    fig.savefig(f"{FIG}/lattice_crosscheck.png", dpi=150)
    plt.close(fig)


def fig_vso_scan():
    """W/W_N vs V_SO at fixed T/Tc (U=95t)."""
    sc = io.load("vso_scan.json")["95"]
    V = np.array(sc["V_SO"])
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for f in ("0.3", "0.5", "0.7", "0.9"):
        y = np.array([x if x is not None else np.nan for x in sc["ratio"][f]], float)
        ax.plot(V, y, "-o", ms=3, lw=1.2, label=f"T/Tc = {f}")
    ax.axhline(1.0, color="0.5", lw=0.7)
    ax.set_xlabel("V_SO [t]")
    ax.set_ylabel("W / W_N")
    ax.set_title("Corrected attenuation vs spin-orbit coupling (U=95t, $\\hbar\\omega$=0.7t)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(f"{FIG}/W_WN_vs_VSO.png", dpi=150)
    plt.close(fig)


def fig_pair_channel():
    """NEW: the pair-creation channel -- integrand over the shell E + E' = hbar*omega,
    and W_p/W_N vs T/Tc with the true threshold 2 E_min = hbar*omega."""
    import sc_loader
    from absorption import C_pair, E_of_ell
    p = sc_loader.sc_params_at(0.4082)
    sc = SCState.from_params(p["VSO"], p["Tc"], p["D0_0"], p["DS_0"], p["Deff_0"])
    st = sc.at(0.95)
    dos = NormalDOS(sc.VSO)
    P = qp_constants(st["D0"], st["DS"], st["mu"])
    emin, kappa = P["Emin"], P["kappa"]
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    # panel a: integrand vs E_b on [Emin, hw/2]
    Eb = np.linspace(emin + 1e-6, OMEGA / 2 - 1e-6, 500)
    Ea = OMEGA - Eb
    tot = np.zeros_like(Eb)
    for (sg, sp), c in (((-1, -1), "tab:blue"), ((+1, -1), "tab:orange"), ((-1, +1), "tab:green")):
        lb = np.sqrt(np.maximum(Eb ** 2 - emin ** 2, 0))
        la = np.sqrt(np.maximum(Ea ** 2 - emin ** 2, 0))
        ebs, xibs, Dbs = branch(lb, sg, P, st["mu"], st["D0"])
        eas, xias, Das = branch(la, sp, P, st["mu"], st["D0"])
        nub = dos.N(ebs) / kappa
        nua = dos.N(eas) / kappa
        Cp = C_pair(Eb, xibs, Dbs, Ea, xias, Das, "a_signed")
        occ = 1.0 - 1.0 / (np.exp(Eb / st["T"]) + 1) - 1.0 / (np.exp(Ea / st["T"]) + 1)
        part = nub * nua * Cp * np.clip(occ, 0, None)
        tot += part
        ax[0].plot(Eb, part, color=c, lw=1.2, label=f"({sg:+d},{sp:+d})")
    ax[0].plot(Eb, tot, "k", lw=1.8, label="sum")
    ax[0].axvline(emin, color="0.5", ls=":", lw=1.0)
    ax[0].axvline(OMEGA - emin, color="0.5", ls=":", lw=1.0)
    ax[0].set_xlabel("$E_b$  [t]")
    ax[0].set_ylabel("pair integrand (arb. units)")
    ax[0].set_title(f"shell $E_b+E_a=\\hbar\\omega$, T/Tc=0.95; window "
                    f"$[E_{{min}},\\hbar\\omega-E_{{min}}]$, $2E_{{min}}$={2*emin:.3f}")
    ax[0].legend(fontsize=7.5)
    # panel b: W_p/W_N vs T/Tc from the production sweep
    r = load_sweep(0.4082)
    f = np.array(r["fracs"])
    pr = np.array([x if x is not None else np.nan for x in r["ratios"]["pair"]], float)
    tem = np.array(r["two_Emin"])
    thr = f[np.argmin(np.abs(tem - OMEGA))]
    ax[1].plot(f, pr, color="tab:red", lw=1.8, label="$W_p/W_N$ (creation channel)")
    ax[1].axvline(thr, color="tab:green", ls=":", lw=1.2,
                  label=f"true threshold $2E_{{min}}=\\hbar\\omega$ (T/Tc={thr:.3f})")
    ax[1].axhline(1.0, color="0.5", lw=0.7)
    ax[1].set_xlabel("T / Tc")
    ax[1].set_ylabel("$W_p/W_N$")
    ax[1].set_title("V_SO = 0.4082: opens exactly at $2E_{min}$, diverges toward Tc")
    ax[1].legend(fontsize=7.5)
    fig.suptitle("Pair creation from the HHM Hamiltonian: one phonon breaks one "
                 "$(\\mathbf{k}s,-\\mathbf{k}s)$ pair into two Bogoliubov quasiparticles")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(f"{FIG}/pair_creation_channel.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    import sys
    os.makedirs(FIG, exist_ok=True)
    which = sys.argv[1:] or ["all"]
    jobs = {"sweep": fig_temperature_sweep, "variants": fig_coherence_variants,
            "branch": fig_branch_structure, "coh": fig_coherence_factors,
            "sumrule": fig_sum_rule_thresholds, "lattice": fig_lattice_crosscheck,
            "vso": fig_vso_scan, "pair": fig_pair_channel}
    for k, fn in jobs.items():
        if "all" in which or k in which:
            fn()
            print("figure:", k, flush=True)
