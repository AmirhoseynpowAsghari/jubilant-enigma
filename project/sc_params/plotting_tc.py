"""
plotting_tc.py -- stage-1 figure gallery (run: python plotting_tc.py).
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io_utils_tc as io

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3, "grid.linewidth": 0.5,
                     "font.size": 9, "legend.fontsize": 8, "figure.facecolor": "white"})


def phi(f):
    return (1.0 - f ** 3.434) ** 0.553 if f < 1 else 0.0


def fig_tc_vs_vso():
    m = io.load("tc_vs_vso.json")["results"]
    u = io.load("u_scans_stage1.json")
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    ax[0].plot(m["V_SO"], m["Tc"], "k-", lw=2.0, label="U=95t (master, 50 pts)")
    ax[1].plot(m["V_SO"], m["Deff_T0"], "k-", lw=2.0)
    for key, c in (("95", "0.0"), ("125", "tab:red"), ("-75", "tab:green")):
        if key == "95":
            c = "tab:blue"
        ax[0].plot(u[key]["V_SO"], u[key]["Tc"], "-o", ms=3, lw=1.0, color=c, label=f"U={key}t")
        ax[1].plot(u[key]["V_SO"], u[key]["Deff_0"], "-o", ms=3, lw=1.0, color=c)
        mu_inside = np.asarray(u[key].get("mu_inside_band", [True] * len(u[key]["V_SO"])), dtype=bool)
        outband = ~mu_inside
        if outband.any():
            ax[0].plot(np.asarray(u[key]["V_SO"])[outband], np.asarray(u[key]["Tc"])[outband],
                       "^", ms=4, mfc="none", mec="0.25", mew=0.8, ls="none",
                       label="mu outside band (relaxed)")
            ax[1].plot(np.asarray(u[key]["V_SO"])[outband], np.asarray(u[key]["Deff_0"])[outband],
                       "^", ms=4, mfc="none", mec="0.25", mew=0.8, ls="none")
    ax[0].set_ylabel("Tc / t")
    ax[1].set_ylabel("$\\Delta_{eff}(0)$ / t")
    for a in ax:
        a.set_xlabel("V_SO [t]")
    ax[0].set_title("critical temperature")
    ax[1].set_title("zero-temperature effective gap")
    ax[0].legend()
    fig.suptitle("Stage-1 V_SO scans (n=1.885, $\\Delta_t$=4.5t); open triangles mark relaxed-policy roots "
                 "with mu outside the bare band (not established stable phases)")
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(f"{FIG}/Tc_vs_VSO.png", dpi=150)
    plt.close(fig)


def fig_gap_vs_T():
    g = io.load("groups_stage1.json")
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.3))
    cols = ["tab:blue", "tab:red", "tab:green", "tab:purple"]
    for grp, panel in (("B", 0), ("C", 1)):
        i = 0
        for p in g:
            if p["set"]["group"] != grp:
                continue
            Tc, Deff = p["Tc"], p["Deff_0"]
            T = np.linspace(0, 1.05 * Tc, 200)
            mu_out = (p.get("mu_window") or {}).get("inside") is False
            status = "; mu outside band" if mu_out else ""
            ax[panel].plot(T, np.array([phi(t / Tc) for t in T]) * Deff, color=cols[i % 4],
                           lw=1.5, label=f"{p['set']['label']} (Tc={Tc:.4f}{status})")
            ax[panel].axvline(Tc, color=cols[i % 4], ls=":", lw=0.7)
            i += 1
        ax[panel].set_xlabel("T [t]")
        ax[panel].set_ylabel("$\\Delta_{eff}(T)$ [t]")
        ax[panel].set_title(f"Group {grp}")
        ax[panel].legend(fontsize=7.5)
    fig.suptitle("Stage-1 gap functions of the paper's parameter sets (repository fit shape)")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{FIG}/gap_vs_T_comparison.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_tc_vs_vso()
    fig_gap_vs_T()
    print("stage-1 figures done")
