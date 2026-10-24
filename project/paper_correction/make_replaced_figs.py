"""
make_paper_figs.py -- regenerated (corrected) figures for arXiv:2607.02812.

* fig_group_{A,B,C,D}.png : the four six-panel composites, same layout/titles as the
  originals, but the two W/W_N panels now show the validated result
  (eta = +1, signed xi xi', scattering channel, hbar*omega = 0.7 t), and the four
  auxiliary panels (mu(T), Delta(T)) are recomputed from the same stage-1 solution.
* fig_scan_{U95,U125,Um75}.png : the three "W/W_N vs spin-orbit coupling" figures.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
# SINGLE SOURCE OF PHYSICS (2026-10-24): the corrected main code, not the vendored copy of an
# earlier stage-2 reconstruction (which still used the 0.02 DOS smoothing and the endpoint-era
# settings).  All W/W_N numbers come from the regenerated JSON; mu(T) uses the same solver.
sys.path.insert(0, os.path.join(HERE, "..", "ultrasonic"))
sys.path.insert(0, os.path.join(HERE, "..", "sc_params"))
import physics_ua as model

G = json.load(open(os.path.join(HERE, "vendor", "results", "groups.json")))
S = json.load(open(os.path.join(HERE, "vendor", "results", "scans.json")))

COL = ["tab:blue", "tab:red", "tab:green", "tab:purple"]
SCANCOL = ["#8B0000", "#E41A1C", "#FF7F00", "#E6AB02"]
plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3, "grid.linewidth": 0.5,
                     "axes.linewidth": 0.8, "font.size": 9, "legend.fontsize": 7.5,
                     "axes.titlesize": 9.5, "figure.facecolor": "white"})

GROUP_TITLES = {
    "A": "Ultrasonic Attenuation — Group A  ($\\Delta_t$ = 0, n = 0.8, V_SO = 0)  [corrected]",
    "B": "Ultrasonic Attenuation — Group B  (U = 75t…125t, $\\Delta_t$ = 4.5, n = 1.8850, V_SO = 0)  [corrected]",
    "C": "Ultrasonic Attenuation — Group C  (U = 75t, n = 1.8850, V_SO = 0.8, $\\Delta_t$ varies)  [corrected]",
    "D": "Ultrasonic Attenuation — Group D  (U = 95t, n = 1.8850, $\\Delta_t$ = 4.5, V_SO varies)  [corrected]",
}


def phi(f):
    """repository gap shape (1-(T/Tc)^p)^q, p=3.434, q=0.553"""
    if f >= 1.0:
        return 0.0
    return (1.0 - f ** 3.434) ** 0.553


NO_SC_NOTE = []
for grp in "ABCD":
    rows = [r for r in G if r["set"]["group"] == grp]
    FS = {"A": (15.45, 14.43), "B": (15.39, 16.29), "C": (15.45, 12.37), "D": (15.95, 10.95)}
    fig, ax = plt.subplots(2, 3, figsize=FS[grp])
    for i, r in enumerate(rows):
        s, st1 = r["set"], r["stage1"]
        if r.get("no_superconductivity") or r.get("a_ratio") is None:
            NO_SC_NOTE.append(f"{s['label']} (group {grp}): no superconducting state in the physical "
                              f"window (chemical potential outside the band); not plotted")
            continue
        fr = np.array(r["fracs"])
        rat = np.array(r["a_ratio"])
        Deff0, Tc, VSO = st1["Deff_0"], st1["Tc"], st1["V_SO"]
        fa = np.linspace(0.0, 1.4, 141)
        pha = np.array([phi(f) for f in fa])
        mu_cur = np.array([model.solve_mu(max(f * Tc, 1e-9), phi(f) * st1["D0_0"],
                                          phi(f) * st1["DS_0"], VSO, n_tgt=s["n"], Nk=200)
                           for f in fa])
        c = COL[i % 4]
        lab = f"{s['label']}  (Tc = {Tc:.4f})"
        ax[0, 0].plot(np.array([phi(f) for f in fr]), rat, color=c, lw=1.6, label=lab)
        ax[0, 1].plot(fr, rat, color=c, lw=1.6, label=lab)
        ax[0, 2].plot(fa, pha, color=c, lw=1.6, label=s["label"])
        ax[1, 0].plot(fa * Tc, mu_cur, color=c, lw=1.4, label=lab)
        ax[1, 0].axvline(Tc, color=c, ls=":", lw=0.7)
        ax[1, 1].plot(fa, mu_cur, color=c, lw=1.4, label=lab)
        ax[1, 1].axvline(1.0, color=c, ls=":", lw=0.7)
        ax[1, 2].plot(fa * Tc, pha * Deff0, color=c, lw=1.4,
                      label=f"{s['label']}  $\\Delta_{{eff}}$ = {Deff0:.4f}")
        ax[1, 2].axvline(Tc, color=c, ls=":", lw=0.7)
    ax[0, 0].set_title("W/W$_N$ vs Normalised Gap")
    ax[0, 1].set_title("W/W$_N$ vs Normalised Temperature")
    ax[0, 2].set_title("Normalised Gap vs T/Tc\n(universal BCS curve)")
    ax[1, 0].set_title("Chemical Potential vs T\n(dotted = each Tc)")
    ax[1, 1].set_title("Chemical Potential vs T/Tc")
    ax[1, 2].set_title("Gap $\\Delta_{eff}$(T) vs T")
    ax[0, 0].set_ylabel("W/W$_N$")
    ax[0, 0].set_xlabel("$\\Delta$(T) / $\\Delta_{eff}$(0)")
    ax[0, 1].set_xlabel("T / Tc")
    ax[0, 2].set_xlabel("T / Tc")
    ax[0, 2].set_ylabel("$\\Delta$(T) / $\\Delta_{eff}$(0)")
    ax[1, 0].set_xlabel("Temperature T  [units of t]")
    ax[1, 0].set_ylabel("$\\mu$(T)  [units of t]")
    ax[1, 1].set_xlabel("T / Tc")
    ax[1, 1].set_ylabel("$\\mu$(T)  [units of t]")
    ax[1, 2].set_xlabel("Temperature T  [units of t]")
    ax[1, 2].set_ylabel("$\\Delta$(T)  [units of t]")
    for a in (ax[0, 0], ax[0, 1]):
        a.set_ylim(0.0, 1.05)
        a.axhline(1.0, color="0.5", lw=0.7)
    ax[0, 1].axvline(1.0, color="0.4", ls="--", lw=0.8)
    ax[0, 2].axvline(1.0, color="0.4", ls="--", lw=0.8)
    ax[0, 0].legend(loc="lower left", title="set")
    ax[0, 1].legend(loc="lower left", title="set")
    ax[0, 2].legend(loc="lower left", title="set" if grp == "A" else None)
    ax[1, 0].legend(loc="best", title="set")
    ax[1, 1].legend(loc="best", title="set")
    ax[1, 2].legend(loc="lower left", fontsize=6.8)
    fig.suptitle(GROUP_TITLES[grp], y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(os.path.join(HERE, "vendor", "figs", f"fig_group_{grp}.png"), dpi=100)
    plt.close(fig)
    print("group", grp, "done", flush=True)

# ---------------------------------------------------------------- V_SO scans
SCAN_FS = {"95": (12.0, 4.07), "125": (12.0, 5.95), "-75": (12.0, 6.1)}
for key, title in (("95", "U=95t"), ("125", "U=125t"), ("-75", "U=-75t")):
    rec = S[key]
    if all(x is None for f in rec["ratio"].values() for x in f):
        NO_SC_NOTE.append(f"U={key}t: no superconducting state at any V_SO in the physical window")
    V = np.array(rec["V_SO"])
    fig, a = plt.subplots(figsize=SCAN_FS[key])
    for j, f in enumerate(("0.3", "0.5", "0.7", "0.9")):
        y = np.array([x if x is not None else np.nan for x in rec["ratio"][f]], float)
        a.plot(V, y, color=SCANCOL[j], lw=1.4, marker="o", ms=3.5, label=f"T/$T_c$ = {f}")
    a.set_ylim(0.0, 1.05)
    a.axhline(1.0, color="0.5", lw=0.7)
    a.set_xlabel("Spin–orbit coupling  V_SO  [units of t]")
    a.set_ylabel("W/W$_N$")
    a.set_title("Ultrasonic Attenuation vs Spin–Orbit Coupling  [corrected]\n"
                f"{title}, n=1.885, omega=0.7t, $\\Delta_t$=4.5t")
    a.legend(loc="lower right", title="fixed T/Tc")
    fig.tight_layout()
    tag = {"95": "U95", "125": "U125", "-75": "Um75"}[key]
    fig.savefig(os.path.join(HERE, "vendor", "figs", f"fig_scan_{tag}.png"), dpi=100)
    plt.close(fig)
    print("scan", key, "done", flush=True)

with open(os.path.join(HERE, "vendor", "figs", "NO_SC_NOTES.txt"), "w") as fh:
    fh.write("\n".join(NO_SC_NOTE) + "\n")
print("notes:", *NO_SC_NOTE, sep="\n  ")
