"""
make_addendum_figs.py -- regenerates the three addendum figures (A1 W/W_N vs T/Tc for the three
master V_SO values and all variants; A2 branch decomposition; A3 limits) from the REGENERATED
production data.  The earlier versions came from a vendored copy of an older code and are not used.
Outputs: vendor/eta1/fig_WWN_vs_T.png, fig_branch_decomposition.png, fig_limits.png
"""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
S2 = os.path.join(ROOT, "project", "ultrasonic", "results")
OUT = os.path.join(HERE, "vendor", "eta1")
os.makedirs(OUT, exist_ok=True)
VS = [("0.0000", "V_SO=0"), ("0.4082", "V_SO=0.4082"), ("0.8000", "V_SO=0.8")]
VARS = [("a_signed", "(a) signed xi xi' (primary)"), ("b_abs", "(b) |xi||xi'|"),
        ("c_paper", "(c) paper's combination"), ("d_std", "(d) no xi xi' (not a BdG element)")]
COL = {"a_signed": "k", "b_abs": "tab:orange", "c_paper": "tab:green", "d_std": "tab:purple"}
LS = {"0.0000": "-", "0.4082": "--", "0.8000": ":"}

# ---- A1: W/W_N vs T/Tc, all variants, three V_SO
fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), sharey=True)
for ax, (v, lab) in zip(axes, VS):
    d = json.load(open(os.path.join(S2, f"temp_sweep_VSO_{v}.json")))
    f = np.array(d["fracs"], float)
    for var, vlab in VARS:
        y = np.array([np.nan if x is None else x for x in d["ratios"][var]], float)
        ax.plot(f, y, color=COL[var], lw=1.5, label=vlab)
    ax.axhline(1.0, color="0.5", lw=0.7)
    ax.set_title(f"U=95t, Dt=4.5t, n=1.885, {lab}  (Tc={d['Tc']:.4f})", fontsize=9)
    ax.set_xlabel("T / Tc"); ax.grid(alpha=0.3)
axes[0].set_ylabel("W / W_N   (scattering channel, hbar*omega=0.7t)")
axes[0].legend(fontsize=7.5, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_WWN_vs_T.png"), dpi=130); plt.close(fig)

# ---- A3: limits (T/Tc -> 0.999) per variant vs V_SO
fig, ax = plt.subplots(figsize=(7.0, 4.4))
vso = [float(v) for v, _ in VS]
for var, vlab in VARS:
    y = []
    for v, _ in VS:
        d = json.load(open(os.path.join(S2, f"temp_sweep_VSO_{v}.json")))
        y.append(d["ratios"][var][-1])
    ax.plot(vso, y, "o-", color=COL[var], label=vlab)
ax.axhline(1.0, color="0.5", lw=0.7)
ax.set_xlabel("V_SO / t"); ax.set_ylabel("lim_{T->Tc} W/W_N  (T/Tc = 0.999)")
ax.set_title("Limits of the corrected ratio, master set"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_limits.png"), dpi=130); plt.close(fig)

# ---- A2: branch decomposition (sigma, sigma') of the scattering sum, master set, V_SO=0.4082
pr = json.load(open(os.path.join(S2, "production_raw.json")))
rows = pr["master"]["0.4082"]["rows"]
keys = sorted(rows[0]["a_signed_parts"].keys())     # file format, e.g. "1,1"
fr = [r["frac"] for r in rows]
tot = [sum(r["a_signed_parts"].values()) for r in rows]
fig, ax = plt.subplots(figsize=(7.4, 4.4))
for k in keys:
    y = [r["a_signed_parts"].get(k, 0.0) / r["S_N"] for r in rows]
    sa, sb = k.split(",")
    ax.plot(fr, y, "-o", ms=3, label=f"(sigma,sigma') = ({int(sa):+d},{int(sb):+d})")
ax.plot(fr, [t / r["S_N"] for t, r in zip(tot, rows)], "k--", lw=1.2, label="sum = W/W_N (a)")
ax.set_xlabel("T / Tc"); ax.set_ylabel("contribution to W/W_N")
ax.set_title("Branch decomposition, master set V_SO=0.4082 (variant a)"); ax.grid(alpha=0.3); ax.legend(fontsize=7.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_branch_decomposition.png"), dpi=130); plt.close(fig)
print("addendum figures regenerated ->", OUT)
