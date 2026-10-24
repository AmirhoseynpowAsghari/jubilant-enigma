"""
plotting_ua_u.py -- U-family figures (groups A and B are U scans).
"""
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import io_utils_ua_u as io

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
plt.rcParams.update({"axes.grid": True, "grid.alpha": 0.3, "font.size": 9,
                     "legend.fontsize": 8, "figure.facecolor": "white"})


def fig_W_vs_U():
    g = io.load("groups_WWN.json")       # U-tree run_groups output (scanning format, groups A and B)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    for p in g:
        s = p["set"]
        lim = p["ratios"]["a_signed"][-1]
        ax.plot(s["U"], lim, "o", ms=7, color="tab:blue" if s["group"] == "A" else "tab:red")
        ax.annotate(f"{s['label']}\n{lim:.3f}", (s["U"], lim), textcoords="offset points",
                    xytext=(6, 6), fontsize=7.5)
    ax.axhline(1.0, color="0.5", lw=0.7)
    ax.set_xlabel("U / t   (attractive group A left, repulsive group B right)")
    ax.set_ylabel("$\\lim_{T\\to T_c^-} W/W_N$")
    ax.set_title("Corrected saturation value vs interaction strength")
    fig.tight_layout()
    fig.savefig(f"{FIG}/W_WN_vs_U.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_W_vs_U()
    print("ultrasonic_U figures done")
