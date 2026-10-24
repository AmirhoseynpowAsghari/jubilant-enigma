"""
make_param_figures.py -- publication figures for the (U, V_SO, Delta_t) peak study.
Reads results/param_features.json (+ param_features_peaks.json) and writes
  figures/param_scan_U.png, param_scan_VSO.png, param_scan_Dt.png  (2x2 panels each)
  figures/param_scan_diagnostics.png
Panels: (a) primary channel a_signed at fixed T/Tc=0.5; (b) a_signed at fixed T=0.1t;
(c) physical TOTAL (scattering+pair) at T=0.1t with pair-open shading and peak stars;
(d) the same total at T=0.2t.  Diagnostics: Tc (dome), E_min(T=0.1) vs hbar*omega/2
(pair threshold), Case-II channel (monotone -- reported honestly), and the primary
channel at fixed T/Tc=0.9.
"""
from __future__ import annotations

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

DATA = json.load(open("results/param_features.json"))
PEAKS = json.load(open("results/param_features_peaks.json"))

COL = {"a": "tab:blue", "tot": "tab:red", "c2": "tab:green"}


def fam(r, which):
    if which == "U":
        return (round(r["VSO"], 4), round(r["Dt"], 4)), r["U"]
    if which == "VSO":
        return (round(r["U"], 4), round(r["Dt"], 4)), r["VSO"]
    return (round(r["U"], 4), round(r["VSO"], 4)), r["Dt"]


def label(key, which):
    if which == "U":
        return f"$V_{{SO}}$={key[0]:.4g}, $\\Delta_t$={key[1]:.4g}"
    if which == "VSO":
        return f"$U$={key[0]:.4g}t, $\\Delta_t$={key[1]:.4g}"
    return f"$U$={key[0]:.4g}t, $V_{{SO}}$={key[1]:.4g}"


ALLOWED = {
    "U":   {(0.0, 4.5), (0.4082, 4.5)},
    "VSO": {(95.0, 4.5), (125.0, 4.5)},
    "Dt":  {(95.0, 0.0), (95.0, 0.4082), (125.0, 0.0), (125.0, 0.4082)},
}


def collect(which):
    fams = {}
    for r in DATA:
        k, x = fam(r, which)
        if k not in ALLOWED[which]:
            continue
        d = fams.setdefault(k, {"x": [], "r": {}})
        d["x"].append(x)
        vals = {"a5": r["fracs"]["0.5"]["a"], "a9": r["fracs"]["0.9"]["a"]}
        nan = float("nan")
        for tk in ("0.10", "0.20"):
            t = r["temps"][tk]
            # F-06: normal-state points (sc=0) carry no SC value -> NaN gap (never plotted as data)
            def _v(key):
                return nan if t.get(key) is None else float(t[key])
            vals.update({f"a{tk}": _v("a"), f"tot{tk}": _v("tot"), f"c2{tk}": _v("c2"),
                         f"po{tk}": t["pair_open"], f"em{tk}": t["Emin"],
                         f"spo{tk}": bool(t["pair_open"] and t["sc"])})
        vals["Tc"] = r["Tc"]
        for kk, v in vals.items():
            d["r"].setdefault(kk, []).append(v)
    for k in fams:
        o = np.argsort(fams[k]["x"])
        fams[k]["x"] = np.asarray(fams[k]["x"])[o]
        for kk in fams[k]["r"]:
            fams[k]["r"][kk] = np.asarray(fams[k]["r"][kk])[o]
    return fams


def mono_tag(fams, chan):
    """Data-driven label: monotone only if EVERY family curve is monotone (NaN gaps skipped)."""
    flags = []
    for f in fams.values():
        y = np.asarray(f["r"][chan], float)
        y = y[np.isfinite(y)]
        if len(y) >= 2:
            d = np.diff(y)
            flags.append(bool(np.all(d >= -1e-12) or np.all(d <= 1e-12)))
    if not flags:
        return "[no SC points]"
    return "[monotone]" if all(flags) else f"[NON-monotone in {sum(not x for x in flags)} of {len(flags)} families]"


def panel(ax, fams, which, chan, title, peaks_for, shade=None):
    cols = plt.get_cmap("tab10").colors
    ymax = max(np.nanmax(f["r"][chan]) for f in fams.values()) * 1.15
    for i, k in enumerate(sorted(fams)):
        x, y = fams[k]["x"], fams[k]["r"][chan]
        if shade is not None:
            ax.fill_between(x, 0, ymax, where=fams[k]["r"][shade], color="0.88", zorder=0)
        ax.plot(x, y, "-o", ms=3, lw=1.4, color=cols[i % 10], label=label(k, which))
    for p in PEAKS:
        if p["scan"] == which and tuple(p["key"]) in fams and p["chan"] == peaks_for:
            ax.plot(p["x_ref"], p["y_ref"], "*", ms=14, color="k", zorder=5)
            ax.annotate(f"{p['x_ref']:.3g}, {p['y_ref']:.2f}",
                        (p["x_ref"], p["y_ref"]), textcoords="offset points",
                        xytext=(6, 6), fontsize=7)
    ax.set_title(title, fontsize=9)
    ax.set_ylim(0, ymax)
    ax.grid(alpha=0.3)


XL = {"U": "$U/t$", "VSO": "$V_{SO}/t$", "Dt": "$\\Delta_t/t$"}


def figure(which, path):
    fams = collect(which)
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0))
    panel(axes[0, 0], fams, which, "a5",
          f"(a) primary: $\\eta=+1$ signed, scattering only, $T/T_c=0.5$  {mono_tag(fams, 'a5')}", None)
    panel(axes[0, 1], fams, which, "a0.10",
          f"(b) primary channel at fixed $T=0.1t$  {mono_tag(fams, 'a0.10')}", None)
    ax = axes[1, 0]
    panel(ax, fams, which, "tot0.10",
          "(c) physical TOTAL (sc+pair), $T=0.1t$; grey = pair-open $2E_{min}<\\hbar\\omega$",
          "total (sc+pair), T=0.1t", shade="spo0.10")
    ax = axes[1, 1]
    panel(ax, fams, which, "tot0.20",
          "(d) physical TOTAL (sc+pair), $T=0.2t$", "total (sc+pair), T=0.2t",
          shade="po0.20")
    for a in axes.ravel():
        a.set_xlabel(XL[which])
        a.set_ylabel("$W/W_N$")
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print("wrote", path)


def diagnostics():
    famsU = collect("U")
    famsD = collect("Dt")
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.0))
    ax = axes[0, 0]
    for k in sorted(famsU):
        ax.plot(famsU[k]["x"], famsU[k]["r"]["Tc"], "-o", ms=3, label=label(k, "U"))
    ax.axhline(0.1, color="k", ls=":", lw=1)
    ax.axhline(0.2, color="k", ls="--", lw=1)
    ax.set_xlabel(XL["U"]); ax.set_ylabel("$T_c/t$")
    ax.set_title("(a) SC dome: $T_c(U)$; lines $T=0.1t, 0.2t$", fontsize=9); ax.grid(alpha=0.3)
    ax.legend(fontsize=7.5)
    ax = axes[0, 1]
    for k in sorted(famsU):
        ax.plot(famsU[k]["x"], famsU[k]["r"]["em0.10"], "-o", ms=3, label=label(k, "U"))
    ax.axhline(0.35, color="tab:red", ls="--", lw=1.2, label="$\\hbar\\omega/2$")
    ax.set_xlabel(XL["U"]); ax.set_ylabel("$E_{min}(T=0.1t)/t$")
    ax.set_title("(b) pair threshold $2E_{min}=\\hbar\\omega$ vs $U$", fontsize=9)
    ax.grid(alpha=0.3); ax.legend(fontsize=7.5)
    ax = axes[1, 0]
    for k in sorted(famsU):
        ax.plot(famsU[k]["x"], famsU[k]["r"]["c20.10"], "-o", ms=3, color="tab:green",
                label=label(k, "U"))
    ax.set_xlabel(XL["U"]); ax.set_ylabel("$W/W_N$ (Case II)")
    ax.set_title("(c) Case-II channel, $T=0.1t$: monotone (no HS resonance in range)",
                 fontsize=9); ax.grid(alpha=0.3); ax.legend(fontsize=7.5)
    ax = axes[1, 1]
    for k in sorted(famsD):
        ax.plot(famsD[k]["x"], famsD[k]["r"]["a9"], "-o", ms=3, label=label(k, "Dt"))
    ax.set_xlabel(XL["Dt"]); ax.set_ylabel("$W/W_N$")
    ax.set_title("(d) primary channel, $T/T_c=0.9$: monotone", fontsize=9)
    ax.grid(alpha=0.3); ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig("figures/param_scan_diagnostics.png", dpi=150)
    plt.close(fig)
    print("wrote figures/param_scan_diagnostics.png")


if __name__ == "__main__":
    figure("U", "figures/param_scan_U.png")
    figure("VSO", "figures/param_scan_VSO.png")
    figure("Dt", "figures/param_scan_Dt.png")
    diagnostics()
