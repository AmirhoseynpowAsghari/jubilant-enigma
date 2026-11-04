#!/usr/bin/env python3
"""Extract temperature-specific W/W_N maxima from the existing W(U) family.

Primary input is the exact relaxed dataset underlying figures/WWN_vs_U_families_VSO.png.
The strict-policy family is additionally included as a policy sensitivity check. This is a
postprocessor only: it does not recompute or overwrite the existing family data or plots.

Two unambiguous reductions are exported:
1) For each (V_SO, T/Tc), maximize over sampled U and retain (U_peak, W_peak).
2) For each (V_SO, U), maximize over the available sampled temperatures (an upper envelope).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
OUT_DEFAULT = HERE / "analysis" / "wwn_temperature_maxima"
VSO_VALUES = (0.0, 0.2, 0.4, 0.6, 0.8)
FRACS = (0.3, 0.5, 0.7, 0.9, 0.999)
FRAC_LABELS = {0.3: "0.3", 0.5: "0.5", 0.7: "0.7", 0.9: "0.9", 0.999: "0.999"}
MARKERS = {0.3: "o", 0.5: "s", 0.7: "^", 0.9: "D", 0.999: "P"}


def _finite(v):
    return isinstance(v, (int, float)) and math.isfinite(v)


def _load_sources(relaxed_path: Path, strict_path: Path):
    relaxed = json.loads(relaxed_path.read_text(encoding="utf-8"))
    strict_root = json.loads(strict_path.read_text(encoding="utf-8"))
    strict = strict_root["uvso_families"]
    return {"relaxed": relaxed, "strict": strict}


def _normalise_records(records, policy):
    """Project family records onto common per-temperature rows with policy diagnostics."""
    out = []
    for r in records:
        try:
            vso, U = float(r["VSO"]), float(r["U"])
            if not any(abs(vso - v) < 1e-10 for v in VSO_VALUES):
                continue
            if policy == "strict":
                if r.get("policy") != "strict":
                    continue
                ratios = r.get("ratio", [])
                band = r.get("mu_window") or {}
                stage2_mu = r.get("stage2_mu", [])
                stage2_inside = r.get("stage2_inside", [])
            else:
                ratios = r.get("ratio", [])
                band = r.get("mu_window") or {}
                stage2_mu = []
                stage2_inside = []
            for j, frac in enumerate(FRACS):
                value = ratios[j] if j < len(ratios) else None
                if not _finite(value):
                    continue
                mu_T = stage2_mu[j] if j < len(stage2_mu) else None
                inside_T = stage2_inside[j] if j < len(stage2_inside) else None
                # A finite strict value should have an accepted finite-T mu; fail loudly if
                # the source schema ever starts including rejected roots in its ratio array.
                if policy == "strict" and inside_T is not True:
                    raise ValueError(f"strict finite ratio lacks accepted mu: U={U}, VSO={vso}, f={frac}")
                out.append(dict(
                    policy=policy, U=U, V_SO=vso, frac=frac, ratio=float(value),
                    Tc=float(r.get("Tc", 0.0)), D0_0=float(r.get("D0_0", 0.0)),
                    DS_0=float(r.get("DS_0", 0.0)), Deff_0=float(r.get("Deff_0", 0.0)),
                    mu0=float(r["mu_0"]) if _finite(r.get("mu_0")) else None,
                    mu0_inside=band.get("inside"), mu_T=float(mu_T) if _finite(mu_T) else None,
                    mu_T_inside=inside_T,
                ))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Malformed {policy} source record: {r}") from exc
    return out


def extract(records):
    peaks, envelope = [], []
    for policy in ("strict", "relaxed"):
        data = [r for r in records if r["policy"] == policy]
        for vso in VSO_VALUES:
            for frac in FRACS:
                group = [r for r in data if r["V_SO"] == vso and r["frac"] == frac]
                if not group:
                    peaks.append(dict(policy=policy, V_SO=vso, frac=frac, U_peak=None,
                                      W_WN_peak=None, n_valid_U=0, U_min=None, U_max=None,
                                      mu0_inside_peak=None, mu_T_inside_peak=None))
                    continue
                # np.argmax is deterministic and returns the first/smallest-U maximum on ties.
                group.sort(key=lambda r: r["U"])
                winner = max(group, key=lambda r: r["ratio"])
                peaks.append(dict(policy=policy, V_SO=vso, frac=frac, U_peak=winner["U"],
                                  W_WN_peak=winner["ratio"], n_valid_U=len(group),
                                  U_min=group[0]["U"], U_max=group[-1]["U"],
                                  mu0_inside_peak=winner["mu0_inside"],
                                  mu_T_inside_peak=winner["mu_T_inside"],
                                  Tc_peak=winner["Tc"], D0_0_peak=winner["D0_0"],
                                  DS_0_peak=winner["DS_0"], Deff_0_peak=winner["Deff_0"],
                                  mu0_peak=winner["mu0"], mu_T_peak=winner["mu_T"]))
            u_values = sorted({r["U"] for r in data if r["V_SO"] == vso})
            for U in u_values:
                group = [r for r in data if r["V_SO"] == vso and r["U"] == U]
                if not group:
                    continue
                winner = max(group, key=lambda r: r["ratio"])
                envelope.append(dict(policy=policy, V_SO=vso, U=U,
                                     W_WN_max_over_temperature=winner["ratio"],
                                     frac_at_max=winner["frac"],
                                     n_valid_temperatures=len(group)))
    return peaks, envelope


def save_csv(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="", encoding="utf-8") as f:
        wr = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)


def make_plots(peaks, envelope, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = {v: plt.get_cmap("viridis")(0.05 + 0.9*i/(len(VSO_VALUES)-1))
              for i, v in enumerate(VSO_VALUES)}
    styles = {"strict": "-", "relaxed": "--"}

    # Main request: one sampled peak (U_peak, W_peak) for each temperature curve and V_SO.
    fig, ax = plt.subplots(figsize=(9.4, 6.2))
    for policy in ("strict", "relaxed"):
        for vso in VSO_VALUES:
            q = [r for r in peaks if r["policy"] == policy and r["V_SO"] == vso
                 and r.get("U_peak") is not None]
            q.sort(key=lambda r: r["frac"])
            ax.plot([r["U_peak"] for r in q], [r["W_WN_peak"] for r in q],
                    color=colors[vso], linestyle=styles[policy], linewidth=1.5,
                    marker="o" if policy == "strict" else "X", markersize=6,
                    label=f"$V_{{SO}}/t={vso:g}$, {policy}")
            for r in q:
                ax.scatter(r["U_peak"], r["W_WN_peak"], marker=MARKERS[r["frac"]],
                           s=52, facecolors="none" if policy == "strict" else colors[vso],
                           edgecolors=colors[vso], linewidths=1.0, zorder=4)
    ax.set(xlabel=r"$U_{peak}/t$ (sampled maximum on $\Delta U=0.5t$ grid)",
           ylabel=r"$(W/W_N)_{peak}$",
           title="Peak W/W_N for each temperature curve, versus its maximizing U")
    ax.grid(alpha=0.25)
    policy_soc = [Line2D([0], [0], color=colors[v], linestyle=styles[p], linewidth=1.7,
                         label=f"$V_{{SO}}/t={v:g}$, {p}")
                  for p in ("strict", "relaxed") for v in VSO_VALUES]
    temp_markers = [Line2D([0], [0], color="0.2", marker=MARKERS[f], linestyle="none",
                           markerfacecolor="0.6", label=rf"$T/T_c={FRAC_LABELS[f]}$")
                    for f in FRACS]
    leg1 = ax.legend(handles=policy_soc, fontsize=8, ncol=2, loc="upper left",
                     title="SOC / policy")
    ax.add_artist(leg1)
    ax.legend(handles=temp_markers, fontsize=8, ncol=1, loc="lower right", title="temperature")
    fig.tight_layout()
    fig.savefig(outdir / "temperature_specific_peak_values_vs_U.png", dpi=180)
    plt.close(fig)

    # Alternative useful summary: maximum over T pointwise at each U, one curve per SOC.
    fig, ax = plt.subplots(figsize=(9.4, 6.2))
    for policy in ("strict", "relaxed"):
        for vso in VSO_VALUES:
            q = [r for r in envelope if r["policy"] == policy and r["V_SO"] == vso]
            q.sort(key=lambda r: r["U"])
            ax.plot([r["U"] for r in q], [r["W_WN_max_over_temperature"] for r in q],
                    color=colors[vso], linestyle=styles[policy], linewidth=1.6,
                    label=f"$V_{{SO}}/t={vso:g}$, {policy}")
    ax.set(xlabel=r"$U/t$", ylabel=r"$\max_{T/T_c}(W/W_N)$",
           title="Pointwise upper envelope over sampled temperatures")
    ax.grid(alpha=0.25); ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(outdir / "temperature_envelope_max_vs_U.png", dpi=180)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--relaxed", type=Path, default=RESULTS / "uvso_families.json")
    ap.add_argument("--strict", type=Path,
                    default=RESULTS / "strict_policy" / "strict_family_data.json")
    ap.add_argument("--output", type=Path, default=OUT_DEFAULT)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    raw = _load_sources(args.relaxed, args.strict)
    records = []
    for policy, data in raw.items():
        records.extend(_normalise_records(data, policy))
    peaks, envelope = extract(records)
    save_csv(args.output / "temperature_specific_maxima.csv", peaks)
    save_csv(args.output / "temperature_upper_envelope.csv", envelope)
    make_plots(peaks, envelope, args.output)
    print(f"Input rows retained: strict={sum(r['policy']=='strict' for r in records)}, "
          f"relaxed={sum(r['policy']=='relaxed' for r in records)}")
    print(f"Wrote {args.output / 'temperature_specific_maxima.csv'}")
    print(f"Wrote {args.output / 'temperature_upper_envelope.csv'}")
    print(f"Wrote {args.output / 'temperature_specific_peak_values_vs_U.png'}")
    print(f"Wrote {args.output / 'temperature_envelope_max_vs_U.png'}")
    print("Peak coordinates (U_peak, W/W_N_peak):")
    for r in peaks:
        if r["U_peak"] is not None:
            print(f"  {r['policy']:7s} VSO={r['V_SO']:.1f} T/Tc={r['frac']:.3f}: "
                  f"U={r['U_peak']:.1f}, W/WN={r['W_WN_peak']:.9g}, "
                  f"mu0_inside={r['mu0_inside_peak']}")


if __name__ == "__main__":
    main()
