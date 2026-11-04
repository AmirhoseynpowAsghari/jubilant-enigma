#!/usr/bin/env python3
"""Publication-style two-panel plot of strict-admissible W/W_N temperature-specific peaks.

Reads the already extracted per-temperature maxima CSV. It does not access the pointwise
upper-envelope CSV, re-solve the physics, interpolate, smooth, or mix in relaxed roots.
Both panels are built from the exact same strict rows written to the new output CSV.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_INPUT = HERE / "analysis" / "wwn_temperature_maxima" / "temperature_specific_maxima.csv"
DEFAULT_OUTPUT = HERE / "analysis" / "wwn_temperature_maxima"
SOC_VALUES = (0.0, 0.2, 0.4, 0.6, 0.8)
TEMPERATURES = (0.3, 0.5, 0.7, 0.9, 0.999)
MARKERS = {0.3: "o", 0.5: "s", 0.7: "^", 0.9: "D", 0.999: "P"}
GRID_STEP = 0.5
BOUNDARY_TOL = 1e-9

FIELDS = ["V_SO_over_t", "T_over_Tc", "U_peak_over_t", "W_WN_peak", "policy",
          "at_admissible_domain_boundary", "admissible_U_min_over_t",
          "admissible_U_max_over_t", "boundary_side", "valid_U_samples"]


def finite_float(value):
    if value is None or value == "":
        return None
    x = float(value)
    return x if math.isfinite(x) else None


def load_strict_peak_rows(input_csv: Path) -> list[dict]:
    """Load exactly one strict root per SOC/temperature; retain missing cells as None."""
    with input_csv.open(newline="", encoding="utf-8") as f:
        source = list(csv.DictReader(f))
    if not source:
        raise ValueError(f"No extracted maxima in {input_csv}")
    strict = [r for r in source if r.get("policy") == "strict"]
    # Duplicate or non-strict records must not silently contaminate a primary-policy plot.
    keyed = {}
    for r in strict:
        key = (float(r["V_SO"]), float(r["frac"]))
        if key in keyed:
            raise ValueError(f"Duplicate strict peak record for {key}")
        keyed[key] = r

    rows = []
    for vso in SOC_VALUES:
        for frac in TEMPERATURES:
            key = (vso, frac)
            r = keyed.get(key)
            if r is None:
                rows.append(dict(V_SO_over_t=vso, T_over_Tc=frac, U_peak_over_t=None,
                                 W_WN_peak=None, policy="strict",
                                 at_admissible_domain_boundary=None,
                                 admissible_U_min_over_t=None, admissible_U_max_over_t=None,
                                 boundary_side=None, valid_U_samples=None))
                continue
            u = finite_float(r.get("U_peak"))
            w = finite_float(r.get("W_WN_peak"))
            umin = finite_float(r.get("U_min"))
            umax = finite_float(r.get("U_max"))
            nvalid = int(r["n_valid_U"]) if r.get("n_valid_U") else None
            side = []
            if u is not None and umin is not None and abs(u - umin) <= BOUNDARY_TOL:
                side.append("lower")
            if u is not None and umax is not None and abs(u - umax) <= BOUNDARY_TOL:
                side.append("upper")
            boundary = bool(side) if u is not None and umin is not None and umax is not None else None
            if u is not None:
                if umin is None or umax is None or nvalid is None:
                    raise ValueError(f"Peak at {key} is missing its admissible-U range")
                # Check that the strict admissible samples form the documented 0.5t grid
                # interval. This makes the boundary flag meaningful rather than inferred
                # from a possibly disconnected set of U samples.
                expected = round((umax - umin) / GRID_STEP) + 1
                if abs((umax - umin) / GRID_STEP - round((umax - umin) / GRID_STEP)) > 1e-8:
                    raise ValueError(f"Strict U endpoints are off the {GRID_STEP}t grid: {vso}, {frac}")
                if nvalid != expected:
                    raise ValueError(f"Strict admissible U samples are not contiguous: "
                                     f"VSO={vso}, frac={frac}, count={nvalid}, expected={expected}")
                if w is None:
                    raise ValueError(f"Finite U peak lacks a finite ratio: VSO={vso}, frac={frac}")
            rows.append(dict(V_SO_over_t=vso, T_over_Tc=frac, U_peak_over_t=u,
                             W_WN_peak=w, policy="strict",
                             at_admissible_domain_boundary=boundary,
                             admissible_U_min_over_t=umin, admissible_U_max_over_t=umax,
                             boundary_side=(";".join(side) if side else ("interior" if boundary is False else None)),
                             valid_U_samples=nvalid))
    return rows


def write_peak_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_two_panel_figure(rows: list[dict], png_path: Path, pdf_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # Tab10 is color-vision friendly enough for five distinct, consistently reused SOC hues.
    cmap = plt.get_cmap("tab10")
    colors = {v: cmap(i) for i, v in enumerate(SOC_VALUES)}
    t_label = {0.3: "0.3", 0.5: "0.5", 0.7: "0.7", 0.9: "0.9", 0.999: "0.999"}
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 11, "axes.titlesize": 11,
        "legend.fontsize": 9, "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.4, 5.4), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.79, bottom=0.14, top=0.82, wspace=0.30)

    # The same row list and exact values are used in both panels. Missing maxima stay gaps.
    for r in rows:
        u, w = r["U_peak_over_t"], r["W_WN_peak"]
        if u is None or w is None:
            continue
        boundary = r["at_admissible_domain_boundary"] is True
        ax_a.scatter(u, w, s=68, marker=MARKERS[r["T_over_Tc"]],
                     facecolors="none" if boundary else colors[r["V_SO_over_t"]],
                     edgecolors="black" if boundary else colors[r["V_SO_over_t"]],
                     linewidths=1.8 if boundary else 0.9, zorder=3)

    # Panel B uses those same strict peak rows, merely ordered along temperature.
    for vso in SOC_VALUES:
        q = [r for r in rows if r["V_SO_over_t"] == vso]
        q.sort(key=lambda r: r["T_over_Tc"])
        x = [r["T_over_Tc"] if r["W_WN_peak"] is not None else float("nan") for r in q]
        y = [r["W_WN_peak"] if r["W_WN_peak"] is not None else float("nan") for r in q]
        ax_b.plot(x, y, color=colors[vso], linewidth=1.8, zorder=2)
        for r in q:
            if r["W_WN_peak"] is None:
                continue
            boundary = r["at_admissible_domain_boundary"] is True
            ax_b.scatter(r["T_over_Tc"], r["W_WN_peak"], s=58,
                         marker=MARKERS[r["T_over_Tc"]],
                         facecolors="none" if boundary else colors[vso],
                         edgecolors="black" if boundary else colors[vso],
                         linewidths=1.8 if boundary else 0.9, zorder=3)

    ax_a.set_xlabel(r"$U_{\mathrm{peak}}/t$")
    ax_a.set_ylabel(r"$(W/W_N)_{\mathrm{peak}}$")
    ax_a.set_title("A  Temperature-specific peaks vs. maximizing U", loc="left", pad=9)
    ax_a.set_xlim(48, 91)
    ax_a.set_ylim(-0.01, 1.02)

    ax_b.set_xlabel(r"$T/T_c$")
    ax_b.set_title("B  Peak magnitude vs. temperature", loc="left", pad=9)
    ax_b.set_xlim(0.27, 1.025)
    ax_b.set_ylim(-0.01, 1.02)

    for ax in (ax_a, ax_b):
        ax.grid(True, color="0.82", linewidth=0.65, alpha=0.7)
        ax.set_axisbelow(True)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    soc_handles = [Line2D([0], [0], color=colors[v], linewidth=2.2,
                          label=rf"$V_{{SO}}/t={v:g}$") for v in SOC_VALUES]
    temp_handles = [Line2D([0], [0], color="0.25", marker=MARKERS[f], linestyle="none",
                           markersize=7, label=rf"$T/T_c={t_label[f]}$") for f in TEMPERATURES]
    boundary_handle = Line2D([0], [0], color="black", marker="o", markerfacecolor="none",
                             linestyle="none", label="admissible-boundary peak")
    fig.legend(handles=soc_handles, loc="upper center", bbox_to_anchor=(0.43, 0.985),
               ncol=5, frameon=False, title=r"SOC (same colors in both panels)",
               title_fontsize=9)
    fig.legend(handles=temp_handles, loc="center left", bbox_to_anchor=(0.81, 0.57),
               ncol=1, frameon=False, title=r"temperature marker", title_fontsize=9)
    if any(r["at_admissible_domain_boundary"] is True for r in rows):
        fig.legend(handles=[boundary_handle], loc="center left", bbox_to_anchor=(0.81, 0.24),
                   frameon=False, fontsize=8)
    fig.text(0.43, 0.035,
             "Strict-admissible roots only; U sampled at 0.5t. No interpolation or smoothing.",
             ha="center", fontsize=8.5, color="0.3")

    fig.savefig(png_path, dpi=360, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = load_strict_peak_rows(args.input)
    csv_path = args.output_dir / "temperature_peak_analysis_two_panel.csv"
    png_path = args.output_dir / "temperature_peak_analysis_two_panel.png"
    pdf_path = args.output_dir / "temperature_peak_analysis_two_panel.pdf"
    write_peak_csv(csv_path, rows)
    make_two_panel_figure(rows, png_path, pdf_path)
    valid = [r for r in rows if r["U_peak_over_t"] is not None]
    boundary_count = sum(r["at_admissible_domain_boundary"] is True for r in valid)
    print(f"Strict peak rows: {len(valid)}/{len(rows)}; admissible-boundary peaks: {boundary_count}")
    print(f"Wrote {csv_path}\nWrote {png_path}\nWrote {pdf_path}")


if __name__ == "__main__":
    main()
