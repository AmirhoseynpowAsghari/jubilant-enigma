#!/usr/bin/env python3
"""Two-panel publication plot of temperature-specific W/W_N peaks.

The selected policy is plotted alone (strict is the default). Values are read from the
existing sampled maxima CSV; no physics is re-solved, and no interpolation/smoothing occurs.
Panel A scatters each temperature curve's (U_peak, W_peak). Panel B connects the same maxima
against temperature. Relaxed mode is explicitly diagnostic and flags sampled-U endpoints.
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
          "at_admissible_domain_boundary", "at_sampled_U_domain_boundary",
          "U_domain_min_over_t", "U_domain_max_over_t", "boundary_side", "valid_U_samples"]


def finite_float(value):
    if value is None or value == "":
        return None
    x = float(value)
    return x if math.isfinite(x) else None


def load_peak_rows(input_csv: Path, policy: str = "strict") -> list[dict]:
    """Load a single policy, retaining absent maxima as missing rather than replacing by 0."""
    if policy not in ("strict", "relaxed"):
        raise ValueError("policy must be 'strict' or 'relaxed'")
    with input_csv.open(newline="", encoding="utf-8") as f:
        source = list(csv.DictReader(f))
    if not source:
        raise ValueError(f"No extracted maxima in {input_csv}")
    selected = [r for r in source if r.get("policy") == policy]
    keyed = {}
    for r in selected:
        key = (float(r["V_SO"]), float(r["frac"]))
        if key in keyed:
            raise ValueError(f"Duplicate {policy} peak record for {key}")
        keyed[key] = r

    rows = []
    for vso in SOC_VALUES:
        for frac in TEMPERATURES:
            key = (vso, frac)
            r = keyed.get(key)
            if r is None:
                rows.append(dict(V_SO_over_t=vso, T_over_Tc=frac, U_peak_over_t=None,
                                 W_WN_peak=None, policy=policy,
                                 at_admissible_domain_boundary=None,
                                 at_sampled_U_domain_boundary=None,
                                 U_domain_min_over_t=None, U_domain_max_over_t=None,
                                 boundary_side=None, valid_U_samples=None))
                continue
            u, w = finite_float(r.get("U_peak")), finite_float(r.get("W_WN_peak"))
            umin, umax = finite_float(r.get("U_min")), finite_float(r.get("U_max"))
            nvalid = int(r["n_valid_U"]) if r.get("n_valid_U") else None
            if u is None:
                # Preserve a missing/nonfinite peak, including its source coverage if present.
                rows.append(dict(V_SO_over_t=vso, T_over_Tc=frac, U_peak_over_t=None,
                                 W_WN_peak=None, policy=policy,
                                 at_admissible_domain_boundary=None,
                                 at_sampled_U_domain_boundary=None,
                                 U_domain_min_over_t=umin, U_domain_max_over_t=umax,
                                 boundary_side=None, valid_U_samples=nvalid))
                continue
            if umin is None or umax is None or nvalid is None or w is None:
                raise ValueError(f"Incomplete finite {policy} peak record at {key}")
            # Existing family files should have a contiguous 0.5t interval of valid points.
            span = (umax - umin) / GRID_STEP
            if abs(span - round(span)) > 1e-8:
                raise ValueError(f"{policy} U endpoints are off the {GRID_STEP}t grid at {key}")
            expected = round(span) + 1
            if nvalid != expected:
                raise ValueError(f"{policy} sampled U values are not contiguous at {key}: "
                                 f"{nvalid} points, expected {expected}")
            if abs((u - umin) / GRID_STEP - round((u - umin) / GRID_STEP)) > 1e-8:
                raise ValueError(f"U_peak is not on the original {GRID_STEP}t grid at {key}")
            lower = abs(u - umin) <= BOUNDARY_TOL
            upper = abs(u - umax) <= BOUNDARY_TOL
            at_boundary = bool(lower or upper)
            side = "lower;upper" if lower and upper else "lower" if lower else "upper" if upper else "interior"
            # Relaxed roots have no band-admissibility cutoff. For them only the scan-domain
            # endpoint is meaningful; keep the admissibility field blank instead of conflating.
            admissible_boundary = at_boundary if policy == "strict" else None
            rows.append(dict(V_SO_over_t=vso, T_over_Tc=frac, U_peak_over_t=u,
                             W_WN_peak=w, policy=policy,
                             at_admissible_domain_boundary=admissible_boundary,
                             at_sampled_U_domain_boundary=at_boundary,
                             U_domain_min_over_t=umin, U_domain_max_over_t=umax,
                             boundary_side=side, valid_U_samples=nvalid))
    return rows


def write_peak_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_two_panel_figure(rows: list[dict], png_path: Path, pdf_path: Path, policy: str) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    cmap = plt.get_cmap("tab10")
    colors = {v: cmap(i) for i, v in enumerate(SOC_VALUES)}
    t_label = {0.3: "0.3", 0.5: "0.5", 0.7: "0.7", 0.9: "0.9", 0.999: "0.999"}
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.labelsize": 11, "axes.titlesize": 11,
        "legend.fontsize": 9, "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11.4, 5.4), sharey=True)
    fig.subplots_adjust(left=0.09, right=0.79, bottom=0.24, top=0.82, wspace=0.30)

    # Panel A: individual maximum points only; do not connect common-U temperatures.
    for r in rows:
        u, w = r["U_peak_over_t"], r["W_WN_peak"]
        if u is None or w is None:
            continue
        boundary = r["at_sampled_U_domain_boundary"] is True
        ax_a.scatter(u, w, s=68, marker=MARKERS[r["T_over_Tc"]],
                     facecolors="none" if boundary else colors[r["V_SO_over_t"]],
                     edgecolors="black" if boundary else colors[r["V_SO_over_t"]],
                     linewidths=1.8 if boundary else 0.9, zorder=3)

    # Panel B: same maxima rows, reorganized by temperature and connected per SOC.
    for vso in SOC_VALUES:
        q = sorted((r for r in rows if r["V_SO_over_t"] == vso), key=lambda r: r["T_over_Tc"])
        x = [r["T_over_Tc"] if r["W_WN_peak"] is not None else float("nan") for r in q]
        y = [r["W_WN_peak"] if r["W_WN_peak"] is not None else float("nan") for r in q]
        ax_b.plot(x, y, color=colors[vso], linewidth=1.8, zorder=2)
        for r in q:
            if r["W_WN_peak"] is None:
                continue
            boundary = r["at_sampled_U_domain_boundary"] is True
            ax_b.scatter(r["T_over_Tc"], r["W_WN_peak"], s=58,
                         marker=MARKERS[r["T_over_Tc"]],
                         facecolors="none" if boundary else colors[vso],
                         edgecolors="black" if boundary else colors[vso],
                         linewidths=1.8 if boundary else 0.9, zorder=3)

    ax_a.set_xlabel(r"$U_{\mathrm{peak}}/t$")
    ax_a.set_ylabel(r"$(W/W_N)_{\mathrm{peak}}$")
    ax_a.set_title("A  Temperature-specific peaks vs. maximizing U", loc="left", pad=9)
    valid_u = [r["U_peak_over_t"] for r in rows if r["U_peak_over_t"] is not None]
    if valid_u:
        span = max(valid_u) - min(valid_u)
        pad = max(1.0, 0.08 * span)
        ax_a.set_xlim(min(valid_u) - pad, max(valid_u) + pad)
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
    fig.legend(handles=soc_handles, loc="upper center", bbox_to_anchor=(0.43, 0.985),
               ncol=5, frameon=False, title=r"SOC (same colors in both panels)", title_fontsize=9)
    fig.legend(handles=temp_handles, loc="center left", bbox_to_anchor=(0.81, 0.57),
               ncol=1, frameon=False, title=r"temperature marker", title_fontsize=9)

    boundary_count = sum(r["at_sampled_U_domain_boundary"] is True for r in rows)
    policy_note = ("Strict-admissible roots only." if policy == "strict" else
                   "Relaxed roots; out-of-band points are diagnostics, not strict physical acceptance.")
    endpoint_note = (f"U sampled at 0.5t; no interpolation/smoothing; {boundary_count} admissible-U boundary peak(s)."
                     if policy == "strict" else
                     f"U sampled at 0.5t; no interpolation/smoothing; {boundary_count} sampled-U endpoint peak(s).")
    fig.text(0.43, 0.12, policy_note, ha="center", fontsize=8.0, color="0.3")
    fig.text(0.43, 0.08, endpoint_note, ha="center", fontsize=8.0, color="0.3")

    fig.savefig(png_path, dpi=360, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--policy", choices=("strict", "relaxed"), default="strict",
                    help="plot only this root policy; relaxed mode is diagnostic")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = load_peak_rows(args.input, policy=args.policy)
    suffix = "" if args.policy == "strict" else "_relaxed"
    stem = f"temperature_peak_analysis_two_panel{suffix}"
    csv_path = args.output_dir / f"{stem}.csv"
    png_path = args.output_dir / f"{stem}.png"
    pdf_path = args.output_dir / f"{stem}.pdf"
    write_peak_csv(csv_path, rows)
    make_two_panel_figure(rows, png_path, pdf_path, args.policy)
    valid = [r for r in rows if r["U_peak_over_t"] is not None]
    boundary_count = sum(r["at_sampled_U_domain_boundary"] is True for r in valid)
    print(f"Policy={args.policy}: peak rows {len(valid)}/{len(rows)}; "
          f"sampled-domain boundary peaks={boundary_count}")
    if args.policy == "strict":
        print(f"Admissible-domain boundary peaks="
              f"{sum(r['at_admissible_domain_boundary'] is True for r in valid)}")
    print(f"Wrote {csv_path}\nWrote {png_path}\nWrote {pdf_path}")


if __name__ == "__main__":
    main()
