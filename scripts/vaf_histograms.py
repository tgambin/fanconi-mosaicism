#!/usr/bin/env python3
"""
VAF histograms — visual support for R1.2 (finer VAF distribution, coloured
by the depth-aware binomial classification from R1.m14).

For each trio × caller, plots the VAF distribution of PASS_all_igv_confirmed
variants, stacked by `binom_class` (mosaic = red, het = blue). Reference
lines at VAF = 0.3741 (old fixed cutoff) and VAF = 0.5 (germline het mode).

Output: outputs/benchmark/vaf_histograms.png (one composite figure,
5 trios × 2 callers) + per-caller standalone PNGs.

Usage:
  python scripts/vaf_histograms.py
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

from _common import (
    OUTPUTS_DIR, TRIOS,
    VAF_HET_HIGH, VAF_MOSAIC_HIGH, VAF_MOSAIC_LOW,
    load_part2_all,
)


COLORS = {
    "mosaic": "#d62728",  # muted red
    "het":    "#1f77b4",  # muted blue
}

# 2.5 % bins across the full mosaic+het range that the pipeline targets
# (0.0955-0.6259). Reviewers specifically asked about the 0.1-0.2 bin.
BIN_EDGES = np.arange(0.0, 0.70 + 1e-9, 0.025)


def _extract(df: pl.DataFrame, cls: str) -> np.ndarray:
    return (
        df.filter(pl.col("PASS_all_igv_confirmed") & (pl.col("binom_class") == cls))
        ["VAF"].drop_nulls().to_numpy()
    )


def _plot_panel(ax, df: pl.DataFrame, title: str):
    """Stacked histogram of binom_class categories, with reference lines.

    Only the two in-range classes (mosaic, het) are plotted — sub-range
    (`low_mosaic`, VAF < 0.0955) and over-range (`high_vaf`, VAF > 0.6259)
    would dominate the Mutect2 y-axis and mask the signal we care about.
    """
    layers = [
        ("het",    _extract(df, "het")),
        ("mosaic", _extract(df, "mosaic")),
    ]
    values = [v for _, v in layers]
    labels = [k for k, _ in layers]
    colors = [COLORS[k] for k in labels]

    # matplotlib stacked histogram
    ax.hist(
        values, bins=BIN_EDGES, stacked=True, color=colors,
        label=labels, edgecolor="white", linewidth=0.3,
    )

    # Reference lines
    ax.axvline(0.3741, color="#666666", linestyle="--", linewidth=0.9,
               label="_nolegend_")
    ax.axvline(0.5, color="#444444", linestyle=":", linewidth=0.9,
               label="_nolegend_")

    # Annotations for thresholds (only one panel per figure gets these labels — see caller)
    ax.set_xlim(0.0, 0.70)
    ax.set_xlabel("VAF")
    ax.set_ylabel("count")
    ax.set_title(title, fontsize=10)
    ax.grid(axis="y", linestyle=":", linewidth=0.4, alpha=0.6)

    # Per-panel counts annotation (top-right)
    n_mosaic = int(len(_extract(df, "mosaic")))
    n_het = int(len(_extract(df, "het")))
    ax.text(
        0.97, 0.93,
        f"mosaic: {n_mosaic}\nhet: {n_het}",
        transform=ax.transAxes,
        ha="right", va="top",
        fontsize=8,
        bbox=dict(facecolor="white", edgecolor="#dddddd", boxstyle="round,pad=0.3", alpha=0.85),
    )


def plot_all(trios: list[str], out_dir):
    """Produce the composite figure + per-caller standalone figures."""
    out_dir.mkdir(parents=True, exist_ok=True)

    callers = [("HC (HaplotypeCaller)", "hc"), ("Mutect2", "mutect2")]

    # Composite: rows=trios, cols=callers
    fig, axes = plt.subplots(
        nrows=len(trios), ncols=len(callers),
        figsize=(10, 2.1 * len(trios)),
        sharex=True,
    )
    if len(trios) == 1:
        axes = np.atleast_2d(axes)

    for i, trio in enumerate(trios):
        for j, (caller_label, caller_key) in enumerate(callers):
            df = load_part2_all(trio, caller_key)
            ax = axes[i, j]
            if df is None or "binom_class" not in df.columns:
                ax.set_title(f"{trio} / {caller_label} — missing", fontsize=10)
                ax.set_axis_off()
                continue
            _plot_panel(ax, df, f"{trio} — {caller_label}")

    # Shared legend (top of figure)
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color=COLORS["mosaic"], label="mosaic (binom_class)"),
        plt.Rectangle((0, 0), 1, 1, color=COLORS["het"],    label="het (binom_class)"),
        plt.Line2D([0], [0], color="#666666", linestyle="--", label="old VAF cutoff (37.41%)"),
        plt.Line2D([0], [0], color="#444444", linestyle=":",  label="germline het (VAF = 0.5)"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.005),
        fontsize=9, frameon=False,
    )
    fig.suptitle(
        "VAF distribution of pipeline- and IGV-confirmed de novo candidates, "
        "coloured by the depth-aware binomial classification (R1.m14)",
        y=1.04, fontsize=11,
    )
    fig.tight_layout()
    composite = out_dir / "vaf_histograms.png"
    fig.savefig(composite, dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Per-caller standalone (one row per caller, 5 trios side by side)
    for caller_label, caller_key in callers:
        fig2, axes2 = plt.subplots(
            nrows=1, ncols=len(trios),
            figsize=(2.6 * len(trios), 3.5),
            sharey=False,
        )
        if len(trios) == 1:
            axes2 = [axes2]
        for i, trio in enumerate(trios):
            df = load_part2_all(trio, caller_key)
            ax = axes2[i]
            if df is None or "binom_class" not in df.columns:
                ax.set_title(f"{trio} — missing", fontsize=10)
                ax.set_axis_off()
                continue
            _plot_panel(ax, df, trio)

        fig2.legend(
            handles=legend_handles,
            loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02),
            fontsize=8, frameon=False,
        )
        fig2.suptitle(f"VAF distribution — {caller_label}", y=1.09, fontsize=11)
        fig2.tight_layout()
        out_path = out_dir / f"vaf_histograms_{caller_key}.png"
        fig2.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig2)

    return composite


def main():
    parser = argparse.ArgumentParser(description="VAF histograms coloured by binom_class")
    parser.add_argument("--trios", default=None,
                        help="Comma-separated trios (default: all 5)")
    args = parser.parse_args()

    trios = args.trios.split(",") if args.trios else TRIOS
    out_dir = OUTPUTS_DIR / "benchmark"
    composite = plot_all(trios, out_dir)
    print(f"Composite figure → {composite}")
    print(f"Per-caller figures → {out_dir}/vaf_histograms_hc.png, "
          f"{out_dir}/vaf_histograms_mutect2.png")


if __name__ == "__main__":
    main()
