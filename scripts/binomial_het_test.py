#!/usr/bin/env python3
"""
Per-variant binomial test: mosaic vs het (Reviewer 1 request R1.m14).

Replaces the fixed VAF cutoff (0.3741) with a depth-aware statistical test.

For each variant:
    H0: VAF = 0.5   (germline heterozygous)
    H1: VAF < 0.5   (subclonal / mosaic)
Test: scipy.stats.binomtest(AD_alt, DP, p=0.5, alternative="less")
BH correction (Benjamini-Hochberg) on the subset of PASS_all variants.

Adds three columns to each outputs/part2/{trio}[_mutect2]_part2_all.tsv:
  - binom_het_pval       raw one-sided p-value, per variant
  - binom_het_pval_BH    BH-adjusted, computed over PASS_all variants only
                         (non-PASS_all rows get 1.0)
  - binom_class          mosaic / het / low_mosaic / high_vaf
                           mosaic      : VAF_MOSAIC_LOW ≤ VAF ≤ VAF_HET_HIGH
                                         AND binom_het_pval_BH < BINOM_HET_ALPHA
                           het         : VAF_MOSAIC_LOW ≤ VAF ≤ VAF_HET_HIGH
                                         AND binom_het_pval_BH ≥ BINOM_HET_ALPHA
                           low_mosaic  : VAF < VAF_MOSAIC_LOW
                           high_vaf    : VAF > VAF_HET_HIGH

The existing VAF-based `VAF_class` column is preserved unchanged for
side-by-side comparison in the report (§2d).

Idempotent. Safe to re-run; drops pre-existing columns before re-adding.

Usage:
  python scripts/binomial_het_test.py                    # all 5 trios
  python scripts/binomial_het_test.py --trios brca1,brca2
"""

import argparse

import numpy as np
import polars as pl
from scipy.stats import binomtest

from _common import (
    BINOM_HET_ALPHA,
    OUTPUTS_DIR, POLARS_INFER_SCHEMA, TRIOS,
    VAF_HET_HIGH, VAF_MOSAIC_LOW,
    load_part2_all, part2_all_path,
)


BINOM_COLS = ("binom_het_pval", "binom_het_pval_BH", "binom_class")


def _binom_p_less(ad_alt: int | None, dp: int | None) -> float:
    """One-sided p-value for H0: VAF=0.5 vs H1: VAF<0.5; 1.0 on bad input."""
    if ad_alt is None or dp is None or dp <= 0 or ad_alt < 0 or ad_alt > dp:
        return 1.0
    return binomtest(int(ad_alt), int(dp), 0.5, alternative="less").pvalue


def _bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment, vectorised (same logic as deepSNV)."""
    n = len(pvals)
    if n == 0:
        return pvals
    order = np.argsort(pvals)
    ranked = pvals[order]
    adj = ranked * n / (np.arange(n) + 1)
    # monotone: take running min from the end
    cummin = np.minimum.accumulate(adj[::-1])[::-1]
    cummin = np.clip(cummin, 0.0, 1.0)
    out = np.empty(n)
    out[order] = cummin
    return out


def apply_to_trio_caller(trio: str, caller: str) -> dict:
    """Mutate outputs/part2/{trio}[_mutect2]_part2_all.tsv in place."""
    df = load_part2_all(trio, caller)
    if df is None:
        return {"trio": trio, "caller": caller, "status": "missing"}

    # idempotency
    df = df.drop([c for c in BINOM_COLS if c in df.columns])

    # Raw p-value for every row (cheap: ~2 µs per binomtest)
    ad_alt = df["AD_alt"].to_list()
    dp = df["DP"].to_list()
    pvals = np.array([_binom_p_less(a, d) for a, d in zip(ad_alt, dp)])

    # BH adjust only on PASS_all subset
    pass_mask = df["PASS_all"].to_numpy()
    pass_idx = np.where(pass_mask)[0]
    bh_full = np.ones(len(df))
    if pass_idx.size > 0:
        bh_full[pass_idx] = _bh_adjust(pvals[pass_idx])

    df = df.with_columns(
        pl.Series("binom_het_pval", pvals),
        pl.Series("binom_het_pval_BH", bh_full),
    )

    # Classification
    in_range = (pl.col("VAF") >= VAF_MOSAIC_LOW) & (pl.col("VAF") <= VAF_HET_HIGH)
    df = df.with_columns(
        pl.when(pl.col("VAF") < VAF_MOSAIC_LOW).then(pl.lit("low_mosaic"))
          .when(pl.col("VAF") > VAF_HET_HIGH).then(pl.lit("high_vaf"))
          .when(in_range & (pl.col("binom_het_pval_BH") < BINOM_HET_ALPHA))
            .then(pl.lit("mosaic"))
          .when(in_range).then(pl.lit("het"))
          .otherwise(pl.lit("unknown"))
          .alias("binom_class")
    )

    df.write_csv(str(part2_all_path(trio, caller)), separator="\t")

    # Summary for console
    pass_df = df.filter(pl.col("PASS_all"))
    old_mosaic = len(pass_df.filter(pl.col("VAF_class") == "mosaic"))
    old_het = len(pass_df.filter(pl.col("VAF_class") == "heterozygous"))
    new_mosaic = len(pass_df.filter(pl.col("binom_class") == "mosaic"))
    new_het = len(pass_df.filter(pl.col("binom_class") == "het"))
    moved_het_to_mosaic = len(pass_df.filter(
        (pl.col("VAF_class") == "heterozygous") & (pl.col("binom_class") == "mosaic")
    ))
    moved_mosaic_to_het = len(pass_df.filter(
        (pl.col("VAF_class") == "mosaic") & (pl.col("binom_class") == "het")
    ))
    return {
        "trio": trio, "caller": caller, "status": "ok",
        "old_mosaic": old_mosaic, "old_het": old_het,
        "new_mosaic": new_mosaic, "new_het": new_het,
        "het_to_mosaic": moved_het_to_mosaic,
        "mosaic_to_het": moved_mosaic_to_het,
    }


def main():
    parser = argparse.ArgumentParser(description="Per-variant binomial het test (R1.m14)")
    parser.add_argument("--trios", default=None,
                        help="Comma-separated trios (default: all 5)")
    args = parser.parse_args()

    trios = args.trios.split(",") if args.trios else TRIOS

    print(f"Binomial test H0: VAF=0.5 vs H1: VAF<0.5  |  BH α = {BINOM_HET_ALPHA}")
    print(f"{'Trio':<8} {'Caller':<8} "
          f"{'VAF mosaic':>10} {'→ binom':>8}  "
          f"{'VAF het':>8} {'→ binom':>8}  "
          f"{'het→mos':>8} {'mos→het':>8}")
    print("-" * 78)
    for trio in trios:
        for caller in ["hc", "mutect2"]:
            r = apply_to_trio_caller(trio, caller)
            if r["status"] != "ok":
                print(f"{trio:<8} {caller:<8}  (skipped — missing _part2_all.tsv)")
                continue
            print(
                f"{r['trio']:<8} {r['caller']:<8} "
                f"{r['old_mosaic']:>10} {r['new_mosaic']:>8}  "
                f"{r['old_het']:>8} {r['new_het']:>8}  "
                f"{r['het_to_mosaic']:>8} {r['mosaic_to_het']:>8}"
            )


if __name__ == "__main__":
    main()
