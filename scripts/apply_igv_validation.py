#!/usr/bin/env python3
"""
Apply manual IGV validation annotations to Part 2 outputs.

Adds four columns to each outputs/part2/{trio}[_mutect2]_part2_all.tsv:
  - IGV_is_real              TRUE / FALSE / "" (not reviewed — outside the
                             master-CSV subset, i.e. already filtered out)
  - IGV_notes                error class p1..p6 (only for FALSE), or ""
  - igvFilter                bool — False only when explicitly FALSE;
                             True for TRUE and for unreviewed variants
                             (neutral, so PASS_all is not blocked where
                             there is no evidence of a false positive).
  - PASS_all_igv_confirmed   PASS_all AND igvFilter

Input (ground truth): outputs/igv_review/all_trios_igv_validated.csv
(which carries per-variant TRUE/FALSE/NA values in `IGV_is_real` and
optional `p1..p6` codes in `IGV_notes`).

NA rows in the master CSV are old-pipeline-only variants (never appeared
in HC or Mutect2) and are therefore absent from any `_part2_all.tsv`
— this step leaves them on the side.

Idempotent: safe to re-run; it drops any pre-existing IGV review columns
(both legacy Haowei_* and current IGV_* names) before re-adding them.

Usage:
  python scripts/apply_igv_validation.py                    # all 5 trios
  python scripts/apply_igv_validation.py --trios brca1,brca2
"""

import argparse

import polars as pl

from _common import (
    IGV_VALIDATION_CSV, OUTPUTS_DIR, POLARS_INFER_SCHEMA, TRIOS,
    normalize_bool_columns, part2_all_path,
)


IGV_REVIEW_COLS = (
    "IGV_is_real", "IGV_notes",
    "igvFilter", "PASS_all_igv_confirmed",
)

# Legacy column names from earlier runs — dropped for idempotency when
# re-running after the Haowei → IGV rename (2026-04-17).
LEGACY_IGV_REVIEW_COLS = (
    "Haowei_is_real", "Haowei_notes",
    "haoweiFilter", "PASS_all_haowei_confirmed",
)


def load_igv_for_trio(validation_df: pl.DataFrame, trio: str) -> pl.DataFrame:
    """Slice the validation CSV for one trio, keep key columns + annotations."""
    return (
        validation_df
        .filter(pl.col("trio") == trio)
        .select(
            pl.col("CHROM"),
            pl.col("POS").cast(pl.Int64),
            pl.col("REF"),
            pl.col("ALT"),
            pl.col("IGV_is_real").cast(pl.Utf8),
            pl.col("IGV_notes").cast(pl.Utf8),
        )
        # strip the "NA" sentinel — treat as missing so we don't join it in
        .with_columns(
            pl.when(pl.col("IGV_is_real") == "NA").then(None)
              .otherwise(pl.col("IGV_is_real")).alias("IGV_is_real"),
            pl.when(pl.col("IGV_notes") == "NA").then(None)
              .otherwise(pl.col("IGV_notes")).alias("IGV_notes"),
        )
    )


def apply_to_trio_caller(trio: str, caller: str, validation_df: pl.DataFrame) -> dict:
    """Mutate outputs/part2/{trio}[_mutect2]_part2_all.tsv in place.

    Returns counts for reporting."""
    path = part2_all_path(trio, caller)
    if not path.exists():
        return {"trio": trio, "caller": caller, "status": "missing"}

    df = pl.read_csv(str(path), separator="\t", infer_schema_length=POLARS_INFER_SCHEMA)
    df = normalize_bool_columns(df)

    # idempotency: drop any pre-existing IGV review columns (both legacy and current)
    stale = [c for c in (*IGV_REVIEW_COLS, *LEGACY_IGV_REVIEW_COLS) if c in df.columns]
    if stale:
        df = df.drop(stale)

    annot = load_igv_for_trio(validation_df, trio)

    df = df.join(
        annot,
        on=["CHROM", "POS", "REF", "ALT"],
        how="left",
    )

    # igvFilter: False only for explicit FALSE; TRUE for everyone else
    # (including unreviewed / not in master CSV — neutral).
    df = df.with_columns(
        (pl.col("IGV_is_real") != "FALSE")
        .fill_null(True)
        .alias("igvFilter")
    )
    df = df.with_columns(
        (pl.col("PASS_all") & pl.col("igvFilter")).alias("PASS_all_igv_confirmed")
    )

    df = df.with_columns(
        pl.col("IGV_is_real").fill_null("").alias("IGV_is_real"),
        pl.col("IGV_notes").fill_null("").alias("IGV_notes"),
    )

    # Normalise the legacy "sub_mosaic" label that comes out of Part2 to the
    # current name "low_mosaic" — avoids a slow full Part2 re-run just to
    # rename a string.  Idempotent.
    if "VAF_class" in df.columns:
        df = df.with_columns(
            pl.when(pl.col("VAF_class") == "sub_mosaic")
              .then(pl.lit("low_mosaic"))
              .otherwise(pl.col("VAF_class"))
              .alias("VAF_class")
        )

    df.write_csv(str(path), separator="\t")

    n_total = len(df)
    n_reviewed = len(df.filter(pl.col("IGV_is_real") != ""))
    n_true = len(df.filter(pl.col("IGV_is_real") == "TRUE"))
    n_false = len(df.filter(pl.col("IGV_is_real") == "FALSE"))
    n_pass_before = len(df.filter(pl.col("PASS_all")))
    n_pass_after = len(df.filter(pl.col("PASS_all_igv_confirmed")))
    return {
        "trio": trio, "caller": caller, "status": "ok",
        "total": n_total, "reviewed": n_reviewed,
        "igv_true": n_true, "igv_false": n_false,
        "pass_before": n_pass_before, "pass_after": n_pass_after,
    }


def main():
    parser = argparse.ArgumentParser(description="Apply manual IGV validation")
    parser.add_argument("--trios", default=None,
                        help="Comma-separated trios (default: all 5)")
    args = parser.parse_args()

    trios = args.trios.split(",") if args.trios else TRIOS

    if not IGV_VALIDATION_CSV.exists():
        raise SystemExit(f"Validation CSV not found: {IGV_VALIDATION_CSV}")

    validation_df = pl.read_csv(
        str(IGV_VALIDATION_CSV), infer_schema_length=POLARS_INFER_SCHEMA
    )
    print(f"Loaded {len(validation_df):,} validated variants "
          f"from {IGV_VALIDATION_CSV.name}")

    rows = []
    for trio in trios:
        for caller in ["hc", "mutect2"]:
            rows.append(apply_to_trio_caller(trio, caller, validation_df))

    print(f"\n{'Trio':<8} {'Caller':<8} {'PASS':>6} → {'IGV':>6}  "
          f"{'Reviewed':>9} {'TRUE':>5} {'FALSE':>5}")
    print("-" * 60)
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['trio']:<8} {r['caller']:<8} {'— skipped (missing file)':>40}")
            continue
        print(
            f"{r['trio']:<8} {r['caller']:<8} "
            f"{r['pass_before']:>6} → {r['pass_after']:>6}  "
            f"{r['reviewed']:>9} {r['igv_true']:>5} {r['igv_false']:>5}"
        )


if __name__ == "__main__":
    main()
