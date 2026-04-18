#!/usr/bin/env python3
"""
MosaicPostProcessPart2 — Local post-processing of Part 1 pileup CSV.

Replicates the R pipeline (fanconi-genome.R) in Python using polars + polars-bio.

Steps:
  1. Load Part 1 CSV (VCF annotations + pileup counts)
  2. Pileup-based de novo validation (parent alt=0, proband alt>=2, depth filters)
  3. Hard filters (MQ, QD, FS, ReadPosRankSum — different thresholds for SNV vs INDEL)
  4. Centromere + Alu/repeat exclusion (polars-bio interval overlap)
  5. gnomAD AF=0 filter (optional, if gnomAD file provided)
  6. VAF calculation and binning (mosaic / heterozygous)
  7. Output: final TSV with all filter columns + classification

Usage:
  python scripts/mosaic_postprocess_part2.py outputs/part1_csv/proband_brca1_150x.part1.csv --sample brca1
  python scripts/mosaic_postprocess_part2.py outputs/part1_csv/proband_brca2_150x.part1.csv --sample brca2
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import polars as pl
import polars_bio as pb
from scipy.stats import betabinom, chi2

from _common import (
    PROJECT_DIR,
    VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH, VAF_HET_HIGH, VAF_MOSAIC_MID,
    MIN_MOTHER_DEPTH, MIN_FATHER_DEPTH, MIN_PROBAND_DEPTH, MIN_VCF_DP,
    MIN_ALT_READS,
    HC_SNV_MIN_MQ, HC_SNV_MIN_READ_POS_RANK_SUM, HC_SNV_MAX_FS, HC_SNV_MIN_QD,
    HC_INDEL_MIN_READ_POS_RANK_SUM, HC_INDEL_MAX_FS, HC_INDEL_MIN_QD,
    DEEPSNV_ESTIMATION_WINDOW, DEEPSNV_MAX_ERROR_ALT, DEEPSNV_MIN_OBSERVATIONS,
    DEEPSNV_DEFAULT_RHO, DEEPSNV_FALLBACK_RHO,
    DEEPSNV_MIN_RHO, DEEPSNV_MAX_RHO, DEEPSNV_EPS,
    BH_SIG_STRONG,
    VEP_DATA, VEP_IMAGE, VEP_FORKS,
    POLARS_INFER_SCHEMA,
)


# --- Constants ---
ANNOTATIONS_DIR = PROJECT_DIR / "data" / "annotations"
CENTROMERES_BED = ANNOTATIONS_DIR / "centromeres.bed"
REPEATS_BED = ANNOTATIONS_DIR / "repeats.bed"  # Pre-filtered Alu + centr


def load_part1_csv(path: str, caller: str = "hc") -> pl.DataFrame:
    """Load Part 1 CSV and parse AD field."""
    df = pl.read_csv(path, infer_schema_length=POLARS_INFER_SCHEMA)
    # Parse AD field: "ref_count,alt_count" → extract alt count
    df = df.with_columns(
        pl.col("AD").str.split(",").list.get(1).cast(pl.Int64).alias("AD_alt"),
        pl.col("AD").str.split(",").list.get(0).cast(pl.Int64).alias("AD_ref"),
    )
    if caller == "hc":
        # Cast string columns that may contain "." for missing values
        for col_name in ["ReadPosRankSum", "BaseQRankSum", "MQRankSum"]:
            if col_name in df.columns and df[col_name].dtype == pl.String:
                df = df.with_columns(
                    pl.when(pl.col(col_name) == ".")
                    .then(pl.lit(None))
                    .otherwise(pl.col(col_name))
                    .cast(pl.Float64)
                    .alias(col_name)
                )
    elif caller == "mutect2":
        # Cast Mutect2-specific string columns that may contain "."
        for col_name in ["TLOD", "GERMQ", "SEQQ", "STRANDQ", "AF"]:
            if col_name in df.columns and df[col_name].dtype == pl.String:
                df = df.with_columns(
                    pl.when(pl.col(col_name) == ".")
                    .then(pl.lit(None))
                    .otherwise(pl.col(col_name))
                    .cast(pl.Float64)
                    .alias(col_name)
                )
    print(f"Loaded {len(df):,} variants from {path} (caller={caller})")
    return df


def classify_variant_type(df: pl.DataFrame) -> pl.DataFrame:
    """Classify as SNV or INDEL based on REF/ALT length."""
    return df.with_columns(
        pl.when(
            (pl.col("REF").str.len_chars() == 1)
            & (pl.col("ALT").str.len_chars() == 1)
            & pl.col("REF").is_in(["A", "C", "G", "T"])
            & pl.col("ALT").is_in(["A", "C", "G", "T"])
        )
        .then(pl.lit("SNV"))
        .otherwise(pl.lit("INDEL"))
        .alias("Type")
    )


def add_pileup_filter(df: pl.DataFrame) -> pl.DataFrame:
    """
    Pileup-based de novo validation — replicates addPileupFilter() from R.

    Filters:
    1. Depth: mother_depth > 20, father_depth > 20, proband_depth > 10
    2. Parent filter: parent ALT base count = 0 (SNV) or parent DEL/INS = 0 (INDEL)
    3. Proband filter: proband ALT base count >= 2 (SNV) or proband DEL/INS >= 2 (INDEL)
    4. VCF-level: AD_alt >= 2, DP > 20
    """
    # Step 1: Depth filter
    df = df.with_columns(
        pl.when(
            (pl.col("mother_depth") > MIN_MOTHER_DEPTH)
            & (pl.col("father_depth") > MIN_FATHER_DEPTH)
            & (pl.col("proband_depth") > MIN_PROBAND_DEPTH)
        )
        .then(pl.lit("PASS"))
        .otherwise(pl.lit("FAIL"))
        .alias("pileupDepthFilterTrio")
    )

    # Step 2: Parent filter — SNV: parent ALT base = 0
    snv_parent_pass = (
        ((pl.col("ALT") == "A") & (pl.col("father_A") == 0) & (pl.col("mother_A") == 0))
        | ((pl.col("ALT") == "C") & (pl.col("father_C") == 0) & (pl.col("mother_C") == 0))
        | ((pl.col("ALT") == "T") & (pl.col("father_T") == 0) & (pl.col("mother_T") == 0))
        | ((pl.col("ALT") == "G") & (pl.col("father_G") == 0) & (pl.col("mother_G") == 0))
    )

    # Parent filter — INDEL: deletion (len(REF)>1) → parent DEL=0; insertion (len(ALT)>1) → parent INS=0
    is_deletion = pl.col("REF").str.len_chars() > pl.col("ALT").str.len_chars()
    is_insertion = pl.col("ALT").str.len_chars() > pl.col("REF").str.len_chars()

    indel_parent_pass = (
        (is_deletion & (pl.col("proband_DEL") > 0) & (pl.col("mother_DEL") == 0) & (pl.col("father_DEL") == 0))
        | (is_insertion & (pl.col("proband_INS") > 0) & (pl.col("mother_INS") == 0) & (pl.col("father_INS") == 0))
    )

    df = df.with_columns(
        pl.when(pl.col("Type") == "SNV")
        .then(pl.when(snv_parent_pass).then(pl.lit("PASS")).otherwise(pl.lit("FAIL")))
        .otherwise(pl.when(indel_parent_pass).then(pl.lit("PASS")).otherwise(pl.lit("FAIL")))
        .alias("pileupParentFilter")
    )

    # Step 3: Proband alt depth filter
    snv_proband_pass = (
        ((pl.col("ALT") == "A") & (pl.col("proband_A") >= MIN_ALT_READS))
        | ((pl.col("ALT") == "C") & (pl.col("proband_C") >= MIN_ALT_READS))
        | ((pl.col("ALT") == "T") & (pl.col("proband_T") >= MIN_ALT_READS))
        | ((pl.col("ALT") == "G") & (pl.col("proband_G") >= MIN_ALT_READS))
    )

    indel_proband_pass = (
        (is_deletion & (pl.col("proband_DEL") >= MIN_ALT_READS))
        | (is_insertion & (pl.col("proband_INS") >= MIN_ALT_READS))
    )

    df = df.with_columns(
        pl.when(pl.col("Type") == "SNV")
        .then(pl.when(snv_proband_pass).then(pl.lit("PASS")).otherwise(pl.lit("FAIL")))
        .otherwise(pl.when(indel_proband_pass).then(pl.lit("PASS")).otherwise(pl.lit("FAIL")))
        .alias("pileupProbandFilter")
    )

    # Step 4: Combined pileup filter (includes VCF AD/DP check from R)
    df = df.with_columns(
        pl.when(
            (pl.col("pileupDepthFilterTrio") == "PASS")
            & (pl.col("pileupParentFilter") == "PASS")
            & (pl.col("pileupProbandFilter") == "PASS")
            & (pl.col("AD_alt") >= MIN_ALT_READS)
            & (pl.col("DP") > MIN_VCF_DP)
        )
        .then(pl.lit("PASS"))
        .otherwise(pl.lit("FAIL"))
        .alias("pileupFullFilter")
    )

    return df


def add_hard_filter(df: pl.DataFrame, caller: str = "hc") -> pl.DataFrame:
    """
    Hard quality filters — replicates addHardFilter() from R.

    HC:      SNV: MQ > 30, ReadPosRankSum > -8, FS < 60, QD > 2
             INDEL: ReadPosRankSum > -20, FS < 200, QD > 2
    Mutect2: All PASS (FilterMutectCalls already applied in Part 1)
    """
    if caller == "mutect2":
        # FilterMutectCalls already applied — all variants here passed FILTER=PASS
        df = df.with_columns(pl.lit("PASS").alias("HardFilter"))
        return df

    snv_pass = (
        (pl.col("MQ") > HC_SNV_MIN_MQ)
        & (pl.col("ReadPosRankSum") > HC_SNV_MIN_READ_POS_RANK_SUM)
        & (pl.col("FS") < HC_SNV_MAX_FS)
        & (pl.col("QD") > HC_SNV_MIN_QD)
    )

    indel_pass = (
        (pl.col("ReadPosRankSum") > HC_INDEL_MIN_READ_POS_RANK_SUM)
        & (pl.col("FS") < HC_INDEL_MAX_FS)
        & (pl.col("QD") > HC_INDEL_MIN_QD)
    )

    df = df.with_columns(
        pl.when(pl.col("Type") == "SNV")
        .then(pl.when(snv_pass).then(pl.lit("PASS")).otherwise(pl.lit("FAIL")))
        .otherwise(pl.when(indel_pass).then(pl.lit("PASS")).otherwise(pl.lit("FAIL")))
        .alias("HardFilter")
    )

    return df


def _get_alt_count_snv(row: dict) -> tuple:
    """Get proband alt count and pooled parent alt count+depth for an SNV."""
    alt = row["ALT"]
    p_alt = row.get(f"proband_{alt}", 0)
    p_depth = row["proband_depth"]
    m_alt = row.get(f"mother_{alt}", 0)
    m_depth = row["mother_depth"]
    f_alt = row.get(f"father_{alt}", 0)
    f_depth = row["father_depth"]
    # Pool parents as control
    c_alt = m_alt + f_alt
    c_depth = m_depth + f_depth
    return int(p_alt), int(p_depth), int(c_alt), int(c_depth)


def _get_alt_count_indel(row: dict) -> tuple:
    """Get proband alt count and pooled parent alt count+depth for an INDEL."""
    ref_len = len(row["REF"])
    alt_len = len(row["ALT"])
    if ref_len > alt_len:  # deletion
        field = "DEL"
    else:  # insertion
        field = "INS"
    p_alt = row.get(f"proband_{field}", 0)
    p_depth = row["proband_depth"]
    m_alt = row.get(f"mother_{field}", 0)
    m_depth = row["mother_depth"]
    f_alt = row.get(f"father_{field}", 0)
    f_depth = row["father_depth"]
    c_alt = m_alt + f_alt
    c_depth = m_depth + f_depth
    return int(p_alt), int(p_depth), int(c_alt), int(c_depth)


def estimate_overdispersion(df: pl.DataFrame) -> float:
    """
    Estimate overdispersion (rho) from error-only positions in parents.

    Uses method-of-moments on per-parent error rates at positions where
    parent alt count <= 2 (likely sequencing error, not real variants).
    """
    snv_positions = df.filter(
        (pl.col("mother_depth") > MIN_MOTHER_DEPTH)
        & (pl.col("father_depth") > MIN_FATHER_DEPTH)
        & (pl.col("Type") == "SNV")
    ).head(DEEPSNV_ESTIMATION_WINDOW)

    if len(snv_positions) < DEEPSNV_MIN_OBSERVATIONS:
        return DEEPSNV_DEFAULT_RHO

    error_rates = []
    for row in snv_positions.iter_rows(named=True):
        ref = row["REF"]
        for who in ["mother", "father"]:
            depth = row[f"{who}_depth"]
            if depth == 0:
                continue
            alt = sum(row.get(f"{who}_{b}", 0) for b in "ACGT" if b != ref)
            # Only use if alt is small (error, not variant)
            if alt <= DEEPSNV_MAX_ERROR_ALT:
                error_rates.append(alt / depth)

    if len(error_rates) < DEEPSNV_MIN_OBSERVATIONS:
        return DEEPSNV_DEFAULT_RHO

    rates = np.array(error_rates)
    mean_p = np.mean(rates)
    var_p = np.var(rates)

    # Method of moments: approximate per-parent depth
    n_mean = (snv_positions["mother_depth"].mean() + snv_positions["father_depth"].mean()) / 2
    binom_var = mean_p * (1 - mean_p) / n_mean if n_mean > 0 else 0

    if var_p <= binom_var or mean_p == 0 or mean_p == 1:
        return DEEPSNV_FALLBACK_RHO

    rho = (var_p - binom_var) / (mean_p * (1 - mean_p) - binom_var)
    rho = max(DEEPSNV_MIN_RHO, min(rho, DEEPSNV_MAX_RHO))
    return rho


def _betabinom_logpmf(k, n, mu, rho):
    """Log-PMF of beta-binomial distribution."""
    if n == 0 or mu <= 0 or mu >= 1:
        return 0.0
    alpha = mu * (1.0 / rho - 1.0)
    beta_param = (1.0 - mu) * (1.0 / rho - 1.0)
    alpha = max(alpha, DEEPSNV_EPS)
    beta_param = max(beta_param, DEEPSNV_EPS)
    return betabinom.logpmf(k, n, alpha, beta_param)


def deepsnv_test_one(k_test, n_test, k_ctrl, n_ctrl, rho):
    """
    Beta-binomial likelihood ratio test (deepSNV-style).

    H0: alt frequency equal in test and control (sequencing error)
    H1: alt frequency higher in test (true variant)

    Returns p-value.
    """
    if n_test == 0 or n_ctrl == 0:
        return 1.0
    if k_test == 0:
        return 1.0

    # Test is one-sided: only interested when test freq > control freq
    freq_test = k_test / n_test
    freq_ctrl = k_ctrl / n_ctrl if n_ctrl > 0 else 0
    if freq_test <= freq_ctrl:
        return 1.0

    # H0: shared frequency (pooled)
    mu0 = (k_test + k_ctrl) / (n_test + n_ctrl)
    mu0 = max(mu0, DEEPSNV_EPS)
    mu0 = min(mu0, 1.0 - DEEPSNV_EPS)

    ll_h0 = (_betabinom_logpmf(k_test, n_test, mu0, rho)
             + _betabinom_logpmf(k_ctrl, n_ctrl, mu0, rho))

    # H1: separate frequencies
    mu_test = max(k_test / n_test, DEEPSNV_EPS)
    mu_test = min(mu_test, 1.0 - DEEPSNV_EPS)
    mu_ctrl = max(k_ctrl / n_ctrl, DEEPSNV_EPS) if n_ctrl > 0 else DEEPSNV_EPS
    mu_ctrl = min(mu_ctrl, 1.0 - DEEPSNV_EPS)

    ll_h1 = (_betabinom_logpmf(k_test, n_test, mu_test, rho)
             + _betabinom_logpmf(k_ctrl, n_ctrl, mu_ctrl, rho))

    # LRT statistic ~ chi2(df=1)
    lrt = -2.0 * (ll_h0 - ll_h1)
    lrt = max(lrt, 0.0)

    return chi2.sf(lrt, df=1)


def add_deepsnv_pvalue(df: pl.DataFrame) -> pl.DataFrame:
    """
    Add deepSNV-style beta-binomial p-value for each variant.

    Uses proband as test, pooled parents as control.
    Estimates overdispersion from the data, then computes per-variant
    p-values. BH correction is applied only to pileup-passing variants
    (true de novo candidates), not all ~178K input variants.
    """
    print("Computing deepSNV beta-binomial p-values...")

    # Estimate overdispersion from data
    rho = estimate_overdispersion(df)
    print(f"  Estimated overdispersion rho = {rho:.4f}")

    # Compute raw p-value for each variant
    pvals = []
    for row in df.iter_rows(named=True):
        if row["Type"] == "SNV":
            k_t, n_t, k_c, n_c = _get_alt_count_snv(row)
        else:
            k_t, n_t, k_c, n_c = _get_alt_count_indel(row)
        pvals.append(deepsnv_test_one(k_t, n_t, k_c, n_c, rho))

    pvals = np.array(pvals)
    df = df.with_columns(pl.Series("deepSNV_pval", pvals))

    # BH correction only on pileup-passing variants (de novo candidates)
    pileup_mask = df["pileupFullFilter"].to_numpy()
    pileup_idx = np.where(pileup_mask)[0]
    n_tested = len(pileup_idx)

    bh_adjusted = np.ones(len(pvals))  # default 1.0 for non-tested
    if n_tested > 0:
        sub_pvals = pvals[pileup_idx]
        sorted_idx = np.argsort(sub_pvals)
        sorted_pvals = sub_pvals[sorted_idx]
        sub_bh = np.zeros(n_tested)
        cummin = 1.0
        for i in range(n_tested - 1, -1, -1):
            adj = sorted_pvals[i] * n_tested / (i + 1)
            cummin = min(cummin, adj)
            sub_bh[sorted_idx[i]] = min(cummin, 1.0)
        bh_adjusted[pileup_idx] = sub_bh

    df = df.with_columns(pl.Series("deepSNV_pval_BH", bh_adjusted))

    n_sig = np.sum(bh_adjusted < BH_SIG_STRONG)
    print(f"  BH correction on {n_tested:,} pileup-passing variants")
    print(f"  Significant at BH-adjusted p < {BH_SIG_STRONG}: {n_sig:,} / {n_tested:,}")

    return df


def add_repeat_filter(df: pl.DataFrame) -> pl.DataFrame:
    """
    Centromere + Alu/repeat exclusion using polars-bio count_overlaps.
    Replicates addRepeatFilter() from R.
    """
    # Load BED files as polars DataFrames (headerless 3-column BED)
    centromere_bed = pl.read_csv(
        str(CENTROMERES_BED), separator="\t", has_header=False,
        new_columns=["chrom", "start", "end"],
    )
    repeats_bed = pl.read_csv(
        str(REPEATS_BED), separator="\t", has_header=False,
        new_columns=["chrom", "start", "end"],
    )

    # Prepare variant positions as intervals (1-based, point intervals)
    variants = df.select(
        pl.col("CHROM").alias("chrom"),
        pl.col("POS").cast(pl.Int64).alias("start"),
        pl.col("POS").cast(pl.Int64).alias("end"),
    ).with_row_index("_idx")

    # --- Centromere filter ---
    centromere_counts = pb.count_overlaps(
        variants, centromere_bed,
        cols1=["chrom", "start", "end"],
        cols2=["chrom", "start", "end"],
        output_type="polars.DataFrame",
    )
    centromere_hits = set(
        centromere_counts.filter(pl.col("count") > 0)["_idx"].to_list()
    )

    # --- Repeat (Alu + centr) filter ---
    repeat_counts = pb.count_overlaps(
        variants, repeats_bed,
        cols1=["chrom", "start", "end"],
        cols2=["chrom", "start", "end"],
        output_type="polars.DataFrame",
    )
    repeat_hits = set(
        repeat_counts.filter(pl.col("count") > 0)["_idx"].to_list()
    )

    # Add filter columns using row index
    idx_series = pl.Series("_idx", list(range(len(df))), dtype=pl.UInt32)
    df = df.with_columns(
        pl.when(idx_series.is_in(list(centromere_hits)))
        .then(pl.lit("FAIL"))
        .otherwise(pl.lit("PASS"))
        .alias("centromereFilter"),
        pl.when(idx_series.is_in(list(repeat_hits)))
        .then(pl.lit("FAIL"))
        .otherwise(pl.lit("PASS"))
        .alias("repeatFilter"),
    )
    df = df.with_columns(
        pl.when(
            (pl.col("centromereFilter") == "PASS") & (pl.col("repeatFilter") == "PASS")
        )
        .then(pl.lit("PASS"))
        .otherwise(pl.lit("FAIL"))
        .alias("repeatFullFilter")
    )

    return df


def add_vaf_and_classify(df: pl.DataFrame) -> pl.DataFrame:
    """Calculate VAF from pileup and classify into mosaic/het bins."""
    # VAF from pileup: alt reads / total depth
    # For SNVs, use the specific ALT base count
    snv_alt_count = (
        pl.when(pl.col("ALT") == "A").then(pl.col("proband_A"))
        .when(pl.col("ALT") == "C").then(pl.col("proband_C"))
        .when(pl.col("ALT") == "T").then(pl.col("proband_T"))
        .when(pl.col("ALT") == "G").then(pl.col("proband_G"))
        .otherwise(pl.lit(0))
    )

    # For INDELs, use DEL or INS count
    is_deletion = pl.col("REF").str.len_chars() > pl.col("ALT").str.len_chars()
    indel_alt_count = (
        pl.when(is_deletion).then(pl.col("proband_DEL"))
        .otherwise(pl.col("proband_INS"))
    )

    df = df.with_columns(
        pl.when(pl.col("Type") == "SNV")
        .then(snv_alt_count)
        .otherwise(indel_alt_count)
        .alias("proband_alt_count")
    )

    # VAF = alt / (alt + ref-like depth) — use AD from VCF as primary
    # R pipeline: VAF = AD / DP (from VCF)
    df = df.with_columns(
        (pl.col("AD_alt").cast(pl.Float64) / pl.col("DP").cast(pl.Float64))
        .round(4)
        .alias("VAF")
    )

    # Classify
    df = df.with_columns(
        pl.when(
            (pl.col("VAF") >= VAF_MOSAIC_LOW) & (pl.col("VAF") < VAF_MOSAIC_HIGH)
        )
        .then(pl.lit("mosaic"))
        .when(
            (pl.col("VAF") >= VAF_MOSAIC_HIGH) & (pl.col("VAF") <= VAF_HET_HIGH)
        )
        .then(pl.lit("heterozygous"))
        .when(pl.col("VAF") < VAF_MOSAIC_LOW)
        .then(pl.lit("low_mosaic"))
        .otherwise(pl.lit("high_vaf"))
        .alias("VAF_class")
    )

    return df


def breakdown_by_vaf(df: pl.DataFrame, label: str) -> None:
    """Print VAF breakdown table like R's breakdownByVAF / getFinalBreakdown."""
    def count_bin(data, low, high):
        return len(data.filter((pl.col("VAF") >= low) & (pl.col("VAF") < high)))

    snv = df.filter(pl.col("Type") == "SNV")
    indel = df.filter(pl.col("Type") == "INDEL")

    bins = [
        ("VAF 9.5-18.7% (mosaic-low)", VAF_MOSAIC_LOW, VAF_MOSAIC_MID),
        ("VAF 18.7-37.4% (mosaic-high)", VAF_MOSAIC_MID, VAF_MOSAIC_HIGH),
        ("VAF 37.4-62.6% (het)", VAF_MOSAIC_HIGH, VAF_HET_HIGH),
    ]

    print(f"\n{'='*60}")
    print(f"  VAF Breakdown: {label}")
    print(f"{'='*60}")
    print(f"{'Bin':<30} {'SNV':>8} {'INDEL':>8} {'All':>8}")
    print(f"{'-'*30} {'-'*8} {'-'*8} {'-'*8}")
    for name, lo, hi in bins:
        s = count_bin(snv, lo, hi)
        i = count_bin(indel, lo, hi)
        a = count_bin(df, lo, hi)
        print(f"{name:<30} {s:>8,} {i:>8,} {a:>8,}")
    print(f"{'-'*30} {'-'*8} {'-'*8} {'-'*8}")
    # Totals for mosaic range
    ms = count_bin(snv, VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    mi = count_bin(indel, VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    ma = count_bin(df, VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    print(f"{'Total mosaic':<30} {ms:>8,} {mi:>8,} {ma:>8,}")
    hs = count_bin(snv, VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    hi_ = count_bin(indel, VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    ha = count_bin(df, VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    print(f"{'Total het':<30} {hs:>8,} {hi_:>8,} {ha:>8,}")


def print_filter_stats(df: pl.DataFrame) -> None:
    """Print filtering statistics at each step."""
    # Filter columns are boolean (True/False) at this point
    total = len(df)
    pileup_pass = len(df.filter(pl.col("pileupFullFilter")))
    hard_pass = len(df.filter(
        pl.col("pileupFullFilter") & pl.col("HardFilter")
    ))
    centromere_pass = len(df.filter(
        pl.col("pileupFullFilter") & pl.col("HardFilter") & pl.col("centromereFilter")
    ))
    repeat_pass = len(df.filter(
        pl.col("pileupFullFilter") & pl.col("HardFilter") & pl.col("repeatFullFilter")
    ))

    print(f"\n--- Filter Statistics ---")
    print(f"Input (from Part 1):        {total:>10,}")
    print(f"After pileup filter:        {pileup_pass:>10,}  ({total - pileup_pass:,} removed)")
    print(f"After hard filter:          {hard_pass:>10,}  ({pileup_pass - hard_pass:,} removed)")
    print(f"After centromere filter:    {centromere_pass:>10,}  ({hard_pass - centromere_pass:,} removed)")
    print(f"After repeat+centr filter:  {repeat_pass:>10,}  ({centromere_pass - repeat_pass:,} removed)")

    # Sub-breakdown of pileup filter
    depth_pass = len(df.filter(pl.col("pileupDepthFilterTrio")))
    parent_pass = len(df.filter(
        pl.col("pileupDepthFilterTrio") & pl.col("pileupParentFilter")
    ))
    print(f"\n  Pileup sub-filters:")
    print(f"    Trio depth filter:      {depth_pass:>10,}  ({total - depth_pass:,} removed)")
    print(f"    + Parent alt=0:         {parent_pass:>10,}  ({depth_pass - parent_pass:,} removed)")

    # Type breakdown
    snv_total = len(df.filter(pl.col("Type") == "SNV"))
    indel_total = len(df.filter(pl.col("Type") == "INDEL"))
    print(f"\n  Variant types: {snv_total:,} SNV, {indel_total:,} INDEL")


def run_vep_gnomad(df: pl.DataFrame, output_dir: Path, sample: str) -> pl.DataFrame:
    """
    Run VEP via Docker to get gnomAD AF, then annotate the DataFrame.

    Exports filtered variants to VCF → runs VEP with --af_gnomadg --pick →
    parses output → joins gnomADg_AF back to DataFrame.
    """
    if not VEP_DATA.exists():
        print("WARNING: VEP cache not found, skipping gnomAD filter")
        return df.with_columns(pl.lit(None).cast(pl.Float64).alias("gnomADg_AF"))

    with tempfile.TemporaryDirectory(prefix="vep_") as tmpdir:
        tmpdir = Path(tmpdir)
        input_vcf = tmpdir / f"{sample}_for_vep.vcf"
        output_tsv = tmpdir / f"{sample}_vep.tsv"

        # Write minimal VCF
        with open(input_vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for row in df.select("CHROM", "POS", "REF", "ALT").iter_rows():
                chrom, pos, ref, alt = row
                vid = f"{chrom}:{pos}_{ref}>{alt}"
                f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t.\t.\t.\n")

        print(f"Running VEP on {len(df):,} variants...")

        # Run VEP in Docker
        cmd = [
            "docker", "run", "--rm",
            "--user", f"{subprocess.getoutput('id -u')}:{subprocess.getoutput('id -g')}",
            "-v", f"{VEP_DATA}:/.vep",
            "-v", f"{tmpdir}:/work",
            VEP_IMAGE, "vep",
            "--cache", "--offline", "--assembly", "GRCh38",
            "--dir_cache", "/.vep",
            "--input_file", f"/work/{input_vcf.name}",
            "--output_file", f"/work/{output_tsv.name}",
            "--tab", "--no_stats", "--quiet", "--force",
            "--fork", str(VEP_FORKS),
            "--af_gnomadg",
            "--pick",
            "--pick_order", "biotype,rank,mane_select,tsl,canonical,appris,ccds,length",
            "--fields", "Uploaded_variation,Location,Allele,gnomADg_AF",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"VEP ERROR: {result.stderr}")
            return df.with_columns(pl.lit(None).cast(pl.Float64).alias("gnomADg_AF"))

        # Parse VEP output
        vep_df = pl.read_csv(
            str(output_tsv),
            separator="\t",
            comment_prefix="##",
            infer_schema_length=5000,
        )
        # VEP output: #Uploaded_variation, Location, Allele, gnomADg_AF
        # Uploaded_variation = "chr1:904493_C>T"
        vep_df = vep_df.rename({vep_df.columns[0]: "variant_key"})

        # Parse gnomADg_AF: "-" means missing/absent
        vep_df = vep_df.with_columns(
            pl.when(pl.col("gnomADg_AF") == "-")
            .then(pl.lit(None))
            .otherwise(pl.col("gnomADg_AF"))
            .cast(pl.Float64)
            .alias("gnomADg_AF")
        )

        # Deduplicate (--pick should give one per variant, but just in case)
        vep_df = vep_df.group_by("variant_key").agg(
            pl.col("gnomADg_AF").first()
        )

        print(f"VEP returned {len(vep_df):,} annotated variants")

    # Create matching key in our DataFrame
    df = df.with_columns(
        (pl.col("CHROM") + ":" + pl.col("POS").cast(pl.String)
         + "_" + pl.col("REF") + ">" + pl.col("ALT"))
        .alias("variant_key")
    )

    # Join
    df = df.join(vep_df.select("variant_key", "gnomADg_AF"), on="variant_key", how="left")
    df = df.drop("variant_key")

    # Stats
    has_af = df.filter(pl.col("gnomADg_AF").is_not_null())
    af_gt0 = has_af.filter(pl.col("gnomADg_AF") > 0)
    print(f"gnomAD: {len(has_af):,} with AF, {len(af_gt0):,} with AF>0 (to be removed)")

    return df


def main():
    parser = argparse.ArgumentParser(description="MosaicPostProcessPart2")
    parser.add_argument("input_csv", help="Part 1 pileup CSV")
    parser.add_argument("--sample", required=True, help="Sample name (brca1, brca2, etc.)")
    parser.add_argument("--output-dir", default=None, help="Output directory (default: outputs/part2/)")
    parser.add_argument("--no-vep", action="store_true", help="Skip VEP/gnomAD annotation")
    parser.add_argument("--caller", default="hc", choices=["hc", "mutect2"],
                        help="Variant caller: hc (HaplotypeCaller) or mutect2 (default: hc)")
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else Path(__file__).parent.parent / "outputs" / "part2"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Load
    df = load_part1_csv(args.input_csv, caller=args.caller)

    # Step 2: Classify SNV vs INDEL
    df = classify_variant_type(df)

    # Step 3: Pileup filters
    df = add_pileup_filter(df)

    # Step 4: Hard filters
    df = add_hard_filter(df, caller=args.caller)

    # Step 5: Repeat/centromere filter
    print("Running repeat/centromere overlap filters...")
    df = add_repeat_filter(df)

    # Step 5b: chrY filter — exclude chrY (haploid in males, HC ploidy=2 unreliable)
    # chrY:11.3Mb region has X-Y homology causing false de novo calls in all trios
    print("Applying sex chromosome filter (exclude chrY)...")
    df = df.with_columns(
        (pl.col("CHROM") != "chrY").alias("sexChrFilter")
    )

    # Step 6: VAF
    df = add_vaf_and_classify(df)

    # Step 7: VEP gnomAD annotation on ALL variants
    if not args.no_vep:
        df = run_vep_gnomad(df, output_dir, args.sample)
        df = df.with_columns(
            (pl.col("gnomADg_AF").is_null() | (pl.col("gnomADg_AF") == 0))
            .alias("gnomAD_filter")
        )
    else:
        print("gnomAD filter: SKIPPED (--no-vep)")
        df = df.with_columns(pl.lit(True).alias("gnomAD_filter"))

    # Convert all PASS/FAIL filter columns to boolean True/False
    filter_cols = [
        "pileupDepthFilterTrio", "pileupParentFilter",
        "pileupProbandFilter", "pileupFullFilter",
        "HardFilter", "centromereFilter", "repeatFilter", "repeatFullFilter",
    ]
    df = df.with_columns([
        (pl.col(c) == "PASS").alias(c) for c in filter_cols
    ])

    # Step 8: deepSNV-style beta-binomial test (proband vs pooled parents)
    # BH correction applied only to pileup-passing variants (true de novo candidates)
    df = add_deepsnv_pvalue(df)

    # Print stats
    print_filter_stats(df)

    # Composite final filter (matches R pipeline logic + chrY exclusion)
    df = df.with_columns(
        (
            pl.col("pileupFullFilter")
            & pl.col("HardFilter")
            & pl.col("centromereFilter")
            & pl.col("repeatFullFilter")
            & pl.col("gnomAD_filter")
            & pl.col("sexChrFilter")
        ).alias("PASS_all")
    )

    final = df.filter(pl.col("PASS_all"))
    print(f"\nFinal candidates (all filters): {len(final):,}")

    # VAF breakdown
    breakdown_by_vaf(final, f"{args.sample} — final candidates")

    # Save outputs
    # Full annotated file — ALL variants with all filter columns (boolean)
    all_output = output_dir / f"{args.sample}_part2_all.tsv"
    df.write_csv(str(all_output), separator="\t")
    print(f"\nAll variants with filters → {all_output}  ({len(df):,} rows)")

    # Final filtered candidates
    final_output = output_dir / f"{args.sample}_part2_final.tsv"
    final.write_csv(str(final_output), separator="\t")
    print(f"Final candidates → {final_output}  ({len(final):,} rows)")

    # Mosaic-only subset
    mosaic = final.filter(pl.col("VAF_class") == "mosaic")
    mosaic_output = output_dir / f"{args.sample}_mosaic.tsv"
    mosaic.write_csv(str(mosaic_output), separator="\t")
    print(f"Mosaic variants → {mosaic_output} ({len(mosaic):,} variants)")

    # Het-only subset
    het = final.filter(pl.col("VAF_class") == "heterozygous")
    het_output = output_dir / f"{args.sample}_het.tsv"
    het.write_csv(str(het_output), separator="\t")
    print(f"Het variants → {het_output} ({len(het):,} variants)")

    # Summary table
    summary = generate_summary(df, final, args.sample)
    summary_output = output_dir / f"{args.sample}_summary.txt"
    with open(summary_output, "w") as f:
        f.write(summary)
    print(f"\nSummary → {summary_output}")
    print(summary)


def generate_summary(df: pl.DataFrame, final: pl.DataFrame, sample: str) -> str:
    """Generate a comprehensive summary table."""
    total = len(df)
    snv_total = len(df.filter(pl.col("Type") == "SNV"))
    indel_total = len(df.filter(pl.col("Type") == "INDEL"))

    # Cumulative filter cascade
    depth = len(df.filter(pl.col("pileupDepthFilterTrio")))
    parent = len(df.filter(pl.col("pileupDepthFilterTrio") & pl.col("pileupParentFilter")))
    pileup = len(df.filter(pl.col("pileupFullFilter")))
    hard = len(df.filter(pl.col("pileupFullFilter") & pl.col("HardFilter")))
    centr = len(df.filter(
        pl.col("pileupFullFilter") & pl.col("HardFilter") & pl.col("centromereFilter")
    ))
    repeat = len(df.filter(
        pl.col("pileupFullFilter") & pl.col("HardFilter") & pl.col("repeatFullFilter")
    ))
    gnomad = len(df.filter(
        pl.col("pileupFullFilter") & pl.col("HardFilter")
        & pl.col("repeatFullFilter") & pl.col("gnomAD_filter")
    ))
    sexchr = len(df.filter(
        pl.col("pileupFullFilter") & pl.col("HardFilter")
        & pl.col("repeatFullFilter") & pl.col("gnomAD_filter")
        & pl.col("sexChrFilter")
    ))
    final_n = len(final)

    # VAF bins on final
    def _vaf(data, lo, hi):
        return len(data.filter((pl.col("VAF") >= lo) & (pl.col("VAF") < hi)))

    snv_f = final.filter(pl.col("Type") == "SNV")
    indel_f = final.filter(pl.col("Type") == "INDEL")

    lines = []
    lines.append(f"\n{'='*70}")
    lines.append(f"  SUMMARY: {sample}")
    lines.append(f"{'='*70}")
    lines.append(f"")
    lines.append(f"  Filter cascade:")
    lines.append(f"  {'Step':<40} {'Count':>10}")
    lines.append(f"  {'-'*40} {'-'*10}")
    lines.append(f"  {'Input (Part 1 de novo candidates)':<40} {total:>10,}")
    lines.append(f"    {'SNV':<38} {snv_total:>10,}")
    lines.append(f"    {'INDEL':<38} {indel_total:>10,}")
    lines.append(f"  {'Pileup: trio depth':<40} {depth:>10,}")
    lines.append(f"  {'Pileup: + parent alt=0':<40} {parent:>10,}")
    lines.append(f"  {'Pileup: full (+ proband alt≥2, AD≥2)':<40} {pileup:>10,}")
    lines.append(f"  {'+ Hard filter (MQ,QD,FS,RPRS)':<40} {hard:>10,}")
    lines.append(f"  {'+ Centromere exclusion':<40} {centr:>10,}")
    lines.append(f"  {'+ Repeat (Alu+centr) exclusion':<40} {repeat:>10,}")
    lines.append(f"  {'+ gnomAD AF=0':<40} {gnomad:>10,}")
    lines.append(f"  {'+ chrY exclusion':<40} {sexchr:>10,}")
    lines.append(f"  {'-'*40} {'-'*10}")
    lines.append(f"  {'FINAL CANDIDATES':<40} {final_n:>10,}")
    lines.append(f"")
    lines.append(f"  VAF breakdown (final):")
    lines.append(f"  {'Bin':<30} {'SNV':>8} {'INDEL':>8} {'All':>8}")
    lines.append(f"  {'-'*30} {'-'*8} {'-'*8} {'-'*8}")
    bins = [
        ("Mosaic 9.5-18.7%", VAF_MOSAIC_LOW, VAF_MOSAIC_MID),
        ("Mosaic 18.7-37.4%", VAF_MOSAIC_MID, VAF_MOSAIC_HIGH),
        ("Het 37.4-62.6%", VAF_MOSAIC_HIGH, VAF_HET_HIGH),
    ]
    for name, lo, hi in bins:
        s = _vaf(snv_f, lo, hi)
        i = _vaf(indel_f, lo, hi)
        a = _vaf(final, lo, hi)
        lines.append(f"  {name:<30} {s:>8,} {i:>8,} {a:>8,}")
    lines.append(f"  {'-'*30} {'-'*8} {'-'*8} {'-'*8}")
    ms = _vaf(snv_f, VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    mi = _vaf(indel_f, VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    ma = _vaf(final, VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    lines.append(f"  {'TOTAL MOSAIC':<30} {ms:>8,} {mi:>8,} {ma:>8,}")
    hs = _vaf(snv_f, VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    hi_ = _vaf(indel_f, VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    ha = _vaf(final, VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    lines.append(f"  {'TOTAL HET':<30} {hs:>8,} {hi_:>8,} {ha:>8,}")
    lines.append(f"{'='*70}")
    lines.append(f"")
    return "\n".join(lines)


if __name__ == "__main__":
    main()
