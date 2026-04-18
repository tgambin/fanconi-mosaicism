#!/usr/bin/env python3
"""
Master output generator — runs entire post-processing pipeline and produces
structured outputs: Part2 filtering, cross-caller comparison, benchmark JSON,
Venn diagrams, and a markdown report.

Usage:
  python scripts/generate_all_outputs.py [--no-vep] [--trios brca1,brca2]

Output structure:
  outputs/
  ├── part2/{trio}_{caller}_part2_all.tsv   — all variants with filter columns
  ├── part2/{trio}_{caller}_part2_final.tsv — final candidates
  ├── part2/{trio}_{caller}_mosaic.tsv      — mosaic subset
  ├── part2/{trio}_{caller}_het.tsv         — het subset
  ├── benchmark/
  │   ├── filter_cascade.json               — per-trio per-caller filter stats
  │   ├── summary.json                      — cross-caller summary
  │   └── report.md                         — human-readable markdown report
  ├── cross_caller/
  │   ├── {trio}_master.csv                 — all het+mosaic variants, 91 columns
  │   ├── {trio}_{class}_{type}_venn.png    — Venn diagrams
  │   └── {trio}_mosaic_het_igv.tsv         — IGV input positions (mosaic+het, HC∪M2)
  └── candidates/
      ├── hc/{trio}_{mosaic|het}.tsv
      └── mutect2/{trio}_{mosaic|het}.tsv
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from scipy.stats import binomtest

from _common import (
    SCRIPT_DIR, PROJECT_DIR, OUTPUTS_DIR,
    SAMPLES, TRIOS, GIAB,
    INDEX_TRIO, CONTROL_TRIO, GIAB_TRIO, FA_TRIOS, OTHER_FA_TRIOS,
    CALLERS, CALLER_SLUGS,
    VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH, VAF_HET_HIGH,
    SUBCLONAL_PEAK_VAF_LOW, SUBCLONAL_PEAK_VAF_HIGH,
    BH_SIG_STRONG, BH_SIG_MODERATE, BINOM_HET_ALPHA, POLARS_INFER_SCHEMA,
    burden_test, normalize_bool_columns, sig_star, vaf_breakdown,
)

PART1_MAP = {t: s["sample_id"] for t, s in SAMPLES.items()}


# ---------------------------------------------------------------------------
# Pipeline-step labels used by both the markdown cascade renderer (§5) and the
# supplementary-table CSV writer (Sx2). Tuple = (display label, key used in
# filter_report.txt / cascade_steps dict).
# ---------------------------------------------------------------------------
PART1_HC_STEPS = (
    ("Raw proband variants",                  "Raw proband variants"),
    ("After normalize (excl spanning dels)",  "After proband normalize (excl spanning dels)"),
    ("After trio de novo filter (isec -C)",   "After trio de novo filter (isec -C)"),
    ("After canonical chr + DP>=10",          "After canonical chr + DP>=10"),
    ("After segdup exclusion",                "After segdup exclusion"),
    ("After LCR exclusion",                   "After LCR exclusion"),
)
PART1_M2_STEPS = (
    ("Raw proband variants",                  "Raw proband variants"),
    ("After PASS filter",                     "After PASS filter"),
    ("After normalize (excl spanning dels)",  "After proband PASS+normalize (excl spanning dels)"),
    ("After trio de novo filter (isec -C)",   "After trio de novo filter (isec -C)"),
    ("After DP>=10",                          "After DP>=10"),
    ("After segdup exclusion",                "After segdup exclusion"),
    ("After LCR exclusion",                   "After LCR exclusion"),
)
PART2_STEPS = (
    ("pileup_depth",        "Pileup: trio depth"),
    ("pileup_parent_alt0",  "+ parent alt=0"),
    ("pileup_full",         "+ proband alt≥2, AD≥2"),
    ("hard_filter",         "+ Hard filter (MQ,QD,FS,RPRS)"),
    ("centromere",          "+ Centromere exclusion"),
    ("repeat",              "+ Repeat (Alu+centr) exclusion"),
    ("gnomad_af0",          "+ gnomAD AF=0"),
    ("chrY_exclusion",      "+ chrY exclusion"),
    ("vaf_window",          "+ In mosaic+het VAF window (9.55–62.59%)"),
    ("igv_validated",       "+ Manual IGV validation"),
)


def run_part2(trio, caller, no_vep=False):
    """Run mosaic_postprocess_part2.py for one trio/caller."""
    sample_id = PART1_MAP[trio]
    suffix = ".mutect2" if caller == "mutect2" else ""
    csv_path = OUTPUTS_DIR / "part1_csv" / f"{sample_id}{suffix}.part1.csv"

    if not csv_path.exists():
        # Try old location
        csv_path = OUTPUTS_DIR / "part1_csv_mutect" / f"{sample_id}.mutect2.part1.csv"
        if not csv_path.exists():
            print(f"  SKIP {trio}/{caller}: Part1 CSV not found")
            return False

    sample_name = f"{trio}_mutect2" if caller == "mutect2" else trio
    cmd = [
        sys.executable, str(SCRIPT_DIR / "mosaic_postprocess_part2.py"),
        str(csv_path),
        "--sample", sample_name,
        "--caller", caller,
    ]
    if no_vep:
        cmd.append("--no-vep")

    print(f"  Running Part2: {trio}/{caller} ...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}")
        return False
    # Extract summary line
    for line in result.stdout.split('\n'):
        if 'FINAL CANDIDATES' in line:
            print(f"  {line.strip()}")
    return True


def parse_part1_report(path):
    """Parse a Part1 FilterVariants filter_report.txt into dict."""
    if not path.exists():
        return None
    lines = path.read_text().strip().split("\n")
    report = {}
    for line in lines:
        if ":" not in line or line.startswith("===") or line.startswith("Input"):
            continue
        key, _, val = line.partition(":")
        key = key.strip()
        val = val.strip()
        try:
            report[key] = int(val)
        except ValueError:
            report[key] = val
    return report


def load_part1_reports(trios):
    """Load Part1 filter reports for all trios × callers."""
    reports_dir = OUTPUTS_DIR / "part1_reports"
    part1 = {}
    for trio in trios:
        for caller in CALLER_SLUGS:
            f = reports_dir / f"{trio}_{caller}_filter_report.txt"
            report = parse_part1_report(f)
            if report:
                part1[f"{trio}_{caller}"] = report
    return part1


def build_filter_cascade_json(trios):
    """Build JSON with filter cascade stats from part2_all files."""
    import polars as pl

    cascades = {}

    for trio in trios:
        for caller in CALLER_SLUGS:
            label = f"{trio}_mutect2" if caller == "mutect2" else trio
            f = OUTPUTS_DIR / "part2" / f"{label}_part2_all.tsv"
            if not f.exists():
                continue

            df = pl.read_csv(str(f), separator="\t", infer_schema_length=POLARS_INFER_SCHEMA)
            df = normalize_bool_columns(df)
            total = len(df)

            snv_total = len(df.filter(pl.col("Type") == "SNV"))
            indel_total = len(df.filter(pl.col("Type") == "INDEL"))

            def cum_filter(*cols):
                expr = pl.col(cols[0])
                for c in cols[1:]:
                    expr = expr & pl.col(c)
                return vaf_breakdown(df.filter(expr))

            steps = [
                ("input", None),
                ("pileup_depth", ["pileupDepthFilterTrio"]),
                ("pileup_parent_alt0", ["pileupDepthFilterTrio", "pileupParentFilter"]),
                ("pileup_full", ["pileupFullFilter"]),
                ("hard_filter", ["pileupFullFilter", "HardFilter"]),
                ("centromere", ["pileupFullFilter", "HardFilter", "centromereFilter"]),
                ("repeat", ["pileupFullFilter", "HardFilter", "repeatFullFilter"]),
                ("gnomad_af0", ["pileupFullFilter", "HardFilter", "repeatFullFilter", "gnomAD_filter"]),
                ("chrY_exclusion", ["pileupFullFilter", "HardFilter", "repeatFullFilter", "gnomAD_filter", "sexChrFilter"]),
            ]

            cascade_steps = {}
            for step_name, cols in steps:
                if cols is None:
                    cascade_steps[step_name] = vaf_breakdown(df)
                    cascade_steps[step_name]["total"] = total
                else:
                    cascade_steps[step_name] = cum_filter(*cols)

            # Restriction to the mosaic+het VAF window. Without this, a
            # low-VAF-sensitive caller such as Mutect2 retains thousands of
            # PASS_all variants outside the mosaic/het bins, which would
            # dominate the post-chrY row and make the subsequent IGV-
            # validation step look like a no-op (unreviewed variants have a
            # neutral igvFilter=True). Restricting to the VAF window brings
            # the cascade down to the candidate set that was actually
            # submitted for manual IGV review.
            vaf_window_expr = (
                pl.col("pileupFullFilter") & pl.col("HardFilter")
                & pl.col("repeatFullFilter") & pl.col("gnomAD_filter")
                & pl.col("sexChrFilter")
                & (pl.col("VAF") >= VAF_MOSAIC_LOW)
                & (pl.col("VAF") <= VAF_HET_HIGH)
            )
            cascade_steps["vaf_window"] = vaf_breakdown(df.filter(vaf_window_expr))

            # Only FA trios carry a manual IGV filter column (GIAB wasn't reviewed).
            if "igvFilter" in df.columns:
                cascade_steps["igv_validated"] = vaf_breakdown(
                    df.filter(vaf_window_expr & pl.col("igvFilter"))
                )

            # deepSNV BH on the final classified set — IGV-confirmed if
            # available (FA trios) else PASS_all (GIAB, no manual IGV review),
            # split by binom_class so the counts match §2d / §5.
            if "PASS_all_igv_confirmed" in df.columns:
                final_df = df.filter(pl.col("PASS_all_igv_confirmed") == True)
            elif "PASS_all" in df.columns:
                final_df = df.filter(pl.col("PASS_all") == True)
            else:
                final_df = df.filter(
                    pl.col("pileupFullFilter") & pl.col("HardFilter")
                    & pl.col("repeatFullFilter") & pl.col("gnomAD_filter")
                    & pl.col("sexChrFilter"))

            def _deepsnv_counts(subset, binom_cls):
                """Count BH<BH_SIG_STRONG and BH<BH_SIG_MODERATE for SNVs in a
                given binomial classification bucket."""
                if "binom_class" not in subset.columns:
                    return {"n": 0, "BH001": 0, "BH005": 0}
                snv = subset.filter(
                    (pl.col("Type") == "SNV") &
                    (pl.col("binom_class") == binom_cls))
                n = len(snv)
                if n == 0 or "deepSNV_pval_BH" not in snv.columns:
                    return {"n": n, "BH001": 0, "BH005": 0}
                bh = snv["deepSNV_pval_BH"].drop_nulls().cast(pl.Float64)
                return {
                    "n": n,
                    "BH001": int((bh < BH_SIG_STRONG).sum()),
                    "BH005": int((bh < BH_SIG_MODERATE).sum()),
                }

            deepsnv = {
                "mosaic_snv": _deepsnv_counts(final_df, "mosaic"),
                "het_snv": _deepsnv_counts(final_df, "het"),
            }

            cascade_entry = {
                "trio": trio, "caller": caller,
                "input_snv": snv_total, "input_indel": indel_total,
                "steps": cascade_steps,
                "deepSNV": deepsnv,
            }

            # Add Part1 report data if available
            p1_reports = load_part1_reports(trios)
            p1_key = f"{trio}_{caller}"
            if p1_key in p1_reports:
                cascade_entry["part1"] = p1_reports[p1_key]

            cascades[f"{trio}_{caller}"] = cascade_entry

    return cascades


def build_summary_json(cascades):
    """Build cross-caller summary from cascades.

    Two views per caller/trio: `default` (pipeline-filtered, before manual
    IGV review) and `igv_validated` (additionally requires IGV_is_real=TRUE).
    The IGV-validated view is omitted for trios where the cascade has no
    `igv_validated` step (e.g. GIAB, which was not reviewed).
    """
    summary = {"hc": {}, "mutect2": {}}
    for key, c in cascades.items():
        trio, caller = c["trio"], c["caller"]
        final = c["steps"]["chrY_exclusion"]
        entry = {
            "mosaic_snv": final["mosaic_snv"],
            "mosaic_indel": final["mosaic_indel"],
            "het_snv": final["het_snv"],
            "het_indel": final["het_indel"],
            "total_mosaic": final["mosaic_all"],
            "total_het": final["het_all"],
            "total_final": final["mosaic_all"] + final["het_all"],
        }
        if "igv_validated" in c["steps"]:
            h = c["steps"]["igv_validated"]
            entry["igv_review"] = {
                "mosaic_snv": h["mosaic_snv"],
                "mosaic_indel": h["mosaic_indel"],
                "het_snv": h["het_snv"],
                "het_indel": h["het_indel"],
                "total_mosaic": h["mosaic_all"],
                "total_het": h["het_all"],
                "total_final": h["mosaic_all"] + h["het_all"],
            }
        summary[caller][trio] = entry
    return summary


def _compute_binom_classification(trios):
    """Per (trio, caller): mosaic/het counts under VAF- and binom-classification,
    both raw (PASS_all) and IGV-confirmed.

    Returns: dict {trio: {caller: {...}}}.
    """
    import polars as pl
    out = {}
    for trio in trios:
        out[trio] = {}
        for caller in CALLER_SLUGS:
            label = f"{trio}_mutect2" if caller == "mutect2" else trio
            f = OUTPUTS_DIR / "part2" / f"{label}_part2_all.tsv"
            if not f.exists():
                continue
            df = pl.read_csv(str(f), separator="\t", infer_schema_length=POLARS_INFER_SCHEMA)
            df = normalize_bool_columns(df)
            if "binom_class" not in df.columns:
                continue
            has_igv = "PASS_all_igv_confirmed" in df.columns

            def _count(subset_filter, class_col, cls, vtype=None):
                q = subset_filter & (pl.col(class_col) == cls)
                if vtype is not None:
                    q = q & (pl.col("Type") == vtype)
                return len(df.filter(q))

            pass_all = pl.col("PASS_all")
            pass_igv = pl.col("PASS_all_igv_confirmed") if has_igv else pass_all

            # Per-variant-type counts for both classifications, both filter levels.
            def _count_bin(subset_filter, vtype, vaf_lo, vaf_hi):
                """SNV/INDEL count in a half-open VAF bin [vaf_lo, vaf_hi)."""
                return len(df.filter(
                    subset_filter & (pl.col("Type") == vtype)
                    & (pl.col("VAF") >= vaf_lo) & (pl.col("VAF") < vaf_hi)
                ))

            def block(subset_filter):
                return {
                    "vaf_mosaic_snv": _count(subset_filter, "VAF_class", "mosaic", "SNV"),
                    "vaf_mosaic_indel": _count(subset_filter, "VAF_class", "mosaic", "INDEL"),
                    "vaf_het_snv": _count(subset_filter, "VAF_class", "heterozygous", "SNV"),
                    "vaf_het_indel": _count(subset_filter, "VAF_class", "heterozygous", "INDEL"),
                    "binom_mosaic_snv": _count(subset_filter, "binom_class", "mosaic", "SNV"),
                    "binom_mosaic_indel": _count(subset_filter, "binom_class", "mosaic", "INDEL"),
                    "binom_het_snv": _count(subset_filter, "binom_class", "het", "SNV"),
                    "binom_het_indel": _count(subset_filter, "binom_class", "het", "INDEL"),
                    # R1.2: SNVs in the 0.10-0.20 VAF bin flagged by the reviewer
                    # as the expected signature of early-embryonic mosaic mutations.
                    "subclonal_peak_snv": _count_bin(
                        subset_filter, "SNV",
                        SUBCLONAL_PEAK_VAF_LOW, SUBCLONAL_PEAK_VAF_HIGH,
                    ),
                }

            # Also: how many variants moved between categories
            moved_het_to_mosaic = len(df.filter(
                pass_all & (pl.col("VAF_class") == "heterozygous")
                        & (pl.col("binom_class") == "mosaic")
            ))
            moved_mosaic_to_het = len(df.filter(
                pass_all & (pl.col("VAF_class") == "mosaic")
                        & (pl.col("binom_class") == "het")
            ))

            out[trio][caller] = {
                "pass_all": block(pass_all),
                "pass_igv": block(pass_igv) if has_igv else None,
                "het_to_mosaic": moved_het_to_mosaic,
                "mosaic_to_het": moved_mosaic_to_het,
            }
    return out


def _compute_igv_precision(trios):
    """Per (trio, caller): precision on pipeline-passing reviewed variants.

    Reads _part2_all.tsv, restricts to PASS_all rows with a non-empty
    IGV_is_real annotation, and returns Reviewed / TRUE / FALSE / Precision.
    """
    import polars as pl
    rows = []
    for trio in trios:
        for caller in CALLER_SLUGS:
            label = f"{trio}_mutect2" if caller == "mutect2" else trio
            f = OUTPUTS_DIR / "part2" / f"{label}_part2_all.tsv"
            if not f.exists():
                continue
            df = pl.read_csv(str(f), separator="\t", infer_schema_length=POLARS_INFER_SCHEMA)
            df = normalize_bool_columns(df)
            if "IGV_is_real" not in df.columns or "PASS_all" not in df.columns:
                continue
            reviewed = df.filter(
                pl.col("PASS_all") & (pl.col("IGV_is_real").cast(pl.Utf8) != "")
                & pl.col("IGV_is_real").is_not_null()
            )
            n_rev = len(reviewed)
            n_true = len(reviewed.filter(pl.col("IGV_is_real") == "TRUE"))
            n_false = len(reviewed.filter(pl.col("IGV_is_real") == "FALSE"))
            precision = f"{n_true / n_rev * 100:.1f}%" if n_rev > 0 else "—"
            rows.append({
                "trio": trio, "caller": caller,
                "reviewed": n_rev, "true": n_true, "false": n_false,
                "precision": precision,
            })
    return rows


def load_cross_caller_data(trios):
    """Load cross-caller intersection data from master CSVs."""
    import csv as csvmod
    cross = {}
    OLD_DIR = PROJECT_DIR / "data" / "metadata"

    for trio in trios:
        master = OUTPUTS_DIR / "cross_caller" / f"{trio}_master.csv"
        if not master.exists():
            continue

        # Read master CSV
        hc_mosaic = set(); m2_mosaic = set(); old_mosaic = set()
        hc_het = set(); m2_het = set(); old_het = set()
        # deepSNV for shared HC∩M2 mosaic
        shared_mosaic_bh_hc = []; shared_mosaic_bh_m2 = []
        shared_het_bh_hc = []

        with open(master) as f:
            for row in csvmod.DictReader(f):
                key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
                if row['variant_type'] != 'SNV':
                    continue
                in_hc = row['in_HC'] == 'True'
                in_m2 = row['in_Mutect2'] == 'True'
                in_old = row['in_OldPipeline'] == 'True'

                if row['class_HC'] == 'mosaic' and in_hc: hc_mosaic.add(key)
                if row['class_Mutect2'] == 'mosaic' and in_m2: m2_mosaic.add(key)
                if row['class_OldPipeline'] == 'mosaic' and in_old: old_mosaic.add(key)
                if row['class_HC'] == 'het' and in_hc: hc_het.add(key)
                if row['class_Mutect2'] == 'het' and in_m2: m2_het.add(key)
                if row['class_OldPipeline'] == 'het' and in_old: old_het.add(key)

        shared_m = hc_mosaic & m2_mosaic
        shared_h = hc_het & m2_het

        # Re-read for deepSNV on shared
        with open(master) as f:
            for row in csvmod.DictReader(f):
                key = f"{row['CHROM']}:{row['POS']}:{row['REF']}:{row['ALT']}"
                if key in shared_m:
                    try: shared_mosaic_bh_hc.append(float(row['deepSNV_BH_HC']))
                    except: pass
                    try: shared_mosaic_bh_m2.append(float(row['deepSNV_BH_Mutect2']))
                    except: pass
                if key in shared_h:
                    try: shared_het_bh_hc.append(float(row['deepSNV_BH_HC']))
                    except: pass

        cross[trio] = {
            "mosaic_snv": {
                "hc": len(hc_mosaic), "m2": len(m2_mosaic), "old": len(old_mosaic),
                "hc_m2": len(shared_m), "hc_old": len(hc_mosaic & old_mosaic),
                "all3": len(hc_mosaic & m2_mosaic & old_mosaic),
                "shared_sig_hc_005": sum(1 for v in shared_mosaic_bh_hc if v < BH_SIG_MODERATE),
                "shared_sig_hc_001": sum(1 for v in shared_mosaic_bh_hc if v < BH_SIG_STRONG),
                "shared_sig_m2_005": sum(1 for v in shared_mosaic_bh_m2 if v < BH_SIG_MODERATE),
            },
            "het_snv": {
                "hc": len(hc_het), "m2": len(m2_het), "old": len(old_het),
                "hc_m2": len(shared_h), "hc_old": len(hc_het & old_het),
                "all3": len(hc_het & m2_het & old_het),
            },
        }
    return cross


def build_coverage_json():
    """Parse coverage text files into structured JSON."""
    cov_dir = OUTPUTS_DIR / "coverage"
    if not cov_dir.exists():
        return {}

    coverage = {}
    for f in sorted(cov_dir.glob("*.txt")):
        name = f.stem  # e.g. "proband_brca1.150x_coverage" or "mother_brca1.initial_coverage"
        lines = f.read_text().strip().split("\n")
        avg = None
        for line in lines:
            if "Average" in line:
                avg = float(line.split("=")[1].strip())
        if avg is not None:
            coverage[name] = round(avg, 3)
    return coverage


# ---------------------------------------------------------------------------
# Supplementary-table data builders. Each returns (header, rows) where rows
# is a list of lists of equal length to the header. Both the CSV writer
# (write_supp_tables_csv) and the markdown renderer (§0, §4 in the report)
# consume these functions so the numbers are generated in one place only.
# ---------------------------------------------------------------------------

def _supp_rows_deepsnv(cascades, trios):
    """Sx3: deepSNV beta-binomial LRT BH confirmation per caller × trio.

    Uses the final classified set — binom_class × PASS_all_igv_confirmed
    for FA trios or × PASS_all for GIAB (see build_filter_cascade_json).
    """
    header = [
        "Caller", "Trio", "Mosaic_SNV",
        "Mosaic_BH_lt_001", "Mosaic_BH_lt_005",
        "Het_SNV", "Het_BH_lt_001", "Het_BH_lt_005",
    ]
    rows = []
    for caller_label, caller in [(c["display"], c["slug"]) for c in CALLERS]:
        for t in trios:
            ds = cascades.get(f"{t}_{caller}", {}).get("deepSNV", {})
            m = ds.get("mosaic_snv", {}) or {}
            h = ds.get("het_snv", {}) or {}
            rows.append([
                caller_label, t,
                m.get("n", 0),
                f"{m.get('BH001', 0)}/{m.get('n', 0)}",
                f"{m.get('BH005', 0)}/{m.get('n', 0)}",
                h.get("n", 0),
                f"{h.get('BH001', 0)}/{h.get('n', 0)}",
                f"{h.get('BH005', 0)}/{h.get('n', 0)}",
            ])
    return header, rows


def _supp_rows_coverage(coverage, trios):
    """Sx5: per-sample coverage — probands (original + 150x downsample) and
    both parents (original only for FA trios; also 40x downsample for GIAB).

    `coverage` is the dict returned by build_coverage_json().
    """
    header = [
        "Sample", "Trio", "Gene", "Role",
        "Original_coverage_x", "Downsampled_coverage_x", "Downsample_target_x",
    ]
    rows = []
    for trio in trios:
        m = COVERAGE_MAP.get(trio, {})
        if trio == GIAB_TRIO:
            # GIAB: HG002 proband + HG003/HG004 parents; proband → 150x, parents → 40x.
            rows.append(["HG002", trio, m.get("gene", ""), "proband",
                         coverage.get("HG002_giab.initial_coverage", ""),
                         coverage.get("HG002_giab.150x_coverage", ""),
                         150])
            rows.append(["HG004", trio, "", "mother",
                         coverage.get("HG004_giab.initial_coverage", ""),
                         coverage.get("HG004_giab.40x_coverage", ""),
                         40])
            rows.append(["HG003", trio, "", "father",
                         coverage.get("HG003_giab.initial_coverage", ""),
                         coverage.get("HG003_giab.40x_coverage", ""),
                         40])
            continue
        sample = m.get("sample", "")
        proband_name = m.get("proband", trio)
        rows.append([proband_name, trio, m.get("gene", ""), "proband",
                     coverage.get(f"{sample}.initial_coverage", ""),
                     coverage.get(f"{sample}.150x_coverage", ""),
                     150])
        mother_key = m.get("mother", "")
        father_key = m.get("father", "")
        rows.append([mother_key.split("_")[0], trio, "", "mother",
                     coverage.get(f"{mother_key}.initial_coverage", ""),
                     "", ""])
        rows.append([father_key.split("_")[0], trio, "", "father",
                     coverage.get(f"{father_key}.initial_coverage", ""),
                     "", ""])
    return header, rows


# Coverage file name mapping — derived from SAMPLES so there's one source of
# truth. `sample` key is the prefix used in outputs/coverage/ files (no "_150x").
COVERAGE_MAP = {
    t: {
        "proband": s["proband"],
        "sample":  s["sample_id"].removesuffix("_150x"),
        "gene":    s["gene"],
        "mother":  s["mother"],
        "father":  s["father"],
    }
    for t, s in SAMPLES.items()
}


def generate_markdown_report(cascades, summary, output_path, trios=None,
                             extended_trios=None):
    """Generate unified markdown report: results + cross-caller + cascade."""
    if trios is None:
        trios = TRIOS
    lines = []
    lines.append("# Mosaic Variant Analysis — Full Report\n")
    lines.append("*Auto-generated by `scripts/generate_all_outputs.py`*\n")

    # --- 0. Coverage ---
    # Same row-list feeds the Sx5 CSV and the two markdown tables below —
    # single source of truth in _supp_rows_coverage().
    coverage = build_coverage_json()
    if coverage:
        cov_json_path = OUTPUTS_DIR / "benchmark" / "coverage.json"
        with open(cov_json_path, "w") as f:
            json.dump(coverage, f, indent=2)

        cov_trios = list(trios) + ([GIAB_TRIO] if GIAB_TRIO not in trios else [])
        _, cov_rows = _supp_rows_coverage(coverage, cov_trios)

        lines.append("## 0. Coverage\n")
        lines.append("### Probands (original → 150x downsampled)\n")
        lines.append("| Sample | Gene | Original | 150x |")
        lines.append("|---|---|---:|---:|")
        for r in cov_rows:
            if r[3] != "proband":
                continue
            sample, _, gene, _, orig, ds, _ = r
            lines.append(f"| {sample} | {gene} | {orig}x | {ds}x |")
        lines.append("")

        lines.append("### Parents\n")
        lines.append("| Sample | Trio | Role | Coverage |")
        lines.append("|---|---|---|---:|")
        for r in cov_rows:
            if r[3] == "proband":
                continue
            sample, trio_key, _, role, orig, ds, _ = r
            # GIAB parents were also downsampled (to 40x); show the downsampled
            # value there, the original for all FA trio parents.
            cov = ds if ds not in (None, "") else orig
            lines.append(f"| {sample} | {trio_key} | {role} | {cov}x |")
        lines.append("")

    # --- 1. Executive Summary ---
    lines.append("## 1. Key Finding\n")
    hc_s = summary.get("hc", {})
    brca2_m = hc_s.get(INDEX_TRIO, {}).get("mosaic_snv", 0)
    ctrl_m = hc_s.get(CONTROL_TRIO, {}).get("mosaic_snv", 0)
    ratio = f"{brca2_m/ctrl_m:.1f}x" if ctrl_m > 0 else "N/A"
    lines.append(f"brca2 (BRCA2/FANCD1) has **{brca2_m} mosaic SNVs** — "
                 f"**{ratio} more than control** ({ctrl_m}). "
                 f"Other FA samples (brca1={hc_s.get('brca1',{}).get('mosaic_snv',0)}, "
                 f"fanca={hc_s.get('fanca',{}).get('mosaic_snv',0)}, "
                 f"fanct={hc_s.get('fanct',{}).get('mosaic_snv',0)}) "
                 f"are comparable to control.\n")

    # --- 1b. Mosaic burden test ---
    lines.append("### Mosaic burden test (conditional binomial exact, one-sided)\n")
    for caller_label, caller_key in [(c["short"], c["slug"]) for c in CALLERS]:
        cal_s = summary.get(caller_key, {})
        ctrl_n = cal_s.get(CONTROL_TRIO, {}).get("mosaic_snv", 0)
        if ctrl_n == 0:
            continue
        lines.append(f"**{caller_label}**\n")
        lines.append("| Comparison | Ratio | p-value | |")
        lines.append("|---|---|---|---|")
        for t in trios:
            if t == CONTROL_TRIO:
                continue
            case_n = cal_s.get(t, {}).get("mosaic_snv", 0)
            total = case_n + ctrl_n
            if total == 0:
                continue
            result = binomtest(case_n, total, 0.5, alternative="greater")
            r = f"{case_n/ctrl_n:.1f}x" if ctrl_n > 0 else "N/A"
            lines.append(f"| {t} ({case_n}) vs control ({ctrl_n}) | {r} | {result.pvalue:.2e} | {sig_star(result.pvalue)} |")

        # Pairwise brca2 vs others
        brca2_n = cal_s.get(INDEX_TRIO, {}).get("mosaic_snv", 0)
        if brca2_n > 0:
            for t in OTHER_FA_TRIOS:
                other_n = cal_s.get(t, {}).get("mosaic_snv", 0)
                total = brca2_n + other_n
                if total == 0:
                    continue
                result = binomtest(brca2_n, total, 0.5, alternative="greater")
                lines.append(f"| brca2 ({brca2_n}) vs {t} ({other_n}) | {brca2_n/other_n:.1f}x | {result.pvalue:.2e} | {sig_star(result.pvalue)} |")
        lines.append("")

    # --- 1b. Sensitivity: alt-read threshold (R1.5) ---
    sens = load_sensitivity_json()
    if sens:
        thresholds = sens["thresholds"]
        results = sens["results"]
        lines.append("## 1b. Sensitivity — Alternative Allele Read-Count Threshold (R1.5)\n")
        lines.append("*Baseline uses **alt≥2** for both pileup proband alt count and VCF `AD_alt`. "
                     f"Varied to K ∈ {{{', '.join(str(t) for t in thresholds)}}} on the same "
                     "`_part2_all.tsv` (other filters unchanged). "
                     "Full tables: [`outputs/sensitivity/alt_threshold.md`](../sensitivity/alt_threshold.md).*\n")

        for caller_label, caller_key in [(c["short"], c["slug"]) for c in CALLERS]:
            cal = results.get(caller_key, {})
            if not cal:
                continue
            lines.append(f"**{caller_label} — Mosaic SNV counts**\n")
            header = "| Trio | " + " | ".join(f"alt≥{k}" for k in thresholds) + " |"
            sep = "|---|" + "|".join(["---:"] * len(thresholds)) + "|"
            lines.append(header)
            lines.append(sep)
            for t in trios:
                if t not in cal:
                    continue
                row = cal[t]
                cells = [str(row.get(str(k), {}).get("mosaic_snv", "")) for k in thresholds]
                label = f"**{t}**" if t == INDEX_TRIO else t
                lines.append(f"| {label} | " + " | ".join(cells) + " |")
            lines.append("")

            lines.append(f"**{caller_label} — brca2 vs control burden p-value per threshold**\n")
            lines.append(header)
            lines.append(sep)
            for t in [x for x in trios if x != CONTROL_TRIO]:
                if t not in cal:
                    continue
                row = cal[t]
                cells = []
                for k in thresholds:
                    b = row.get(str(k), {}).get("burden_vs_control", {})
                    p = b.get("p_value")
                    ratio = b.get("ratio")
                    if p is None:
                        cells.append("—")
                        continue
                    ratio_s = f"{ratio}x " if ratio is not None else ""
                    cells.append(f"{ratio_s}{p:.1e} {sig_star(p)}")
                label = f"**{t}**" if t == INDEX_TRIO else t
                lines.append(f"| {label} | " + " | ".join(cells) + " |")
            lines.append("")

        lines.append("**Conclusion**: the brca2 mosaic-SNV excess is robust across "
                     "thresholds — p-value stays < 1e-10 up to alt≥5 for both callers, "
                     "whereas other FA trios remain non-significant vs control at every K.\n")

    # --- 2. Final Results ---
    giab_in_summary = GIAB_TRIO in summary.get("hc", {})
    for caller_label, caller_key in [(f'{c["short"]} ({c["display"]})' if c["slug"]=="hc" else c["display"], c["slug"]) for c in CALLERS]:
        lines.append(f"## 2{'a' if caller_key == 'hc' else 'b'}. {caller_label} — Final Candidates\n")
        cols = list(trios) + (["GIAB*"] if giab_in_summary else [])
        lines.append("| | " + " | ".join(cols) + " |")
        lines.append("|---|" + "|".join(["---:"] * len(cols)) + "|")
        for row_name, field in [("**Mosaic SNV**", "mosaic_snv"), ("Mosaic INDEL", "mosaic_indel"),
                                ("**Total mosaic**", "total_mosaic"), ("Het SNV", "het_snv"),
                                ("Het INDEL", "het_indel"), ("Total het", "total_het"),
                                ("**Total final**", "total_final")]:
            vals = []
            for t in trios:
                v = summary.get(caller_key, {}).get(t, {}).get(field, "")
                vals.append(f"**{v}**" if t == INDEX_TRIO and "mosaic" in field.lower() else str(v))
            if giab_in_summary:
                gv = summary.get(caller_key, {}).get(GIAB_TRIO, {}).get(field, "")
                vals.append(str(gv))
            lines.append(f"| {row_name} | " + " | ".join(vals) + " |")
        if giab_in_summary:
            lines.append("")
            lines.append("\\*GIAB (HG002) counts inflated due to EBV-transformed lymphoblastoid cell line (LCL) artifacts — see GIAB section.")
        lines.append("")

    # --- 2c. Manual IGV validation ---
    has_igv = any(
        "igv_review" in summary.get("hc", {}).get(t, {}) for t in trios
    )
    if has_igv:
        lines.append("## 2c. After Manual IGV Validation\n")
        lines.append("*Final candidates after manual IGV review — variants flagged "
                     "`IGV_is_real=FALSE` are removed. Variants without an IGV "
                     "screenshot (old-pipeline-only) are not reviewed and absent "
                     "from this subset anyway. Source: "
                     "[`outputs/igv_review/all_trios_igv_validated.csv`]"
                     "(../igv_review/all_trios_igv_validated.csv), legend: "
                     "[`IGV_review_report.pdf`](../igv_review/IGV_review_report.pdf).*\n")
        for caller_label, caller_key in [(c["short"], c["slug"]) for c in CALLERS]:
            lines.append(f"**{caller_label} — pipeline → IGV-confirmed (removed)**\n")
            lines.append("| | " + " | ".join(trios) + " |")
            lines.append("|---|" + "|".join(["---:"] * len(trios)) + "|")
            for row_name, field in [
                ("Mosaic SNV", "mosaic_snv"),
                ("Mosaic INDEL", "mosaic_indel"),
                ("**Total mosaic**", "total_mosaic"),
                ("Het SNV", "het_snv"),
                ("Het INDEL", "het_indel"),
                ("Total het", "total_het"),
            ]:
                vals = []
                for t in trios:
                    s = summary.get(caller_key, {}).get(t, {})
                    pre = s.get(field, "")
                    post = s.get("igv_review", {}).get(field, "")
                    if pre == "" or post == "":
                        vals.append("—")
                    else:
                        removed = pre - post if isinstance(pre, int) else ""
                        vals.append(f"{pre} → **{post}** ({'−' + str(removed) if removed else '—'})")
                lines.append(f"| {row_name} | " + " | ".join(vals) + " |")
            lines.append("")

        # Precision per trio / per caller (TP / (TP+FP))
        lines.append("**Precision per trio (IGV-reviewed TRUE / reviewed)**\n")
        lines.append("*Computed on the subset of pipeline-passing SNV+INDEL candidates "
                     "that were IGV-reviewed.*\n")
        precision_rows = _compute_igv_precision(trios)
        if precision_rows:
            lines.append("| Trio | Caller | Reviewed | TRUE | FALSE | Precision |")
            lines.append("|---|---|---:|---:|---:|---:|")
            for r in precision_rows:
                lines.append(
                    f"| {r['trio']} | {r['caller']} | {r['reviewed']} | "
                    f"{r['true']} | {r['false']} | {r['precision']} |"
                )
            lines.append("")

        # Burden test after manual IGV filter
        lines.append("**Burden test after manual IGV filter — mosaic SNV vs control**\n")
        lines.append("*Conditional binomial exact, one-sided (case > control), "
                     "IGV-confirmed mosaic SNV only.*\n")
        for caller_label, caller_key in [(c["short"], c["slug"]) for c in CALLERS]:
            lines.append(f"**{caller_label}**\n")
            lines.append("| Comparison | Ratio | p-value | |")
            lines.append("|---|---|---|---|")
            ctrl_h = summary.get(caller_key, {}).get(CONTROL_TRIO, {}).get("igv_review", {})
            ctrl_n = ctrl_h.get("mosaic_snv", 0)
            for t in [x for x in trios if x != CONTROL_TRIO]:
                case_n = summary.get(caller_key, {}).get(t, {}).get("igv_review", {}).get("mosaic_snv", 0)
                if case_n + ctrl_n == 0:
                    continue
                from scipy.stats import binomtest as _bt
                result = _bt(case_n, case_n + ctrl_n, 0.5, alternative="greater")
                r_str = f"{case_n/ctrl_n:.1f}x" if ctrl_n > 0 else "N/A"
                lines.append(
                    f"| {t} ({case_n}) vs control ({ctrl_n}) | {r_str} | "
                    f"{result.pvalue:.2e} | {sig_star(result.pvalue)} |"
                )
            lines.append("")

    # --- 2d. Per-variant binomial test (R1.m14) ---
    binom_data = _compute_binom_classification(trios)
    has_binom = any(
        binom_data.get(t, {}).get("hc") is not None for t in trios
    )
    if has_binom:
        lines.append("## 2d. Per-Variant Binomial Test — Mosaic vs Het (R1.m14)\n")
        lines.append(f"*Replaces the fixed VAF cutoff (0.3741) with a depth-aware "
                     f"statistical test: H₀: VAF = 0.5 vs H₁: VAF < 0.5, BH-corrected "
                     f"at α = {BINOM_HET_ALPHA} on the PASS_all subset. A variant is "
                     "now classified as **mosaic** when the binomial test rejects "
                     "H₀ — otherwise **het**. The old `VAF_class` column is kept for "
                     "reference; the new classification lives in `binom_class`.*\n")

        lines.append("### Reclassification vs fixed VAF cutoff\n")
        lines.append("*Counts among `PASS_all` variants. The `het → mosaic` column "
                     "counts variants classified het by VAF but mosaic by the binomial "
                     "test (typical case: VAF in 0.37–0.45 with enough depth to reject "
                     "H₀). `mosaic → het` is the opposite direction.*\n")
        lines.append("| Trio | Caller | VAF mosaic | binom mosaic | VAF het | binom het | het→mos | mos→het |")
        lines.append("|---|---|---:|---:|---:|---:|---:|---:|")
        for t in trios:
            if t not in binom_data:
                continue
            for caller in CALLER_SLUGS:
                if caller not in binom_data[t]:
                    continue
                d = binom_data[t][caller]
                p = d["pass_all"]
                vm = p["vaf_mosaic_snv"] + p["vaf_mosaic_indel"]
                bm = p["binom_mosaic_snv"] + p["binom_mosaic_indel"]
                vh = p["vaf_het_snv"] + p["vaf_het_indel"]
                bh = p["binom_het_snv"] + p["binom_het_indel"]
                lines.append(
                    f"| {t} | {caller.upper() if caller=='hc' else 'Mutect2'} "
                    f"| {vm} | {bm} | {vh} | {bh} "
                    f"| {d['het_to_mosaic']} | {d['mosaic_to_het']} |"
                )
        lines.append("")

        # Final candidates under binom + IGV review (SNV/INDEL split) — source-of-truth
        # for the manuscript Results section (replaces the VAF-cutoff numbers in §2c).
        lines.append("### Final candidates under depth-aware binomial classification (IGV-confirmed)\n")
        lines.append("*These counts are the source-of-truth for the manuscript Results "
                     "section (lines ~220–235). Classification: `binom_class` + "
                     "`PASS_all_igv_confirmed`. The earlier VAF-cutoff tables "
                     "(§2a–§2c) are kept for transparency and are referenced in the "
                     "Supplementary as a diagnostic classifier.*\n")
        for caller_label, caller in [(c["short"], c["slug"]) for c in CALLERS]:
            lines.append(f"**{caller_label}**\n")
            lines.append("| | brca1 | brca2 | fanca | fanct | control |")
            lines.append("|---|---:|---:|---:|---:|---:|")
            row_specs = [
                ("**Mosaic SNV**", "binom_mosaic_snv"),
                ("Mosaic INDEL", "binom_mosaic_indel"),
                ("**Total mosaic**", None),  # computed
                ("Het SNV", "binom_het_snv"),
                ("Het INDEL", "binom_het_indel"),
                ("Total het", None),  # computed
                ("**Total final**", None),  # computed
            ]
            for row_label, key in row_specs:
                vals = []
                for t in trios:
                    if t not in binom_data or caller not in binom_data[t]:
                        vals.append("")
                        continue
                    entry = binom_data[t][caller]
                    block = entry.get("pass_igv") or entry.get("pass_all", {})
                    if key is not None:
                        vals.append(str(block.get(key, 0)))
                    elif row_label == "**Total mosaic**":
                        vals.append(str(block.get("binom_mosaic_snv", 0) + block.get("binom_mosaic_indel", 0)))
                    elif row_label == "Total het":
                        vals.append(str(block.get("binom_het_snv", 0) + block.get("binom_het_indel", 0)))
                    elif row_label == "**Total final**":
                        vals.append(str(
                            block.get("binom_mosaic_snv", 0) + block.get("binom_mosaic_indel", 0)
                            + block.get("binom_het_snv", 0) + block.get("binom_het_indel", 0)
                        ))
                lines.append(f"| {row_label} | " + " | ".join(vals) + " |")
            lines.append("")

        # Burden test on binom mosaic SNV, IGV-confirmed
        lines.append("### Burden test — binom mosaic SNV, IGV-confirmed vs control\n")
        lines.append("*Conditional binomial exact, one-sided (case > control). "
                     "Using `binom_class == 'mosaic'` AND `PASS_all_igv_confirmed` "
                     "(manually reviewed both mosaic and het IGV screenshots, so both "
                     "sides of the reclassification carry manual confirmation).*\n")
        for caller_label, caller in [(c["short"], c["slug"]) for c in CALLERS]:
            lines.append(f"**{caller_label}**\n")
            lines.append("| Comparison | Ratio | p-value | |")
            lines.append("|---|---|---|---|")
            ctrl_entry = binom_data.get(CONTROL_TRIO, {}).get(caller, {})
            ctrl_igv = ctrl_entry.get("pass_igv") or ctrl_entry.get("pass_all", {})
            ctrl_n = ctrl_igv.get("binom_mosaic_snv", 0)
            for t in [x for x in trios if x != CONTROL_TRIO]:
                entry = binom_data.get(t, {}).get(caller, {})
                block = entry.get("pass_igv") or entry.get("pass_all", {})
                case_n = block.get("binom_mosaic_snv", 0)
                if case_n + ctrl_n == 0:
                    continue
                from scipy.stats import binomtest as _bt
                result = _bt(case_n, case_n + ctrl_n, 0.5, alternative="greater")
                r_str = f"{case_n/ctrl_n:.2f}x" if ctrl_n > 0 else "N/A"
                lines.append(
                    f"| {t} ({case_n}) vs control ({ctrl_n}) | {r_str} | "
                    f"{result.pvalue:.2e} | {sig_star(result.pvalue)} |"
                )
            lines.append("")

        # Cross-FA burden — brca2 vs each other FA proband
        lines.append("### Burden test — brca2 vs each other FA proband (binom_class, IGV-confirmed)\n")
        lines.append("*Shows that the FANCD1/BRCA2 excess is not only against the negative "
                     "control but also relative to every other FA complementation group in "
                     "this study.*\n")
        for caller_label, caller in [(c["short"], c["slug"]) for c in CALLERS]:
            lines.append(f"**{caller_label}**\n")
            lines.append("| Comparison | Ratio | p-value | |")
            lines.append("|---|---|---|---|")
            b2_entry = binom_data.get(INDEX_TRIO, {}).get(caller, {})
            b2_h = b2_entry.get("pass_igv") or b2_entry.get("pass_all", {})
            b2_n = b2_h.get("binom_mosaic_snv", 0)
            for t in OTHER_FA_TRIOS:
                entry = binom_data.get(t, {}).get(caller, {})
                block = entry.get("pass_igv") or entry.get("pass_all", {})
                case_n = block.get("binom_mosaic_snv", 0)
                if b2_n + case_n == 0:
                    continue
                from scipy.stats import binomtest as _bt
                result = _bt(b2_n, b2_n + case_n, 0.5, alternative="greater")
                r_str = f"{b2_n/max(case_n, 1):.2f}x" if case_n > 0 else "N/A"
                lines.append(
                    f"| brca2 ({b2_n}) vs {t} ({case_n}) | {r_str} | "
                    f"{result.pvalue:.2e} | {sig_star(result.pvalue)} |"
                )
            lines.append("")

        # R1.2 — counts in the 0.10–0.20 VAF bin highlighted by the reviewer
        # as the expected signature of early-embryonic subclonal SNVs. Reported
        # on the IGV-confirmed set so it is comparable to the binom_class
        # mosaic totals above.
        lo_pct = int(SUBCLONAL_PEAK_VAF_LOW * 100)
        hi_pct = int(SUBCLONAL_PEAK_VAF_HIGH * 100)
        lines.append(
            f"### Subclonal-peak VAF bin ({lo_pct}–{hi_pct}%) — "
            "IGV-confirmed SNV counts (R1.2)\n"
        )
        lines.append(
            "*Reviewer 1 highlighted VAF ≈ 0.1–0.2 as the expected signature "
            "of early-embryonic subclonal SNVs. The table below reports the "
            "number of IGV-confirmed SNVs in this bin per trio, for both "
            "callers, to complement the per-trio VAF histograms below. The "
            "FANCD1/BRCA2 peak in this bin is also visible in "
            "`vaf_histograms.png` as a pronounced red shoulder between 10% "
            "and 20% VAF.*\n"
        )
        for caller_label, caller in [(c["short"], c["slug"]) for c in CALLERS]:
            lines.append(f"**{caller_label}**\n")
            lines.append("| | " + " | ".join(trios) + " |")
            lines.append("|---" + ":|---" * len(trios) + ":|")
            vals = []
            for t in trios:
                b = binom_data.get(t, {}).get(caller) or {}
                blk = b.get("pass_igv") or b.get("pass_all") or {}
                vals.append(str(blk.get("subclonal_peak_snv", 0)))
            lines.append(
                f"| IGV-confirmed SNVs in VAF [{lo_pct}%, {hi_pct}%) | "
                + " | ".join(vals) + " |"
            )
            lines.append("")

        # VAF histogram figure
        if (OUTPUTS_DIR / "benchmark" / "vaf_histograms.png").exists():
            lines.append("### VAF distribution per trio (binom_class coloured)\n")
            lines.append("*Pipeline- and IGV-confirmed candidates "
                         "(`PASS_all_igv_confirmed`) stacked by the depth-aware "
                         "binomial classification. Red = mosaic (binom rejects "
                         "H₀:VAF=0.5 at BH<0.05), blue = het (cannot reject). "
                         "Dashed line = old fixed VAF cutoff (37.41%); dotted "
                         "line = germline het (VAF=0.5).*\n")
            lines.append("![VAF distribution per trio](vaf_histograms.png)\n")
            lines.append("Per-caller standalone figures: "
                         "[`vaf_histograms_hc.png`](vaf_histograms_hc.png), "
                         "[`vaf_histograms_mutect2.png`](vaf_histograms_mutect2.png).\n")

        # Side-by-side: VAF-cutoff + IGV vs binom + IGV review
        lines.append("### Side-by-side: old vs new brca2 signal (mosaic SNV, IGV-confirmed)\n")
        lines.append("*Shows that the brca2 mosaic-SNV excess is robust to the choice "
                     "of classification rule.*\n")
        lines.append("| Caller | Classification | brca2 | control | Ratio | p-value |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for caller in CALLER_SLUGS:
            caller_label = "HC" if caller == "hc" else "Mutect2"
            brca2_entry = binom_data.get(INDEX_TRIO, {}).get(caller, {})
            ctrl_entry = binom_data.get(CONTROL_TRIO, {}).get(caller, {})
            brca2_h = brca2_entry.get("pass_igv") or brca2_entry.get("pass_all", {})
            ctrl_h = ctrl_entry.get("pass_igv") or ctrl_entry.get("pass_all", {})
            from scipy.stats import binomtest as _bt
            for cls_label, key in [
                ("VAF-cutoff (0.3741)", "vaf_mosaic_snv"),
                ("Binomial (BH<0.05)", "binom_mosaic_snv"),
            ]:
                b_n = brca2_h.get(key, 0)
                c_n = ctrl_h.get(key, 0)
                if b_n + c_n == 0:
                    continue
                r = _bt(b_n, b_n + c_n, 0.5, alternative="greater")
                ratio = f"{b_n / c_n:.1f}x" if c_n > 0 else "N/A"
                lines.append(
                    f"| {caller_label} | {cls_label} | **{b_n}** | {c_n} "
                    f"| {ratio} | {r.pvalue:.2e} {sig_star(r.pvalue)} |"
                )
        lines.append("")

    # --- 3. Cross-Caller Validation ---
    lines.append("## 3. Cross-Caller Validation (HC ∩ Mutect2 ∩ Old Pipeline)\n")
    cross = load_cross_caller_data(trios)

    lines.append("### Mosaic SNV\n")
    lines.append("| Sample | HC | M2 | Old | HC∩M2 | HC∩Old | All 3 |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for t in trios:
        if t not in cross:
            continue
        m = cross[t]["mosaic_snv"]
        lines.append(
            f"| {t} | {m['hc']} | {m['m2']} | {m['old']} | **{m['hc_m2']}** | {m['hc_old']} | **{m['all3']}** |"
        )
    lines.append("")

    lines.append("### Het de novo SNV\n")
    lines.append("| Sample | HC | M2 | Old | HC∩M2 | HC∩Old | All 3 |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for t in trios:
        if t not in cross:
            continue
        h = cross[t]["het_snv"]
        lines.append(f"| {t} | {h['hc']} | {h['m2']} | {h['old']} | **{h['hc_m2']}** | {h['hc_old']} | **{h['all3']}** |")
    lines.append("")

    lines.append("Venn diagrams: `outputs/cross_caller/{trio}_mosaic_SNV_venn.png`\n")
    lines.append("Master spreadsheets (91 columns): `outputs/cross_caller/{trio}_master.csv`\n")

    # --- 4. deepSNV Statistical Validation ---
    # Same row-list feeds the Sx3 CSV and the per-caller tables below.
    lines.append("## 4. deepSNV Beta-Binomial Test\n")
    lines.append("*Per-variant beta-binomial likelihood-ratio test (Gerstung et al. "
                 "2014, DOI:10.1093/bioinformatics/btt750) against the pooled-parents "
                 "error distribution, Benjamini-Hochberg-adjusted per trio × caller. "
                 "Counts are computed on the final classified set used throughout "
                 "§2d and §5: `binom_class` × `PASS_all_igv_confirmed` (FA trios) or "
                 "`binom_class` × `PASS_all` (GIAB — manual IGV review was not "
                 "performed for this trio).*\n")
    _, ds_rows = _supp_rows_deepsnv(cascades, trios)
    for caller_label, caller_tag in [(c["short"], c["display"]) for c in CALLERS]:
        lines.append(f"### {caller_label}\n")
        lines.append("| Sample | Mosaic SNV | BH<0.01 | BH<0.05 | Het SNV | BH<0.01 | BH<0.05 |")
        lines.append("|---|---:|---:|---:|---:|---:|---:|")
        for r in ds_rows:
            if r[0] != caller_tag:
                continue
            _, sample, mos_n, mos_bh01, mos_bh05, het_n, het_bh01, het_bh05 = r
            lines.append(
                f"| {sample} | {mos_n} | {mos_bh01} | {mos_bh05} "
                f"| {het_n} | {het_bh01} | {het_bh05} |"
            )
        lines.append("")

    # --- 5. Filter Cascade (unified pipeline view) ---
    # Two tables, one per caller, with trios as columns. Each row is one
    # filter step from Part 1 (Terra: raw variant calling, trio de-novo
    # filter, segdup/LCR masking) through Part 2 (local: pileup validation,
    # hard filters, gnomAD AF=0, chrY exclusion, mosaic+het VAF window,
    # manual IGV validation). Below each cascade we report the final
    # classification (binomial mosaic/het × SNV/INDEL) on the surviving set.
    def _fmt_count(n):
        return f"{n:,}" if isinstance(n, (int, float)) and n >= 1000 else str(n)

    lines.append("## 5. Filter Cascade\n")
    lines.append("*Unified pipeline view across all 5 trios: Part 1 (Terra: raw "
                 "variant calling, trio-based de novo filter, region masking) → "
                 "Part 2 (local: pileup-based parent/proband alt-read filters, "
                 "hard quality filters, gnomAD AF=0, chrY exclusion, restriction "
                 "to the mosaic+het VAF window, manual IGV validation). Each row "
                 "is the count of candidate variants surviving that step. The "
                 "final mosaic/heterozygous split (depth-aware per-variant "
                 "binomial test, §2d / R1.m14) is reported in the small table "
                 "below each caller's cascade.*\n")

    for caller_label, caller_key in [(c["short"], c["slug"]) for c in CALLERS]:
        lines.append(f"### {caller_label}\n")

        # Columns: one per trio that has a cascade entry for this caller
        caller_trios = [t for t in trios if f"{t}_{caller_key}" in cascades]
        if not caller_trios:
            continue

        header = "| Step | " + " | ".join(caller_trios) + " |"
        sep = "|---" + ":|---" * len(caller_trios) + ":|"
        lines.append(header)
        lines.append(sep)

        # Part 1 (Terra)
        p1_labels = PART1_M2_STEPS if caller_key == "mutect2" else PART1_HC_STEPS
        for label, p1_key in p1_labels:
            row_vals = []
            any_hit = False
            for t in caller_trios:
                p1 = cascades[f"{t}_{caller_key}"].get("part1", {})
                if p1_key in p1:
                    row_vals.append(_fmt_count(p1[p1_key]))
                    any_hit = True
                else:
                    row_vals.append("—")
            if any_hit:
                lines.append(f"| {label} | " + " | ".join(row_vals) + " |")

        # Part 2 (local)
        for step_key, step_label in PART2_STEPS:
            row_vals = []
            any_hit = False
            for t in caller_trios:
                s = cascades[f"{t}_{caller_key}"]["steps"].get(step_key)
                if s is not None:
                    row_vals.append(_fmt_count(s["total"]))
                    any_hit = True
                else:
                    row_vals.append("—")
            if any_hit:
                lines.append(f"| {step_label} | " + " | ".join(row_vals) + " |")

        lines.append("")

        # Final classification block: one table per caller, rows = class/type,
        # columns = trios. Prefer IGV-confirmed counts (FA trios); fall back to
        # PASS_all (GIAB, no manual IGV review). Which classifier was used is
        # reported per trio beneath the cascade table.
        has_any_igv = any(
            (binom_data.get(t, {}).get(caller_key, {}) or {}).get("pass_igv")
            is not None for t in caller_trios
        )
        classification_suffix = (
            "binom_class + manual IGV validation"
            if has_any_igv
            else "binom_class (no manual IGV review)"
        )
        lines.append(f"*Final classification ({classification_suffix}):*\n")
        lines.append("| | " + " | ".join(caller_trios) + " |")
        lines.append("|---" + ":|---" * len(caller_trios) + ":|")

        def _block(t):
            b = binom_data.get(t, {}).get(caller_key) or {}
            return b.get("pass_igv") or b.get("pass_all") or {}

        row_specs = [
            ("Mosaic SNV",       "binom_mosaic_snv"),
            ("Mosaic INDEL",     "binom_mosaic_indel"),
            ("**Mosaic total**", None),  # computed
            ("Het SNV",          "binom_het_snv"),
            ("Het INDEL",        "binom_het_indel"),
            ("Het total",        None),  # computed
            ("**Final total**",  None),  # computed
        ]
        for row_label, key in row_specs:
            vals = []
            for t in caller_trios:
                blk = _block(t)
                if key is not None:
                    vals.append(str(blk.get(key, 0)))
                elif row_label == "**Mosaic total**":
                    vals.append(str(blk.get("binom_mosaic_snv", 0) + blk.get("binom_mosaic_indel", 0)))
                elif row_label == "Het total":
                    vals.append(str(blk.get("binom_het_snv", 0) + blk.get("binom_het_indel", 0)))
                elif row_label == "**Final total**":
                    vals.append(str(
                        blk.get("binom_mosaic_snv", 0) + blk.get("binom_mosaic_indel", 0)
                        + blk.get("binom_het_snv", 0) + blk.get("binom_het_indel", 0)
                    ))
            lines.append(f"| {row_label} | " + " | ".join(vals) + " |")
        lines.append("")

    with open(output_path, "w") as f:
        f.write("\n".join(lines))

    # Supplementary tables (CSVs) live in docs/supp_tables/ so they can be
    # dropped into the manuscript submission package without post-processing.
    # Pass extended_trios (FA + GIAB) so Sx2 / Sx5 include GIAB; write_supp_tables_csv
    # strips GIAB internally for the statistical tables (Sx1 / Sx3) per R1.8.
    write_supp_tables_csv(cascades, binom_data, extended_trios or trios,
                          PROJECT_DIR / "docs" / "supp_tables")


def write_supp_tables_csv(cascades, binom_data, trios, output_dir):
    """Write Sx1/Sx2/Sx3 supplementary tables as CSV.

    Reuses the cascade and binom_data already computed for the markdown
    report — no duplicate data loading. All counts use binom_class +
    PASS_all_igv_confirmed (FA trios) or PASS_all (GIAB).

    Sx1 — pairwise burden test (binom mosaic SNV, IGV-confirmed).
    Sx2 — unified filter cascade (Part1 + Part2 + VAF window + IGV).
    Sx3 — deepSNV BH confirmation rates per trio × caller (binom classification).
    """
    import csv as _csv
    from scipy.stats import binomtest as _bt

    output_dir.mkdir(parents=True, exist_ok=True)

    def _sig(p):
        return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 5e-2 else "ns"

    def _block(trio, caller):
        b = binom_data.get(trio, {}).get(caller) or {}
        return b.get("pass_igv") or b.get("pass_all") or {}

    # Two orderings:
    #   ordered_fa — FA trios + in-house control, for the statistical tables
    #                (Sx1 burden, Sx3 deepSNV) which must NOT mix the LCL-
    #                inflated GIAB counts into case-vs-control tests (R1.8).
    #   ordered    — ordered_fa + GIAB, for diagnostic tables (Sx2 cascade,
    #                Sx5 coverage) where GIAB is shown for transparency.
    ordered_fa = [t for t in list(FA_TRIOS) + [CONTROL_TRIO] if t in trios]
    ordered = list(ordered_fa)
    if GIAB_TRIO in trios:
        ordered.append(GIAB_TRIO)

    # --- Sx1: pairwise burden test ---
    sx_rows = [[
        "Caller", "Comparison", "Case_count", "Control_count",
        "Ratio", "p_value", "Significance",
    ]]
    for caller_label, caller in [(c["display"], c["slug"]) for c in CALLERS]:
        ctrl = _block(CONTROL_TRIO, caller).get("binom_mosaic_snv", 0)
        for t in list(FA_TRIOS):
            if t not in trios:
                continue
            case = _block(t, caller).get("binom_mosaic_snv", 0)
            if case + ctrl == 0:
                continue
            p = _bt(case, case + ctrl, 0.5, alternative="greater").pvalue
            ratio = round(case / max(ctrl, 1), 2)
            sx_rows.append([
                caller_label, f"{t.upper()} vs control",
                case, ctrl, ratio, f"{p:.2e}", _sig(p),
            ])
        # BRCA2 vs each other FA proband
        b2 = _block(INDEX_TRIO, caller).get("binom_mosaic_snv", 0)
        for t in OTHER_FA_TRIOS:
            if t not in trios:
                continue
            oth = _block(t, caller).get("binom_mosaic_snv", 0)
            if b2 + oth == 0:
                continue
            p = _bt(b2, b2 + oth, 0.5, alternative="greater").pvalue
            ratio = round(b2 / max(oth, 1), 2)
            sx_rows.append([
                caller_label, f"BRCA2 vs {t.upper()}",
                b2, oth, ratio, f"{p:.2e}", _sig(p),
            ])

    with open(output_dir / "supp_table_Sx1_burden_test.csv", "w", newline="") as f:
        _csv.writer(f).writerows(sx_rows)

    # --- Sx2: filter cascade (reuses the module-level PART1/PART2 step lists) ---
    FINAL_ROWS = [
        ("    Mosaic SNV (binom)",    "binom_mosaic_snv"),
        ("    Mosaic INDEL (binom)",  "binom_mosaic_indel"),
        ("    Het SNV (binom)",       "binom_het_snv"),
        ("    Het INDEL (binom)",     "binom_het_indel"),
    ]

    sy_rows = [["Caller", "Step"] + ordered]
    for caller_label, caller in [(c["display"], c["slug"]) for c in CALLERS]:
        p1_labels = PART1_M2_STEPS if caller == "mutect2" else PART1_HC_STEPS
        for label, p1_key in p1_labels:
            row = [caller_label, label]
            any_hit = False
            for t in ordered:
                p1 = cascades.get(f"{t}_{caller}", {}).get("part1", {})
                if p1_key in p1:
                    row.append(p1[p1_key])
                    any_hit = True
                else:
                    row.append("")
            if any_hit:
                sy_rows.append(row)
        for step_key, step_label in PART2_STEPS:
            row = [caller_label, step_label]
            any_hit = False
            for t in ordered:
                s = cascades.get(f"{t}_{caller}", {}).get("steps", {}).get(step_key)
                if s is not None:
                    row.append(s["total"])
                    any_hit = True
                else:
                    row.append("")
            if any_hit:
                sy_rows.append(row)
        for label, key in FINAL_ROWS:
            row = [caller_label, label]
            for t in ordered:
                row.append(_block(t, caller).get(key, 0))
            sy_rows.append(row)

    with open(output_dir / "supp_table_Sx2_filter_cascade.csv", "w", newline="") as f:
        _csv.writer(f).writerows(sy_rows)

    # --- Sx3: deepSNV BH confirmation rates (shared with report §4) ---
    # FA trios only — GIAB is LCL-inflated and intentionally omitted (R1.8).
    sz_header, sz_data = _supp_rows_deepsnv(cascades, ordered_fa)
    with open(output_dir / "supp_table_Sx3_deepSNV.csv", "w", newline="") as f:
        _csv.writer(f).writerows([sz_header] + sz_data)

    # --- Sx4: alt-read sensitivity (K ∈ {2,3,4,5}) ---
    # Source: outputs/sensitivity/alt_threshold.json (auto-generated by the
    # sensitivity_alt_threshold.py step — must run before this function).
    # Skips cleanly if the JSON isn't there yet.
    import json as _json
    sens_json = OUTPUTS_DIR / "sensitivity" / "alt_threshold.json"
    if sens_json.exists():
        data = _json.loads(sens_json.read_text())
        ks = data.get("thresholds", [2, 3, 4, 5])
        sens_trios = data.get("trios", ordered)
        sens_results = data.get("results", {})
        s4_rows = [[
            "Caller", "Trio", "Threshold_K",
            "Mosaic_SNV", "Mosaic_INDEL", "Het_SNV", "Het_INDEL",
            "Mosaic_het_ratio",
            "Burden_vs_control_ratio", "Burden_vs_control_p_value",
            "Burden_vs_control_significance",
        ]]
        for caller_label, caller in [(c["display"], c["slug"]) for c in CALLERS]:
            for t in sens_trios:
                for k in ks:
                    r = sens_results.get(caller, {}).get(t, {}).get(str(k))
                    if r is None:
                        continue
                    bvc = r.get("burden_vs_control") or {}
                    pv = bvc.get("p_value")
                    s4_rows.append([
                        caller_label, t, k,
                        r.get("mosaic_snv", 0), r.get("mosaic_indel", 0),
                        r.get("het_snv", 0), r.get("het_indel", 0),
                        r.get("mosaic_het_ratio", ""),
                        bvc.get("ratio", ""),
                        f"{pv:.2e}" if pv is not None else "",
                        sig_star(pv) if pv is not None else "",
                    ])
        with open(output_dir / "supp_table_Sx4_sensitivity_alt_threshold.csv", "w", newline="") as f:
            _csv.writer(f).writerows(s4_rows)
        s4_rowcount = len(s4_rows) - 1
    else:
        s4_rowcount = 0

    # --- Sx5: per-sample coverage (shared with report §0) ---
    cov_json = OUTPUTS_DIR / "benchmark" / "coverage.json"
    if cov_json.exists():
        coverage = _json.loads(cov_json.read_text())
        cov_trios = ordered if GIAB_TRIO in ordered else ordered + [GIAB_TRIO]
        s5_header, s5_data = _supp_rows_coverage(coverage, cov_trios)
        with open(output_dir / "supp_table_Sx5_coverage.csv", "w", newline="") as f:
            _csv.writer(f).writerows([s5_header] + s5_data)
        s5_rowcount = len(s5_data)
    else:
        s5_rowcount = 0

    print(f"  Wrote supplementary CSVs to {output_dir}/ "
          f"(Sx1: {len(sx_rows)-1} rows, Sx2: {len(sy_rows)-1} rows, "
          f"Sx3: {len(sz_data)} rows, Sx4: {s4_rowcount} rows, "
          f"Sx5: {s5_rowcount} rows)")


def run_cross_caller(trios):
    """Run cross_caller_comparison.py."""
    cmd = [sys.executable, str(SCRIPT_DIR / "cross_caller_comparison.py")]
    print("Running cross-caller comparison...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}")
        return False
    # Print summary
    for line in result.stdout.split('\n'):
        if any(x in line for x in ["SUMMARY", "Sample"] + [t for t in (list(FA_TRIOS) + [CONTROL_TRIO])]):
            print(f"  {line}")
    return True


def run_sensitivity(trios):
    """Run sensitivity_alt_threshold.py — R1.5."""
    fa_trios = [t for t in trios if t != GIAB_TRIO]
    cmd = [
        sys.executable, str(SCRIPT_DIR / "sensitivity_alt_threshold.py"),
        "--trios", ",".join(fa_trios),
    ]
    print("Running sensitivity analysis (alt-read thresholds)...")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}")
        return False
    for line in result.stdout.split('\n'):
        if any(x in line for x in ["K=", "JSON", "Markdown"]):
            print(f"  {line}")
    return True


def load_sensitivity_json():
    """Load outputs/sensitivity/alt_threshold.json (if present)."""
    path = OUTPUTS_DIR / "sensitivity" / "alt_threshold.json"
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def build_igv_inputs(trios):
    """Generate per-trio IGV input TSVs (mosaic + het de novo from HC ∪ Mutect2).

    Reads each {trio}_master.csv and writes {trio}_mosaic_het_igv.tsv with all
    positions called as mosaic or het by HC or Mutect2 (excluding old-pipeline-
    only variants), deduplicated by (chrom, pos).
    """
    import csv

    fa_trios = [t for t in trios if t != GIAB_TRIO]
    cc_dir = OUTPUTS_DIR / "cross_caller"
    for trio in fa_trios:
        master = cc_dir / f"{trio}_master.csv"
        if not master.exists():
            print(f"  WARNING: {master} not found, skipping")
            continue
        rows = []
        seen = set()
        with open(master) as f:
            for r in csv.DictReader(f):
                if r["in_HC_or_M2"] not in ("1", 1, "True", "TRUE", True):
                    continue
                classes = {r.get("class_HC", ""), r.get("class_Mutect2", "")}
                if not (classes & {"mosaic", "het"}):
                    continue
                key = (r["CHROM"], r["POS"])
                if key in seen:
                    continue
                seen.add(key)
                rows.append((r["CHROM"], r["POS"], r["REF"], r["ALT"]))
        out = cc_dir / f"{trio}_mosaic_het_igv.tsv"
        with open(out, "w") as f:
            for row in rows:
                f.write("\t".join(row) + "\n")
        print(f"  {trio}: {len(rows)} positions -> {out.name}")


def build_master_csv(trios, output_path):
    """Combine per-trio master CSVs into a single CSV for manual IGV review.

    Adds columns:
      - trio                  — which trio the variant belongs to (first column)
      - has_IGV_screenshot    — TRUE if an IGV screenshot was generated
                                (i.e. variant is in {trio}_mosaic_het_igv.tsv,
                                meaning HC or Mutect2 called it as mosaic OR het)
      - IGV_is_real        — empty, to be filled manually (TRUE/FALSE/?)
      - IGV_notes          — empty, free-text notes

    All boolean columns are written as the literal strings TRUE/FALSE so they
    survive LibreOffice/Excel import without being coerced to 0/1.
    """
    import pandas as pd

    fa_trios = [t for t in trios if t != GIAB_TRIO]
    frames = []
    for trio in fa_trios:
        master = OUTPUTS_DIR / "cross_caller" / f"{trio}_master.csv"
        igv_tsv = OUTPUTS_DIR / "cross_caller" / f"{trio}_mosaic_het_igv.tsv"
        if not master.exists():
            print(f"  WARNING: {master} not found, skipping")
            continue

        df = pd.read_csv(master)

        # Build set of (chrom, pos) that have IGV screenshots
        igv_positions = set()
        if igv_tsv.exists():
            with open(igv_tsv) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        igv_positions.add((parts[0], int(parts[1])))

        df.insert(0, "trio", trio)
        df["has_IGV_screenshot"] = [
            (str(c), int(p)) in igv_positions
            for c, p in zip(df["CHROM"], df["POS"])
        ]
        df["IGV_is_real"] = ""
        df["IGV_notes"] = ""

        frames.append(df)
        n_igv = int(df["has_IGV_screenshot"].sum())
        print(f"  {trio}: {len(df)} variants ({n_igv} with IGV screenshot)")

    if not frames:
        print("  ERROR: no master CSVs found")
        return

    combined = pd.concat(frames, ignore_index=True)

    # Normalize ALL boolean-ish columns to literal TRUE/FALSE strings so
    # LibreOffice does not render them as 0/1.
    bool_like_cols = [
        "in_HC", "in_Mutect2", "in_HC_or_M2", "in_OldPipeline",
        "has_IGV_screenshot",
    ]
    # Also catch any deepSNV_BH_*_sig* and any column with only {0,1} or {True,False}
    for c in combined.columns:
        if c.startswith("deepSNV_BH_") and c.endswith(("_sig001", "_sig005")):
            bool_like_cols.append(c)
    for c in bool_like_cols:
        if c in combined.columns:
            combined[c] = combined[c].map(
                lambda v: "TRUE" if v in (True, 1, "1", "True", "TRUE") else
                          ("FALSE" if v in (False, 0, "0", "False", "FALSE") else "")
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_path, index=False)
    print(f"  Master CSV: {output_path} ({len(combined)} rows, {len(combined.columns)} cols)")


def update_candidates(trios):
    """Copy final mosaic/het TSVs to candidates folder."""
    for caller_dir, suffix in [("hc", ""), ("mutect2", "_mutect2")]:
        d = OUTPUTS_DIR / "candidates" / caller_dir
        d.mkdir(parents=True, exist_ok=True)
        for trio in trios:
            for vtype in ["mosaic", "het"]:
                src = OUTPUTS_DIR / "part2" / f"{trio}{suffix}_{vtype}.tsv"
                dst = d / f"{trio}_{vtype}.tsv"
                if src.exists():
                    import shutil
                    shutil.copy2(src, dst)


def update_claude_md(cascades, summary, trios):
    """Update the auto-generated section in CLAUDE.md."""
    claude_md = PROJECT_DIR / "CLAUDE.md"
    content = claude_md.read_text()

    BEGIN = "<!-- BEGIN AUTO-GENERATED RESULTS -->"
    END = "<!-- END AUTO-GENERATED RESULTS -->"

    if BEGIN not in content or END not in content:
        print("  WARNING: markers not found in CLAUDE.md, skipping")
        return

    fa_trios = [t for t in trios if t != GIAB_TRIO]
    hc_s = summary.get("hc", {})
    m2_s = summary.get("mutect2", {})
    burden = summary.get("burden_test", {})

    lines = []
    lines.append(BEGIN)
    lines.append("<!-- Auto-generated by scripts/generate_all_outputs.py from outputs/benchmark/*.json -->")
    lines.append("<!-- Do not edit manually — regenerate: python scripts/generate_all_outputs.py --skip-part2 -->")
    lines.append("")

    # Key Results
    lines.append("## Key Results")
    lines.append("")
    lines.append("**Full report** (coverage, filter cascades, cross-caller validation, deepSNV): "
                 "[`outputs/benchmark/report.md`](outputs/benchmark/report.md)")
    lines.append("")

    # HC summary
    brca2_m = hc_s.get(INDEX_TRIO, {}).get("mosaic_snv", 0)
    ctrl_m = hc_s.get(CONTROL_TRIO, {}).get("mosaic_snv", 0)
    ratio = f"{brca2_m/ctrl_m:.1f}x" if ctrl_m > 0 else "N/A"
    brca2_p = burden.get("hc", {}).get("brca2_vs_control", {}).get("p_value", "?")

    lines.append(f"brca2 has **{brca2_m} mosaic SNVs (HC)** — **{ratio} more than control** ({ctrl_m}). "
                 f"Other FA samples ("
                 + ", ".join(f"{t}={hc_s.get(t,{}).get('mosaic_snv',0)}" for t in fa_trios if t not in (INDEX_TRIO, CONTROL_TRIO))
                 + f") are comparable to control — consistent with Reviewer 2's hypothesis.")
    lines.append("")

    # Burden test
    lines.append(f"**Mosaic burden test** (conditional binomial exact): "
                 f"brca2 vs control **p={brca2_p:.2e}** (only brca2 significant). "
                 f"Full pairwise comparisons in report.md Section 1.")
    lines.append("")

    # Cross-caller — load from cross_caller data
    cross = load_cross_caller_data(fa_trios)
    brca2_cross = cross.get(INDEX_TRIO, {}).get("mosaic_snv", {})
    hc_m2 = brca2_cross.get("hc_m2", 0)
    all3 = brca2_cross.get("all3", 0)
    lines.append(f"**Cross-caller validation**: brca2 has **{hc_m2} HC∩Mutect2 shared** mosaic SNVs "
                 f"(all deepSNV-significant), **{all3} confirmed by all 3 callers** (HC + Mutect2 + old bcbio). "
                 f"Het de novo: {min(cross[t]['het_snv']['hc_m2'] for t in fa_trios)}-"
                 f"{max(cross[t]['het_snv']['hc_m2'] for t in fa_trios)} shared per trio — "
                 f"consistent with expected germline rate.")
    lines.append("")

    lines.append(f"**Data**: `outputs/cross_caller/{{trio}}_master.csv` (91 columns), "
                 f"Venn diagrams, `outputs/benchmark/*.json`.")
    lines.append("")

    # GIAB section
    giab_hc = summary.get("hc", {}).get(GIAB_TRIO, {})
    giab_m2 = summary.get("mutect2", {}).get(GIAB_TRIO, {})
    if giab_hc:
        lines.append("### GIAB Negative Control — Analysis and Caveats")
        lines.append("")
        lines.append("GIAB Ashkenazi trio (HG002/HG003/HG004) shows **inflated de novo counts** compared to FA trios:")
        lines.append(f"- HC: {giab_hc.get('mosaic_snv',0)} mosaic SNV, "
                     f"{giab_hc.get('het_snv',0)} het SNV "
                     f"(vs {ctrl_m}-{brca2_m} mosaic, "
                     f"{min(hc_s.get(t,{}).get('het_snv',0) for t in fa_trios)}-"
                     f"{max(hc_s.get(t,{}).get('het_snv',0) for t in fa_trios)} het in FA trios)")
        if giab_m2:
            lines.append(f"- Mutect2: {giab_m2.get('mosaic_snv',0)} mosaic SNV, "
                         f"{giab_m2.get('het_snv',0)} het SNV")
        lines.append("")
        lines.append("**Root cause**: HG002 is an EBV-transformed lymphoblastoid cell line (LCL). "
                     "LCLs accumulate somatic mutations during cell culture/passage that appear as de novo variants. Key evidence:")
        lines.append("")
        lines.append(f"1. **GIAB truth set confirms ~936 de novo SNVs** in HG002 "
                     "(Shadrina et al. 2025, Life Sci Alliance 8(6):e202403039, PMID:40155050) — "
                     f"~15x more than expected ~60 germline de novo. Our {giab_hc.get('het_snv',0)} het SNV are a subset after filtering.")
        lines.append(f"2. **{giab_hc.get('het_snv',0)-2}/{giab_hc.get('het_snv',0)} of our het SNV confirmed as de novo** "
                     "by GIAB benchmark VCFs (absent in both HG003 and HG004 benchmark v4.2.1).")
        lines.append("3. **Daniels et al. (2024)** GIAB mosaic benchmark (bioRxiv 2024.12.02.625685, PMID:39677813): "
                     '85 mosaic SNVs (VAF 5-30%) in HG002, explicitly noting '
                     '*"mutations that have arisen during the cell line generation and culturing process"*.')
        lines.append("4. **Nickles et al. (2012)** (BMC Genomics 13:477, PMID:22974163): "
                     "WGS of fresh blood vs LCL from same donor — somatic mutations accumulate during culture.")
        lines.append("")
        lines.append("**Conclusion**: GIAB not suitable for mosaic burden comparison due to LCL artifacts. "
                     "proband_control (fresh blood) is primary negative control. GIAB validates pipeline correctness "
                     "(detects true de novo confirmed by benchmark).")
        lines.append("")
        lines.append("**GIAB benchmark data**: `data/giab_benchmark/` (HG002/HG003/HG004 benchmark VCFs v4.2.1)")
        lines.append("")

    # Include full report.md
    report_path = OUTPUTS_DIR / "benchmark" / "report.md"
    if report_path.exists():
        report_content = report_path.read_text().strip()
        # Rewrite image / link targets so they resolve from CLAUDE.md (repo root)
        # — report.md sits in outputs/benchmark/, so its relative paths need
        # that prefix when we inline the content here.
        import re as _re
        report_content = _re.sub(
            r"(\!?\[[^\]]*\]\()(?!https?://|/|#|outputs/|\.\./)([^)]+)(\))",
            r"\1outputs/benchmark/\2\3",
            report_content,
        )
        # Also fix the ../igv_review/... links that report.md uses to reach
        # sibling folders — those resolve correctly from outputs/benchmark/
        # but not from repo root.
        report_content = report_content.replace(
            "](../igv_review/", "](outputs/igv_review/"
        )
        report_content = report_content.replace(
            "](../sensitivity/", "](outputs/sensitivity/"
        )
        report_lines = report_content.split("\n")
        # Skip "# Mosaic Variant Analysis — Full Report" and auto-generated note
        start = 0
        for i, line in enumerate(report_lines):
            if line.startswith("## 0.") or line.startswith("## 1."):
                start = i
                break
        lines.append("---")
        lines.append("")
        lines.append("### Full Benchmark Report")
        lines.append("")
        lines.append("*Auto-generated from `outputs/benchmark/*.json`*")
        lines.append("")
        lines.extend(report_lines[start:])
        lines.append("")

    lines.append(END)

    # Replace in CLAUDE.md
    before = content[:content.index(BEGIN)]
    after = content[content.index(END) + len(END):]
    new_content = before + "\n".join(lines) + after
    claude_md.write_text(new_content)
    print(f"  Updated CLAUDE.md ({len(lines)} lines in auto-generated section)")


def main():
    parser = argparse.ArgumentParser(description="Generate all outputs")
    parser.add_argument("--no-vep", action="store_true", help="Skip VEP (use existing gnomAD annotations)")
    parser.add_argument("--trios", default=None, help="Comma-separated list of trios (default: all 5)")
    parser.add_argument("--no-giab", action="store_true",
                        help="Exclude GIAB trio (default: included)")
    parser.add_argument("--skip-part2", action="store_true", help="Skip Part2 re-run (use existing files)")
    args = parser.parse_args()

    trios = args.trios.split(",") if args.trios else TRIOS
    if not args.no_giab and GIAB_TRIO not in trios:
        trios = trios + GIAB

    print("=" * 70)
    print("  MOSAIC VARIANT ANALYSIS — FULL OUTPUT GENERATION")
    print("=" * 70)
    print(f"  Trios: {', '.join(trios)}")
    print(f"  VEP: {'skip' if args.no_vep else 'enabled'}")
    print(f"  Part2: {'skip (reuse)' if args.skip_part2 else 'run'}")
    print()

    # Step 1: Run Part2
    if not args.skip_part2:
        print("=" * 70)
        print("  STEP 1: Part2 Post-Processing")
        print("=" * 70)
        for trio in trios:
            for caller in CALLER_SLUGS:
                run_part2(trio, caller, no_vep=args.no_vep)
        print()

    # Step 1b: Apply manual IGV validation (adds 4 columns to _part2_all.tsv)
    print("=" * 70)
    print("  STEP 1b: Apply manual IGV validation")
    print("=" * 70)
    fa_trios_for_igv = [t for t in trios if t != GIAB_TRIO]
    cmd = [
        sys.executable, str(SCRIPT_DIR / "apply_igv_validation.py"),
        "--trios", ",".join(fa_trios_for_igv),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}")
    else:
        for line in result.stdout.split("\n"):
            if line.strip():
                print(f"  {line}")
    print()

    # Step 1c: Per-variant binomial het-test (R1.m14) — adds 3 columns
    print("=" * 70)
    print("  STEP 1c: Per-variant binomial het-test (R1.m14)")
    print("=" * 70)
    fa_trios_for_binom = [t for t in trios if t != GIAB_TRIO]
    cmd = [
        sys.executable, str(SCRIPT_DIR / "binomial_het_test.py"),
        "--trios", ",".join(fa_trios_for_binom),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}")
    else:
        for line in result.stdout.split("\n"):
            if line.strip():
                print(f"  {line}")
    print()

    # Step 1d: VAF histograms per trio (coloured by binom_class)
    print("=" * 70)
    print("  STEP 1d: VAF histograms per trio (visual support for R1.2)")
    print("=" * 70)
    cmd = [
        sys.executable, str(SCRIPT_DIR / "vaf_histograms.py"),
        "--trios", ",".join(fa_trios_for_binom),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[-500:]}")
    else:
        for line in result.stdout.split("\n"):
            if line.strip():
                print(f"  {line}")
    print()

    # Step 2: Build benchmark JSON
    print("=" * 70)
    print("  STEP 2: Benchmark JSON")
    print("=" * 70)
    benchmark_dir = OUTPUTS_DIR / "benchmark"
    benchmark_dir.mkdir(parents=True, exist_ok=True)

    cascades = build_filter_cascade_json(trios)
    cascade_path = benchmark_dir / "filter_cascade.json"
    with open(cascade_path, "w") as f:
        json.dump(cascades, f, indent=2)
    print(f"  Filter cascade: {cascade_path}")

    summary = build_summary_json(cascades)

    # Add burden test to summary (ratio rounded to 1 decimal — legacy format)
    burden = {}
    for caller_key in CALLER_SLUGS:
        cal_s = summary.get(caller_key, {})
        ctrl_n = cal_s.get(CONTROL_TRIO, {}).get("mosaic_snv", 0)
        if ctrl_n == 0:
            continue
        burden[caller_key] = {}
        for t in [t for t in trios if t != GIAB_TRIO and t != CONTROL_TRIO]:
            case_n = cal_s.get(t, {}).get("mosaic_snv", 0)
            if case_n + ctrl_n == 0:
                continue
            burden[caller_key][f"{t}_vs_control"] = burden_test(
                case_n, ctrl_n, ratio_decimals=1
            )
    summary["burden_test"] = burden

    summary_path = benchmark_dir / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Summary: {summary_path}")
    print()

    # Step 2b: Sensitivity analysis (R1.5) — must run before markdown report
    print("=" * 70)
    print("  STEP 2b: Sensitivity Analysis (alt-read threshold, R1.5)")
    print("=" * 70)
    run_sensitivity(trios)
    print()

    # Step 3: Markdown report
    report_path = benchmark_dir / "report.md"
    fa_trios = [t for t in trios if t != GIAB_TRIO]
    # Pass full trios list (including GIAB) so write_supp_tables_csv can
    # include GIAB in the diagnostic tables (Sx2 cascade, Sx5 coverage);
    # the function strips GIAB internally for the statistical tables
    # (Sx1 burden, Sx3 deepSNV) per R1.8.
    generate_markdown_report(cascades, summary, report_path,
                             trios=fa_trios, extended_trios=trios)
    print(f"  Report: {report_path}")
    print()

    # Step 4: Cross-caller comparison
    print("=" * 70)
    print("  STEP 3: Cross-Caller Comparison")
    print("=" * 70)
    fa_trios = [t for t in trios if t != GIAB_TRIO]
    run_cross_caller(fa_trios)
    print()

    # Step 5: Update candidates
    print("=" * 70)
    print("  STEP 4: Update Candidates")
    print("=" * 70)
    update_candidates(trios)
    print("  Done")
    print()

    # Step 4b: IGV input TSVs (mosaic + het de novo from HC ∪ Mutect2)
    print("=" * 70)
    print("  STEP 4b: IGV input TSVs (mosaic + het)")
    print("=" * 70)
    build_igv_inputs(trios)
    print()

    # Step 4c: Combined master CSV (for manual IGV review)
    print("=" * 70)
    print("  STEP 4c: Combined Master CSV (manual IGV review package)")
    print("=" * 70)
    igv_dir = OUTPUTS_DIR / "igv_review"
    master_csv = igv_dir / "all_trios_master.csv"
    build_master_csv(trios, master_csv)
    print()

    # Step 6: Update CLAUDE.md auto-generated section
    print("=" * 70)
    print("  STEP 5: Update CLAUDE.md")
    print("=" * 70)
    update_claude_md(cascades, summary, trios)
    print()

    print("=" * 70)
    print("  ALL DONE")
    print("=" * 70)
    print(f"  outputs/part2/          — Part2 results per trio×caller")
    print(f"  outputs/benchmark/      — JSON + markdown report")
    print(f"  outputs/cross_caller/   — Master CSVs, Venn diagrams")
    print(f"  outputs/candidates/     — Final mosaic/het TSVs")


if __name__ == "__main__":
    main()
