"""Shared constants and helpers for the mosaic post-processing scripts.

Kept deliberately small — no filtering or statistical logic lives here, only
values and light helpers that were duplicated verbatim across scripts.
"""

import os
from pathlib import Path

import polars as pl
from scipy.stats import binomtest


# --- Paths ---------------------------------------------------------------
SCRIPT_DIR = Path(__file__).parent
PROJECT_DIR = SCRIPT_DIR.parent
OUTPUTS_DIR = PROJECT_DIR / "outputs"


# --- VAF bin boundaries (from the R pipeline; keep in sync) --------------
VAF_MOSAIC_LOW = 0.0955
VAF_MOSAIC_HIGH = 0.3741
VAF_HET_HIGH = 0.6259
# Mid-point inside the mosaic range, used only for the cosmetic 4-bin
# breakdown (mosaic-low 9.5-18.7 / mosaic-high 18.7-37.4 / het 37.4-62.6).
VAF_MOSAIC_MID = 0.1869

# Reviewer 1 (R1.2) highlighted the 0.1–0.2 VAF range as the expected
# signature of early-embryonic subclonal SNVs ("we may hope to see more
# mutations in a VAF bin of ~0.1-0.2"). We report per-trio counts in this
# bin in outputs/benchmark/report.md §2d.
SUBCLONAL_PEAK_VAF_LOW = 0.10
SUBCLONAL_PEAK_VAF_HIGH = 0.20


# --- Pileup depth thresholds (trio-level de novo validation) -------------
MIN_MOTHER_DEPTH = 20   # pileup mother_depth >
MIN_FATHER_DEPTH = 20   # pileup father_depth >
MIN_PROBAND_DEPTH = 10  # pileup proband_depth > (proband covered deeper than parents)
MIN_VCF_DP = 20         # proband VCF DP >


# --- Alt-read thresholds -------------------------------------------------
# Minimum proband alt-read count (both pileup base count and VCF AD_alt).
# The sensitivity analysis (R1.5) varies this across SENSITIVITY_THRESHOLDS.
MIN_ALT_READS = 2
SENSITIVITY_THRESHOLDS = [2, 3, 4, 5]


# --- HC hard filters (GATK VariantFiltration thresholds) -----------------
# SNV
HC_SNV_MIN_MQ = 30
HC_SNV_MIN_READ_POS_RANK_SUM = -8
HC_SNV_MAX_FS = 60
HC_SNV_MIN_QD = 2
# INDEL
HC_INDEL_MIN_READ_POS_RANK_SUM = -20
HC_INDEL_MAX_FS = 200
HC_INDEL_MIN_QD = 2


# --- deepSNV overdispersion estimation (beta-binomial MoM) ---------------
DEEPSNV_ESTIMATION_WINDOW = 20000   # first N SNV positions used
DEEPSNV_MAX_ERROR_ALT = 2           # per-parent alt ≤ this is treated as error
DEEPSNV_MIN_OBSERVATIONS = 100      # need ≥ this many error observations
DEEPSNV_DEFAULT_RHO = 0.01          # used when not enough data
DEEPSNV_FALLBACK_RHO = 0.001        # used when MoM variance ≤ binomial variance
DEEPSNV_MIN_RHO = 0.001             # clamp lower bound for rho
DEEPSNV_MAX_RHO = 0.05              # clamp upper bound for rho
DEEPSNV_EPS = 1e-10                 # numerical floor for mu / alpha / beta


# --- BH significance levels (used by all scripts) ------------------------
BH_SIG_STRONG = 0.01    # stringent
BH_SIG_MODERATE = 0.05  # moderate

# Reviewer 1 (R1.m14) — per-variant binomial test H0: VAF=0.5 vs H1: VAF<0.5,
# BH-corrected on PASS_all candidates. Used to classify mosaic vs het in a
# depth-aware way (replaces the fixed VAF=0.3741 cutoff in `binom_class`).
BINOM_HET_ALPHA = BH_SIG_MODERATE


# --- Manual IGV validation ------------------------------------------------
IGV_VALIDATION_CSV = (
    OUTPUTS_DIR / "igv_review" / "all_trios_igv_validated.csv"
)


# --- VEP config ----------------------------------------------------------
# Path to the local VEP cache. Override with the VEP_DATA environment
# variable if your cache lives somewhere other than the default below.
# The cache must contain the `homo_sapiens/<release>_GRCh38` subtree (see
# https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html).
VEP_DATA = Path(os.environ.get("VEP_DATA", "/path/to/vep_data"))
VEP_IMAGE = "ensemblorg/ensembl-vep:release_115.1"
VEP_FORKS = 8


# --- Polars CSV parsing --------------------------------------------------
POLARS_INFER_SCHEMA = 10000  # rows used to infer dtypes on VCF-derived CSVs


# --- Sample / trio map ---------------------------------------------------
# Trio → sample metadata. Proband / parent identifiers are deliberately
# replaced with generic placeholders in this public release; the GIAB trio
# keeps its publicly available accessions (HG002/HG003/HG004). To reproduce
# the analysis on your own trios, edit this dict (or override from a site-
# local `samples.json`) so that the `sample_id`, `mother` and `father`
# fields match the file-name stems used by your BAM / VCF outputs under
# `outputs/part1_csv/` and `outputs/coverage/`.
SAMPLES: dict[str, dict[str, str]] = {
    "brca1":   {"proband": "proband_brca1",   "sample_id": "proband_brca1_150x",
                "gene": "BRCA1",
                "mother": "mother_brca1",
                "father": "father_brca1"},
    "brca2":   {"proband": "proband_brca2",   "sample_id": "proband_brca2_150x",
                "gene": "BRCA2",
                "mother": "mother_brca2",
                "father": "father_brca2"},
    "fanca":   {"proband": "proband_fanca",   "sample_id": "proband_fanca_150x",
                "gene": "FANCA",
                "mother": "mother_fanca",
                "father": "father_fanca"},
    "fanct":   {"proband": "proband_fanct",   "sample_id": "proband_fanct_150x",
                "gene": "UBE2T",
                "mother": "mother_fanct",
                "father": "father_fanct"},
    "control": {"proband": "proband_control", "sample_id": "proband_control_150x",
                "gene": "(none)",
                "mother": "mother_control",
                "father": "father_control"},
    "giab":    {"proband": "HG002",           "sample_id": "HG002_giab_150x",
                "gene": "GIAB control",
                "mother": "HG004_giab",
                "father": "HG003_giab"},
}

TRIOS = ["brca1", "brca2", "fanca", "fanct", "control"]
GIAB = ["giab"]

# Named groupings used by downstream analyses and report generators so that
# trio membership is never hard-coded by a string literal.
INDEX_TRIO = "brca2"                              # FANCD1/BRCA2 — primary case
CONTROL_TRIO = "control"                          # in-house blood-derived negative control
GIAB_TRIO = "giab"                                # GIAB HG002 — external negative control (R1.8)
OTHER_FA_TRIOS = ("brca1", "fanca", "fanct")      # FA probands other than INDEX_TRIO
FA_TRIOS = (INDEX_TRIO,) + OTHER_FA_TRIOS         # all FA probands (case set)
# TRIOS (list above) = FA_TRIOS + (CONTROL_TRIO,) — kept as a list for call-site
# compatibility. GIAB_TRIO is appended only when --include-giab is requested.


# --- Variant callers ------------------------------------------------------
# Each caller is represented by three identifiers used throughout the codebase:
#   display — long name used in supplementary-table "Caller" columns
#   slug    — short key used in dict lookups, file-name suffixes, and cascade keys
#   short   — label used in section titles and short text ("HC", "Mutect2")
# Tuple ordering is stable; iterate CALLERS instead of inlining the pair.
CALLERS = (
    {"display": "HaplotypeCaller", "slug": "hc",      "short": "HC"},
    {"display": "Mutect2",         "slug": "mutect2", "short": "Mutect2"},
)
CALLER_SLUGS = tuple(c["slug"] for c in CALLERS)


# --- VAF classification --------------------------------------------------
def classify_vaf(vaf: float) -> str:
    """Classify a VAF into mosaic/het/other (cross-caller convention)."""
    if VAF_MOSAIC_LOW <= vaf < VAF_MOSAIC_HIGH:
        return "mosaic"
    if VAF_MOSAIC_HIGH <= vaf <= VAF_HET_HIGH:
        return "het"
    return "other"


# --- Helpers for reading _part2_all.tsv ----------------------------------
PART2_BOOL_COLS = (
    "pileupDepthFilterTrio", "pileupParentFilter", "pileupProbandFilter",
    "pileupFullFilter", "HardFilter", "centromereFilter", "repeatFilter",
    "repeatFullFilter", "gnomAD_filter", "sexChrFilter", "PASS_all",
)


def normalize_bool_columns(df: pl.DataFrame) -> pl.DataFrame:
    """Cast Part2 boolean-ish columns written as "true"/"false" strings to bool.

    polars re-reads TSVs with bool columns as Utf8 when any value is missing.
    Keeps columns that are already bool untouched.
    """
    for c in PART2_BOOL_COLS:
        if c in df.columns and df[c].dtype == pl.Utf8:
            df = df.with_columns(
                (pl.col(c).cast(pl.Utf8).str.to_lowercase() == "true").alias(c)
            )
    return df


def part2_all_path(trio: str, caller: str) -> Path:
    """Path to outputs/part2/{trio}[_mutect2]_part2_all.tsv."""
    suffix = "_mutect2" if caller == "mutect2" else ""
    return OUTPUTS_DIR / "part2" / f"{trio}{suffix}_part2_all.tsv"


def load_part2_all(trio: str, caller: str) -> pl.DataFrame | None:
    """Load _part2_all.tsv with boolean columns normalised. None if missing."""
    path = part2_all_path(trio, caller)
    if not path.exists():
        return None
    df = pl.read_csv(str(path), separator="\t", infer_schema_length=20000)
    return normalize_bool_columns(df)


# --- Stats ---------------------------------------------------------------
def burden_test(case_n: int, ctrl_n: int, ratio_decimals: int = 2) -> dict:
    """Conditional binomial exact test, one-sided (case > control).

    Treats the two counts as Poisson (equal exposure). Returns a dict suitable
    for serialisation with ratio, p-value and a boolean significance flag.

    `ratio_decimals` controls rounding of the ratio only — different call sites
    historically use different precision (summary.json uses 1, sensitivity uses 2).
    """
    total = case_n + ctrl_n
    if total == 0:
        return {"case": case_n, "control": ctrl_n,
                "ratio": None, "p_value": None, "significant": False}
    result = binomtest(case_n, total, 0.5, alternative="greater")
    return {
        "case": case_n,
        "control": ctrl_n,
        "ratio": round(case_n / ctrl_n, ratio_decimals) if ctrl_n > 0 else None,
        "p_value": float(f"{result.pvalue:.4e}"),
        "significant": bool(result.pvalue < 0.05),
    }


def sig_star(p: float | None) -> str:
    """Convention: *** p<.001, ** p<.01, * p<.05, ns otherwise."""
    if p is None:
        return "N/A"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


# --- VAF breakdown -------------------------------------------------------
def vaf_breakdown(df: pl.DataFrame) -> dict:
    """Count mosaic/het × SNV/INDEL for a filtered subset of variants.

    Expects columns: Type ∈ {SNV, INDEL}, VAF (float).
    Returns dict with mosaic_{snv,indel,all}, het_{snv,indel,all}, total.
    """
    def bin_count(vtype: str, lo: float, hi: float) -> int:
        return len(df.filter(
            (pl.col("Type") == vtype) & (pl.col("VAF") >= lo) & (pl.col("VAF") < hi)
        ))

    m_snv = bin_count("SNV", VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    m_ind = bin_count("INDEL", VAF_MOSAIC_LOW, VAF_MOSAIC_HIGH)
    h_snv = bin_count("SNV", VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    h_ind = bin_count("INDEL", VAF_MOSAIC_HIGH, VAF_HET_HIGH)
    return {
        "total": len(df),
        "mosaic_snv": m_snv, "mosaic_indel": m_ind, "mosaic_all": m_snv + m_ind,
        "het_snv": h_snv, "het_indel": h_ind, "het_all": h_snv + h_ind,
    }
