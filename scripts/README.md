# Scripts

Analysis scripts for the Fanconi Anemia somatic mosaicism study.

## Master Pipeline

```bash
# Regenerate all outputs from Part1 CSVs (skip VEP if gnomAD already annotated)
python scripts/generate_all_outputs.py --skip-part2

# Full re-run including Part2 (with VEP gnomAD annotation)
python scripts/generate_all_outputs.py --include-giab
```

## Core Scripts

| Script | Purpose |
|---|---|
| **`_common.py`** | Shared constants and helpers (no CLI). One source of truth for VAF bins, pileup depth / alt-read / HC hard-filter thresholds, deepSNV MoM parameters, BH sig levels, sample metadata (`SAMPLES`), VEP config, and small helpers (`burden_test`, `vaf_breakdown`, `classify_vaf`, `sig_star`, `normalize_bool_columns`, `load_part2_all`). **Edit here if you want to retune a threshold.** |
| **`generate_all_outputs.py`** | Master orchestrator. Runs Part2 → manual IGV validation → benchmark JSONs → cross-caller comparison → IGV inputs → CLAUDE.md update. Single entry point for reproducing all `outputs/`. |
| **`mosaic_postprocess_part2.py`** | Part2 post-processing pipeline. Takes Part1 CSV (from Terra), applies pileup validation, hard filters, centromere/repeat exclusion, chrY exclusion, gnomAD AF=0 (VEP), VAF binning, deepSNV beta-binomial test. Produces `_part2_all.tsv`, `_final.tsv`, `_mosaic.tsv`, `_het.tsv`, `_summary.txt`. |
| **`apply_igv_validation.py`** | Adds four columns (`IGV_is_real`, `IGV_notes`, `igvFilter`, `PASS_all_igv_confirmed`) to each `_part2_all.tsv` by joining with manual IGV review (`outputs/igv_review/all_trios_igv_validated.csv`). Idempotent. |
| **`binomial_het_test.py`** | Reviewer 1 (R1.m14) per-variant test. Adds three columns to each `_part2_all.tsv`: `binom_het_pval` (binomtest H₀:VAF=0.5, one-sided), `binom_het_pval_BH` (BH-adjusted on PASS_all subset), `binom_class` — depth-aware replacement for the fixed VAF=0.3741 cutoff. Keeps the old `VAF_class` alongside for comparison. Idempotent. |
| **`vaf_histograms.py`** | Plots VAF distribution per trio (2.5% bins, 0-70%) coloured by `binom_class` (mosaic = red, het = blue). Candidate set: `PASS_all_igv_confirmed`. Reference lines at old VAF cutoff (0.3741) and germline het (0.5). Outputs `outputs/benchmark/vaf_histograms.png` (composite 5×2 grid) + per-caller standalone figures. Embedded in `report.md` §2d. |
| **`cross_caller_comparison.py`** | 3-way comparison: new HC vs Mutect2 vs old pipeline (bcbio). Produces master CSVs (95 columns with per-caller annotations + manual IGV reviewer columns), Venn diagrams, intersection tables. |
| **`sensitivity_alt_threshold.py`** | Reviewer 1 (R1.5) sensitivity analysis. Re-classifies existing `_part2_all.tsv` at alt-read thresholds K ∈ {2,3,4,5} without re-running pileup/VEP. Outputs `outputs/sensitivity/alt_threshold.{json,md}` with mosaic/het counts per threshold + burden test vs control per K. |

