# Fanconi Anemia — Mosaic Variant Analysis

WDL workflows (Terra/Cromwell) and Python post-processing scripts for the
detection of somatic mosaic and germline *de novo* variants from trio WGS
data, as described in the accompanying manuscript ("Non-invasive
assessment of early embryonic somatic mosaicism and hypermutation
phenomenon in patients with Fanconi anemia", Gambin, Zhou et al.).

This repository is the Supplemental Code accompanying the manuscript. Raw
sequencing data and per-variant result tables are not included here; the
in-house FA trios are deposited in the European Nucleotide Archive under
the accessions reported in the manuscript, and the GIAB Ashkenazi trio
(HG002/HG003/HG004) is publicly available.

## Quick start

```bash
# Python dependencies
pip install -r requirements.txt

# VEP cache for gnomAD annotation — download once and point VEP_DATA to it
export VEP_DATA=/path/to/vep_data
docker pull ensemblorg/ensembl-vep:release_115.1

# Post-processing (runs Part 2 through report generation; requires Part 1 CSVs
# from Terra under outputs/part1_csv/ — see wdl/README.md for how to produce them)
python scripts/generate_all_outputs.py
```

## Pipeline overview

The analysis is organised in two stages. **Stage 1 (cloud)** runs on
Terra/Cromwell and performs the compute-heavy steps (variant calling,
per-position pileup validation). **Stage 2 (local)** runs in Python on
the small Stage-1 CSV outputs and performs classification, statistical
tests, manual IGV review integration and report generation.

```
                                                                            ┌─── STAGE 1 (Terra/Cromwell, WDL) ───────────────────────────┐
  raw proband BAM  ──► DownsampleProband.wdl ──► 150x proband BAM ──┬──► HaplotypeCallerScatter.wdl ───┐
                                                                   │                                     │
  parental BAMs  ────────────────────────────────────────────────────┼──► HaplotypeCallerScatter.wdl ───┤
                                                                   │    (parents, ploidy=2)             │
                                                                   │                                     │
                                                                   └──► Mutect2TumorOnlyScatter.wdl     │
                                                                        (proband, tumour-only)          │
                                                                                                         │
                                                                        ┌────────────────────────────────┘
                                                                        │
           MosaicPostProcessPart1.wdl (HC)  ◄─────────────────────────┤   trio VCFs +
           MutectPostProcessPart1.wdl  (Mutect2) ◄─────────────────────┘   slice BAMs
              │  trio de novo filter (bcftools isec -C)
              │  canonical chr + DP pre-filter
              │  segdup / LCR region masking
              │  per-position pileup (samtools mpileup)
              ▼
           Part 1 CSV ──────────────────────────────────────────────┐
                                                                      │
                            ┌──── STAGE 2 (local, Python) ────────────┘
                            │
                            ▼
        mosaic_postprocess_part2.py
          │  hard filters (MQ, QD, FS, ReadPosRankSum)
          │  centromere / repeat / chrY exclusion
          │  gnomAD AF=0 annotation via VEP
          │  VAF classification + deepSNV LRT (Gerstung 2014)
          ▼
        binomial_het_test.py        ──►  depth-aware mosaic / het classifier
          │                               (binomtest H₀:VAF=0.5, BH < 0.05)
          ▼
        IGVTrio.wdl  ──►  trio screenshots  ──►  manual TRUE/FALSE scoring
                                                         │
                                                         ▼
        apply_igv_validation.py      ──►  PASS_all_igv_confirmed
          │
          ▼
        generate_all_outputs.py      ──►  per-trio burden test vs control
                                          ──►  sensitivity_alt_threshold.py
                                          ──►  cross_caller_comparison.py
                                          ──►  vaf_histograms.py
                                          ──►  final report + CSV tables
```

## Analysis steps

1. **Proband downsampling to a shared target depth** (`wdl/mosaic/DownsampleProband.wdl`).
   Probands are sequenced at ≥150x; parents at ~40x. Each proband BAM is
   downsampled to a shared 150x target with `samtools view -s` (fixed
   random seed), and post-downsample depth is measured with
   `samtools depth`. Matched-depth BAMs are the input for every step that
   follows.

2. **Variant calling on the depth-matched BAMs** —
   `HaplotypeCallerScatter.wdl` (ploidy=2, scattered 100–300-way) for
   both probands and parents; `Mutect2TumorOnlyScatter.wdl` (100-way
   scatter) for the probands in parallel.

3. **Trio-based *de novo* filtering and per-position pileup validation**
   (Part 1, `wdl/mosaic/{Mosaic,Mutect}PostProcessPart1.wdl`). Proband
   sites absent from both parents' VCFs are kept (`bcftools isec -C -w 1`),
   restricted to canonical chromosomes, masked against segmental
   duplications and low-complexity regions, and re-examined with
   per-position pileup (`samtools mpileup -r CHROM:POS-POS` across the
   proband and both parents, parallelised with a Python
   `ThreadPoolExecutor`). Output: one Part 1 CSV per trio × caller with
   VCF annotations plus per-position pileup counts.

4. **Local Part 2 post-processing** (`scripts/mosaic_postprocess_part2.py`,
   orchestrated by `scripts/generate_all_outputs.py`). Applies the
   pileup-based parental-absence filter, the proband alt-read threshold,
   the GATK-style hard filters (SNV: `INFO/MQ` > 30, `QD` > 2, `FS` < 60,
   `ReadPosRankSum` > −8; INDEL: `QD` > 2, `FS` < 200, `ReadPosRankSum`
   > −20), centromere / repeat / chrY exclusion, and the gnomAD AF = 0
   filter (via Ensembl VEP). Produces a `_part2_all.tsv` per trio ×
   caller with all filter columns retained.

5. **Depth-aware mosaic / heterozygous classification**
   (`scripts/binomial_het_test.py`). For each surviving candidate, a
   per-variant one-sided exact binomial test of `AD_alt` against `DP`
   with H₀: VAF = 0.5 and H₁: VAF < 0.5 (Benjamini-Hochberg adjusted over
   the PASS-filter candidates per trio × caller). A variant is called
   `mosaic` when H₀ is rejected at BH < 0.05 and `het` otherwise. The
   classification is stored as the `binom_class` column.

6. **Per-variant beta-binomial LRT against sequencing noise** (Gerstung
   et al., 2014), fitted on the pooled parental alt-read distribution
   (method-of-moments error rate and overdispersion), BH-corrected per
   trio × caller. Raw and adjusted p-values are retained as
   `deepSNV_pval` and `deepSNV_pval_BH`.

7. **Manual IGV review** (`wdl/IGVTrio.wdl` +
   `scripts/apply_igv_validation.py`). Trio IGV screenshots
   (proband + mother + father, same locus, same window) are generated in
   batch for every pipeline-passing candidate and scored manually as
   TRUE (real call) or FALSE (artefact), with an error-class code p1–p6
   assigned to each FALSE call. The scoring CSV is joined back into the
   `_part2_all.tsv` and a `PASS_all_igv_confirmed` column is emitted.
   All counts reported in the manuscript are computed on the
   IGV-confirmed subset.

8. **Statistical tests and robustness analyses** (all in
   `scripts/generate_all_outputs.py`):
   - **Burden test** — pairwise conditional binomial exact test
     (one-sided, case > control) for the per-proband mosaic-SNV count of
     each FA trio vs the negative control.
   - **Alt-read sensitivity** (`scripts/sensitivity_alt_threshold.py`) —
     re-classification of the full candidate set at K ∈ {3, 4, 5}
     confirms that the main finding is not sensitive to the alt-read
     cut-off.
   - **Cross-caller intersection** (`scripts/cross_caller_comparison.py`)
     — HaplotypeCaller ∩ Mutect2 ∩ previous-pipeline per-variant
     intersection with Venn diagrams.
   - **VAF distributions** (`scripts/vaf_histograms.py`) — per-trio
     2.5%-bin VAF histograms stacked by the binomial classifier.

9. **Report generation** (`scripts/generate_all_outputs.py`). Produces
   the final summary tables (coverage, filter cascade, deepSNV
   confirmation rates, burden p-values, cross-caller intersections) and
   writes them as both markdown and CSV alongside each analysis output.

## Reproducing the analysis

Stages 1 and 2 are independent; Stage 2 operates on the Stage-1 CSVs and
can be re-run without repeating the expensive cloud steps.

```bash
# --- Stage 1 (Terra / Cromwell) ------------------------------------
# Submit each WDL with the inputs JSON format documented in wdl/README.md.
# The output of Part 1 is one CSV per trio × caller, which needs to be
# staged under outputs/part1_csv/ for Stage 2.

# --- Stage 2 (local) -----------------------------------------------
export VEP_DATA=/path/to/vep_data
docker pull ensemblorg/ensembl-vep:release_115.1

python scripts/generate_all_outputs.py                 # full Part 2 + report
python scripts/generate_all_outputs.py --skip-part2    # re-render report only
python scripts/generate_all_outputs.py --no-giab       # FA trios only
```

## Repository layout

```
wdl/
├── mosaic/        Active WDL workflows (HaplotypeCaller, Mutect2, Part 1
│                  post-processing, downsampling, read-group fix, BAM fetch)
├── legacy/        Earlier WDLs, superseded (kept for reference)
└── IGVTrio.wdl    IGV trio screenshots (used for manual review)

scripts/
├── _common.py                       Shared constants and helpers
├── mosaic_postprocess_part2.py      Local Part 2 post-processing
├── apply_igv_validation.py          Apply manual IGV review annotations
├── binomial_het_test.py             Depth-aware mosaic/het classifier
├── sensitivity_alt_threshold.py     Alt-read-threshold sensitivity analysis
├── vaf_histograms.py                Per-trio VAF histograms
├── cross_caller_comparison.py       HaplotypeCaller vs Mutect2 vs prior-pipeline
└── generate_all_outputs.py          Master orchestrator (calls the above)
```

See `wdl/README.md` and `scripts/README.md` for per-file details.

## Prerequisites

- Python 3.12+ (`pip install -r requirements.txt`)
- Docker (for VEP via `ensemblorg/ensembl-vep:release_115.1`)
- Cromwell (for local WDL execution) or Terra access (for cloud execution)
- Ensembl VEP GRCh38 cache (~25 GB) — path set via `VEP_DATA` env var

## Tech stack

- **WDL/Cromwell** on Terra (Google Cloud) for variant calling and Part 1
  post-processing (HaplotypeCaller ploidy=2, Mutect2 tumour-only).
- **Python** (polars, polars-bio, scipy) for local Part 2 analysis,
  statistical classification, cross-caller comparison and report
  generation.
- **GATK** HaplotypeCaller + Mutect2, `bcftools`, `samtools`, Ensembl VEP.
- Reference: GRCh38 / hg38.

## Sample identifiers

The `SAMPLES` dict in `scripts/_common.py` carries generic placeholder
identifiers for the in-house probands and parents. Edit the dict so that
`sample_id`, `mother` and `father` fields match the file-name stems used
by your own BAM / VCF outputs under `outputs/part1_csv/` and
`outputs/coverage/`. The GIAB row keeps the publicly available HG002 /
HG003 / HG004 accessions.

## License and citation

If you use this pipeline or any of the workflows, please cite the
accompanying manuscript. Authorship, acknowledgments and funding are
listed in the manuscript.
