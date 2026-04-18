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
