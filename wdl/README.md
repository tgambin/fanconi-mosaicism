# WDL Workflows

Terra/Cromwell workflows for the mosaic variant analysis pipeline.

## Active Workflows (`wdl/mosaic/`)

| WDL | Purpose | Inputs |
|---|---|---|
| **HaplotypeCallerScatter.wdl** | HC scatter (ploidy=2, 100-300 shards) | Proband or parent BAM |
| **Mutect2TumorOnlyScatter.wdl** | Mutect2 tumor-only scatter (100 shards) | Proband BAM (150x) |
| **MosaicPostProcessPart1.wdl** | HC Part1: trio VCF filter + pileup | Trio VCFs + BAMs |
| **MutectPostProcessPart1.wdl** | Mutect2 Part1: same as above for M2 | Mutect2 VCF + parent HC VCFs + BAMs |
| **DownsampleProband.wdl** | Downsample BAMs + coverage stats | BAM + target coverages |
| **AddReadGroups.wdl** | Add @RG header to BAMs (GIAB fix) | BAM without @RG |
| **DownloadBAM.wdl** | Download BAMs from URL to GCS | URL |
| **MosdepthCoverage.wdl** | Coverage statistics via mosdepth | BAM |

## IGV Screenshots

| WDL | Purpose |
|---|---|
| **IGVTrio.wdl** | BAM slicing + IGV headless screenshots for a trio (self-contained; all tasks inlined). |

## Legacy (`wdl/legacy/`)

Superseded workflows — kept for reference.

| WDL | Replaced by |
|---|---|
| MosaicVariantCalling.wdl | HaplotypeCallerScatter + Mutect2TumorOnlyScatter |
| MosaicPostProcessing.wdl | MosaicPostProcessPart1 + Part2 (local Python) |
| MosaicDownsampling.wdl | DownsampleProband.wdl |
| Mutect2Standalone.wdl | Mutect2TumorOnlyScatter.wdl |
| TrioPileup.wdl | MosaicPostProcessPart1 (NIO per-position pileup) |
