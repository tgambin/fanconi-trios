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

## IGV Screenshots (`wdl/`)

| WDL | Purpose |
|---|---|
| **IGVTrio.wdl** | BAM slicing + IGV screenshots for trio |

`IGVTrio.wdl` is self-contained; its tasks are defined in the same file.

All GATK tasks are pinned to `broadinstitute/gatk:4.6.2.0`. Samtools tasks
use `staphb/samtools:1.19`. The root README describes workflow order, required
data, and rendering input JSON templates before Terra upload.
