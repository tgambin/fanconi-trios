# Mutation Spectrum Analysis (7-bin + SBS96)

An R script that classifies confirmed true-positive single-nucleotide variants (SNVs)
into a **7-bin mutation spectrum** and a full **96-channel SBS profile** using
[MutationalPatterns](https://bioconductor.org/packages/MutationalPatterns/)
and the human reference genome (hg38).

The script computes the spectrum two independent ways — by collapsing the
MutationalPatterns 96-channel matrix, and by direct trinucleotide classification —
and cross-checks that they agree, so the counts are self-validating.

---

## What it does

Given a master table of manually reviewed variants, the script:

1. **Loads and filters** the master CSV, keeping only manually confirmed
   true-positive SNVs with single-base REF/ALT alleles.
2. **Builds `GRanges`** objects with UCSC-style sequence names.
3. **Extracts the trinucleotide context** of every variant from hg38.
4. **Classifies** each SNV into one of 7 bins (pyrimidine-normalized; `C>T`
   split by CpG context — see below).
5. **Generates the 96-channel SBS profile** via MutationalPatterns (writes
   per-sample VCFs, reads them back as `GRanges`, and builds the mutation matrix).
6. **Collapses** the 96 channels down to the 7 bins.
7. **Independently counts** the 7 bins directly from the trinucleotide
   classification.
8. **Plots** per-sample 7-bin bar charts for both routes.
9. **Computes 96-channel cosine similarity** between samples and draws a heatmap.
10. **Cross-checks** the MutationalPatterns-collapsed counts against the manual
    counts and reports any mismatch.

### The 7 bins

| Bin | Description |
|-----|-------------|
| `C>A` | C→A (and G→T on the complementary strand) |
| `C>G` | C→G / G→C |
| `C>T@CpG` | C→T where the 3′ base is G (CpG dinucleotide) |
| `C>T_nonCpG` | C→T outside a CpG context |
| `T>A` | T→A / A→T |
| `T>C` | T→C / A→G |
| `T>G` | T→G / A→C |

All mutations are represented from the pyrimidine (C or T) strand — variants with
a purine (A/G) reference are reverse-complemented before binning. Splitting `C>T`
by CpG context isolates the elevated deamination signal at methylated CpG sites.

---

## Requirements

- **R** (≥ 4.1, for the native `|>` pipe)
- **Bioconductor** packages:
  ```r
  if (!require("BiocManager")) install.packages("BiocManager")
  BiocManager::install(c("MutationalPatterns", "BSgenome.Hsapiens.UCSC.hg38"))
  ```
- **CRAN** packages:
  ```r
  install.packages(c("tidyverse", "ggplot2"))
  ```
  (`GenomicRanges` is pulled in as a MutationalPatterns dependency.)

---

## Input format

A single CSV of reviewed variants. The following columns are required:

| Column | Type | Notes |
|--------|------|-------|
| `CHROM` | string | Chromosome (any style; converted to UCSC `chr…` internally) |
| `POS` | integer | 1-based genomic position on **hg38** |
| `REF` | string | Reference allele — single base for SNVs |
| `ALT` | string | Alternate allele — single base for SNVs |
| `variant_type` | string | Variant class; only rows equal to `"SNV"` are kept |
| `trio` | string | Sample/group code — must match the `trio_map` keys (see below) |
| `is_true_positive` | logical | `TRUE` for manually confirmed true positives; only these are kept |

Positions and alleles must be on the **hg38** build to match the reference genome.

### Sample labels

Sample display names come from a small lookup defined near the top of the script:

```r
trio_map <- tibble(
  trio        = c("brca1","brca2","fanca","fanct","control"),
  sample_name = c("BRCA1","BRCA2","FANCA","FANCT","Control")
)
```

Edit `trio_map` so that the `trio` values match those in your CSV and
`sample_name` holds the labels you want on the plots.

---

## Usage

1. Open the script and set the paths in the **Configuration** block:
   ```r
   input_csv  <- "all_trios_master.csv"   # your master table
   output_dir <- "mutpatterns_output"     # where results are written
   ```
   Relative paths resolve against the directory you launch the script from;
   absolute paths also work.

2. Run it:
   ```bash
   Rscript mutation_spectrum_mutpatterns.R
   ```
   or, from an interactive R session, `source("mutation_spectrum_mutpatterns.R")`.

Progress is printed to the console, and all results are written to `output_dir/`.

---

## Outputs

All files are written under `output_dir/` (default `mutpatterns_output/`):

| File | Description |
|------|-------------|
| `sbs96_matrix_mutpatterns.csv` | 96-channel SBS count matrix (96 types × samples) |
| `sbs96_profile_mutpatterns.pdf` | 96-channel SBS profile bar plot |
| `seven_bin_spectrum_mutpatterns.csv` | 7-bin counts, from collapsing the 96-channel matrix |
| `seven_bin_spectrum_manual.csv` | 7-bin counts, from direct trinucleotide classification |
| `spectrum_7bin_mutpatterns.pdf` | 7-bin bar plot (MutationalPatterns route) |
| `spectrum_7bin_manual.pdf` | 7-bin bar plot (manual route) |
| `cosine_similarity_96.csv` | Pairwise 96-channel cosine-similarity matrix |
| `cosine_similarity_heatmap.pdf` | Heatmap of the cosine-similarity matrix |
| `comparison_mp_vs_manual.csv` | Per-bin, per-sample MP vs manual counts with match/difference |
| `vcf_files/` | Per-sample VCF intermediates generated for MutationalPatterns |

---

## Methodology notes

- **Reference genome:** hg38 (`BSgenome.Hsapiens.UCSC.hg38`). Trinucleotide
  contexts are read directly from this build.
- **DBS/MBS handling:** `read_vcfs_as_granges(..., predefined_dbs_mbs = TRUE)`
  prevents MutationalPatterns from auto-merging adjacent SNVs into doublet/multi-base
  substitutions, so every manually reviewed SNV is preserved in the 96-channel matrix.
- **Self-validation:** the 7-bin spectrum is derived twice — once by collapsing the
  MutationalPatterns 96-channel matrix and once by classifying trinucleotides
  directly. `comparison_mp_vs_manual.csv` reports any discrepancy; a clean run shows
  all bins matching.

---
