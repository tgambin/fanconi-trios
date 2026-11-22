# Mutation Spectrum Analysis (7-bin + SBS96, somatic vs germline)

An R script that classifies confirmed true-positive single-nucleotide variants (SNVs)
into a **7-bin mutation spectrum** and a full **96-channel SBS profile** using
[MutationalPatterns](https://bioconductor.org/packages/MutationalPatterns/)
and the human reference genome (hg38).

The script computes the spectrum two independent ways — by collapsing the
MutationalPatterns 96-channel matrix, and by direct trinucleotide classification —
and cross-checks that they agree, so the counts are self-validating. It then
compares the 7-bin spectrum of the somatic (mosaic) SNVs of each sample with a
pooled germline reference (Table 2 and Supplemental Figs. S2–S3 of the paper).

---

## What it does

Given a master table of manually reviewed variants, the script:

1. **Loads and filters** the master CSV, keeping SNVs with single-base REF/ALT
   alleles (the table already contains only confirmed true positives).
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
11. **Builds the somatic and germline sets.** Somatic: the mosaic SNVs of each of
    the five trios (`consensus_class == "mosaic"`) and of the published RAD51C
    case. Germline reference: the heterozygous de novo SNVs (`consensus_class ==
    "het"`) of the five trios and of the RAD51C case, pooled (n = 239). Sites
    whose REF does not match hg38 are dropped and logged.
12. **Counts the 7 bins** for every somatic sample and for the pooled germline
    reference.
13. **Goodness-of-fit against the germline reference:** Pearson chi-square test
    with a Monte Carlo p-value (B = 100,000, fixed seed; the smallest attainable
    p is about 1/(B+1)), Benjamini–Hochberg adjustment across samples, cosine
    similarity and Jensen–Shannon divergence (in nats) to the germline spectrum.
14. **Pairwise homogeneity** between somatic samples: Fisher–Freeman–Halton test
    on each 2 × 7 table (Monte Carlo p-value), BH-adjusted across all pairs.

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

Three CSV files. The master table of reviewed trio variants needs these columns:

| Column | Type | Notes |
|--------|------|-------|
| `CHROM` | string | Chromosome (any style; converted to UCSC `chr…` internally) |
| `POS` | integer | 1-based genomic position on **hg38** |
| `REF` | string | Reference allele — single base for SNVs |
| `ALT` | string | Alternate allele — single base for SNVs |
| `variant_type` | string | Variant class; only rows equal to `"SNV"` are kept |
| `trio` | string | Sample/group code — must match the `trio_map` keys (see below) |
| `consensus_class` | string | `mosaic` or `het` (heterozygous de novo); used from Step 11 |

The table must contain only manually confirmed true-positive variants.

The two RAD51C tables (mosaic and heterozygous SNVs of the published case) need
the columns `Chr`, `POS`, `Ref` and `Alt`.

Positions and alleles must be on the **hg38** build to match the reference genome.

### Data availability

The input tables are not distributed with this repository because they contain
participant-level variants. The trio variants are listed in Supplemental Tables S2
(mosaic) and S3 (heterozygous de novo) of the paper. The RAD51C variants are from
the published case: Zemet R, Du H, Gambin T, Lupski JR, Liu P, Stankiewicz P. 2023.
SNV/indel hypermutator phenotype in biallelic RAD51C variant: Fanconi anemia.
*Hum Genet* 142: 721–733. doi:10.1007/s00439-023-02550-4.

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
   rad51c_mosaic_csv <- "mosaic-vars-rad51c_truelist.csv"
   rad51c_het_csv    <- "het-vars-rad51c_truelist.csv"
   n_monte_carlo <- 1e5                   # Monte Carlo replicates
   rng_seed      <- 1                     # fixed seed: reproducible p-values
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
| `somatic_germline/germline_reference_composition.csv` | Germline reference: SNVs per source and total (239) |
| `somatic_germline/ref_check_log.csv` | Sites dropped because REF does not match hg38 |
| `somatic_germline/seven_bin_counts.csv` | 7-bin counts, somatic samples and pooled germline |
| `somatic_germline/gof_vs_germline.csv` | Goodness-of-fit p-value, BH-adjusted p-value, cosine similarity and JSD per sample (Table 2) |
| `somatic_germline/pairwise_homogeneity_qvals.csv` | BH-adjusted q-values of the pairwise homogeneity tests |

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
