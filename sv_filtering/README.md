# SV filtering

Study-specific filtering of the structural variants called with ARC-SV v0.9.6
(`arcsv/`), as described in the Methods ("Structural variants") and in Fig. 4.
The filtering was done with `bedtools intersect` (BEDTools) and manual value
selection; there is no separate script. This folder lists the steps and contains
the annotation files that were used.

## Annotation files (`annotations/`, hg38)

| File | Content | Filter |
|---|---|---|
| `hg38_str_sim_lc_sat_merged.bed` | simple repeats, low-complexity regions, satellite repeats and simple tandem repeats (merged) | 1 |
| `hg38_segdups.bedpe` | segmental duplications, as pairs; every pair is listed in both orientations | 2 |
| `hg38_svblacklist.bed` | SV blacklist regions (Methods: https://cf.10xgenomics.com/supp/genome/GRCh38/sv_blacklist.bed) | 4 |
| `hg38_gap.bed` | reference assembly gaps (UCSC `gap` table; identical to `arcsv/resources/hg38_gap.bed`) | 5 |

Verify the files with `(cd annotations && shasum -a 256 -c SHA256SUMS)`.

## Exclusion filters

Input: ARC-SV calls merged across chromosomes (`arcsv_out.tab`). Its first three
columns (`chrom`, `minbp`, `maxbp`) are 0-based and BED-compatible. Each filter
removes calls independently, so the order does not matter.

```bash
A=annotations
tail -n +2 arcsv_out.tab > calls.bed      # drop the header line

# (1) >85% overlap with simple repeat, low-complexity, satellite and simple tandem repeat annotations
bedtools intersect -f 0.85 -v -a calls.bed -b $A/hg38_str_sim_lc_sat_merged.bed > f1.bed

# (2) >50% overlap with segmental duplications (first end of each pair covers all segdups)
cut -f1-3 $A/hg38_segdups.bedpe > segdups.bed
bedtools intersect -f 0.5 -v -a f1.bed -b segdups.bed > f2.bed

# (3) span < 50 bp
awk -F'\t' '$3 - $2 >= 50' f2.bed > f3.bed

# (4) overlap with SV blacklist regions
bedtools intersect -f 0.5 -v -a f3.bed -b $A/hg38_svblacklist.bed > f4.bed

# (5) overlap with reference gaps
bedtools intersect -v -a f4.bed -b $A/hg38_gap.bed > f5.bed
```

The commands for filters 1 and 4 are the ones used in the study. Filters 2, 3
and 5 follow the thresholds given in the Methods.

## Selection of putative somatic SVs (manual)

From the remaining calls (deep-coverage alignments), keep SVs that are:

1. 500 bp to 10 kb in size, with an estimated VAF of 5–50%;
2. not overlapping (≥50% reciprocal overlap, e.g. `bedtools intersect -f 0.5 -r -v`)
   with SVs previously detected in the human population. The population SVs are
   those reported in Zhou et al. 2024 (PMID 39353437); the file is not included;
3. unique to the sample in which they were detected, i.e., not detected in any
   other sample of this study.
