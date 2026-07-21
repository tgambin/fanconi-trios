# Scripts

Analysis scripts for the Fanconi Anemia somatic mosaicism study.

## Master Pipeline

The master pipeline expects Terra Part1 CSVs under `outputs/part1_csv/` using
the `<sample_id>.part1.csv` and `<sample_id>.mutect2.part1.csv` names defined
in `_common.py`. Study-trio runs additionally require the controlled-access
`outputs/igv_review/all_trios_igv_validated.csv` file.
Its required columns are `trio`, `CHROM`, `POS`, `REF`, `ALT`, `IGV_is_real`,
and `IGV_notes`; `IGV_is_real` accepts `TRUE`, `FALSE`, or `NA`.

```bash
# Full local run after supplying the controlled-access Part1 CSVs
# and all_trios_igv_validated.csv
python scripts/generate_all_outputs.py

# First pass: create IGV input TSVs and a blank manual-review CSV, then stop
python scripts/generate_all_outputs.py --prepare-igv --no-giab

# Regenerate downstream reports from Part2 tables created by a previous run
python scripts/generate_all_outputs.py --skip-part2

# Omit GIAB when its two Part1 files were not supplied
python scripts/generate_all_outputs.py --no-giab
```

## Core Scripts

| Script | Purpose |
|---|---|
| **`_common.py`** | Shared constants and helpers (no CLI). One source of truth for VAF bins, pileup depth / alt-read / HC hard-filter thresholds, deepSNV MoM parameters, BH sig levels, sample metadata (`SAMPLES`), VEP config, and small helpers (`burden_test`, `vaf_breakdown`, `classify_vaf`, `sig_star`, `normalize_bool_columns`, `load_part2_all`). **Edit here if you want to retune a threshold.** |
| **`generate_all_outputs.py`** | Master orchestrator. Its `--prepare-igv` pass runs Part2 and creates the cross-caller/IGV review package. After manual review, the normal or `--skip-part2` pass applies IGV decisions, runs depth-aware classification, and creates reports and final candidate tables. |
| **`mosaic_postprocess_part2.py`** | Part2 post-processing pipeline. Takes Part1 CSV (from Terra), applies pileup validation, hard filters, centromere/repeat exclusion, chrY exclusion, gnomAD AF=0 (VEP), VAF binning, deepSNV beta-binomial test. Produces `_part2_all.tsv`, `_final.tsv`, `_mosaic.tsv`, `_het.tsv`, `_summary.txt`. |
| **`apply_igv_validation.py`** | Adds four columns (`IGV_is_real`, `IGV_notes`, `igvFilter`, `PASS_all_igv_confirmed`) to each `_part2_all.tsv` by joining with manual IGV review (`outputs/igv_review/all_trios_igv_validated.csv`). Idempotent. |
| **`binomial_het_test.py`** | Per-variant binomial test. Adds three columns to each `_part2_all.tsv`: `binom_het_pval` (binomtest H₀:VAF=0.5, one-sided), `binom_het_pval_BH` (BH-adjusted on PASS_all subset), `binom_class` — depth-aware replacement for the fixed VAF=0.3741 cutoff. Keeps `VAF_class` as a diagnostic reference. Idempotent. |
| **`vaf_histograms.py`** | Plots VAF distribution per trio (2.5% bins, 0-70%) coloured by `binom_class` (mosaic = red, het = blue). Candidate set: `PASS_all_igv_confirmed`. Reference lines at the fixed VAF cutoff (0.3741) and germline het (0.5). Outputs `outputs/benchmark/vaf_histograms.png` (composite 5×2 grid) + per-caller standalone figures. Embedded in `report.md` §2d. |
| **`cross_caller_comparison.py`** | Two-caller comparison: HaplotypeCaller vs Mutect2. Produces master CSVs, Venn diagrams, and intersection tables. |
| **`sensitivity_alt_threshold.py`** | Sensitivity analysis. Re-classifies existing `_part2_all.tsv` at alt-read thresholds K ∈ {2,3,4,5} without re-running pileup/VEP. Outputs `outputs/sensitivity/alt_threshold.{json,md}` with mosaic/het counts per threshold + burden test vs control per K. |

Use the workflow manifests and the scripts listed above.
