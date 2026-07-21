# Fanconi Anemia — Mosaic Variant Analysis

WDL workflows for Terra/Cromwell and Python post-processing scripts used to
identify mosaic de novo variants in Fanconi anemia trios. The repository is a
supplemental-code package: it contains code, workflow manifests, small reference
annotations, and public-data download instructions, but no participant variants,
BAM/VCF files, manual IGV decisions, or generated results.

## What can be reproduced

- The WDL calling and Part1 post-processing methods can be run on appropriately
  authorized trio BAMs or on the public GIAB Ashkenazi trio.
- The local Part2 filtering and report generation can reproduce the study
  analysis after the controlled-access Part1 CSVs and manual IGV review table
  have been supplied.
- Public GIAB inputs and benchmark download instructions provide a non-study
  dataset for exercising the workflow. They do not reproduce the study cohort
  results.

## Repository layout

```text
wdl/mosaic/          Terra/Cromwell calling and Part1 workflows
wdl/IGVTrio.wdl      BAM slicing and headless IGV screenshots
inputs/              Terra and local input-manifest templates
scripts/             Local Part2 filtering, statistics, and reports
data/annotations/    Small GRCh38 filtering annotations
data/metadata/       Non-variant test intervals only
data/giab_benchmark/ Checksums and GIAB download instructions
```

Generated files are written under `outputs/`, which is intentionally ignored
and must not be committed.

## Software setup

Python 3.13 and the dependency versions pinned in `requirements.txt` were
validated for this package. Docker is required for local VEP annotation and
for building the IGV image used by Terra.
Rendering the supplied JSON
templates uses `envsubst` (from GNU gettext); uploading reference annotations
from the command line uses `gsutil` from the Google Cloud CLI.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The gnomAD filtering step uses VEP 115.1 with its matching GRCh38 cache:

```bash
mkdir -p vep_data
docker pull ensemblorg/ensembl-vep:release_115.1
docker run --rm -it \
  -v "$PWD/vep_data:/data" \
  ensemblorg/ensembl-vep:release_115.1 \
  INSTALL.pl -a cf -s homo_sapiens -y GRCh38
export VEP_DATA="$PWD/vep_data"
```

`--no-vep` is available for software testing only. It disables the gnomAD
filter and therefore does not reproduce the reported analysis.

For IGV screenshots, build the versioned image from the included Dockerfile,
push it to a registry accessible by Terra, and use that full image tag when
rendering the IGV input manifests:

```bash
export IGV_DOCKER_IMAGE=ghcr.io/your-account/igv-headless:2.16.2
docker build -t "$IGV_DOCKER_IMAGE" .
docker push "$IGV_DOCKER_IMAGE"
```

## Data availability

The repository does not contain participant-level sequence or variant data.
To reproduce the study analysis, obtain the authorized BAMs through the data
access process described in the associated manuscript and supply the generated
Part1 and IGV-review files locally.

Public GRCh38 references in the Terra manifests point to the Broad public GCS
reference bucket. GIAB BAM and benchmark acquisition is described in
[`data/README.md`](data/README.md). The annotations consumed directly by the
Python Part2 scripts are included in `data/annotations/`.

## Running the Terra workflows

1. Import the required WDL from `wdl/mosaic/` into Terra.
2. Choose the corresponding JSON from `inputs/<WorkflowName>/`.
3. Render `${TERRA_BUCKET}` before uploading the JSON. The value is the bucket
   name without the `gs://` prefix. IGV manifests additionally require
   `${IGV_DOCKER_IMAGE}`:

   ```bash
   export TERRA_BUCKET=my-terra-workspace-bucket
   envsubst '$TERRA_BUCKET' \
     < inputs/HaplotypeCallerScatter/HaplotypeCallerScatter.giab_proband.terra.json \
     > /tmp/HaplotypeCallerScatter.giab_proband.json
   ```

4. Replace any cohort-specific object paths with locations available to your
   Terra workspace. `${TERRA_BUCKET}` is a template placeholder, not a Terra
   runtime variable.
5. Run the workflows in this order:

   ```text
   DownloadBAM (public GIAB only)
       → DownsampleProband
       → AddReadGroups, if required
       → HaplotypeCallerScatter for proband and both parents
       → Mutect2TumorOnlyScatter for the proband
       → MosaicPostProcessPart1 and MutectPostProcessPart1
       → download both Part1 CSVs and their filter reports
   ```

After each Terra run, either put the returned object at the path expected by
the next manifest or replace that manifest value with the actual output URI.
The template paths are not an automatic workflow linkage.

The production workflow expects the full GRCh38 segmental-duplication BED to be
available as `${TERRA_BUCKET}/annotations/GRCh38GenomicSuperDup.bed`. Instructions
for producing it are in [`data/README.md`](data/README.md).

## Running local Part2 and reports

Copy Terra outputs into the following local layout. The exact sample IDs used
by the study are defined once in `scripts/_common.py`.

```text
outputs/
├── part1_csv/
│   ├── <sample_id>.part1.csv
│   └── <sample_id>.mutect2.part1.csv
├── part1_reports/                         # optional report inputs
│   ├── <trio>_hc_filter_report.txt
│   └── <trio>_mutect2_filter_report.txt
└── igv_review/
    └── all_trios_igv_validated.csv        # controlled-access study trios
```

For the complete authorized study dataset:

```bash
export VEP_DATA=/absolute/path/to/vep_data
python scripts/generate_all_outputs.py
```

That single command assumes the completed controlled-access IGV table is
already present. For a first run, prepare the manual-review package instead:

```bash
python scripts/generate_all_outputs.py --prepare-igv --no-giab
gsutil -m cp outputs/cross_caller/*_mosaic_het_igv.tsv \
  "gs://${TERRA_BUCKET}/igv_inputs/"
```

Render the five manifests in `inputs/IGVTrio/` with both `TERRA_BUCKET` and
`IGV_DOCKER_IMAGE`, for example:

```bash
envsubst '$TERRA_BUCKET $IGV_DOCKER_IMAGE' \
  < inputs/IGVTrio/IGVTrio.brca1_mosaic_het.terra.json \
  > /tmp/IGVTrio.brca1.json
```

Run `wdl/IGVTrio.wdl` in Terra and inspect the returned screenshot archives.
Start with
`outputs/igv_review/all_trios_master.csv`, fill `IGV_is_real` with
`TRUE`, `FALSE`, or `NA` (and optionally `IGV_notes`), then save the controlled
table as `outputs/igv_review/all_trios_igv_validated.csv`. Resume without
re-running Part2:

```bash
python scripts/generate_all_outputs.py --skip-part2 --no-giab
```

Both IGV CSVs contain per-variant study data. They remain under the ignored
`outputs/` directory and must not be committed.

Useful variants:

```bash
# Exclude the GIAB control when its Part1 CSVs were not supplied
python scripts/generate_all_outputs.py --no-giab

# Rebuild reports from existing outputs/part2 tables
python scripts/generate_all_outputs.py --skip-part2
```

The generator stops with a concrete list of missing inputs. Study-trio runs
require `all_trios_igv_validated.csv`; this private table is intentionally not
distributed. Cross-caller analysis compares HaplotypeCaller with Mutect2 only.

For a public GIAB-only Part2 test, place the two downloaded Terra outputs at
`outputs/part1_csv/HG002_giab_150x.part1.csv` and
`outputs/part1_csv/HG002_giab_150x.mutect2.part1.csv`, then run:

```bash
python scripts/mosaic_postprocess_part2.py \
  outputs/part1_csv/HG002_giab_150x.part1.csv \
  --sample giab --caller hc
python scripts/mosaic_postprocess_part2.py \
  outputs/part1_csv/HG002_giab_150x.mutect2.part1.csv \
  --sample giab_mutect2 --caller mutect2
python scripts/binomial_het_test.py --trios giab
```

This path exercises Part2 without the private IGV table. It still requires the
matching VEP cache unless `--no-vep` is used for a non-equivalent smoke test.

Directory-specific details are documented in `data/README.md`,
`inputs/README.md`, `scripts/README.md`, and `wdl/README.md`.
