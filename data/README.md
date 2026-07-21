# Data Directory

Reference annotations and instructions for optional benchmarks. Variant result tables are intentionally excluded.

## Structure

```
data/
├── annotations/          [IN GIT]
│   ├── centromeres.bed               Centromere regions (GRCh38)
│   ├── repeats.bed                   Alu + centromere repeat regions
│   ├── LCR-hs38.bed.gz              Low Complexity Regions (Heng Li)
│   ├── GRCh38GenomicSuperDup.chr14.bed  Segmental duplications (chr14 test)
│   └── SHA256SUMS                    Checksums for included annotations
│
├── metadata/             [IN GIT]
│   └── intervals.list                WGS calling intervals (chr14 test)
│
├── giab_benchmark/       [BENCHMARK VCFs DOWNLOADED LOCALLY]
│   ├── HG002_benchmark.vcf.gz       GIAB v4.2.1 benchmark — HG002 (son)
│   ├── HG003_benchmark.vcf.gz       GIAB v4.2.1 benchmark — HG003 (father)
│   ├── HG004_benchmark.vcf.gz       GIAB v4.2.1 benchmark — HG004 (mother)
│   ├── giab_denovo_truth.vcf.gz      Derived locally; not tracked
│   └── SHA256SUMS                     Checksums for downloaded VCFs and indexes [IN GIT]
```

## Notes

- `centromeres.bed` and `repeats.bed` are used by local Part2; the LCR and
  segmental-duplication annotations are used by the Part1 WDL workflows.
- Variant result tables and truth positions are private and must be supplied separately when reproducing analyses.

No participant-level variants, BAMs, VCFs, or generated result tables are
distributed in this directory.

## Segmental-duplication BED for Terra

The repository contains a chr14 subset for local testing. Production Part1
manifests expect the complete GRCh38 UCSC `genomicSuperDups` table converted to
three-column BED and uploaded to the user's Terra bucket:

```bash
curl -L \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz \
  -o genomicSuperDups.txt.gz
gzip -cd genomicSuperDups.txt.gz | cut -f2-4 \
  > GRCh38GenomicSuperDup.bed
gsutil cp GRCh38GenomicSuperDup.bed \
  "gs://${TERRA_BUCKET}/annotations/GRCh38GenomicSuperDup.bed"
```

`TERRA_BUCKET` is the bucket name without the `gs://` prefix. The included LCR,
centromere, and repeat files are the exact annotations consumed across Part1
and Part2.
Verify the included annotation files with:

```bash
(cd data/annotations && sha256sum -c SHA256SUMS)
```

## Obtaining the public GIAB trio BAMs

The repository does not redistribute the large GIAB alignments. The three
`inputs/DownloadBAM/DownloadBAM.giab_*.terra.json` manifests contain the
official NCBI URLs for the GRCh38 300x Ashkenazi trio BAMs and indexes. In
Terra, run the included workflows in this order:

1. `DownloadBAM.wdl` with each of the three `DownloadBAM` manifests.
2. `DownsampleProband.wdl` with the corresponding GIAB manifest: HG002 to
   150x and HG003/HG004 to 40x.
3. `AddReadGroups.wdl` with the corresponding GIAB manifest.
4. Use the resulting `HG002.rg.bam`, `HG003.rg.bam`, and `HG004.rg.bam` in
   the GIAB HaplotypeCaller, Mutect2, and Part1 manifests.

The GCS paths in the manifests describe the expected layout. After each run,
copy the output there or update the next manifest to the actual Terra output
URI.

## Downloading the GIAB benchmark files

The GIAB benchmark VCFs are public reference data and are intentionally not
stored in this repository. Download the NIST Genome in a Bottle v4.2.1 small-
variant benchmarks for the Ashkenazi trio, GRCh38, chromosomes 1-22, from the
official NCBI-hosted GIAB release:

```bash
mkdir -p data/giab_benchmark

base=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio

curl -L "$base/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
  -o data/giab_benchmark/HG002_benchmark.vcf.gz
curl -L "$base/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
  -o data/giab_benchmark/HG002_benchmark.vcf.gz.tbi

curl -L "$base/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
  -o data/giab_benchmark/HG003_benchmark.vcf.gz
curl -L "$base/HG003_NA24149_father/NISTv4.2.1/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
  -o data/giab_benchmark/HG003_benchmark.vcf.gz.tbi

curl -L "$base/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
  -o data/giab_benchmark/HG004_benchmark.vcf.gz
curl -L "$base/HG004_NA24143_mother/NISTv4.2.1/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
  -o data/giab_benchmark/HG004_benchmark.vcf.gz.tbi
```

Verify the downloads from `data/giab_benchmark/` on Linux:

```bash
(cd data/giab_benchmark && sha256sum -c SHA256SUMS)
```

Create the proband-private truth set used by the project (requires
`bcftools`):

```bash
cd data/giab_benchmark
bcftools isec -C -w 1 \
  HG002_benchmark.vcf.gz \
  HG003_benchmark.vcf.gz \
  HG004_benchmark.vcf.gz \
  -O z -o giab_denovo_truth.vcf.gz
bcftools index -t giab_denovo_truth.vcf.gz
```

See the [NIST Genome in a Bottle project](https://www.nist.gov/programs-projects/genome-bottle)
and cite the benchmark release as described by GIAB when reusing these data.
