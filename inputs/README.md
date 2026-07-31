# Input JSONs

Input configurations for WDL workflows on Terra. Organized by workflow name.

## Active Workflows

| Directory | WDL | Inputs | Description |
|---|---|---:|---|
| `HaplotypeCallerScatter/` | HC scatter | 20 | 5 probands (150x, 300 shards) + 10 parents + 3 GIAB + test |
| `Mutect2TumorOnlyScatter/` | Mutect2 scatter | 8 | 5 probands + GIAB + test/local |
| `MosaicPostProcessPart1/` | HC PostProcess | 8 | 5 trios + GIAB + test/local |
| `MutectPostProcessPart1/` | M2 PostProcess | 6 | 5 trios + GIAB |
| `DownsampleProband/` | Downsampling | 19 | 5 probands + 10 parents (coverage only) + 3 GIAB + test |
| `IGVTrio/` | IGV screenshots | 5 | 5 trios, mosaic + het de novo positions (HC ∪ M2) |
| `AddReadGroups/` | @RG fix | 3 | GIAB HG002/HG003/HG004 |
| `DownloadBAM/` | BAM download | 3 | GIAB trio from URL |

## Naming Convention

```
{WorkflowName}/{WorkflowName}.{sample_or_trio}[.terra].json
```

- `.terra.json` — for Terra execution (GCS paths)
- `.local.json` — for local Cromwell (local paths)
- `.json` — general / either

Superseded input manifests are not included. The JSON files present here are
the supported Terra/local configurations for the workflows listed above.

Terra manifests use the `${TERRA_BUCKET}` placeholder. Set `TERRA_BUCKET` to
the bucket name (without `gs://`) and render a manifest with `envsubst` before
uploading. Local manifests use `${LOCAL_TERRA_DIR}` and should be rendered in
the same way. IGV manifests also use `${IGV_DOCKER_IMAGE}` for the versioned
container built from the root Dockerfile. These placeholders are not expanded
by Terra or Cromwell.

The GIAB manifests reference public NCBI/Broad data and can be used as the
public test path. Study-cohort manifests use role-based placeholders such as
`BRCA1_PROBAND`, `BRCA1_MOTHER`, and `BRCA1_FATHER`; they do not contain the
original study sample labels. Their GCS paths are templates: authorized users
must replace them with objects in their own workspace. No credentials, tokens,
participant variants, or original Terra submission identifiers are stored in
these files.
