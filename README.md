# IgA nephropathy molecular subtyping pipeline

This repository contains a publication-oriented, end-to-end R pipeline to identify molecular subtypes in the single allowed bulk dataset, **GSE104948**, and validate subtype signatures in one public human kidney biopsy single-cell dataset, **GSE171314**.

## Outputs

All generated assets are written under:

- `output/figures`
- `output/tables`
- `output/models`

## Quick start

```bash
Rscript run_pipeline.R
```

The pipeline automatically creates missing directories, downloads GEO data, performs QC, subtype discovery, differential expression, GSVA/ssGSEA immune interpretation, LASSO modeling, Seurat analysis, and single-cell subtype signature scoring.


## GitHub Actions execution

Because some Codex containers block apt access, the full pipeline is configured to run in GitHub Actions via `.github/workflows/run-iga-pipeline.yml`. Trigger **Run IgA nephropathy pipeline** manually from the Actions tab; the workflow installs R and all CRAN/Bioconductor dependencies, runs `Rscript run_pipeline.R`, creates `output.zip`, and uploads `iga_pipeline_output`.

## Reproducibility

- Configuration is centralized in `config/pipeline.yml`.
- Random seeds are set from the configuration.
- `renv.lock` can be created by running `Rscript scripts/init_renv.R` after package installation if strict environment pinning is required by a journal or HPC environment.

## Dataset design rules

- Bulk RNA-seq/transcriptome dataset: `GSE104948` only.
- Single-cell validation dataset: `GSE171314`, human kidney biopsy scRNA-seq from IgA nephropathy and control subjects.
