#!/usr/bin/env Rscript

# Lightweight dependency check for GitHub Actions.
# Core Bioconductor packages should be installed by apt in the workflow,
# so this script must NOT call BiocManager online.

dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)

options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"),
  timeout = 1200
)

cran_core <- c(
  "yaml",
  "ggplot2",
  "pheatmap",
  "matrixStats",
  "reshape2",
  "glmnet",
  "data.table"
)

bioc_core <- c(
  "GEOquery",
  "Biobase",
  "limma",
  "ConsensusClusterPlus"
)

optional <- c(
  "sva",
  "GSVA",
  "clusterProfiler",
  "org.Hs.eg.db",
  "msigdbr"
)

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
    try(install.packages(missing, dependencies = TRUE), silent = FALSE)
  }
}

install_cran_if_missing(cran_core)

required <- c(cran_core, bioc_core)
missing_required <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_required)) {
  stop(
    "Required packages are still missing: ",
    paste(missing_required, collapse = ", "),
    "\nThese should be installed by apt-get in .github/workflows/run-iga-pipeline.yml."
  )
}

missing_optional <- optional[!vapply(optional, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_optional)) {
  message("Optional packages missing and will be skipped: ", paste(missing_optional, collapse = ", "))
}

writeLines(capture.output(sessionInfo()), "output/logs/r_session_info.txt")
