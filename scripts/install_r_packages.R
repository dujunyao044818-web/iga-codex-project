#!/usr/bin/env Rscript
# Robust lightweight installer for GitHub Actions.
# Heavy single-cell / enrichment packages are optional and must not stop the workflow.
dir.create("output/logs", recursive = TRUE, showWarnings = FALSE)
options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/latest"),
  Ncpus = max(1L, parallel::detectCores() - 1L),
  timeout = 1200
)
Sys.setenv(
  BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "FALSE",
  R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true"
)

install_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing CRAN packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
}

install_bioc <- function(pkgs, required = TRUE) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))
  message("Installing Bioconductor packages: ", paste(missing, collapse = ", "))
  ok <- tryCatch({
    BiocManager::install(missing, ask = FALSE, update = FALSE, force = FALSE)
    TRUE
  }, error = function(e) {
    message("Bioconductor install warning/error: ", conditionMessage(e))
    FALSE
  })
  still_missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (required && length(still_missing)) {
    stop("Required packages still missing after installation: ", paste(still_missing, collapse = ", "))
  }
  if (!required && length(still_missing)) {
    message("Optional packages still missing and will be skipped: ", paste(still_missing, collapse = ", "))
  }
  invisible(ok)
}

required_cran <- c("yaml", "ggplot2", "pheatmap", "matrixStats", "reshape2", "glmnet", "data.table", "BiocManager")
required_bioc <- c("GEOquery", "Biobase", "limma", "ConsensusClusterPlus")
optional_bioc <- c("sva", "GSVA", "clusterProfiler", "org.Hs.eg.db")
optional_cran <- c("msigdbr")
# Seurat/tidyverse are intentionally not required for this bulk-first run.

install_cran(required_cran)
install_bioc(required_bioc, required = TRUE)
try(install_cran(optional_cran), silent = TRUE)
try(install_bioc(optional_bioc, required = FALSE), silent = TRUE)

writeLines(capture.output(sessionInfo()), "output/logs/r_session_info.txt")
