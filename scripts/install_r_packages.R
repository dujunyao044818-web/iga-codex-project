#!/usr/bin/env Rscript

# Lightweight dependency check for GitHub Actions.
# Core Bioconductor packages are expected from apt where possible.
# Optional enrichment dependencies are installed best-effort and never block the main pipeline.

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

optional_cran_enrichment <- c(
  "msigdbr"
)

optional_bioc_enrichment <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "enrichplot",
  "DOSE",
  "AnnotationDbi"
)

optional_other <- c(
  "sva",
  "GSVA"
)

install_cran_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
    try(install.packages(missing, dependencies = TRUE), silent = FALSE)
  }
}

install_optional_missing <- function(pkgs, installer, source) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(TRUE))
  message("Installing optional ", source, " packages: ", paste(missing, collapse = ", "))
  tryCatch(
    installer(missing),
    error = function(e) warning(
      "Optional ", source, " package installation failed for ",
      paste(missing, collapse = ", "), ": ", conditionMessage(e),
      call. = FALSE
    )
  )
  invisible(TRUE)
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

# Optional enrichment dependencies: best-effort only.
install_optional_missing(optional_cran_enrichment, function(x) install.packages(x, dependencies = TRUE), "CRAN enrichment")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  tryCatch(
    install.packages("BiocManager", dependencies = TRUE),
    error = function(e) warning("BiocManager installation failed: ", conditionMessage(e), call. = FALSE)
  )
}

if (requireNamespace("BiocManager", quietly = TRUE)) {
  install_optional_missing(
    optional_bioc_enrichment,
    function(x) BiocManager::install(x, ask = FALSE, update = FALSE),
    "Bioconductor enrichment"
  )
} else {
  warning("Optional Bioconductor enrichment packages were skipped because BiocManager is unavailable.", call. = FALSE)
}

optional_all <- c(optional_other, optional_cran_enrichment, optional_bioc_enrichment)
optional_status <- vapply(optional_all, requireNamespace, logical(1), quietly = TRUE)
message("Optional package availability:")
for (pkg in optional_all) {
  message("  - ", pkg, ": ", ifelse(optional_status[[pkg]], "available", "unavailable"))
}

writeLines(capture.output(sessionInfo()), "output/logs/r_session_info.txt")
