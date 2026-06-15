#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1L, parallel::detectCores() - 1L))
cran <- c("yaml", "ggplot2", "pheatmap", "matrixStats", "reshape2", "glmnet", "Seurat", "tidyverse", "msigdbr", "data.table", "BiocManager")
bioc <- c("GEOquery", "Biobase", "limma", "sva", "ConsensusClusterPlus", "GSVA", "clusterProfiler", "org.Hs.eg.db")
install_missing <- function(pkgs, installer) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) installer(missing)
}
install_missing(cran, function(x) install.packages(x, dependencies = TRUE))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install_missing(bioc, function(x) BiocManager::install(x, ask = FALSE, update = FALSE))
all_pkgs <- c(cran, bioc)
missing <- all_pkgs[!vapply(all_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Packages still missing after installation: ", paste(missing, collapse = ", "))
writeLines(capture.output(sessionInfo()), "output/logs/r_session_info.txt")
