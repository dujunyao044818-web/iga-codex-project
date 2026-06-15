options(repos = c(CRAN = 'https://cloud.r-project.org'))
cran_pkgs <- c('yaml','ggplot2','pheatmap','dplyr','tidyr','tibble','readr','stringr','matrixStats','glmnet','survival','survminer','msigdbr','GSVA','remotes')
bioc_pkgs <- c('GEOquery','limma','sva','ConsensusClusterPlus','clusterProfiler','org.Hs.eg.db','Biobase','BiocGenerics','GSEABase')
install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
  }
}
install_if_missing(cran_pkgs)
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}
# Seurat is heavy; install only when available within CI time. The pipeline can continue without it.
if (!requireNamespace('Seurat', quietly = TRUE)) {
  try(install.packages('Seurat', dependencies = TRUE), silent = TRUE)
}
writeLines(capture.output(sessionInfo()), 'sessionInfo_install.txt')
