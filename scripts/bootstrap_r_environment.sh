#!/usr/bin/env bash
set -euo pipefail

LOG_DIR="output/logs"
mkdir -p "$LOG_DIR"
exec > >(tee "$LOG_DIR/bootstrap_r_environment.log") 2>&1

export DEBIAN_FRONTEND=noninteractive

echo "[bootstrap] Checking for Rscript"
if ! command -v Rscript >/dev/null 2>&1; then
  echo "[bootstrap] Rscript not found; installing R and binary packages with apt-get"
  apt-get update
  apt-get install -y \
    r-base r-base-dev \
    r-cran-yaml r-cran-ggplot2 r-cran-pheatmap r-cran-matrixstats \
    r-cran-reshape2 r-cran-glmnet r-cran-seurat r-cran-tidyverse \
    r-bioc-geoquery r-bioc-biobase r-bioc-limma r-bioc-sva \
    r-bioc-consensusclusterplus r-bioc-gsva r-cran-msigdbr zip
else
  echo "[bootstrap] Found $(Rscript --version 2>&1)"
fi

Rscript scripts/install_r_packages.R
