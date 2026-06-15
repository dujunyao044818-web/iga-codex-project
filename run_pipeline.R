#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(yaml)
})
source("R/utils.R")
source("R/bulk_analysis.R")
source("R/functional_modeling.R")
source("R/single_cell_analysis.R")

cfg <- yaml::read_yaml("config/pipeline.yml")
set.seed(cfg$seed)
make_dirs(cfg)
dir.create(file.path(cfg$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(cfg$output_dir, "logs", "pipeline.log")
sink(log_file, split = TRUE)
sink(log_file, type = "message", append = TRUE)
on.exit({ sink(type = "message"); sink() }, add = TRUE)

log_step <- function(label, expr, required = TRUE) {
  message("[", Sys.time(), "] START: ", label)
  out <- tryCatch(force(expr), error = function(e) {
    message("[", Sys.time(), "] ERROR in ", label, ": ", conditionMessage(e))
    if (required) stop(e) else NULL
  })
  message("[", Sys.time(), "] END: ", label)
  out
}

ensure_packages(c("GEOquery", "Biobase", "limma", "matrixStats", "ConsensusClusterPlus", "pheatmap", "reshape2", "glmnet", "yaml", "data.table"))

bulk <- log_step("Download and preprocess GSE104948", load_bulk_gse(cfg))
group <- log_step("Detect sample groups from phenotype metadata", infer_bulk_groups(bulk$pheno, cfg), required = FALSE)
if (is.null(group)) group <- factor(rep("Unknown", ncol(bulk$expr)))
log_step("Annotated PCA and sample annotation export", bulk_qc(bulk$expr, bulk$pheno, group, cfg))
subtype <- log_step("Consensus clustering k=2..4", run_consensus(bulk$expr, cfg))
deg <- log_step("limma subtype differential expression", run_limma(bulk$expr, subtype, cfg))
markers <- log_step("Subtype marker export", get_subtype_markers(deg, cfg))
log_step("Hallmark GSVA and immune signatures", run_functional(bulk$expr, subtype, cfg), required = FALSE)
log_step("LASSO subtype model", run_lasso(bulk$expr, subtype, markers, cfg), required = FALSE)
log_step("Single-cell validation", run_single_cell(markers, cfg), required = FALSE)

message("Pipeline complete. Outputs written to ", normalizePath(cfg$output_dir))
