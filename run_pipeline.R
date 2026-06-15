#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(yaml)
})
source("R/utils.R")
source("R/bulk_analysis.R")
source("R/functional_modeling.R")
source("R/single_cell_analysis.R")
source("R/subtype_interpretation.R")

cfg <- yaml::read_yaml("config/pipeline.yml")
set.seed(cfg$seed)
make_dirs(cfg)
dir.create(file.path(cfg$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(cfg$output_dir, "logs", "pipeline.log")
log_con <- file(log_file, open = "at")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

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

igan_samples <- intersect(names(group)[as.character(group) == "IgAN"], colnames(bulk$expr))
if (length(igan_samples) < 10) {
  stop("Too few IgAN samples detected for subtype analysis: ", length(igan_samples))
}
igan_expr <- bulk$expr[, igan_samples, drop = FALSE]
save_analysis_sample_report(bulk$expr, group, igan_expr, cfg)
message("Using IgAN-only subset for molecular subtyping: n=", ncol(igan_expr))

subtype <- log_step("IgAN-only consensus clustering k=2..4", run_consensus(igan_expr, cfg))
deg <- log_step("IgAN subtype limma differential expression", run_limma(igan_expr, subtype, cfg))
markers <- log_step("IgAN subtype marker export", get_subtype_markers(deg, cfg))
log_step("IgAN subtype interpretation and QC", run_subtype_interpretation_qc(cfg), required = FALSE)
log_step("IgAN Hallmark GSVA and immune signatures", run_functional(igan_expr, subtype, cfg), required = FALSE)
log_step("IgAN LASSO subtype model", run_lasso(igan_expr, subtype, markers, cfg), required = FALSE)
log_step("Single-cell validation", run_single_cell(markers, cfg), required = FALSE)

message("Pipeline complete. Outputs written to ", normalizePath(cfg$output_dir))
