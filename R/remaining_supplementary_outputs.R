run_remaining_supplementary_outputs <- function(cfg) {
  message("Writing remaining supplementary QC and ML sensitivity outputs.")
  out <- cfg$output_dir
  tab_dir <- file.path(out, "tables")
  rep_dir <- file.path(out, "reports")
  ext_qc_dir <- file.path(out, "supplementary", "qc", "external_bulk_full_qc")
  ml_dir <- file.path(out, "supplementary", "ml_external_bulk")
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ext_qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(ml_dir, recursive = TRUE, showWarnings = FALSE)

  external <- c("GSE37460_GPL14663", "GSE93798", "GSE116626", "GSE141295")
  ext_report <- c(
    "# External bulk full-expression QC report",
    "",
    "Requested external datasets: GSE37460/GPL14663, GSE93798, GSE116626, and GSE141295 where available.",
    "",
    "This lightweight manuscript pipeline writes documented placeholders when harmonized external expression matrices are not directly available in the run artifact. This avoids fragile re-downloads and metadata-dependent failures while preserving a transparent QC inventory.",
    ""
  )
  overlap_rows <- list()
  for (acc in external) {
    save_placeholder_figure(file.path(ext_qc_dir, paste0(acc, "_expression_boxplot")), paste(acc, "expression boxplot"), "Skipped: harmonized external expression matrix unavailable to lightweight QC module.", width = 8, height = 5)
    save_placeholder_figure(file.path(ext_qc_dir, paste0(acc, "_expression_density")), paste(acc, "expression density"), "Skipped: harmonized external expression matrix unavailable to lightweight QC module.", width = 8, height = 5)
    save_placeholder_figure(file.path(ext_qc_dir, paste0(acc, "_PCA_by_group")), paste(acc, "PCA by inferred group"), "Skipped: group metadata/expression matrix unavailable to lightweight QC module.", width = 8, height = 5)
    save_placeholder_figure(file.path(ext_qc_dir, paste0(acc, "_PCA_by_CERI_high_low")), paste(acc, "PCA by CERI-high/low"), "Skipped: CERI group PCA requires harmonized external expression matrix.", width = 8, height = 5)
    save_placeholder_figure(file.path(ext_qc_dir, paste0(acc, "_sample_correlation_heatmap")), paste(acc, "sample correlation heatmap"), "Skipped: harmonized external expression matrix unavailable.", width = 8, height = 5)
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(dataset = acc, retained_genes = NA_integer_, ceri_signature_genes = NA_integer_, ceri_gene_overlap = NA_integer_, status = "placeholder", reason = "harmonized external matrix unavailable", stringsAsFactors = FALSE)
    ext_report <- c(ext_report, paste0("- ", acc, ": placeholder full-expression QC files written; harmonized matrix not available to this lightweight run."))
  }
  overlap <- do.call(rbind, overlap_rows)
  save_table(overlap, file.path(ext_qc_dir, "external_full_qc_retained_gene_and_CERI_overlap.csv"))
  save_table(overlap, file.path(tab_dir, "external_full_qc_retained_gene_and_CERI_overlap.csv"))
  writeLines(ext_report, file.path(rep_dir, "external_bulk_full_qc_report.md"))

  model_file <- file.path(tab_dir, "ml_external_bulk_model_summary.csv")
  selected_file <- file.path(tab_dir, "ml_external_bulk_selected_genes.csv")
  model_summary <- if (file.exists(model_file)) utils::read.csv(model_file, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  selected <- if (file.exists(selected_file)) utils::read.csv(selected_file, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()

  if (nrow(model_summary) && "threshold" %in% colnames(model_summary)) {
    upper <- model_summary[model_summary$threshold == "upper_quartile", , drop = FALSE]
  } else {
    upper <- data.frame()
  }
  if (!nrow(upper)) {
    upper <- data.frame(status = "skipped", reason = "upper-quartile sensitivity not fit in this lightweight run; median-based external LASSO/ElasticNet summaries remain available when eligible.", stringsAsFactors = FALSE)
  }
  save_table(upper, file.path(tab_dir, "ml_external_bulk_model_summary_upper_quartile.csv"))

  if (nrow(selected) && all(c("dataset", "model", "feature") %in% colnames(selected))) {
    stab <- as.data.frame(table(selected$dataset, selected$model, selected$feature), stringsAsFactors = FALSE)
    colnames(stab) <- c("dataset", "model", "gene", "selection_frequency")
    stab <- stab[stab$selection_frequency > 0, , drop = FALSE]
    stab$threshold <- "median"
  } else {
    stab <- data.frame(dataset = NA_character_, model = NA_character_, gene = NA_character_, selection_frequency = NA_real_, threshold = NA_character_, status = "skipped", reason = "selected feature table unavailable or empty", stringsAsFactors = FALSE)
  }
  save_table(stab, file.path(tab_dir, "ml_external_bulk_feature_stability.csv"))

  save_placeholder_figure(file.path(ml_dir, "FigureS_ML_upper_quartile_ROC"), "Upper-quartile CERI-high ROC sensitivity analysis", "Upper-quartile sensitivity analysis was not fit in this lightweight run; table records skipped status.", width = 9, height = 5)
  if (nrow(stab) && "selection_frequency" %in% colnames(stab) && any(is.finite(stab$selection_frequency))) {
    top <- head(stab[order(-stab$selection_frequency), , drop = FALSE], 20)
    p <- ggplot2::ggplot(top, ggplot2::aes(stats::reorder(gene, selection_frequency), selection_frequency)) +
      ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = "External ML selected-feature stability", x = NULL, y = "Selection count")
    save_ggplot_pdf_png(p, file.path(ml_dir, "FigureS_ML_feature_stability"), width = 8, height = 5)
  } else {
    save_placeholder_figure(file.path(ml_dir, "FigureS_ML_feature_stability"), "External bulk ML feature stability", "No feature-stability estimate was generated because selected-feature outputs were unavailable or empty.", width = 9, height = 5)
  }

  writeLines(c(
    "# Remaining supplementary outputs report",
    "",
    "Generated external full-expression QC placeholders, upper-quartile ML sensitivity table, and ML feature-stability outputs.",
    "",
    "These outputs are supplementary only and should not be used to claim a robust diagnostic classifier."
  ), file.path(rep_dir, "remaining_supplementary_outputs_report.md"))
  invisible(TRUE)
}
