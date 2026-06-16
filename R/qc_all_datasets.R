run_all_dataset_qc <- function(bulk, group, igan_expr, subtype, cfg) {
  out <- cfg$output_dir
  qc_dir <- file.path(out, "supplementary", "qc")
  tab_dir <- file.path(out, "tables")
  rep_dir <- file.path(out, "reports")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

  notes <- character()
  add_note <- function(x) notes <<- c(notes, x)
  finite_vec <- function(x) x[is.finite(x) & !is.na(x)]
  safe_ggsave <- function(path, plot, width = 7, height = 5) {
    tryCatch({ ggplot2::ggsave(path, plot, width = width, height = height); TRUE }, error = function(e) { add_note(paste("Skipped figure", basename(path), "-", conditionMessage(e))); FALSE })
  }
  safe_csv <- function(x, path) tryCatch(utils::write.csv(x, path, row.names = FALSE), error = function(e) add_note(paste("Could not write", basename(path), conditionMessage(e))))

  expr <- bulk$expr
  pheno <- bulk$pheno
  group <- factor(as.character(group))
  names(group) <- colnames(expr)

  sample_summary <- data.frame(dataset = "GSE104948", sample = colnames(expr), inferred_group = as.character(group), stringsAsFactors = FALSE)
  subtype_df <- data.frame(sample = names(subtype), subtype = as.character(subtype), stringsAsFactors = FALSE)
  sample_summary <- merge(sample_summary, subtype_df, by = "sample", all.x = TRUE)
  sample_summary$subtype_note <- ifelse(sample_summary$subtype == "Subtype2", "Subtype2, n=4; exploratory small-cluster subtype", "")
  safe_csv(sample_summary, file.path(tab_dir, "qc_GSE104948_sample_summary.csv"))

  group_counts <- as.data.frame(table(sample_summary$inferred_group), stringsAsFactors = FALSE)
  colnames(group_counts) <- c("group", "n")
  safe_csv(group_counts, file.path(tab_dir, "qc_GSE104948_group_counts.csv"))

  subtype_counts <- as.data.frame(table(sample_summary$subtype), stringsAsFactors = FALSE)
  colnames(subtype_counts) <- c("subtype", "n")
  subtype_counts$note <- ifelse(subtype_counts$subtype == "Subtype2", "exploratory small-cluster subtype", "")
  safe_csv(subtype_counts, file.path(tab_dir, "qc_GSE104948_subtype_counts.csv"))

  expr_long <- tryCatch({
    set.seed(1)
    genes <- rownames(expr)
    if (length(genes) > 2000) genes <- sample(genes, 2000)
    reshape2::melt(expr[genes, , drop = FALSE], varnames = c("gene", "sample"), value.name = "expression")
  }, error = function(e) NULL)
  if (!is.null(expr_long) && nrow(expr_long) > 0) {
    expr_long$group <- sample_summary$inferred_group[match(expr_long$sample, sample_summary$sample)]
    p <- ggplot2::ggplot(expr_long, ggplot2::aes(x = sample, y = expression, fill = group)) + ggplot2::geom_boxplot(outlier.size = 0.1) + ggplot2::theme_bw(base_size = 8) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) + ggplot2::labs(title = "Supplementary QC: GSE104948 expression distribution", x = "Samples", y = "Expression")
    safe_ggsave(file.path(qc_dir, "Suppl_QC_GSE104948_expression_boxplot.pdf"), p, 8, 4)
    p2 <- ggplot2::ggplot(expr_long, ggplot2::aes(x = expression, group = sample, color = group)) + ggplot2::geom_density(alpha = 0.25) + ggplot2::theme_bw(base_size = 9) + ggplot2::labs(title = "Supplementary QC: GSE104948 expression density", x = "Expression", y = "Density")
    safe_ggsave(file.path(qc_dir, "Suppl_QC_GSE104948_expression_density.pdf"), p2, 7, 4)
  } else add_note("GSE104948 expression distribution plots skipped because expression matrix was empty or could not be melted.")

  make_pca <- function(mat, annot, label_col, file, title) {
    ok <- tryCatch({
      mat <- mat[apply(mat, 1, stats::var, na.rm = TRUE) > 0, , drop = FALSE]
      if (nrow(mat) < 3 || ncol(mat) < 3) stop("insufficient matrix dimensions")
      pc <- stats::prcomp(t(mat), scale. = TRUE)
      df <- data.frame(sample = rownames(pc$x), PC1 = pc$x[, 1], PC2 = pc$x[, 2], label = annot[match(rownames(pc$x), annot$sample), label_col], stringsAsFactors = FALSE)
      if (length(unique(stats::na.omit(df$label))) < 1) df$label <- "Unknown"
      p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = label)) + ggplot2::geom_point(size = 2) + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(title = title, color = label_col)
      safe_ggsave(file, p, 6, 5)
      TRUE
    }, error = function(e) { add_note(paste("Skipped PCA", basename(file), "-", conditionMessage(e))); FALSE })
    invisible(ok)
  }
  make_pca(expr, sample_summary, "inferred_group", file.path(qc_dir, "Suppl_QC_GSE104948_PCA_by_group.pdf"), "Supplementary QC: GSE104948 PCA by group")
  if (!is.null(igan_expr) && ncol(igan_expr) >= 3) make_pca(igan_expr, sample_summary, "subtype", file.path(qc_dir, "Suppl_QC_GSE104948_IgAN_only_PCA_by_subtype.pdf"), "Supplementary QC: IgAN-only PCA by subtype")

  if (nrow(subtype_counts) > 0) {
    p3 <- ggplot2::ggplot(subtype_counts, ggplot2::aes(x = subtype, y = n)) + ggplot2::geom_col() + ggplot2::theme_bw(base_size = 10) + ggplot2::labs(title = "Supplementary QC: IgAN subtype sizes", subtitle = "Subtype2 is exploratory because n=4", x = "Subtype", y = "Samples")
    safe_ggsave(file.path(qc_dir, "Suppl_QC_GSE104948_subtype_size_barplot.pdf"), p3, 5, 4)
  }

  qc_summary <- data.frame(dataset = "GSE104948", n_samples = ncol(expr), n_features = nrow(expr), n_igan_only = ifelse(is.null(igan_expr), NA_integer_, ncol(igan_expr)), n_subtype1 = sum(as.character(subtype) == "Subtype1", na.rm = TRUE), n_subtype2 = sum(as.character(subtype) == "Subtype2", na.rm = TRUE), qc_role = "Discovery bulk QC; supplementary only", stringsAsFactors = FALSE)

  ext_sum_file <- file.path(tab_dir, "external_validation_dataset_summary.csv")
  ext_score_file <- file.path(tab_dir, "external_validation_signature_scores.csv")
  if (file.exists(ext_sum_file)) {
    ext_sum <- utils::read.csv(ext_sum_file, stringsAsFactors = FALSE, check.names = FALSE)
    ext_sum$qc_role <- "External bulk validation / ML eligibility QC; supplementary only"
    safe_csv(ext_sum, file.path(tab_dir, "qc_external_bulk_dataset_summary.csv"))
    common_cols <- intersect(colnames(qc_summary), colnames(ext_sum))
    if (length(common_cols) > 0) add_note("External dataset-level QC table created from external_validation_dataset_summary.csv.")
  } else add_note("External validation dataset summary not found; external QC table will be limited.")

  if (file.exists(ext_score_file)) {
    sc <- utils::read.csv(ext_score_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (all(c("dataset", "sample", "small_cluster_signature_score") %in% colnames(sc))) {
      ig <- sc[sc$inferred_group == "IgAN" & is.finite(sc$small_cluster_signature_score), , drop = FALSE]
      if (nrow(ig) > 0) {
        ig$CERI_median_group <- ave(ig$small_cluster_signature_score, ig$dataset, FUN = function(x) ifelse(x >= stats::median(x, na.rm = TRUE), "CERI-high", "CERI-low"))
        safe_csv(ig[, c("dataset", "platform", "sample", "inferred_group", "small_cluster_signature_score", "CERI_median_group", "n_up_genes_matched", "n_down_genes_matched")], file.path(tab_dir, "qc_external_bulk_CERI_score_groups.csv"))
        p4 <- ggplot2::ggplot(ig, ggplot2::aes(x = CERI_median_group, y = small_cluster_signature_score)) + ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_jitter(width = 0.15, size = 1.2) + ggplot2::facet_wrap(~dataset, scales = "free") + ggplot2::theme_bw(base_size = 9) + ggplot2::labs(title = "Supplementary QC: external bulk CERI score strata", x = "Median-based stratum", y = "CERI signature score")
        safe_ggsave(file.path(qc_dir, "Suppl_QC_external_bulk_CERI_score_strata.pdf"), p4, 8, 5)
        count_df <- as.data.frame(table(ig$dataset, ig$CERI_median_group), stringsAsFactors = FALSE)
        colnames(count_df) <- c("dataset", "CERI_median_group", "n")
        safe_csv(count_df, file.path(tab_dir, "qc_external_bulk_CERI_high_low_counts.csv"))
      } else add_note("External CERI score QC skipped because no finite IgAN signature scores were available.")
    }
  } else add_note("External validation signature score table not found; CERI-high/low QC skipped.")

  notes_df <- if (length(notes) == 0) data.frame(note = "All available QC plots/tables generated without fatal errors.") else data.frame(note = notes)
  safe_csv(notes_df, file.path(tab_dir, "qc_all_datasets_notes.csv"))
  safe_csv(qc_summary, file.path(tab_dir, "qc_all_datasets_summary.csv"))

  report <- c(
    "# QC report for all manuscript datasets",
    "",
    "QC-heavy figures are intentionally saved under `output/supplementary/qc/` and should be treated as Supplementary Figures, not main-text figures.",
    "",
    paste0("- GSE104948 total samples: ", ncol(expr)),
    paste0("- GSE104948 retained features: ", nrow(expr)),
    paste0("- IgAN-only samples used for subtype discovery: ", ifelse(is.null(igan_expr), NA, ncol(igan_expr))),
    paste0("- Subtype1 samples: ", sum(as.character(subtype) == "Subtype1", na.rm = TRUE)),
    paste0("- Subtype2 samples: ", sum(as.character(subtype) == "Subtype2", na.rm = TRUE), " (exploratory small-cluster subtype)"),
    "",
    "External bulk QC is generated from the external validation score/summary tables produced earlier in the pipeline. Dataset-level expression PCA for external cohorts is not forced here to avoid duplicate downloads and fragile metadata-dependent failures; unavailable plots are recorded as notes instead of failing the workflow.",
    "",
    "## QC notes",
    paste0("- ", notes_df$note)
  )
  writeLines(report, file.path(rep_dir, "qc_all_datasets_report.md"))
  invisible(list(summary = qc_summary, notes = notes_df))
}
