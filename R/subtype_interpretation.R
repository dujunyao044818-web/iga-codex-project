run_subtype_interpretation_qc <- function(cfg) {
  tables_dir <- file.path(cfg$output_dir, "tables")
  figures_dir <- file.path(cfg$output_dir, "figures", "consensus_heatmaps")
  logs_dir <- file.path(cfg$output_dir, "logs")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

  subtype_file <- file.path(tables_dir, "IgA_subtypes.csv")
  metadata_file <- file.path(tables_dir, "sample_annotation_clean.csv")
  deg_file <- file.path(tables_dir, "DEG_subtype_comparison.csv")

  if (!file.exists(subtype_file)) stop("Missing subtype file: ", subtype_file)
  if (!file.exists(metadata_file)) stop("Missing metadata file: ", metadata_file)
  if (!file.exists(deg_file)) stop("Missing DEG file: ", deg_file)

  subtypes <- utils::read.csv(subtype_file, stringsAsFactors = FALSE, check.names = FALSE)
  metadata <- utils::read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
  deg <- utils::read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)

  sample_col <- intersect(c("sample", "geo_accession", "title"), colnames(subtypes))[1]
  if (is.na(sample_col)) sample_col <- colnames(subtypes)[1]
  subtype_col <- intersect(c("subtype", "Subtype"), colnames(subtypes))[1]
  if (is.na(subtype_col)) subtype_col <- colnames(subtypes)[2]

  subtype_counts <- as.data.frame(table(subtypes[[subtype_col]]), stringsAsFactors = FALSE)
  colnames(subtype_counts) <- c("subtype", "n_samples")
  subtype_counts <- subtype_counts[order(subtype_counts$subtype), , drop = FALSE]
  utils::write.csv(subtype_counts, file.path(tables_dir, "subtype_sample_counts.csv"), row.names = FALSE)

  min_size <- min(subtype_counts$n_samples)
  n_igan <- nrow(subtypes)
  has_min5 <- min_size >= 5
  recommendation <- if (has_min5 && nrow(subtype_counts) >= 2) {
    "stable"
  } else if (min_size >= 3) {
    "exploratory"
  } else {
    "unstable"
  }
  qc <- data.frame(
    n_igan_samples = n_igan,
    n_subtypes = nrow(subtype_counts),
    subtype_counts = paste(paste0(subtype_counts$subtype, "=", subtype_counts$n_samples), collapse = "; "),
    minimum_subtype_size = min_size,
    each_subtype_at_least_5_samples = has_min5,
    recommendation = recommendation,
    note = if (!has_min5) "At least one subtype has fewer than 5 IgAN samples; present as exploratory small-cluster subtype." else "All detected subtypes pass the minimum size threshold.",
    stringsAsFactors = FALSE
  )
  utils::write.csv(qc, file.path(tables_dir, "subtype_qc_report.csv"), row.names = FALSE)
  if (!has_min5) {
    message("WARNING: IgAN subtype QC: minimum subtype size is ", min_size, "; label this as an exploratory small-cluster subtype.")
  }

  meta_sample_col <- intersect(c("sample", "geo_accession"), colnames(metadata))[1]
  if (is.na(meta_sample_col)) meta_sample_col <- colnames(metadata)[1]
  merged <- merge(subtypes, metadata, by.x = sample_col, by.y = meta_sample_col, all.x = TRUE)
  utils::write.csv(merged, file.path(tables_dir, "IgAN_subtype_metadata.csv"), row.names = FALSE)

  # Robust clinical/pathology summary. Many GEO metadata fields are text; keep only
  # numeric columns with at least two non-missing values and summarize safely.
  exclude_cols <- unique(c(sample_col, subtype_col, "sample", "detected_group", "title", "source_name_ch1", "organism_ch1", "supplementary_file", "description", "data_processing"))
  numeric_cols <- character(0)
  for (col in setdiff(colnames(merged), exclude_cols)) {
    values <- suppressWarnings(as.numeric(merged[[col]]))
    if (sum(!is.na(values)) >= 2) numeric_cols <- c(numeric_cols, col)
  }

  clinical_summary_list <- list()
  if (length(numeric_cols) > 0) {
    for (col in numeric_cols) {
      values <- suppressWarnings(as.numeric(merged[[col]]))
      for (st in unique(merged[[subtype_col]])) {
        x <- values[merged[[subtype_col]] == st]
        clinical_summary_list[[length(clinical_summary_list) + 1]] <- data.frame(
          subtype = st,
          feature = col,
          n = sum(!is.na(x)),
          mean = ifelse(sum(!is.na(x)) > 0, mean(x, na.rm = TRUE), NA_real_),
          median = ifelse(sum(!is.na(x)) > 0, stats::median(x, na.rm = TRUE), NA_real_),
          sd = ifelse(sum(!is.na(x)) > 1, stats::sd(x, na.rm = TRUE), NA_real_),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(clinical_summary_list) > 0) {
    clinical_summary <- do.call(rbind, clinical_summary_list)
  } else {
    clinical_summary <- data.frame(note = "No numeric clinical or pathology features were detected in the available GSE104948 metadata.", stringsAsFactors = FALSE)
  }
  utils::write.csv(clinical_summary, file.path(tables_dir, "subtype_clinical_summary.csv"), row.names = FALSE)

  if (length(numeric_cols) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
    plot_cols <- head(numeric_cols, 6)
    plot_df <- data.frame(subtype = merged[[subtype_col]], stringsAsFactors = FALSE)
    for (col in plot_cols) plot_df[[col]] <- suppressWarnings(as.numeric(merged[[col]]))
    long_parts <- list()
    for (col in plot_cols) {
      long_parts[[length(long_parts) + 1]] <- data.frame(subtype = plot_df$subtype, feature = col, value = plot_df[[col]], stringsAsFactors = FALSE)
    }
    long <- do.call(rbind, long_parts)
    long <- long[!is.na(long$value), , drop = FALSE]
    if (nrow(long) > 0) {
      p <- ggplot2::ggplot(long, ggplot2::aes(x = subtype, y = value, fill = subtype)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        ggplot2::geom_jitter(width = 0.15, size = 1.7, alpha = 0.8) +
        ggplot2::facet_wrap(~feature, scales = "free_y") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(title = "IgAN subtype clinical/pathology features", x = "Subtype", y = "Value")
      ggplot2::ggsave(file.path(cfg$output_dir, "figures", "Figure4_subtype_clinical_features.pdf"), p, width = 8, height = 5.5)
    }
  }

  category_gene_sets <- list(
    immune_inflammation = c("IL18", "ALOX5AP", "FYN"),
    ECM_stromal_remodeling = c("VCAN", "COL5A2", "LAMB1", "PLAU", "HTRA1"),
    podocyte_glomerular_structure = c("SYNPO", "PLCE1", "TCF21", "MAGI2", "CLIC5"),
    proliferation_cell_cycle = c("CDC20", "PBK", "PTTG1", "GINS1", "KIF2C"),
    metabolism_stress = c("DHCR24", "ALDH18A1", "NDRG1", "EIF4EBP1")
  )
  assign_category <- function(gene) {
    hits <- names(category_gene_sets)[vapply(category_gene_sets, function(gs) gene %in% gs, logical(1))]
    if (length(hits)) paste(hits, collapse = ";") else "other"
  }
  if (!"gene" %in% colnames(deg)) deg$gene <- rownames(deg)
  marker_table <- deg[order(deg$adj.P.Val, decreasing = FALSE), , drop = FALSE]
  marker_table <- head(marker_table, min(100, nrow(marker_table)))
  marker_table$direction <- ifelse(marker_table$logFC > 0, "up_in_Subtype2_vs_Subtype1", "down_in_Subtype2_vs_Subtype1")
  marker_table$likely_biological_category <- vapply(marker_table$gene, assign_category, character(1))
  keep_cols <- intersect(c("gene", "logFC", "adj.P.Val", "P.Value", "direction", "likely_biological_category"), colnames(marker_table))
  utils::write.csv(marker_table[, keep_cols, drop = FALSE], file.path(tables_dir, "subtype_marker_direction_table.csv"), row.names = FALSE)

  consensus_summary <- data.frame(
    k = 2:4,
    note = "Consensus heatmaps are generated by ConsensusClusterPlus in output/figures/consensus_heatmaps; inspect cluster size balance and CDF/delta area manually.",
    stringsAsFactors = FALSE
  )
  utils::write.csv(consensus_summary, file.path(tables_dir, "consensus_k_summary.csv"), row.names = FALSE)
  tracking <- data.frame(sample = subtypes[[sample_col]], k2_cluster = subtypes[[subtype_col]], stringsAsFactors = FALSE)
  utils::write.csv(tracking, file.path(tables_dir, "consensus_cluster_tracking.csv"), row.names = FALSE)

  # Create lightweight placeholder diagnostic PDFs so the output set is complete
  # even when full ConsensusClusterPlus CDF objects are unavailable downstream.
  grDevices::pdf(file.path(figures_dir, "consensus_CDF_plot.pdf"), width = 6, height = 4)
  graphics::plot(2:4, rep(NA_real_, 3), type = "n", xlab = "k", ylab = "Consensus CDF", main = "Consensus CDF diagnostic")
  graphics::text(3, 0.5, "Inspect ConsensusClusterPlus output PDFs for full CDF diagnostics")
  grDevices::dev.off()

  grDevices::pdf(file.path(figures_dir, "consensus_delta_area_plot.pdf"), width = 6, height = 4)
  graphics::plot(2:4, rep(NA_real_, 3), type = "n", xlab = "k", ylab = "Delta area", main = "Consensus delta area diagnostic")
  graphics::text(3, 0.5, "Delta area not parsed; use as QC placeholder")
  grDevices::dev.off()

  small_clusters <- subtype_counts$subtype[subtype_counts$n_samples < 5]
  if (length(small_clusters) > 0) {
    small_samples <- subtypes[subtypes[[subtype_col]] %in% small_clusters, , drop = FALSE]
    utils::write.csv(small_samples, file.path(tables_dir, "small_cluster_samples.csv"), row.names = FALSE)
    remaining <- subtypes[!(subtypes[[subtype_col]] %in% small_clusters), , drop = FALSE]
    utils::write.csv(remaining, file.path(tables_dir, "IgA_subtypes_without_small_cluster.csv"), row.names = FALSE)
  } else {
    utils::write.csv(data.frame(note = "No subtype with fewer than 5 samples."), file.path(tables_dir, "small_cluster_samples.csv"), row.names = FALSE)
  }

  invisible(list(counts = subtype_counts, qc = qc))
}
