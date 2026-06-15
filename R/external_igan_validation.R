run_external_igan_validation <- function(cfg) {
  tables_dir <- file.path(cfg$output_dir, "tables")
  figures_dir <- file.path(cfg$output_dir, "figures")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  ext_cfg <- cfg$external_validation
  if (is.null(ext_cfg) || isFALSE(ext_cfg$enabled)) {
    message("External IgAN validation is disabled in config.")
    return(invisible(NULL))
  }
  datasets <- ext_cfg$datasets
  if (is.null(datasets) || length(datasets) == 0) {
    message("No external validation datasets configured.")
    utils::write.csv(data.frame(note = "No external validation datasets configured."), file.path(tables_dir, "external_validation_dataset_summary.csv"), row.names = FALSE)
    return(invisible(NULL))
  }

  marker_file <- file.path(tables_dir, "subtype_marker_direction_table.csv")
  fallback_file <- file.path(tables_dir, "DEG_subtype_comparison.csv")
  if (file.exists(marker_file)) {
    markers <- utils::read.csv(marker_file, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (file.exists(fallback_file)) {
    markers <- utils::read.csv(fallback_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
    markers$direction <- ifelse(markers$logFC > 0, "up_in_Subtype2_vs_Subtype1", "down_in_Subtype2_vs_Subtype1")
  } else {
    stop("Missing marker/DEG table for external validation.")
  }
  if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
  markers$gene <- unique(as.character(markers$gene))

  up_genes <- unique(markers$gene[grepl("up_in_Subtype2", markers$direction, ignore.case = TRUE)])
  down_genes <- unique(markers$gene[grepl("down_in_Subtype2", markers$direction, ignore.case = TRUE)])
  if (length(up_genes) < 3) up_genes <- head(unique(markers$gene), min(25, length(unique(markers$gene))))
  if (length(down_genes) < 3 && "logFC" %in% colnames(markers)) down_genes <- unique(markers$gene[markers$logFC < 0])

  map_external_symbols <- function(eset) {
    expr <- Biobase::exprs(eset)
    fdat <- Biobase::fData(eset)
    gene_col <- NA_character_
    if (ncol(fdat) > 0) {
      candidates <- grep("symbol|gene.symbol|gene_assignment|gene.?title|gene", colnames(fdat), ignore.case = TRUE, value = TRUE)
      if (length(candidates) > 0) gene_col <- candidates[1]
    }
    if (!is.na(gene_col)) {
      sym <- as.character(fdat[[gene_col]])
      # Clean common platform formats: multiple symbols separated by ///, //, ;, comma, or spaces.
      sym <- gsub("///.*$", "", sym)
      sym <- gsub("//.*$", "", sym)
      sym <- gsub(";.*$", "", sym)
      sym <- gsub(",.*$", "", sym)
      sym <- trimws(sym)
      good <- !is.na(sym) & sym != "" & !grepl("^[0-9]+$", sym)
      if (sum(good) >= 100) {
        expr <- expr[good, , drop = FALSE]
        rownames(expr) <- sym[good]
      }
    }
    # Collapse duplicate gene symbols by mean expression.
    if (any(duplicated(rownames(expr)))) {
      expr <- rowsum(expr, group = rownames(expr), reorder = FALSE) / as.numeric(table(rownames(expr))[rownames(rowsum(expr, group = rownames(expr), reorder = FALSE))])
    }
    expr
  }

  infer_igan_labels <- function(pheno) {
    combined <- apply(pheno, 1, function(x) paste(x, collapse = " | "))
    igan <- grepl("IgA nephropathy|IgAN|iga neph", combined, ignore.case = TRUE)
    non_igan <- grepl("control|healthy|normal|living donor|tumou?r nephrectomy", combined, ignore.case = TRUE)
    label <- ifelse(igan, "IgAN", ifelse(non_igan, "Control/Other", "Other/Unknown"))
    factor(label, levels = c("IgAN", "Control/Other", "Other/Unknown"))
  }

  zscore_rows <- function(mat) {
    m <- rowMeans(mat, na.rm = TRUE)
    s <- matrixStats::rowSds(mat, na.rm = TRUE)
    s[is.na(s) | s == 0] <- 1
    sweep(sweep(mat, 1, m, "-"), 1, s, "/")
  }

  all_scores <- list()
  summary_rows <- list()

  for (gse in datasets) {
    message("External validation: processing ", gse)
    res <- tryCatch({
      gsets <- GEOquery::getGEO(gse, GSEMatrix = TRUE)
      if (!is.list(gsets)) gsets <- list(gsets)
      dataset_scores <- list()
      dataset_summaries <- list()
      for (i in seq_along(gsets)) {
        eset <- gsets[[i]]
        platform <- Biobase::annotation(eset)
        expr <- map_external_symbols(eset)
        pheno <- Biobase::pData(eset)
        labels <- infer_igan_labels(pheno)
        common_up <- intersect(up_genes, rownames(expr))
        common_down <- intersect(down_genes, rownames(expr))
        common_all <- unique(c(common_up, common_down))
        status <- "usable"
        if (length(common_up) < 3 && length(common_down) < 3) status <- "insufficient_signature_overlap"
        if (sum(labels == "IgAN", na.rm = TRUE) < ext_cfg$min_igan_samples) status <- "too_few_igan_samples"

        if (status == "usable") {
          z <- zscore_rows(expr[common_all, , drop = FALSE])
          up_score <- if (length(common_up) > 0) colMeans(z[common_up, , drop = FALSE], na.rm = TRUE) else rep(NA_real_, ncol(z))
          down_score <- if (length(common_down) > 0) colMeans(z[common_down, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(z))
          small_cluster_score <- up_score - down_score
          score_df <- data.frame(
            dataset = gse,
            platform = platform,
            sample = colnames(expr),
            inferred_group = as.character(labels),
            small_cluster_signature_score = as.numeric(small_cluster_score),
            up_signature_score = as.numeric(up_score),
            down_signature_score = as.numeric(down_score),
            n_up_genes_matched = length(common_up),
            n_down_genes_matched = length(common_down),
            stringsAsFactors = FALSE
          )
          dataset_scores[[length(dataset_scores) + 1]] <- score_df
        }
        dataset_summaries[[length(dataset_summaries) + 1]] <- data.frame(
          dataset = gse,
          platform = platform,
          n_samples = ncol(expr),
          n_features = nrow(expr),
          n_igan = sum(labels == "IgAN", na.rm = TRUE),
          n_control_or_other = sum(labels == "Control/Other", na.rm = TRUE),
          n_unknown = sum(labels == "Other/Unknown", na.rm = TRUE),
          n_up_genes_matched = length(common_up),
          n_down_genes_matched = length(common_down),
          status = status,
          stringsAsFactors = FALSE
        )
      }
      list(scores = dataset_scores, summaries = dataset_summaries)
    }, error = function(e) {
      message("External validation failed for ", gse, ": ", conditionMessage(e))
      list(scores = list(), summaries = list(data.frame(
        dataset = gse,
        platform = NA_character_,
        n_samples = NA_integer_,
        n_features = NA_integer_,
        n_igan = NA_integer_,
        n_control_or_other = NA_integer_,
        n_unknown = NA_integer_,
        n_up_genes_matched = NA_integer_,
        n_down_genes_matched = NA_integer_,
        status = paste0("download_or_parse_failed: ", conditionMessage(e)),
        stringsAsFactors = FALSE
      )))
    })
    all_scores <- c(all_scores, res$scores)
    summary_rows <- c(summary_rows, res$summaries)
  }

  summary_df <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame(note = "No external datasets processed.")
  utils::write.csv(summary_df, file.path(tables_dir, "external_validation_dataset_summary.csv"), row.names = FALSE)

  if (length(all_scores) > 0) {
    score_df <- do.call(rbind, all_scores)
    utils::write.csv(score_df, file.path(tables_dir, "external_validation_signature_scores.csv"), row.names = FALSE)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p <- ggplot2::ggplot(score_df, ggplot2::aes(x = inferred_group, y = small_cluster_signature_score)) +
        ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.6) +
        ggplot2::geom_jitter(width = 0.18, size = 1.5, alpha = 0.75) +
        ggplot2::facet_wrap(~dataset + platform, scales = "free_x") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
        ggplot2::labs(title = "External IgAN small-cluster signature validation", x = "Inferred group", y = "Small-cluster signature score")
      ggplot2::ggsave(file.path(figures_dir, "Figure6_external_signature_validation.pdf"), p, width = 9, height = 5.5)
    }
  } else {
    utils::write.csv(data.frame(note = "No usable external validation scores. Check external_validation_dataset_summary.csv."), file.path(tables_dir, "external_validation_signature_scores.csv"), row.names = FALSE)
  }

  invisible(list(summary = summary_df))
}
