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

  curated_gene_sets <- list(
    T_cell = c("CD3D", "CD3E", "CD2", "TRAC", "IL7R", "LCK", "CD247"),
    CD8_T_cell = c("CD8A", "CD8B", "GZMB", "GZMK", "PRF1", "NKG7"),
    B_cell = c("MS4A1", "CD79A", "CD79B", "CD74", "BANK1", "CD19"),
    Plasma_cell = c("MZB1", "XBP1", "JCHAIN", "SDC1", "IGHG1", "IGKC"),
    Macrophage_monocyte = c("CD68", "LST1", "AIF1", "C1QA", "C1QB", "C1QC", "TYROBP", "FCGR3A"),
    Neutrophil = c("S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "MPO"),
    Complement = c("C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFB", "CFD"),
    Interferon_inflammation = c("ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT3", "MX1", "OAS1", "STAT1", "CXCL10"),
    ECM_stromal_remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A2", "FN1", "VCAN", "LAMB1", "MMP2", "MMP9", "TIMP1", "PLAU", "HTRA1"),
    Podocyte_glomerular_structure = c("NPHS1", "NPHS2", "SYNPO", "PODXL", "WT1", "PLCE1", "TCF21", "MAGI2", "CLIC5"),
    Endothelial = c("PECAM1", "VWF", "KDR", "FLT1", "EMCN", "CDH5", "SEMA3G"),
    Tubular_injury_stress = c("HAVCR1", "LCN2", "KRT8", "KRT18", "VCAM1", "NDRG1", "HMOX1", "ATF3"),
    Cell_cycle_proliferation = c("MKI67", "TOP2A", "CDC20", "PBK", "PTTG1", "GINS1", "KIF2C", "PCNA"),
    Metabolism_stress = c("DHCR24", "ALDH18A1", "EIF4EBP1", "NDRG1", "DDIT3", "HSPA5")
  )

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
    if (any(duplicated(rownames(expr)))) {
      summed <- rowsum(expr, group = rownames(expr), reorder = FALSE)
      counts <- as.numeric(table(rownames(expr))[rownames(summed)])
      expr <- summed / counts
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

  score_gene_sets <- function(z, gene_sets) {
    out <- list()
    for (nm in names(gene_sets)) {
      genes <- intersect(gene_sets[[nm]], rownames(z))
      if (length(genes) >= 2) {
        score <- colMeans(z[genes, , drop = FALSE], na.rm = TRUE)
      } else {
        score <- rep(NA_real_, ncol(z))
      }
      out[[nm]] <- data.frame(sample = colnames(z), signature = nm, score = as.numeric(score), n_genes_matched = length(genes), stringsAsFactors = FALSE)
    }
    do.call(rbind, out)
  }

  compare_high_low <- function(pathway_df, score_df) {
    igan_scores <- score_df[score_df$inferred_group == "IgAN" & !is.na(score_df$small_cluster_signature_score), , drop = FALSE]
    if (nrow(igan_scores) < 8) return(data.frame())
    q25 <- stats::quantile(igan_scores$small_cluster_signature_score, 0.25, na.rm = TRUE)
    q75 <- stats::quantile(igan_scores$small_cluster_signature_score, 0.75, na.rm = TRUE)
    igan_scores$signature_stratum <- ifelse(igan_scores$small_cluster_signature_score >= q75, "signature_high", ifelse(igan_scores$small_cluster_signature_score <= q25, "signature_low", "middle"))
    keep <- igan_scores[igan_scores$signature_stratum %in% c("signature_high", "signature_low"), c("sample", "signature_stratum"), drop = FALSE]
    df <- merge(pathway_df, keep, by = "sample")
    rows <- list()
    for (sig in unique(df$signature)) {
      x <- df[df$signature == sig & df$signature_stratum == "signature_high", "score"]
      y <- df[df$signature == sig & df$signature_stratum == "signature_low", "score"]
      x <- x[!is.na(x)]; y <- y[!is.na(y)]
      if (length(x) >= 2 && length(y) >= 2) {
        p <- tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
        rows[[length(rows) + 1]] <- data.frame(
          signature = sig,
          n_high = length(x),
          n_low = length(y),
          mean_high = mean(x),
          mean_low = mean(y),
          mean_difference_high_minus_low = mean(x) - mean(y),
          wilcox_p = p,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(rows) == 0) return(data.frame())
    out <- do.call(rbind, rows)
    out$wilcox_fdr <- stats::p.adjust(out$wilcox_p, method = "BH")
    out[order(out$wilcox_p), , drop = FALSE]
  }

  all_scores <- list()
  all_pathway_scores <- list()
  all_high_low <- list()
  summary_rows <- list()

  for (gse in datasets) {
    message("External validation: processing ", gse)
    res <- tryCatch({
      gsets <- GEOquery::getGEO(gse, GSEMatrix = TRUE)
      if (!is.list(gsets)) gsets <- list(gsets)
      dataset_scores <- list()
      dataset_pathways <- list()
      dataset_high_low <- list()
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
          z <- zscore_rows(expr[unique(c(common_all, unlist(curated_gene_sets))), , drop = FALSE])
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

          pathway_df <- score_gene_sets(z, curated_gene_sets)
          pathway_df$dataset <- gse
          pathway_df$platform <- platform
          pathway_df <- merge(pathway_df, score_df[, c("sample", "inferred_group", "small_cluster_signature_score")], by = "sample", all.x = TRUE)
          dataset_pathways[[length(dataset_pathways) + 1]] <- pathway_df

          high_low_df <- compare_high_low(pathway_df, score_df)
          if (nrow(high_low_df) > 0) {
            high_low_df$dataset <- gse
            high_low_df$platform <- platform
            dataset_high_low[[length(dataset_high_low) + 1]] <- high_low_df
          }
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
      list(scores = dataset_scores, pathways = dataset_pathways, high_low = dataset_high_low, summaries = dataset_summaries)
    }, error = function(e) {
      message("External validation failed for ", gse, ": ", conditionMessage(e))
      list(scores = list(), pathways = list(), high_low = list(), summaries = list(data.frame(
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
    all_pathway_scores <- c(all_pathway_scores, res$pathways)
    all_high_low <- c(all_high_low, res$high_low)
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

  if (length(all_pathway_scores) > 0) {
    pathway_scores <- do.call(rbind, all_pathway_scores)
    utils::write.csv(pathway_scores, file.path(tables_dir, "external_validation_pathway_signature_scores.csv"), row.names = FALSE)
  } else {
    pathway_scores <- data.frame()
    utils::write.csv(data.frame(note = "No pathway/immune signature scores computed."), file.path(tables_dir, "external_validation_pathway_signature_scores.csv"), row.names = FALSE)
  }

  if (length(all_high_low) > 0) {
    high_low <- do.call(rbind, all_high_low)
    high_low <- high_low[order(high_low$dataset, high_low$wilcox_p), , drop = FALSE]
    utils::write.csv(high_low, file.path(tables_dir, "external_validation_high_low_pathway_comparison.csv"), row.names = FALSE)
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      p2 <- ggplot2::ggplot(high_low, ggplot2::aes(x = signature, y = dataset, fill = mean_difference_high_minus_low)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(wilcox_fdr), "", ifelse(wilcox_fdr < 0.05, "*", ""))), size = 4) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::labs(title = "Pathway/immune differences in external IgAN signature-high vs signature-low samples", x = "Signature", y = "External dataset", fill = "High - Low")
      ggplot2::ggsave(file.path(figures_dir, "Figure7_external_high_low_pathway_heatmap.pdf"), p2, width = 11, height = 5.5)
    }
  } else {
    utils::write.csv(data.frame(note = "No high-vs-low pathway comparison was possible."), file.path(tables_dir, "external_validation_high_low_pathway_comparison.csv"), row.names = FALSE)
  }

  invisible(list(summary = summary_df))
}
