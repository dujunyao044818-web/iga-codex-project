ceri_modules <- function() {
  list(
    CERI = c("C3", "C4A", "C4B", "C1QA", "C1QB", "C1QC", "CFB", "VCAN", "COL5A2", "LAMB1", "PLAU", "HTRA1", "IL18", "ALOX5AP", "FYN"),
    Complement = c("C3", "C4A", "C4B", "C1QA", "C1QB", "C1QC", "CFB"),
    ECM_remodeling = c("VCAN", "COL5A2", "LAMB1", "PLAU", "HTRA1", "COL1A1", "COL1A2", "COL3A1", "FN1"),
    Podocyte_glomerular_structure = c("SYNPO", "PLCE1", "TCF21", "MAGI2", "CLIC5", "NPHS1", "NPHS2", "WT1", "PODXL"),
    Macrophage_monocyte = c("CD68", "LST1", "AIF1", "C1QA", "C1QB", "C1QC", "FCGR3A"),
    Tubular_injury = c("HAVCR1", "LCN2", "VCAM1", "NDRG1", "ALDH18A1", "EIF4EBP1"),
    B_cell = c("MS4A1", "CD79A", "CD79B", "BANK1"),
    T_cell = c("CD3D", "CD3E", "TRAC", "IL7R", "CD2"),
    Plasma_cell = c("MZB1", "JCHAIN", "XBP1", "IGHG1", "IGKC"),
    Cell_cycle_repair = c("CDC20", "PBK", "PTTG1", "GINS1", "KIF2C")
  )
}

module_scores_from_expr <- function(expr, modules) {
  scores <- lapply(names(modules), function(module) {
    genes <- intersect(modules[[module]], rownames(expr))
    if (!length(genes)) return(rep(NA_real_, ncol(expr)))
    colMeans(expr[genes, , drop = FALSE], na.rm = TRUE)
  })
  out <- as.data.frame(scores, check.names = FALSE)
  colnames(out) <- names(modules)
  rownames(out) <- colnames(expr)
  out
}

run_ceri_module_correlation_network <- function(bulk, cfg) {
  message("Running CERI module correlation network.")
  modules <- ceri_modules()
  scores <- module_scores_from_expr(bulk$expr, modules)
  valid <- vapply(scores, function(x) sum(is.finite(x)) >= 3 && stats::sd(x, na.rm = TRUE) > 0, logical(1))
  scores_valid <- scores[, valid, drop = FALSE]
  if (ncol(scores_valid) >= 2) {
    cor_mat <- stats::cor(scores_valid, use = "pairwise.complete.obs", method = "spearman")
  } else {
    cor_mat <- matrix(NA_real_, nrow = length(modules), ncol = length(modules), dimnames = list(names(modules), names(modules)))
  }
  save_table(as.data.frame(cor_mat), file.path(cfg$output_dir, "tables", "CERI_module_correlation_matrix.csv"))
  if ("CERI" %in% rownames(cor_mat)) {
    mods <- setdiff(colnames(cor_mat), "CERI")
    consensus <- data.frame(
      module = mods,
      discovery_spearman_with_CERI = as.numeric(cor_mat["CERI", mods]),
      consensus_direction = ifelse(as.numeric(cor_mat["CERI", mods]) >= 0, "positive", "negative"),
      datasets_available = cfg$bulk$geo_accession,
      stringsAsFactors = FALSE
    )
  } else {
    consensus <- data.frame(module = names(modules), discovery_spearman_with_CERI = NA_real_, consensus_direction = NA_character_, datasets_available = cfg$bulk$geo_accession, stringsAsFactors = FALSE)
  }
  save_table(consensus, file.path(cfg$output_dir, "tables", "CERI_module_correlation_consensus.csv"))
  supp_dir <- file.path(cfg$output_dir, "supplementary_figures")
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  if (all(is.na(cor_mat))) {
    save_placeholder_figure(file.path(supp_dir, "FigureS_CERI_module_correlation_heatmap"), "CERI module correlation heatmap", "Insufficient finite module scores.", width = 7, height = 6)
    save_placeholder_figure(file.path(supp_dir, "FigureS_CERI_module_correlation_network"), "CERI module correlation network", "Insufficient finite module scores.", width = 7, height = 6)
  } else {
    pheatmap::pheatmap(cor_mat, color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
                       breaks = seq(-1, 1, length.out = 102), filename = file.path(supp_dir, "FigureS_CERI_module_correlation_heatmap.pdf"),
                       width = 7, height = 6)
    grDevices::png(file.path(supp_dir, "FigureS_CERI_module_correlation_heatmap.png"), width = 2100, height = 1800, res = 300)
    pheatmap::pheatmap(cor_mat, color = grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(101), breaks = seq(-1, 1, length.out = 102))
    grDevices::dev.off()
    grDevices::pdf(file.path(supp_dir, "FigureS_CERI_module_correlation_network.pdf"), width = 8, height = 6)
    par(mar = c(1, 1, 3, 1))
    plot.new(); plot.window(xlim = c(-1.2, 1.2), ylim = c(-1.05, 1.05), asp = 1)
    title("CERI module correlation network")
    mods <- setdiff(colnames(cor_mat), "CERI")
    theta <- seq(0, 2 * pi, length.out = length(mods) + 1)[-1]
    xy <- data.frame(module = mods, x = cos(theta), y = sin(theta))
    points(0, 0, pch = 21, bg = "#D55E00", cex = 4)
    text(0, 0, "CERI", col = "white", font = 2)
    for (i in seq_len(nrow(xy))) {
      r <- cor_mat["CERI", xy$module[i]]
      col <- ifelse(is.na(r), "grey80", ifelse(r >= 0, "#B2182B", "#2166AC"))
      lwd <- ifelse(is.na(r), 1, 1 + 4 * abs(r))
      segments(0, 0, xy$x[i], xy$y[i], col = col, lwd = lwd)
      points(xy$x[i], xy$y[i], pch = 21, bg = "white", cex = 3)
      text(xy$x[i], xy$y[i], xy$module[i], cex = 0.75)
    }
    grDevices::dev.off()
  }
  invisible(consensus)
}

run_immune_signature_visualizations <- function(bulk, cfg) {
  message("Running immune signature CERI-high/low visualizations.")
  modules <- ceri_modules()
  scores <- module_scores_from_expr(bulk$expr, modules)
  immune <- intersect(c("Macrophage_monocyte", "T_cell", "B_cell", "Plasma_cell"), colnames(scores))
  if (!"CERI" %in% colnames(scores) || !length(immune) || all(!is.finite(scores$CERI))) {
    save_table(data.frame(status = "skipped", reason = "No finite CERI/immune scores"), file.path(cfg$output_dir, "tables", "immune_signature_CERI_high_low_comparison.csv"))
    save_placeholder_figure(file.path(cfg$output_dir, "supplementary_figures", "FigureS_immune_signature_heatmap"), "Immune signature heatmap", "No finite CERI/immune scores.", width = 7, height = 5)
    save_placeholder_figure(file.path(cfg$output_dir, "supplementary_figures", "FigureS_immune_signature_boxplot"), "Immune signature boxplot", "No finite CERI/immune scores.", width = 7, height = 5)
    return(invisible(NULL))
  }
  ceri_group <- ifelse(scores$CERI >= stats::median(scores$CERI, na.rm = TRUE), "CERI-high", "CERI-low")
  comp <- do.call(rbind, lapply(immune, function(sig) {
    x <- scores[[sig]]
    p <- tryCatch(stats::wilcox.test(x ~ ceri_group)$p.value, error = function(e) NA_real_)
    data.frame(signature = sig, mean_CERI_high = mean(x[ceri_group == "CERI-high"], na.rm = TRUE),
               mean_CERI_low = mean(x[ceri_group == "CERI-low"], na.rm = TRUE),
               difference = mean(x[ceri_group == "CERI-high"], na.rm = TRUE) - mean(x[ceri_group == "CERI-low"], na.rm = TRUE),
               p_value = p, stringsAsFactors = FALSE)
  }))
  comp$FDR <- p.adjust(comp$p_value, method = "BH")
  save_table(comp, file.path(cfg$output_dir, "tables", "immune_signature_CERI_high_low_comparison.csv"))
  plot_mat <- t(scale(scores[, immune, drop = FALSE]))
  ann <- data.frame(CERI_group = ceri_group)
  rownames(ann) <- rownames(scores)
  pheatmap::pheatmap(plot_mat, annotation_col = ann, show_colnames = FALSE,
                     filename = file.path(cfg$output_dir, "supplementary_figures", "FigureS_immune_signature_heatmap.pdf"), width = 7, height = 4.5)
  df <- reshape2::melt(data.frame(sample = rownames(scores), CERI_group = ceri_group, scores[, immune, drop = FALSE]), id.vars = c("sample", "CERI_group"))
  p <- ggplot2::ggplot(df, ggplot2::aes(CERI_group, value, fill = CERI_group)) + ggplot2::geom_boxplot(outlier.size = 0.6) +
    ggplot2::facet_wrap(~variable, scales = "free_y") + ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Immune signatures by CERI-high/low status", x = NULL, y = "Module score") +
    ggplot2::theme(legend.position = "none")
  save_ggplot_pdf_png(p, file.path(cfg$output_dir, "supplementary_figures", "FigureS_immune_signature_boxplot"), width = 8, height = 5)
  invisible(comp)
}

run_ppi_hub_gene_package <- function(cfg) {
  message("Running PPI/hub gene fallback package.")
  deg_path <- file.path(cfg$output_dir, "tables", "DEG_subtype_comparison.csv")
  deg <- if (file.exists(deg_path)) utils::read.csv(deg_path, stringsAsFactors = FALSE, check.names = FALSE) else data.frame(gene = character(), adj.P.Val = numeric(), logFC = numeric())
  if (!"gene" %in% colnames(deg)) deg$gene <- rownames(deg)
  if (!"adj.P.Val" %in% colnames(deg)) deg$adj.P.Val <- seq_len(nrow(deg))
  modules <- ceri_modules()
  genes <- unique(c(head(deg$gene[order(deg$adj.P.Val)], 80), unlist(modules)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  membership <- utils::stack(modules); colnames(membership) <- c("gene", "module")
  membership <- membership[membership$gene %in% genes, ]
  edge_list <- do.call(rbind, lapply(split(membership$gene, membership$module), function(g) {
    if (length(g) < 2) return(NULL)
    cmb <- utils::combn(unique(g), 2)
    data.frame(from = cmb[1, ], to = cmb[2, ], evidence = "curated_CERI_module_co_membership", stringsAsFactors = FALSE)
  }))
  if (is.null(edge_list)) edge_list <- data.frame(from = character(), to = character(), evidence = character())
  deg_count <- table(c(edge_list$from, edge_list$to))
  ranking <- data.frame(gene = genes, degree = as.integer(deg_count[genes]), stringsAsFactors = FALSE)
  ranking$degree[is.na(ranking$degree)] <- 0L
  ranking$betweenness_fallback <- ranking$degree / pmax(max(ranking$degree), 1)
  ranking$curated_module_overlap <- vapply(ranking$gene, function(g) sum(vapply(modules, function(m) g %in% m, logical(1))), integer(1))
  ranking <- ranking[order(-ranking$degree, -ranking$curated_module_overlap), ]
  save_table(ranking, file.path(cfg$output_dir, "tables", "PPI_hub_gene_ranking.csv"))
  save_table(edge_list, file.path(cfg$output_dir, "tables", "PPI_STRING_compatible_edge_table.csv"))
  grDevices::pdf(file.path(cfg$output_dir, "supplementary_figures", "FigureS_PPI_hub_gene_network.pdf"), width = 8, height = 6)
  par(mar = c(1, 1, 3, 1)); plot.new(); plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), asp = 1)
  title("Curated CERI PPI/hub-gene fallback network")
  top <- head(ranking$gene, 25)
  if (length(top)) {
    theta <- seq(0, 2*pi, length.out = length(top) + 1)[-1]
    xy <- data.frame(gene = top, x = cos(theta), y = sin(theta))
    for (i in seq_len(nrow(edge_list))) {
      if (edge_list$from[i] %in% top && edge_list$to[i] %in% top) {
        a <- xy[xy$gene == edge_list$from[i], ]; b <- xy[xy$gene == edge_list$to[i], ]
        segments(a$x, a$y, b$x, b$y, col = "grey80")
      }
    }
    points(xy$x, xy$y, pch = 21, bg = "#E69F00", cex = 2)
    text(xy$x, xy$y, xy$gene, cex = 0.75, pos = 3)
  }
  grDevices::dev.off()
  writeLines(c("# PPI / hub gene report", "", "STRING online access is not required. This package writes a STRING-compatible edge table and a documented fallback network based on curated CERI module co-membership.", "", "Hub ranking uses degree, a fallback betweenness-like scaled degree, and overlap with curated CERI modules."), file.path(cfg$output_dir, "reports", "PPI_hub_gene_report.md"))
  invisible(ranking)
}

run_wgcna_exploratory_or_report <- function(bulk, cfg) {
  n_available <- ncol(bulk$expr)
  report <- c("# WGCNA exploratory report", "", paste0("Discovery matrix sample size available to this module: ", n_available), "", "WGCNA is not forced as a main analysis because the IgAN-only discovery sample size and the exploratory small-cluster subtype structure are not ideal for stable module discovery.", "", "Recommendation: run WGCNA only in larger external bulk cohorts after harmonized preprocessing and subject-level metadata curation.")
  writeLines(report, file.path(cfg$output_dir, "reports", "WGCNA_exploratory_report.md"))
  save_table(data.frame(status = "not_run", reason = "underpowered discovery cohort for robust WGCNA main-story claim", stringsAsFactors = FALSE), file.path(cfg$output_dir, "tables", "WGCNA_module_trait_correlation.csv"))
  save_placeholder_figure(file.path(cfg$output_dir, "supplementary_figures", "FigureS_WGCNA_module_trait_heatmap"), "Exploratory WGCNA", "Not run: discovery cohort is underpowered for robust WGCNA claims.", width = 7, height = 5)
  invisible(TRUE)
}

run_ceri_supplementary_extensions <- function(bulk, cfg) {
  run_ceri_module_correlation_network(bulk, cfg)
  run_immune_signature_visualizations(bulk, cfg)
  run_ppi_hub_gene_package(cfg)
  run_wgcna_exploratory_or_report(bulk, cfg)
  invisible(TRUE)
}
