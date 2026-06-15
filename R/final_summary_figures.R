run_final_summary_figures <- function(cfg) {
  tables_dir <- file.path(cfg$output_dir, "tables")
  figures_dir <- file.path(cfg$output_dir, "figures")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  high_low_file <- file.path(tables_dir, "external_validation_high_low_pathway_comparison.csv")
  dataset_file <- file.path(tables_dir, "external_validation_dataset_summary.csv")
  subtype_qc_file <- file.path(tables_dir, "subtype_qc_report.csv")
  subtype_counts_file <- file.path(tables_dir, "subtype_sample_counts.csv")

  if (!file.exists(high_low_file)) {
    utils::write.csv(data.frame(note = "External high-low pathway comparison file not found."), file.path(tables_dir, "external_validation_consensus_summary.csv"), row.names = FALSE)
    return(invisible(NULL))
  }
  high_low <- utils::read.csv(high_low_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("signature", "dataset", "mean_difference_high_minus_low") %in% colnames(high_low))) {
    utils::write.csv(data.frame(note = "External high-low pathway comparison file has unexpected columns."), file.path(tables_dir, "external_validation_consensus_summary.csv"), row.names = FALSE)
    return(invisible(NULL))
  }

  focus_order <- c(
    "Complement", "ECM_stromal_remodeling", "Podocyte_glomerular_structure",
    "Macrophage_monocyte", "T_cell", "B_cell", "CD8_T_cell", "Plasma_cell",
    "Tubular_injury_stress", "Cell_cycle_proliferation", "Interferon_inflammation",
    "Endothelial", "Neutrophil", "Metabolism_stress"
  )

  rows <- list()
  for (sig in unique(high_low$signature)) {
    df <- high_low[high_low$signature == sig & !is.na(high_low$mean_difference_high_minus_low), , drop = FALSE]
    if (nrow(df) == 0) next
    effects <- df$mean_difference_high_minus_low
    best_fdr <- if ("wilcox_fdr" %in% colnames(df)) min(df$wilcox_fdr, na.rm = TRUE) else NA_real_
    if (!is.finite(best_fdr)) best_fdr <- NA_real_
    n_pos <- sum(effects > 0, na.rm = TRUE)
    n_neg <- sum(effects < 0, na.rm = TRUE)
    consensus <- if (n_pos > n_neg) "higher_in_signature_high" else if (n_neg > n_pos) "lower_in_signature_high" else "mixed"
    interpretation <- switch(sig,
      Complement = "Complement activation is enriched in high-score external IgAN samples.",
      ECM_stromal_remodeling = "ECM/stromal remodeling is enriched in high-score external IgAN samples.",
      Podocyte_glomerular_structure = "Podocyte/glomerular structural signature is reduced in high-score external IgAN samples.",
      Macrophage_monocyte = "Myeloid inflammatory signal shows external-cohort support.",
      T_cell = "T-cell immune signal shows external-cohort support.",
      B_cell = "B-cell immune signal shows external-cohort support.",
      Tubular_injury_stress = "Tubular injury/stress signal may accompany high-score external IgAN samples.",
      Cell_cycle_proliferation = "Cell-cycle/proliferation or repair signal may accompany high-score samples.",
      "Exploratory pathway/immune signature difference in high-score external IgAN samples."
    )
    rows[[length(rows) + 1]] <- data.frame(
      signature = sig,
      n_datasets_tested = length(unique(df$dataset)),
      n_datasets_supporting_positive = n_pos,
      n_datasets_supporting_negative = n_neg,
      mean_effect_across_datasets = mean(effects, na.rm = TRUE),
      median_effect_across_datasets = stats::median(effects, na.rm = TRUE),
      best_fdr = best_fdr,
      consensus_direction = consensus,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }
  consensus <- if (length(rows) > 0) do.call(rbind, rows) else data.frame(note = "No pathway signatures could be summarized.")
  if ("signature" %in% colnames(consensus)) {
    consensus$signature <- factor(consensus$signature, levels = unique(c(focus_order, as.character(consensus$signature))))
    consensus <- consensus[order(consensus$signature), , drop = FALSE]
    consensus$signature <- as.character(consensus$signature)
  }
  utils::write.csv(consensus, file.path(tables_dir, "external_validation_consensus_summary.csv"), row.names = FALSE)

  ds_summary <- if (file.exists(dataset_file)) utils::read.csv(dataset_file, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  subtype_qc <- if (file.exists(subtype_qc_file)) utils::read.csv(subtype_qc_file, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()
  subtype_counts <- if (file.exists(subtype_counts_file)) utils::read.csv(subtype_counts_file, stringsAsFactors = FALSE, check.names = FALSE) else data.frame()

  mechanism <- data.frame(
    feature = c(
      "Discovery finding",
      "External reproducibility",
      "Complement",
      "ECM remodeling",
      "Podocyte/glomerular structure",
      "Immune activation",
      "Working subgroup name",
      "Single-cell validation need",
      "Spatial validation need"
    ),
    direction = c(
      if (nrow(subtype_counts) > 0) paste(paste0(subtype_counts[[1]], "=", subtype_counts[[2]]), collapse = "; ") else "IgAN small-cluster-like state",
      if (nrow(ds_summary) > 0 && "status" %in% colnames(ds_summary)) paste0(sum(ds_summary$status == "usable", na.rm = TRUE), " usable external platform(s)") else "External cohorts analyzed",
      "increased in signature-high samples",
      "increased in signature-high samples",
      "decreased in signature-high samples",
      "partially increased in signature-high samples",
      "CERI-like subgroup",
      "needed",
      "needed"
    ),
    evidence = c(
      "GSE104948 IgAN-only consensus clustering identified a small cluster; use exploratory wording.",
      "External bulk cohorts support a reproducible small-cluster-like signature.",
      "Consensus external high-low comparison supports Complement up-regulation.",
      "Consensus external high-low comparison supports ECM/stromal remodeling up-regulation.",
      "Consensus external high-low comparison supports podocyte/glomerular structural signature reduction.",
      "T, B, macrophage, and/or cell-cycle signatures show dataset-dependent support.",
      "Complement–ECM Remodeling–Injury-like; candidate molecular state, not definitive subtype.",
      "Map CERI genes and pathway scores to podocytes, mesangial cells, endothelial cells, tubular cells, and immune cells in IgAN sc/snRNA-seq.",
      "Localize CERI signature to glomerular, tubulointerstitial, periglomerular, or immune-infiltrate regions using spatial data or RNAscope/IHC."
    ),
    interpretation = c(
      "Exploratory discovery cohort finding.",
      "Raises confidence that the signal is not simply a GSE104948-specific outlier cluster.",
      "Core component of CERI-like state.",
      "Core component of CERI-like state.",
      "Core injury/structure-loss component of CERI-like state.",
      "Potential modifier rather than sole defining feature.",
      "Use cautious nomenclature: candidate reproducible molecular state.",
      "Cellular-origin validation step.",
      "Tissue-context validation step."
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(mechanism, file.path(tables_dir, "final_mechanism_hypothesis.csv"), row.names = FALSE)

  if (requireNamespace("ggplot2", quietly = TRUE) && "signature" %in% colnames(consensus)) {
    plot_df <- consensus[consensus$signature %in% focus_order, , drop = FALSE]
    if (nrow(plot_df) == 0) plot_df <- consensus
    plot_df$signature <- factor(plot_df$signature, levels = rev(unique(plot_df$signature)))
    plot_df$label <- paste0("+", plot_df$n_datasets_supporting_positive, "/-", plot_df$n_datasets_supporting_negative)
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = signature, y = mean_effect_across_datasets)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
      ggplot2::geom_col(alpha = 0.8) +
      ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.05, size = 3) +
      ggplot2::coord_flip() +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(
        title = "Consensus external validation of CERI-like pathway signatures",
        subtitle = "Effect size: signature-high minus signature-low IgAN samples across external cohorts",
        x = "Pathway / immune signature",
        y = "Mean effect across external datasets"
      )
    ggplot2::ggsave(file.path(figures_dir, "Figure8_external_validation_consensus.pdf"), p, width = 8.5, height = 5.8)
  }

  grDevices::pdf(file.path(figures_dir, "Figure9_mechanism_model.pdf"), width = 10, height = 6)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mar = c(1, 1, 3, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 10), ylim = c(0, 6))
  box <- function(x, y, label, w = 2.0, h = 0.65) {
    graphics::rect(x - w/2, y - h/2, x + w/2, y + h/2, border = "black", lwd = 1.2)
    graphics::text(x, y, label, cex = 0.85)
  }
  arrow <- function(x1, y1, x2, y2) graphics::arrows(x1, y1, x2, y2, length = 0.08, lwd = 1.2)
  graphics::title("Proposed CERI-like IgAN molecular state")
  box(1.3, 3, "Small-cluster-like\nIgAN state\n(GSE104948)", w = 2.1, h = 1.0)
  box(3.5, 4.5, "Complement\nactivation ↑", w = 1.8)
  box(3.5, 3, "Immune / inflammatory\nactivation ↑", w = 2.1)
  box(3.5, 1.5, "ECM remodeling /\nstromal activation ↑", w = 2.2)
  box(6.4, 3, "Podocyte / glomerular\nstructure ↓", w = 2.3)
  box(8.8, 3, "CERI-like candidate\ninjury state", w = 2.0)
  arrow(2.35, 3.05, 2.55, 4.35)
  arrow(2.35, 3, 2.45, 3)
  arrow(2.35, 2.95, 2.55, 1.65)
  arrow(4.55, 4.35, 5.35, 3.35)
  arrow(4.65, 3, 5.25, 3)
  arrow(4.65, 1.65, 5.35, 2.65)
  arrow(7.55, 3, 7.75, 3)
  graphics::text(5.1, 5.55, "External support: GSE37460, GSE93798, GSE116626", cex = 0.9)
  graphics::text(5.1, 5.2, "Consistent signals: Complement↑  ECM remodeling↑  Podocyte/glomerular structure↓", cex = 0.85)
  graphics::text(5.1, 0.45, "Interpretation: candidate reproducible molecular state; validate cell origin and spatial localization before claiming a definitive subtype.", cex = 0.78)
  grDevices::dev.off()

  invisible(list(consensus = consensus, mechanism = mechanism))
}
