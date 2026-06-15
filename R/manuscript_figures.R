run_manuscript_figures <- function(cfg) {
  figures_dir <- file.path(cfg$output_dir, "figures")
  tables_dir <- file.path(cfg$output_dir, "tables")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

  has_gg <- requireNamespace("ggplot2", quietly = TRUE)

  save_plot <- function(plot, filename, width = 7, height = 5, dpi = 300) {
    if (!has_gg) return(invisible(FALSE))
    pdf_file <- file.path(figures_dir, paste0(filename, ".pdf"))
    png_file <- file.path(figures_dir, paste0(filename, ".png"))
    ggplot2::ggsave(pdf_file, plot, width = width, height = height, device = grDevices::cairo_pdf)
    ggplot2::ggsave(png_file, plot, width = width, height = height, dpi = dpi)
    invisible(TRUE)
  }

  placeholder_plot <- function(title, subtitle) {
    if (!has_gg) return(NULL)
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.6, label = title, fontface = "bold", size = 5) +
      ggplot2::annotate("text", x = 0.5, y = 0.42, label = subtitle, size = 3.8) +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
  }

  read_csv_safe <- function(path) {
    if (!file.exists(path)) return(NULL)
    tryCatch(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  }

  if (!has_gg) {
    message("ggplot2 not available; manuscript figures skipped.")
    return(invisible(NULL))
  }

  # Figure 5: consensus QC summary, figure-first version of clustering diagnostics.
  consensus <- read_csv_safe(file.path(tables_dir, "consensus_k_summary.csv"))
  if (!is.null(consensus) && nrow(consensus) > 0) {
    k_col <- intersect(c("k", "K", "cluster_k", "n_clusters"), colnames(consensus))[1]
    if (is.na(k_col)) k_col <- colnames(consensus)[1]
    numeric_cols <- colnames(consensus)[vapply(consensus, is.numeric, logical(1))]
    metric_cols <- setdiff(numeric_cols, k_col)
    if (length(metric_cols) > 0) {
      metric <- metric_cols[1]
      p <- ggplot2::ggplot(consensus, ggplot2::aes(x = .data[[k_col]], y = .data[[metric]])) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::geom_point(size = 2.3) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold")) +
        ggplot2::labs(
          title = "Consensus clustering QC summary",
          subtitle = "Numerical clustering diagnostics are retained in supplementary tables",
          x = "Candidate cluster number", y = metric
        )
    } else {
      p <- placeholder_plot("Consensus clustering QC", "No numeric consensus metric was available for plotting; see supplementary table.")
    }
  } else {
    p <- placeholder_plot("Consensus clustering QC", "Consensus summary table was not available; see pipeline QC outputs.")
  }
  save_plot(p, "Figure5_consensus_QC", width = 6.5, height = 4.8)

  # Figure 6: CERI marker mechanism summary.
  marker_tbl <- read_csv_safe(file.path(tables_dir, "subtype_marker_direction_table.csv"))
  if (!is.null(marker_tbl) && nrow(marker_tbl) > 0) {
    gene_col <- intersect(c("gene", "Gene", "symbol", "Symbol", "marker_gene"), colnames(marker_tbl))[1]
    effect_col <- intersect(c("logFC", "log2FC", "effect", "direction_score", "estimate"), colnames(marker_tbl))[1]
    if (!is.na(gene_col) && !is.na(effect_col)) {
      plot_tbl <- marker_tbl[is.finite(marker_tbl[[effect_col]]), , drop = FALSE]
      plot_tbl <- plot_tbl[order(abs(plot_tbl[[effect_col]]), decreasing = TRUE), , drop = FALSE]
      plot_tbl <- head(plot_tbl, 20)
      plot_tbl[[gene_col]] <- factor(plot_tbl[[gene_col]], levels = rev(plot_tbl[[gene_col]]))
      p <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = .data[[effect_col]], y = .data[[gene_col]])) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold")) +
        ggplot2::labs(
          title = "CERI-like marker program",
          subtitle = "Top marker genes highlight inflammatory, ECM-remodeling and glomerular-structure programs",
          x = effect_col, y = "Marker gene"
        )
    } else {
      p <- placeholder_plot("CERI-like marker program", "Marker table did not contain expected gene/effect columns; see supplementary table.")
    }
  } else {
    p <- placeholder_plot("CERI-like marker program", "Marker direction table was not available; see DEG and marker outputs.")
  }
  save_plot(p, "Figure6_CERI_marker_mechanisms", width = 7, height = 5.2)

  manifest <- data.frame(
    figure = c("Figure5_consensus_QC", "Figure6_CERI_marker_mechanisms"),
    pdf = file.path("output", "figures", c("Figure5_consensus_QC.pdf", "Figure6_CERI_marker_mechanisms.pdf")),
    png = file.path("output", "figures", c("Figure5_consensus_QC.png", "Figure6_CERI_marker_mechanisms.png")),
    role = c("Figure-first consensus clustering QC summary", "Figure-first CERI marker mechanism summary"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(tables_dir, "manuscript_figure_manifest.csv"), row.names = FALSE)
  message("Figure-first manuscript figure package completed.")
  invisible(manifest)
}
