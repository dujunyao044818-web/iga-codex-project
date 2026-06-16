suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
})

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "),
         ". Install them before running the pipeline.", call. = FALSE)
  }
}

make_dirs <- function(cfg) {
  dirs <- file.path(cfg$output_dir, c(
    "figures", "tables", "models", "reports", "main_figures", "supplementary_figures",
    "supplementary/qc", "supplementary/ml_external_bulk", "supplementary/functional_enrichment"
  ))
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  invisible(lapply(c("data/raw", "data/processed"), dir.create, recursive = TRUE, showWarnings = FALSE))
}

save_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, path, row.names = TRUE)
}

clean_gene_symbols <- function(x) {
  x <- sub("///.*$", "", x)
  x <- trimws(x)
  make.unique(x)
}

save_ggplot_pdf_png <- function(plot, stem, width = 7, height = 5, dpi = 300) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi)
  invisible(TRUE)
}

save_placeholder_figure <- function(stem, title, message, width = 7, height = 5) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  p <- ggplot2::ggplot(data.frame(x = 0, y = 0), ggplot2::aes(x, y)) +
    ggplot2::geom_blank() +
    ggplot2::annotate("text", x = 0, y = 0.2, label = title, fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = 0, y = -0.1, label = message, size = 3.5) +
    ggplot2::theme_void() +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1)
  save_ggplot_pdf_png(p, stem, width = width, height = height)
  invisible(TRUE)
}
