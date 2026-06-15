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
  dirs <- file.path(cfg$output_dir, c("figures", "tables", "models"))
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
