run_functional <- function(expr, subtype, cfg) {
  out <- list()
  immune_sets <- list(
    T_cell = c("CD3D", "CD3E", "CD2", "TRAC"),
    CD8_T = c("CD8A", "CD8B", "GZMB", "NKG7"),
    B_cell = c("MS4A1", "CD79A", "CD79B", "CD19"),
    Macrophage = c("CD68", "C1QA", "C1QB", "CSF1R"),
    Neutrophil = c("S100A8", "S100A9", "FCGR3B", "CXCR2")
  )
  score_one <- function(genes) {
    g <- intersect(genes, rownames(expr))
    if (!length(g)) return(rep(NA_real_, ncol(expr)))
    colMeans(expr[g, , drop = FALSE], na.rm = TRUE)
  }
  immune <- t(vapply(immune_sets, score_one, numeric(ncol(expr))))
  colnames(immune) <- colnames(expr)
  save_table(immune, file.path(cfg$output_dir, "tables", "immune_signature_scores.csv"))
  out$immune <- immune

  if (requireNamespace("GSVA", quietly = TRUE) && requireNamespace("msigdbr", quietly = TRUE)) {
    gs <- tryCatch({
      hall <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
      gsets <- split(hall$gene_symbol, hall$gs_name)
      GSVA::gsva(as.matrix(expr), gsets, verbose = FALSE)
    }, error = function(e) {
      message("GSVA skipped: ", conditionMessage(e))
      NULL
    })
    if (!is.null(gs)) {
      save_table(gs, file.path(cfg$output_dir, "tables", "GSVA_hallmark_scores.csv"))
      ann <- data.frame(Subtype = subtype)
      rownames(ann) <- colnames(expr)
      pheatmap::pheatmap(gs, scale = "row", annotation_col = ann, show_colnames = FALSE, filename = file.path(cfg$output_dir, "figures", "Figure4_GSVA_heatmap.pdf"))
      out$gsva <- gs
    }
  } else {
    message("GSVA/msigdbr not available; skipping Hallmark GSVA.")
  }
  invisible(out)
}

run_lasso <- function(expr, subtype, markers, cfg) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    message("glmnet unavailable; skipping LASSO.")
    return(NULL)
  }
  subtype <- factor(subtype)
  if (nlevels(subtype) < 2) {
    message("Only one subtype; skipping LASSO.")
    return(NULL)
  }
  genes <- intersect(markers, rownames(expr))
  if (length(genes) < 2) {
    vars <- sort(apply(expr, 1, stats::var), decreasing = TRUE)
    genes <- names(vars)[seq_len(min(100, length(vars)))]
  }
  x <- t(expr[genes, , drop = FALSE])
  y <- as.numeric(subtype == levels(subtype)[2])
  cv <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = cfg$lasso$alpha, nfolds = min(cfg$lasso$nfolds, length(y)))
  saveRDS(cv, file.path(cfg$output_dir, "models", "lasso_model.rds"))
  grDevices::pdf(file.path(cfg$output_dir, "figures", "Figure5_LASSO_cv.pdf"), width = 7, height = 5)
  plot(cv)
  grDevices::dev.off()
  coef_df <- as.data.frame(as.matrix(stats::coef(cv, s = "lambda.min")))
  colnames(coef_df) <- "coefficient"
  coef_df$gene <- rownames(coef_df)
  utils::write.csv(coef_df, file.path(cfg$output_dir, "tables", "LASSO_coefficients.csv"), row.names = FALSE)
  cv
}
