run_gsva_scores <- function(expr, gene_sets, method = "gsva") {
  if (method == "ssgsea" && exists("ssgseaParam", where = asNamespace("GSVA"), mode = "function")) {
    param <- GSVA::ssgseaParam(as.matrix(expr), gene_sets)
    return(GSVA::gsva(param, verbose = FALSE))
  }
  if (exists("gsvaParam", where = asNamespace("GSVA"), mode = "function")) {
    param <- GSVA::gsvaParam(as.matrix(expr), gene_sets, kcdf = "Gaussian")
    return(GSVA::gsva(param, verbose = FALSE))
  }
  GSVA::gsva(as.matrix(expr), gene_sets, method = method, kcdf = "Gaussian", verbose = FALSE)
}

run_functional <- function(expr, subtype, cfg) {
  if (!requireNamespace("GSVA", quietly = TRUE) || !requireNamespace("msigdbr", quietly = TRUE)) {
    message("GSVA or msigdbr is unavailable; skipping Hallmark GSVA.")
    return(invisible(NULL))
  }
  tryCatch({
    ensure_packages(c("pheatmap", "reshape2"))
    hall <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    hsets <- split(hall$gene_symbol, hall$gs_name)
    imm <- hsets[grep("IMMUNE|INFLAM|INTERFERON|COMPLEMENT|IL6|TNFA|ALLOGRAFT", names(hsets))]
    gsva_scores <- run_gsva_scores(expr, hsets, method = "gsva")
    ssgsea_scores <- tryCatch(run_gsva_scores(expr, imm, method = "ssgsea"), error = function(e) {
      message("Immune ssGSEA failed: ", conditionMessage(e)); NULL
    })
    ann <- data.frame(subtype = subtype); rownames(ann) <- names(subtype)
    pheatmap::pheatmap(gsva_scores, annotation_col = ann, scale = "row",
                       filename = file.path(cfg$output_dir, "figures", "Figure4_GSVA_heatmap.pdf"),
                       width = 9, height = 10)
    save_table(gsva_scores, file.path(cfg$output_dir, "tables", "GSVA_hallmark_scores.csv"))
    if (!is.null(ssgsea_scores)) {
      pheatmap::pheatmap(ssgsea_scores, annotation_col = ann, scale = "row", filename = file.path(cfg$output_dir, "figures", "immune_ssgsea_heatmap.pdf"), width = 8, height = 7)
      save_table(ssgsea_scores, file.path(cfg$output_dir, "tables", "immune_ssgsea_scores.csv"))
      long <- reshape2::melt(ssgsea_scores); colnames(long) <- c("signature", "sample", "score"); long$subtype <- subtype[long$sample]
      ggplot(long, aes(subtype, score, fill = subtype)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = .15, size = .8) + facet_wrap(~signature, scales = "free_y") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      ggsave(file.path(cfg$output_dir, "figures", "immune_infiltration_boxplots.pdf"), width = 12, height = 8)
    }
    invisible(list(gsva = gsva_scores, immune = ssgsea_scores))
  }, error = function(e) {
    message("Hallmark GSVA failed and will be skipped: ", conditionMessage(e))
    invisible(NULL)
  })
}

run_lasso <- function(expr, subtype, markers, cfg) {
  ensure_packages(c("glmnet"))
  genes <- intersect(unique(unlist(markers)), rownames(expr))
  x <- t(expr[genes, , drop = FALSE]); y <- as.integer(subtype) - 1
  cv <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = cfg$lasso$alpha, nfolds = min(cfg$lasso$nfolds, length(y)))
  co <- as.matrix(coef(cv, s = "lambda.min")); co <- data.frame(gene = rownames(co), coefficient = co[,1]); co <- co[co$coefficient != 0,]
  risk <- as.numeric(predict(cv, newx = x, s = "lambda.min", type = "link"))
  saveRDS(cv, file.path(cfg$output_dir, "models", "lasso_subtype_model.rds"))
  save_table(co, file.path(cfg$output_dir, "tables", "lasso_coefficients.csv"))
  save_table(data.frame(sample = rownames(x), subtype = subtype, risk_score = risk), file.path(cfg$output_dir, "tables", "lasso_risk_scores.csv"))
  ggplot(co[co$gene != "(Intercept)",], aes(reorder(gene, coefficient), coefficient)) + geom_col() + coord_flip() + theme_bw()
  ggsave(file.path(cfg$output_dir, "figures", "lasso_coefficients.pdf"), width = 7, height = 6)
  cv
}
