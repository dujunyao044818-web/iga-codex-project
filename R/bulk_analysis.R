load_bulk_gse <- function(cfg) {
  ensure_packages(c("GEOquery", "Biobase", "limma"))
  geo_id <- cfg$bulk$geo_accession
  message("Downloading ", geo_id, " from GEO")
  gse <- GEOquery::getGEO(geo_id, GSEMatrix = TRUE, getGPL = TRUE)
  eset <- if (length(gse) > 1) gse[[which.max(vapply(gse, ncol, integer(1)))]] else gse[[1]]
  expr <- Biobase::exprs(eset)
  pheno <- Biobase::pData(eset)
  fdata <- Biobase::fData(eset)

  message(geo_id, " phenotype columns: ", paste(colnames(pheno), collapse = ", "))
  save_table(data.frame(column = colnames(pheno)), file.path(cfg$output_dir, "tables", "GSE104948_pData_columns.csv"))

  gene_cols <- grep("symbol|gene", colnames(fdata), ignore.case = TRUE, value = TRUE)
  if (length(gene_cols) > 0) {
    gene_col <- gene_cols[1]
    genes <- clean_gene_symbols(as.character(fdata[[gene_col]]))
    keep <- !is.na(genes) & genes != "" & genes != "---"
    expr <- expr[keep, , drop = FALSE]
    genes <- genes[keep]
    expr <- limma::avereps(expr, ID = genes)
  }

  expr <- as.matrix(expr)
  mode(expr) <- "numeric"
  if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
  expr <- limma::normalizeBetweenArrays(expr, method = "quantile")
  expr <- expr[stats::complete.cases(expr), , drop = FALSE]
  expr <- expr[apply(expr, 1, stats::var, na.rm = TRUE) > 0, , drop = FALSE]

  save_table(expr, file.path(cfg$output_dir, "tables", "bulk_expression_matrix.csv"))
  save_table(pheno, file.path(cfg$output_dir, "tables", "bulk_sample_metadata_raw.csv"))
  list(expr = expr, pheno = pheno)
}

metadata_text <- function(pheno, columns) {
  present <- intersect(columns, colnames(pheno))
  if (!length(present)) present <- colnames(pheno)
  apply(pheno[, present, drop = FALSE], 1, paste, collapse = " ")
}

infer_bulk_groups <- function(pheno, cfg) {
  candidate_cols <- c("title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1", "characteristics_ch1.2", "characteristics_ch1.3", "diagnosis", "disease", "tissue", "group")
  txt <- tolower(metadata_text(pheno, candidate_cols))
  disease_pat <- paste(tolower(cfg$bulk$disease_terms), collapse = "|")
  control_pat <- paste(tolower(cfg$bulk$control_terms), collapse = "|")
  group <- rep("Unknown", length(txt))
  group[grepl(control_pat, txt)] <- "Control"
  group[grepl(disease_pat, txt)] <- "IgAN"
  factor(group)
}

bulk_qc <- function(expr, pheno, group, cfg) {
  pca <- stats::prcomp(t(expr), scale. = TRUE)
  pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
  df <- data.frame(sample = colnames(expr), PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, stringsAsFactors = FALSE)
  p <- ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = group)) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "GSE104948 annotated PCA", x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"), color = "Detected group")
  ggplot2::ggsave(file.path(cfg$output_dir, "figures", "Figure1_annotated_PCA.pdf"), p, width = 7, height = 5)
  ann <- data.frame(sample = rownames(pheno), detected_group = as.character(group), pheno, check.names = FALSE)
  save_table(ann, file.path(cfg$output_dir, "tables", "sample_annotation_clean.csv"))
  save_table(df, file.path(cfg$output_dir, "tables", "bulk_pca_coordinates.csv"))
  invisible(df)
}

run_consensus <- function(expr, cfg) {
  ensure_packages(c("ConsensusClusterPlus", "pheatmap"))
  vars <- apply(expr, 1, stats::var, na.rm = TRUE)
  n_genes <- min(cfg$bulk$consensus$max_genes, length(vars))
  top <- names(sort(vars, decreasing = TRUE))[seq_len(n_genes)]
  out_dir <- file.path(cfg$output_dir, "figures", "consensus_heatmaps")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cc <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(expr[top, ]), maxK = cfg$bulk$consensus$max_k, reps = cfg$bulk$consensus$reps, pItem = cfg$bulk$consensus$p_item, pFeature = cfg$bulk$consensus$p_feature, title = out_dir, plot = "pdf", seed = cfg$seed)
  k <- min(2, length(cc))
  subtype <- factor(paste0("Subtype", cc[[k]]$consensusClass))
  names(subtype) <- colnames(expr)
  tab <- data.frame(sample = names(subtype), subtype = as.character(subtype), stringsAsFactors = FALSE)
  utils::write.csv(tab, file.path(cfg$output_dir, "tables", "IgA_subtypes.csv"), row.names = FALSE)
  ann <- data.frame(Subtype = subtype)
  rownames(ann) <- colnames(expr)
  pheatmap::pheatmap(expr[top[seq_len(min(100, length(top)))], , drop = FALSE], scale = "row", annotation_col = ann, show_colnames = FALSE, filename = file.path(cfg$output_dir, "figures", "Figure2_subtype_heatmap.pdf"))
  subtype
}

run_limma <- function(expr, subtype, cfg) {
  ensure_packages("limma")
  subtype <- factor(subtype)
  design <- stats::model.matrix(~0 + subtype)
  colnames(design) <- levels(subtype)
  fit <- limma::lmFit(expr, design)
  if (nlevels(subtype) >= 2) {
    contrast <- paste(levels(subtype)[2], levels(subtype)[1], sep = "-")
    cont <- limma::makeContrasts(contrasts = contrast, levels = design)
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
    deg <- limma::topTable(fit2, number = Inf, sort.by = "P")
  } else {
    deg <- data.frame(logFC = rep(0, nrow(expr)), P.Value = rep(1, nrow(expr)), adj.P.Val = rep(1, nrow(expr)), row.names = rownames(expr))
  }
  deg$gene <- rownames(deg)
  utils::write.csv(deg, file.path(cfg$output_dir, "tables", "DEG_subtype_comparison.csv"), row.names = FALSE)
  deg$neglog10P <- -log10(pmax(deg$P.Value, 1e-300))
  p <- ggplot2::ggplot(deg, ggplot2::aes(logFC, neglog10P)) + ggplot2::geom_point(alpha = 0.6, size = 1.2) + ggplot2::theme_bw(base_size = 12) + ggplot2::labs(title = "Subtype differential expression", x = "log2 fold-change", y = "-log10(P)")
  ggplot2::ggsave(file.path(cfg$output_dir, "figures", "Figure3_volcano.pdf"), p, width = 7, height = 5)
  deg
}

get_subtype_markers <- function(deg, cfg) {
  n <- cfg$bulk$markers$top_n_per_subtype
  genes <- head(deg$gene[order(deg$P.Value, decreasing = FALSE)], n)
  writeLines(genes, file.path(cfg$output_dir, "tables", "subtype_marker_genes.txt"))
  genes
}
