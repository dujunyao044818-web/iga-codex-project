load_bulk_gse <- function(cfg) {
  ensure_packages(c("GEOquery", "Biobase", "limma"))
  message("Downloading ", cfg$bulk$geo_accession, " from GEO")
  gse <- GEOquery::getGEO(cfg$bulk$geo_accession, GSEMatrix = TRUE, getGPL = TRUE)
  eset <- if (length(gse) > 1) gse[[which.max(vapply(gse, ncol, integer(1)))]] else gse[[1]]
  expr <- Biobase::exprs(eset)
  pd <- Biobase::pData(eset)
  fd <- Biobase::fData(eset)
  message("GSE104948 phenotype columns: ", paste(colnames(pd), collapse = ", "))
  gene_col <- intersect(c("Gene Symbol", "gene_assignment", "GENE_SYMBOL", "Symbol", "gene_symbol"), colnames(fd))[1]
  if (!is.na(gene_col)) rownames(expr) <- clean_gene_symbols(as.character(fd[[gene_col]]))
  keep <- nzchar(rownames(expr)) & !is.na(rownames(expr))
  expr <- expr[keep, , drop = FALSE]
  expr <- limma::avereps(expr, ID = rownames(expr))
  if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
  expr <- limma::normalizeBetweenArrays(expr, method = "quantile")
  if (requireNamespace("sva", quietly = TRUE)) {
    batch_cols <- intersect(c("batch", "source_name_ch1", "characteristics_ch1.1"), colnames(pd))
    if (length(batch_cols) && length(unique(pd[[batch_cols[1]]])) > 1) {
      batch <- factor(pd[[batch_cols[1]]])
      mod <- model.matrix(~1, data = pd)
      expr <- sva::ComBat(dat = expr, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
    }
  }
  list(expr = expr, pheno = pd)
}

metadata_text <- function(pheno, columns) {
  present <- intersect(columns, colnames(pheno))
  if (!length(present)) present <- colnames(pheno)
  apply(pheno[, present, drop = FALSE], 1, paste, collapse = " ")
}

infer_bulk_groups <- function(pheno, cfg) {
  candidate_cols <- c("title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1",
                      "characteristics_ch1.2", "characteristics_ch1.3", "diagnosis", "disease",
                      "tissue", "group")
  txt <- tolower(metadata_text(pheno, candidate_cols))
  disease_pat <- paste(cfg$bulk$disease_terms, collapse = "|")
  control_pat <- paste(cfg$bulk$control_terms, collapse = "|")
  group <- rep(NA_character_, length(txt))
  group[grepl(control_pat, txt)] <- "Control"
  group[grepl(disease_pat, txt)] <- "IgAN"
  group <- ifelse(is.na(group), "Unknown", group)
  factor(group, levels = unique(c("Control", "IgAN", "Unknown", sort(unique(group)))))
}

export_sample_annotation <- function(pheno, group, cfg) {
  candidate_cols <- c("title", "source_name_ch1", "characteristics_ch1", "characteristics_ch1.1",
                      "characteristics_ch1.2", "characteristics_ch1.3", "diagnosis", "disease",
                      "tissue", "group")
  present <- intersect(candidate_cols, colnames(pheno))
  ann <- data.frame(sample = rownames(pheno), detected_group = as.character(group), stringsAsFactors = FALSE)
  ann <- cbind(ann, pheno[, present, drop = FALSE])
  save_table(ann, file.path(cfg$output_dir, "tables", "sample_annotation_clean.csv"))
  save_table(data.frame(column = colnames(pheno)), file.path(cfg$output_dir, "tables", "GSE104948_pData_columns.csv"))
  ann
}

bulk_qc <- function(expr, pheno, group, cfg) {
  export_sample_annotation(pheno, group, cfg)
  pca <- prcomp(t(expr), scale. = TRUE)
  var_pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)
  pc <- data.frame(sample = colnames(expr), group = group, PC1 = pca$x[,1], PC2 = pca$x[,2])
  p <- ggplot(pc, aes(PC1, PC2, color = group)) +
    geom_point(size = 3.5, alpha = 0.9) +
    labs(title = "GSE104948 IgA nephropathy bulk RNA-seq PCA",
         x = paste0("PC1 (", var_pct[1], "% variance)"),
         y = paste0("PC2 (", var_pct[2], "% variance)"),
         color = "Detected group") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "right")
  ggsave(file.path(cfg$output_dir, "figures", "Figure1_annotated_PCA.pdf"), p, width = 7, height = 5.5)
  save_table(pc, file.path(cfg$output_dir, "tables", "bulk_pca_coordinates.csv"))
}

select_variable_genes <- function(expr, n = 2000) {
  vars <- matrixStats::rowVars(expr)
  rownames(expr)[order(vars, decreasing = TRUE)[seq_len(min(n, length(vars)))]]
}

run_consensus <- function(expr, cfg) {
  ensure_packages(c("ConsensusClusterPlus", "matrixStats", "pheatmap"))
  top_genes <- select_variable_genes(expr, n = 2000)
  message("Running consensus clustering on top ", length(top_genes), " variable genes, k=2..", cfg$bulk$consensus$max_k)
  consensus_dir <- file.path(cfg$output_dir, "figures", "consensus_heatmaps")
  dir.create(consensus_dir, recursive = TRUE, showWarnings = FALSE)
  res <- ConsensusClusterPlus::ConsensusClusterPlus(as.matrix(expr[top_genes, ]), maxK = cfg$bulk$consensus$max_k,
    reps = cfg$bulk$consensus$reps, pItem = cfg$bulk$consensus$p_item, pFeature = cfg$bulk$consensus$p_feature,
    clusterAlg = "hc", distance = "pearson", seed = cfg$seed, plot = "pdf", title = consensus_dir)
  k <- 2
  clusters <- factor(res[[k]]$consensusClass, labels = paste0("Subtype", seq_len(k)))
  names(clusters) <- colnames(expr)
  save_table(data.frame(sample = names(clusters), subtype = clusters), file.path(cfg$output_dir, "tables", "IgA_subtypes.csv"))
  plot_subtype_heatmap(expr, clusters, cfg, top_genes)
  clusters
}

plot_subtype_heatmap <- function(expr, subtype, cfg, top_genes = NULL) {
  ensure_packages("pheatmap")
  if (is.null(top_genes)) top_genes <- select_variable_genes(expr, n = 2000)
  heat_genes <- head(top_genes, min(100, length(top_genes)))
  ann <- data.frame(subtype = subtype)
  rownames(ann) <- names(subtype)
  pheatmap::pheatmap(expr[heat_genes, names(subtype), drop = FALSE], scale = "row", show_rownames = FALSE,
                     annotation_col = ann, clustering_distance_cols = "correlation",
                     filename = file.path(cfg$output_dir, "figures", "Figure2_subtype_heatmap.pdf"),
                     width = 8, height = 9)
}

run_limma <- function(expr, subtype, cfg) {
  ensure_packages("limma")
  if (nlevels(subtype) < 2) stop("Need at least two consensus subtypes for limma comparison")
  design <- model.matrix(~0 + subtype)
  colnames(design) <- levels(subtype)
  fit <- limma::lmFit(expr, design)
  cont <- limma::makeContrasts(Subtype2_vs_Subtype1 = Subtype2 - Subtype1, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))
  deg <- limma::topTable(fit2, number = Inf, sort.by = "P")
  deg$gene <- rownames(deg)
  deg$significant <- deg$adj.P.Val < cfg$bulk$differential_expression$fdr_cutoff & abs(deg$logFC) > cfg$bulk$differential_expression$logfc_cutoff
  save_table(deg, file.path(cfg$output_dir, "tables", "DEG_subtype_comparison.csv"))
  p <- ggplot(deg, aes(logFC, -log10(adj.P.Val), color = significant)) +
    geom_point(alpha = .75, size = 1.5) +
    geom_vline(xintercept = c(-cfg$bulk$differential_expression$logfc_cutoff, cfg$bulk$differential_expression$logfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(cfg$bulk$differential_expression$fdr_cutoff), linetype = "dashed") +
    scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#D55E00")) +
    labs(title = "Subtype differential expression", x = "log2 fold-change (Subtype2 vs Subtype1)", y = "-log10 adjusted P value", color = "Significant") +
    theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  ggsave(file.path(cfg$output_dir, "figures", "Figure3_volcano.pdf"), p, width = 7, height = 5.5)
  deg
}

get_subtype_markers <- function(deg, cfg) {
  n <- cfg$bulk$markers$top_n_per_subtype
  up <- head(deg$gene[order(-deg$logFC)], n)
  dn <- head(deg$gene[order(deg$logFC)], n)
  markers <- list(Subtype1 = dn, Subtype2 = up)
  saveRDS(markers, file.path(cfg$output_dir, "models", "bulk_subtype_signatures.rds"))
  utils::write.csv(stack(markers), file.path(cfg$output_dir, "tables", "subtype_marker_genes.csv"), row.names = FALSE)
  markers
}
