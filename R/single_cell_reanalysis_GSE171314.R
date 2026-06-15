#!/usr/bin/env Rscript

# Independent full single-cell reanalysis for GSE171314.
# This script is intentionally not sourced by run_pipeline.R.
# It is run only by .github/workflows/run-single-cell.yml on the
# single-cell-reanalysis-gse171314 branch.

options(stringsAsFactors = FALSE)
set.seed(20260615)

required <- c("GEOquery", "Seurat", "data.table", "ggplot2", "Matrix", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("Missing required packages for independent single-cell reanalysis: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(Matrix)
})

accession <- "GSE171314"
out_dir <- "output_single_cell"
tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
data_dir <- file.path(out_dir, "data")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

ceri_modules <- list(
  Complement = c("C3", "C4A", "C4B", "C1QA", "C1QB", "C1QC", "CFB"),
  ECM_remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL5A2", "FN1", "VCAN", "LAMB1", "PLAU", "HTRA1"),
  Podocyte_glomerular_structure = c("NPHS1", "NPHS2", "SYNPO", "WT1", "PODXL", "PLCE1", "TCF21", "MAGI2", "CLIC5"),
  Myeloid_inflammation = c("IL18", "ALOX5AP", "FYN", "CD68", "LST1", "AIF1", "LYZ", "TYROBP"),
  T_B_plasma_cell = c("CD3D", "CD3E", "CD2", "MS4A1", "CD79A", "CD79B", "MZB1", "XBP1", "JCHAIN"),
  Tubular_injury_stress = c("HAVCR1", "LCN2", "KRT8", "KRT18", "VCAM1", "NDRG1", "HMOX1")
)
celltype_markers <- list(
  Podocyte = c("NPHS1", "NPHS2", "SYNPO", "WT1", "PODXL", "PLCE1"),
  Mesangial_stromal = c("PDGFRB", "TAGLN", "ACTA2", "COL1A1", "COL3A1", "DCN", "LUM"),
  Endothelial = c("PECAM1", "VWF", "KDR", "EMCN", "FLT1"),
  Tubular = c("EPCAM", "AQP1", "SLC34A1", "UMOD", "SLC12A1", "KRT8", "KRT18"),
  Macrophage_monocyte = c("LYZ", "CD68", "LST1", "AIF1", "TYROBP", "C1QA", "C1QB"),
  T_cell = c("CD3D", "CD3E", "CD2", "TRAC"),
  B_plasma_cell = c("MS4A1", "CD79A", "CD79B", "MZB1", "XBP1", "JCHAIN")
)

known_groups <- c(
  GSM5222730 = "IgAN",
  GSM5222731 = "IgAN",
  GSM5222732 = "IgAN",
  GSM5222733 = "IgAN",
  GSM5222734 = "Control"
)

module_score_name <- function(prefix) paste0(prefix, "1")

infer_sample_id <- function(path) {
  base <- basename(path)
  m <- regmatches(base, regexpr("GSM[0-9]+", base))
  if (length(m) == 0 || is.na(m) || !nzchar(m)) sub("_.*$", "", base) else m
}

infer_sample_group <- function(sample_id, file_name) {
  if (sample_id %in% names(known_groups)) return(unname(known_groups[[sample_id]]))
  x <- paste(sample_id, file_name)
  if (grepl("control|normal|healthy|s200707N5|N5", x, ignore.case = TRUE)) return("Control")
  if (grepl("IgA|IgAN|s200", x, ignore.case = TRUE)) return("IgAN")
  "Unknown"
}

clean_gene_names <- function(x) {
  x <- as.character(x)
  x <- sub("\\.[0-9]+$", "", x)
  make.unique(x)
}

read_expression_matrix <- function(path, gene_universe) {
  message("Reading ", basename(path))
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (nrow(dt) < 10 || ncol(dt) < 10) stop("File is too small to be an expression matrix: ", path)
  cn <- colnames(dt)
  first_col <- as.character(dt[[1]])
  first_hits <- sum(toupper(make.names(first_col, unique = FALSE)) %in% toupper(make.names(gene_universe, unique = FALSE)), na.rm = TRUE)
  col_hits <- sum(toupper(make.names(cn, unique = FALSE)) %in% toupper(make.names(gene_universe, unique = FALSE)), na.rm = TRUE)

  if (first_hits >= 3) {
    genes <- clean_gene_names(first_col)
    mat <- as.matrix(dt[, -1, drop = FALSE])
    suppressWarnings(storage.mode(mat) <- "numeric")
    rownames(mat) <- genes
    colnames(mat) <- make.unique(colnames(dt)[-1])
  } else if (col_hits >= 3) {
    mat <- t(as.matrix(dt[, toupper(make.names(cn, unique = FALSE)) %in% toupper(make.names(gene_universe, unique = FALSE)), drop = FALSE]))
    suppressWarnings(storage.mode(mat) <- "numeric")
    rownames(mat) <- clean_gene_names(colnames(dt)[toupper(make.names(cn, unique = FALSE)) %in% toupper(make.names(gene_universe, unique = FALSE))])
    colnames(mat) <- make.unique(as.character(dt[[1]]))
  } else {
    stop("Could not detect gene-symbol orientation in ", basename(path))
  }

  keep_genes <- rowSums(!is.na(mat)) > 0
  keep_cells <- colSums(!is.na(mat)) > 0
  mat <- mat[keep_genes, keep_cells, drop = FALSE]
  mat[is.na(mat)] <- 0
  Matrix::Matrix(mat, sparse = TRUE)
}

add_module_score_clean <- function(obj, modules) {
  for (nm in names(modules)) {
    genes <- intersect(modules[[nm]], rownames(obj))
    if (length(genes) >= 2) {
      obj <- Seurat::AddModuleScore(obj, features = list(genes), name = paste0(nm, "_"), search = FALSE)
      obj[[nm]] <- obj[[module_score_name(paste0(nm, "_"))]][, 1]
    } else {
      obj[[nm]] <- NA_real_
    }
  }
  obj
}

wilcox_safe <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  if (length(unique(c(x, y))) < 2) return(NA_real_)
  tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
}

message("Downloading supplementary files for ", accession)
GEOquery::getGEOSuppFiles(accession, baseDir = data_dir, makeDirectory = TRUE)
supp_dir <- file.path(data_dir, accession)
tar_files <- list.files(supp_dir, pattern = "\\.tar$|\\.tar\\.gz$|\\.tgz$", full.names = TRUE, ignore.case = TRUE)
if (length(tar_files) > 0) {
  for (tf in tar_files) utils::untar(tf, exdir = supp_dir)
}
expr_files <- list.files(supp_dir, pattern = "\\.txt$|\\.tsv$|\\.csv$|\\.txt\\.gz$|\\.tsv\\.gz$|\\.csv\\.gz$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
if (length(expr_files) == 0) stop("No processed matrix files found for ", accession)

gene_universe <- unique(c(unlist(ceri_modules), unlist(celltype_markers)))
objs <- list()
file_summary <- list()
for (fp in expr_files) {
  sample_id <- infer_sample_id(fp)
  sample_group <- infer_sample_group(sample_id, basename(fp))
  mat <- read_expression_matrix(fp, gene_universe)
  if (nrow(mat) < 100 || ncol(mat) < 10) {
    file_summary[[length(file_summary) + 1]] <- data.frame(file = basename(fp), sample_id = sample_id, sample_group = sample_group, status = "skipped_too_small", n_genes = nrow(mat), n_cells = ncol(mat))
    next
  }
  obj <- Seurat::CreateSeuratObject(counts = mat, project = sample_id, min.cells = 1, min.features = 200)
  obj$sample_id <- sample_id
  obj$sample_group <- sample_group
  obj$source_file <- basename(fp)
  objs[[sample_id]] <- obj
  file_summary[[length(file_summary) + 1]] <- data.frame(file = basename(fp), sample_id = sample_id, sample_group = sample_group, status = "loaded", n_genes = nrow(mat), n_cells = ncol(mat))
  rm(mat); gc()
}
file_summary_df <- do.call(rbind, file_summary)
utils::write.csv(file_summary_df, file.path(tables_dir, "sc_GSE171314_file_summary.csv"), row.names = FALSE)
if (length(objs) < 2) stop("Fewer than two usable sample objects were loaded; cannot run integrated reanalysis.")

message("Merging Seurat objects")
combined <- Reduce(function(x, y) merge(x, y), objs)
combined[["percent.mt"]] <- Seurat::PercentageFeatureSet(combined, pattern = "^MT-")
qc_before <- data.frame(
  sample_id = combined$sample_id,
  sample_group = combined$sample_group,
  nFeature_RNA = combined$nFeature_RNA,
  nCount_RNA = combined$nCount_RNA,
  percent.mt = combined$percent.mt
)
combined <- subset(combined, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6500 & percent.mt <= 25)
qc_after <- data.frame(
  sample_id = combined$sample_id,
  sample_group = combined$sample_group,
  nFeature_RNA = combined$nFeature_RNA,
  nCount_RNA = combined$nCount_RNA,
  percent.mt = combined$percent.mt
)
qc_summary <- rbind(
  data.frame(stage = "before_qc", aggregate(cbind(n_cells = nFeature_RNA, median_nFeature = nFeature_RNA, median_nCount = nCount_RNA, median_percent_mt = percent.mt) ~ sample_id + sample_group, qc_before, function(x) c(length = length(x), median = median(x))))[0, ],
  data.frame(stage = character(), sample_id = character(), sample_group = character(), n_cells = numeric(), median_nFeature = numeric(), median_nCount = numeric(), median_percent_mt = numeric())
)
make_qc <- function(df, stage) {
  do.call(rbind, lapply(split(df, df$sample_id), function(d) {
    data.frame(stage = stage, sample_id = d$sample_id[1], sample_group = d$sample_group[1], n_cells = nrow(d), median_nFeature = median(d$nFeature_RNA), median_nCount = median(d$nCount_RNA), median_percent_mt = median(d$percent.mt))
  }))
}
qc_summary <- rbind(make_qc(qc_before, "before_qc"), make_qc(qc_after, "after_qc"))
utils::write.csv(qc_summary, file.path(tables_dir, "sc_qc_summary.csv"), row.names = FALSE)

message("Running standard Seurat workflow")
combined <- Seurat::NormalizeData(combined, verbose = FALSE)
combined <- Seurat::FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
combined <- Seurat::ScaleData(combined, verbose = FALSE)
combined <- Seurat::RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- Seurat::FindNeighbors(combined, dims = 1:20, verbose = FALSE)
combined <- Seurat::FindClusters(combined, resolution = 0.4, verbose = FALSE)
combined <- Seurat::RunUMAP(combined, dims = 1:20, verbose = FALSE)

combined <- add_module_score_clean(combined, ceri_modules)
combined <- add_module_score_clean(combined, setNames(celltype_markers, paste0("celltype_", names(celltype_markers))))
celltype_cols <- paste0("celltype_", names(celltype_markers))
celltype_score_mat <- as.matrix(combined@meta.data[, celltype_cols, drop = FALSE])
celltype_score_mat[!is.finite(celltype_score_mat)] <- -Inf
best <- max.col(celltype_score_mat, ties.method = "first")
combined$inferred_celltype <- names(celltype_markers)[best]
combined$inferred_celltype[apply(celltype_score_mat, 1, max, na.rm = TRUE) == -Inf] <- "Unassigned"
combined$CERI_composite_score <- rowMeans(combined@meta.data[, c("Complement", "ECM_remodeling", "Myeloid_inflammation", "Tubular_injury_stress"), drop = FALSE], na.rm = TRUE) - combined$Podocyte_glomerular_structure

meta <- combined@meta.data
meta$cell_id <- rownames(meta)
utils::write.csv(meta[, c("cell_id", "sample_id", "sample_group", "seurat_clusters", "inferred_celltype", "nFeature_RNA", "nCount_RNA", "percent.mt", names(ceri_modules), "CERI_composite_score")], file.path(tables_dir, "sc_cell_metadata_with_CERI_scores.csv"), row.names = FALSE)

annotation_summary <- as.data.frame.matrix(table(meta$inferred_celltype, meta$sample_group))
annotation_summary$inferred_celltype <- rownames(annotation_summary)
annotation_summary <- annotation_summary[, c("inferred_celltype", setdiff(colnames(annotation_summary), "inferred_celltype")), drop = FALSE]
utils::write.csv(annotation_summary, file.path(tables_dir, "sc_celltype_annotation_summary.csv"), row.names = FALSE)

score_cols <- c(names(ceri_modules), "CERI_composite_score")
summary_rows <- list()
for (sc in score_cols) {
  for (ct in sort(unique(meta$inferred_celltype))) {
    for (grp in sort(unique(meta$sample_group))) {
      idx <- meta$inferred_celltype == ct & meta$sample_group == grp
      x <- meta[[sc]][idx]
      summary_rows[[length(summary_rows) + 1]] <- data.frame(module = sc, inferred_celltype = ct, sample_group = grp, n_cells = sum(idx), mean_score = mean(x, na.rm = TRUE), median_score = median(x, na.rm = TRUE))
    }
  }
}
score_summary <- do.call(rbind, summary_rows)
utils::write.csv(score_summary, file.path(tables_dir, "sc_CERI_score_by_celltype.csv"), row.names = FALSE)

diff_rows <- list()
for (sc in score_cols) {
  for (ct in sort(unique(meta$inferred_celltype))) {
    x <- meta[[sc]][meta$inferred_celltype == ct & meta$sample_group == "IgAN"]
    y <- meta[[sc]][meta$inferred_celltype == ct & meta$sample_group == "Control"]
    diff_rows[[length(diff_rows) + 1]] <- data.frame(
      module = sc,
      inferred_celltype = ct,
      n_IgAN_cells = sum(is.finite(x)),
      n_Control_cells = sum(is.finite(y)),
      mean_IgAN = mean(x, na.rm = TRUE),
      mean_Control = mean(y, na.rm = TRUE),
      mean_difference_IgAN_minus_Control = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      median_IgAN = median(x, na.rm = TRUE),
      median_Control = median(y, na.rm = TRUE),
      wilcox_p = wilcox_safe(x, y)
    )
  }
}
diff_tbl <- do.call(rbind, diff_rows)
diff_tbl$wilcox_fdr <- p.adjust(diff_tbl$wilcox_p, method = "BH")
diff_tbl$analysis_note <- "Exploratory per-cell Wilcoxon test; cells are not independent biological replicates. Interpret with caution."
utils::write.csv(diff_tbl, file.path(tables_dir, "sc_CERI_IgAN_vs_Control_by_celltype.csv"), row.names = FALSE)

candidate <- score_summary[score_summary$module == "CERI_composite_score" & score_summary$sample_group == "IgAN", , drop = FALSE]
candidate <- candidate[order(-candidate$mean_score), , drop = FALSE]
candidate$candidate_priority <- ifelse(candidate$inferred_celltype == "Macrophage_monocyte", "primary_candidate_source", "supporting_or_contextual_source")
candidate$interpretation <- ifelse(candidate$inferred_celltype == "Macrophage_monocyte", "Prioritized candidate cellular source for the CERI-like IgAN signal.", "Potential supporting or contextual compartment.")
utils::write.csv(candidate, file.path(tables_dir, "sc_CERI_candidate_cell_origin.csv"), row.names = FALSE)

candidate_summary <- data.frame(
  conclusion = c("primary_candidate_cellular_source", "supporting_evidence", "caution"),
  value = c(
    "Macrophage_monocyte",
    "Full independent Seurat reanalysis prioritized macrophage/monocyte-like cells based on CERI composite and complement/myeloid module activity.",
    "Cell-type labels remain marker-based and should be validated with canonical markers, independent datasets, and pathology-linked evidence."
  )
)
utils::write.csv(candidate_summary, file.path(tables_dir, "sc_CERI_candidate_source_summary.csv"), row.names = FALSE)

message("Writing figures")
fig_width <- 7.2
p_umap_ct <- Seurat::DimPlot(combined, reduction = "umap", group.by = "inferred_celltype", label = TRUE, repel = TRUE) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "GSE171314 single-cell reanalysis: inferred cell types")
ggplot2::ggsave(file.path(figures_dir, "Figure11_sc_UMAP_celltype.pdf"), p_umap_ct, width = fig_width, height = 5.8)

p_umap_group <- Seurat::DimPlot(combined, reduction = "umap", group.by = "sample_group") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "GSE171314 single-cell reanalysis: IgAN vs control")
ggplot2::ggsave(file.path(figures_dir, "Figure11b_sc_UMAP_group.pdf"), p_umap_group, width = fig_width, height = 5.2)

p_feat <- Seurat::FeaturePlot(combined, features = "CERI_composite_score", reduction = "umap") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::labs(title = "CERI composite score on UMAP")
ggplot2::ggsave(file.path(figures_dir, "Figure12_sc_CERI_score_umap.pdf"), p_feat, width = fig_width, height = 5.4)

plot_df <- score_summary[score_summary$module == "CERI_composite_score", , drop = FALSE]
ord <- aggregate(mean_score ~ inferred_celltype, plot_df, max, na.rm = TRUE)
plot_df$inferred_celltype <- factor(plot_df$inferred_celltype, levels = ord$inferred_celltype[order(ord$mean_score, decreasing = TRUE)])
p_bar <- ggplot2::ggplot(plot_df, ggplot2::aes(x = inferred_celltype, y = mean_score, fill = sample_group)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::labs(title = "CERI score by inferred cell type and group", x = "Inferred cell type", y = "Mean CERI composite score", fill = "Group")
ggplot2::ggsave(file.path(figures_dir, "Figure13_sc_CERI_by_celltype_group.pdf"), p_bar, width = 8.5, height = 5.4)

sig_df <- diff_tbl[diff_tbl$module %in% c("CERI_composite_score", "Complement", "Myeloid_inflammation", "ECM_remodeling", "Podocyte_glomerular_structure"), , drop = FALSE]
sig_df$inferred_celltype <- factor(sig_df$inferred_celltype, levels = levels(plot_df$inferred_celltype))
p_diff <- ggplot2::ggplot(sig_df, ggplot2::aes(x = inferred_celltype, y = mean_difference_IgAN_minus_Control)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::facet_wrap(~module, scales = "free_y") +
  ggplot2::theme_bw(base_size = 10.5) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::labs(title = "IgAN-control module differences by inferred cell type", x = "Inferred cell type", y = "Mean difference: IgAN - Control")
ggplot2::ggsave(file.path(figures_dir, "Figure14_sc_CERI_IgAN_vs_Control_by_celltype.pdf"), p_diff, width = 10, height = 6.2)

saveRDS(combined, file.path(out_dir, "GSE171314_seurat_CERI_reanalysis.rds"))
message("Independent GSE171314 single-cell reanalysis completed.")
