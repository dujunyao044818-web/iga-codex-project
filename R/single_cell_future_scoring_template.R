# Future-ready CERI single-cell scoring template
#
# This file is intentionally NOT sourced by run_pipeline.R and is NOT used by the
# GitHub Actions pipeline. It is a template for future validation once a processed
# IgAN sc/snRNA-seq object is available.
#
# Requires Seurat; not used in GitHub Actions pipeline.

# compute_ceri_module_scores_seurat <- function(seurat_obj) {
#   if (!requireNamespace("Seurat", quietly = TRUE)) {
#     stop("Seurat is required for this future template but is not required by the main pipeline.")
#   }
#
#   gene_sets <- list(
#     Complement = c("C3", "C4A", "C4B", "C1QA", "C1QB", "C1QC", "CFB"),
#     ECM_remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL5A2", "FN1", "VCAN", "LAMB1", "PLAU", "HTRA1"),
#     Podocyte_glomerular_structure = c("NPHS1", "NPHS2", "SYNPO", "WT1", "PODXL", "PLCE1", "TCF21", "MAGI2", "CLIC5"),
#     Myeloid_inflammation = c("CD68", "LST1", "AIF1", "C1QA", "C1QB", "TYROBP", "IL18", "ALOX5AP"),
#     T_B_plasma_cell = c("CD3D", "CD3E", "CD2", "MS4A1", "CD79A", "CD79B", "MZB1", "XBP1", "JCHAIN"),
#     Tubular_injury = c("HAVCR1", "LCN2", "KRT8", "KRT18", "VCAM1", "NDRG1", "HMOX1")
#   )
#
#   for (nm in names(gene_sets)) {
#     genes <- intersect(gene_sets[[nm]], rownames(seurat_obj))
#     if (length(genes) >= 2) {
#       seurat_obj <- Seurat::AddModuleScore(
#         object = seurat_obj,
#         features = list(genes),
#         name = paste0("CERI_", nm, "_")
#       )
#     }
#   }
#   seurat_obj
# }

# summarize_ceri_by_celltype <- function(seurat_obj, celltype_col, disease_col) {
#   meta <- seurat_obj@meta.data
#   score_cols <- grep("^CERI_", colnames(meta), value = TRUE)
#   if (length(score_cols) == 0) stop("No CERI module score columns found.")
#
#   out <- list()
#   for (score_col in score_cols) {
#     for (ct in unique(meta[[celltype_col]])) {
#       for (grp in unique(meta[[disease_col]])) {
#         idx <- meta[[celltype_col]] == ct & meta[[disease_col]] == grp
#         x <- meta[[score_col]][idx]
#         out[[length(out) + 1]] <- data.frame(
#           module = score_col,
#           celltype = ct,
#           group = grp,
#           n_cells = sum(idx, na.rm = TRUE),
#           mean_score = mean(x, na.rm = TRUE),
#           median_score = stats::median(x, na.rm = TRUE),
#           stringsAsFactors = FALSE
#         )
#       }
#     }
#   }
#   do.call(rbind, out)
# }

# plot_ceri_by_celltype <- function(seurat_obj, output_file, celltype_col = "celltype") {
#   if (!requireNamespace("Seurat", quietly = TRUE)) stop("Seurat is required.")
#   if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
#
#   meta <- seurat_obj@meta.data
#   score_cols <- grep("^CERI_", colnames(meta), value = TRUE)
#   if (length(score_cols) == 0) stop("No CERI module score columns found.")
#
#   long <- do.call(rbind, lapply(score_cols, function(sc) {
#     data.frame(
#       celltype = meta[[celltype_col]],
#       module = sc,
#       score = meta[[sc]],
#       stringsAsFactors = FALSE
#     )
#   }))
#
#   p <- ggplot2::ggplot(long, ggplot2::aes(x = celltype, y = score)) +
#     ggplot2::geom_boxplot(outlier.shape = NA) +
#     ggplot2::facet_wrap(~module, scales = "free_y") +
#     ggplot2::theme_bw(base_size = 11) +
#     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
#     ggplot2::labs(x = "Cell type", y = "CERI module score")
#   ggplot2::ggsave(output_file, p, width = 10, height = 6)
# }
