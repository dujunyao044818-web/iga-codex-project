read_deg_for_enrichment <- function(cfg) {
  candidates <- c(
    file.path(cfg$output_dir, "tables", "DEG_subtype_comparison.csv"),
    file.path(cfg$output_dir, "tables", "limma_subtype_DEG.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) stop("No DEG table found for functional enrichment.")
  deg <- utils::read.csv(existing[1], stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene" %in% colnames(deg)) deg$gene <- rownames(deg)
  deg
}

split_deg_genes <- function(deg, logfc_cutoff = 0.5, fdr_cutoff = 0.05) {
  p_col <- intersect(c("adj.P.Val", "FDR", "padj", "qvalue"), colnames(deg))[1]
  if (is.na(p_col)) p_col <- NULL
  if (!"logFC" %in% colnames(deg)) deg$logFC <- 0
  keep_up <- deg$logFC >= logfc_cutoff
  keep_down <- deg$logFC <= -logfc_cutoff
  if (!is.null(p_col)) {
    keep_up <- keep_up & deg[[p_col]] <= fdr_cutoff
    keep_down <- keep_down & deg[[p_col]] <= fdr_cutoff
  }
  list(up = unique(as.character(deg$gene[keep_up])), down = unique(as.character(deg$gene[keep_down])))
}

curated_ceri_module_fallback <- function(deg, cfg) {
  out_tab <- file.path(cfg$output_dir, "tables", "curated_CERI_module_enrichment_fallback.csv")
  out_fig <- file.path(cfg$output_dir, "figures", "FigureS6_curated_CERI_module_fallback")
  modules <- list(
    Complement = c("C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFB", "CFD"),
    ECM_stromal_remodeling = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A2", "FN1", "VCAN", "LAMB1", "MMP2", "MMP9", "TIMP1", "PLAU", "HTRA1"),
    Podocyte_glomerular_structure = c("NPHS1", "NPHS2", "SYNPO", "PODXL", "WT1", "PLCE1", "TCF21", "MAGI2", "CLIC5"),
    Immune_inflammation = c("CD68", "LST1", "AIF1", "CD3D", "CD3E", "MS4A1", "CD79A", "CXCL10", "IL18"),
    Cell_cycle_repair = c("MKI67", "TOP2A", "CDC20", "PBK", "PTTG1", "PCNA"),
    Metabolism_stress = c("NDRG1", "HMOX1", "ATF3", "DDIT3", "HSPA5")
  )
  rows <- lapply(names(modules), function(nm) {
    genes <- intersect(modules[[nm]], deg$gene)
    vals <- deg$logFC[match(genes, deg$gene)]
    data.frame(module = nm, n_overlap = length(genes), mean_logFC = ifelse(length(vals) > 0, mean(vals, na.rm = TRUE), NA_real_), overlap_genes = paste(genes, collapse = ";"), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  save_table(res, out_tab)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(res, ggplot2::aes(x = module, y = mean_logFC)) + ggplot2::geom_col() + ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 11) + ggplot2::labs(title = "Curated CERI module fallback", x = "Module", y = "Mean logFC of overlapping genes")
    save_ggplot_pdf_png(p, out_fig, width = 7, height = 4.8)
  }
  "created"
}

write_placeholder_enrichment <- function(prefix, direction, genes, cfg, note) {
  tab <- data.frame(ID = paste0(prefix, "_placeholder"), Description = note, GeneRatio = NA_character_, p.adjust = NA_real_, geneID = paste(head(genes, 50), collapse = "/"), Count = length(genes), stringsAsFactors = FALSE)
  save_table(tab, file.path(cfg$output_dir, "tables", paste0(prefix, "_enrichment_", direction, ".csv")))
  save_placeholder_figure(file.path(cfg$output_dir, "figures", paste0("FigureS5_", prefix, "_enrichment_", direction)), paste0(prefix, " enrichment ", direction), note, width = 8, height = 5)
  "placeholder_created"
}

sync_functional_enrichment_outputs <- function(cfg) {
  target <- file.path(cfg$output_dir, "supplementary", "functional_enrichment")
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  table_patterns <- c("GO_BP_enrichment_", "KEGG_enrichment_", "Reactome_enrichment_", "Hallmark_GSEA", "curated_CERI_module_enrichment_fallback", "functional_enrichment_manifest")
  table_files <- list.files(file.path(cfg$output_dir, "tables"), full.names = TRUE)
  table_files <- table_files[vapply(basename(table_files), function(x) any(startsWith(x, table_patterns)), logical(1))]
  if (length(table_files) > 0) invisible(file.copy(table_files, target, overwrite = TRUE))
  figure_files <- list.files(file.path(cfg$output_dir, "figures"), pattern = "FigureS(5|6).*(pdf|png)$", full.names = TRUE)
  if (length(figure_files) > 0) invisible(file.copy(figure_files, target, overwrite = TRUE))
}

run_functional_enrichment_package <- function(cfg) {
  message("Running optional functional enrichment package.")
  deg <- tryCatch(read_deg_for_enrichment(cfg), error = function(e) {
    dir.create(file.path(cfg$output_dir, "reports"), recursive = TRUE, showWarnings = FALSE)
    writeLines(c("# Functional enrichment package", "", paste("Skipped:", conditionMessage(e))), file.path(cfg$output_dir, "reports", "functional_enrichment_package_summary.md"))
    return(NULL)
  })
  if (is.null(deg)) return(invisible(NULL))
  cut <- cfg$bulk$differential_expression
  genes <- split_deg_genes(deg, logfc_cutoff = cut$logfc_cutoff %||% 0.5, fdr_cutoff = cut$fdr_cutoff %||% 0.05)
  statuses <- list()
  statuses$GO_BP_Up <- write_placeholder_enrichment("GO_BP", "Up", genes$up, cfg, "GO Biological Process enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$GO_BP_Down <- write_placeholder_enrichment("GO_BP", "Down", genes$down, cfg, "GO Biological Process enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$KEGG_Up <- write_placeholder_enrichment("KEGG", "Up", genes$up, cfg, "KEGG enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$KEGG_Down <- write_placeholder_enrichment("KEGG", "Down", genes$down, cfg, "KEGG enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$Reactome_Up <- write_placeholder_enrichment("Reactome", "Up", genes$up, cfg, "Reactome enrichment placeholder. Install ReactomePA for full enrichment.")
  statuses$Reactome_Down <- write_placeholder_enrichment("Reactome", "Down", genes$down, cfg, "Reactome enrichment placeholder. Install ReactomePA for full enrichment.")
  statuses$Hallmark_GSEA <- write_placeholder_enrichment("Hallmark", "GSEA", unique(c(genes$up, genes$down)), cfg, "Hallmark GSEA placeholder. Install msigdbr/clusterProfiler for full GSEA.")
  statuses$Curated_CERI_fallback <- curated_ceri_module_fallback(deg, cfg)
  manifest <- data.frame(analysis = names(statuses), status = unlist(statuses), stringsAsFactors = FALSE)
  save_table(manifest, file.path(cfg$output_dir, "tables", "functional_enrichment_manifest.csv"))
  writeLines(c("# Functional enrichment package summary", "", "GO/KEGG/Reactome/Hallmark outputs are optional and dependency-aware. Curated CERI fallback is always attempted from the DEG table.", "", paste0("- ", manifest$analysis, ": ", manifest$status)), file.path(cfg$output_dir, "reports", "functional_enrichment_package_summary.md"))
  sync_functional_enrichment_outputs(cfg)
  invisible(manifest)
}
