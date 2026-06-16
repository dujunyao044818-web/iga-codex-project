read_deg_for_enrichment <- function(cfg) {
  candidates <- c(
    file.path(cfg$output_dir, "tables", "DEG_subtype_comparison.csv"),
    file.path(cfg$output_dir, "tables", "limma_subtype_DEG.csv")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) stop("No DEG table found for functional enrichment.")
  deg <- utils::read.csv(existing[1], stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene" %in% colnames(deg)) deg$gene <- rownames(deg)
  if (!"logFC" %in% colnames(deg)) deg$logFC <- 0
  if (!"adj.P.Val" %in% colnames(deg)) {
    p_col <- intersect(c("FDR", "padj", "qvalue", "P.Value", "pvalue"), colnames(deg))[1]
    deg$adj.P.Val <- if (!is.na(p_col)) deg[[p_col]] else 1
  }
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

plot_enrichment_dotbar <- function(df, stem, title) {
  if (is.null(df) || nrow(df) == 0 || !"Description" %in% colnames(df)) {
    save_placeholder_figure(stem, title, "No enriched terms passed filtering.", width = 8.5, height = 5.8)
    return(invisible(FALSE))
  }
  df <- head(df[order(df$p.adjust), , drop = FALSE], 20)
  df$minus_log10_FDR <- -log10(pmax(df$p.adjust, .Machine$double.xmin))
  if (!"Count" %in% colnames(df)) df$Count <- seq_len(nrow(df))
  p <- ggplot2::ggplot(df, ggplot2::aes(stats::reorder(Description, minus_log10_FDR), minus_log10_FDR, size = Count, color = p.adjust)) +
    ggplot2::geom_point(alpha = 0.9) + ggplot2::coord_flip() + ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = title, x = NULL, y = "-log10(FDR)", color = "FDR", size = "Count")
  save_ggplot_pdf_png(p, stem, width = 8.5, height = 5.8)
  invisible(TRUE)
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

run_real_enrichment_outputs <- function(deg, cfg) {
  outdir <- file.path(cfg$output_dir, "supplementary", "functional_enrichment")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  statuses <- list()
  write_real_enrichment_skip <- function(stem, table_name, title, reason) {
    save_table(data.frame(status = "skipped", reason = reason, stringsAsFactors = FALSE), file.path(cfg$output_dir, "tables", table_name))
    save_placeholder_figure(file.path(outdir, stem), title, reason, width = 8.5, height = 5.8)
  }
  sig <- unique(deg$gene[deg$adj.P.Val <= 0.05 & is.finite(deg$logFC)])
  if (length(sig) < 5) {
    for (x in c("GO_BP", "KEGG", "Reactome")) write_real_enrichment_skip(paste0(x, "_dotplot"), paste0(x, "_enrichment_real.csv"), paste(x, "enrichment"), "Skipped: fewer than 5 significant genes.")
    statuses$GO_BP_real <- "skipped: fewer than 5 significant genes"
    statuses$KEGG_real <- "skipped: fewer than 5 significant genes"
    statuses$Reactome_real <- "skipped: fewer than 5 significant genes"
  } else if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    write_real_enrichment_skip("GO_BP_dotplot", "GO_BP_enrichment_real.csv", "GO biological process enrichment", "Skipped: clusterProfiler/org.Hs.eg.db unavailable.")
    write_real_enrichment_skip("KEGG_dotplot", "KEGG_enrichment_real.csv", "KEGG enrichment", "Skipped: clusterProfiler/org.Hs.eg.db unavailable.")
    write_real_enrichment_skip("Reactome_dotplot", "Reactome_enrichment_real.csv", "Reactome pathway enrichment", "Skipped: clusterProfiler/org.Hs.eg.db unavailable.")
    statuses$GO_BP_real <- "skipped: clusterProfiler/org.Hs.eg.db unavailable"
    statuses$KEGG_real <- "skipped: clusterProfiler/org.Hs.eg.db unavailable"
    statuses$Reactome_real <- "skipped: clusterProfiler/org.Hs.eg.db unavailable"
  } else {
    entrez <- tryCatch(clusterProfiler::bitr(sig, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db), error = function(e) NULL)
    if (is.null(entrez) || !nrow(entrez)) {
      write_real_enrichment_skip("GO_BP_dotplot", "GO_BP_enrichment_real.csv", "GO biological process enrichment", "Skipped: SYMBOL to ENTREZ mapping failed.")
      write_real_enrichment_skip("KEGG_dotplot", "KEGG_enrichment_real.csv", "KEGG enrichment", "Skipped: SYMBOL to ENTREZ mapping failed.")
      write_real_enrichment_skip("Reactome_dotplot", "Reactome_enrichment_real.csv", "Reactome pathway enrichment", "Skipped: SYMBOL to ENTREZ mapping failed.")
      statuses$GO_BP_real <- "skipped: SYMBOL to ENTREZ mapping failed"
      statuses$KEGG_real <- "skipped: SYMBOL to ENTREZ mapping failed"
      statuses$Reactome_real <- "skipped: SYMBOL to ENTREZ mapping failed"
    } else {
      ids <- unique(entrez$ENTREZID)
      go <- tryCatch(clusterProfiler::enrichGO(gene = ids, OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", readable = TRUE), error = function(e) NULL)
      go_df <- tryCatch(as.data.frame(go), error = function(e) data.frame())
      save_table(go_df, file.path(cfg$output_dir, "tables", "GO_BP_enrichment_real.csv"))
      plot_enrichment_dotbar(go_df, file.path(outdir, "GO_BP_dotplot"), "GO biological process enrichment")
      statuses$GO_BP_real <- if (nrow(go_df)) "created" else "created empty table; no terms passed filtering"

      kegg <- tryCatch(clusterProfiler::enrichKEGG(gene = ids, organism = "hsa", pAdjustMethod = "BH"), error = function(e) NULL)
      kegg_df <- tryCatch(as.data.frame(kegg), error = function(e) data.frame())
      save_table(kegg_df, file.path(cfg$output_dir, "tables", "KEGG_enrichment_real.csv"))
      plot_enrichment_dotbar(kegg_df, file.path(outdir, "KEGG_dotplot"), "KEGG enrichment")
      statuses$KEGG_real <- if (nrow(kegg_df)) "created" else "created empty table; no terms passed filtering"

      if (requireNamespace("ReactomePA", quietly = TRUE)) {
        react <- tryCatch(ReactomePA::enrichPathway(gene = ids, organism = "human", pAdjustMethod = "BH", readable = TRUE), error = function(e) NULL)
        react_df <- tryCatch(as.data.frame(react), error = function(e) data.frame())
        save_table(react_df, file.path(cfg$output_dir, "tables", "Reactome_enrichment_real.csv"))
        plot_enrichment_dotbar(react_df, file.path(outdir, "Reactome_dotplot"), "Reactome pathway enrichment")
        statuses$Reactome_real <- if (nrow(react_df)) "created" else "created empty table; no terms passed filtering"
      } else {
        save_table(data.frame(status = "skipped", reason = "ReactomePA unavailable"), file.path(cfg$output_dir, "tables", "Reactome_enrichment_real.csv"))
        save_placeholder_figure(file.path(outdir, "Reactome_dotplot"), "Reactome pathway enrichment", "Skipped because ReactomePA is unavailable.", width = 8.5, height = 5.8)
        statuses$Reactome_real <- "skipped: ReactomePA unavailable"
      }
    }
  }
  statuses$Hallmark_GSEA_real <- run_hallmark_gsea_real(deg, cfg)
  statuses
}

run_hallmark_gsea_real <- function(deg, cfg) {
  outdir <- file.path(cfg$output_dir, "supplementary", "functional_enrichment")
  if (!requireNamespace("msigdbr", quietly = TRUE) || !requireNamespace("clusterProfiler", quietly = TRUE)) {
    save_table(data.frame(status = "skipped", reason = "msigdbr or clusterProfiler unavailable"), file.path(cfg$output_dir, "tables", "Hallmark_GSEA_real.csv"))
    save_placeholder_figure(file.path(outdir, "Hallmark_GSEA_dotplot"), "Hallmark GSEA", "Skipped because msigdbr or clusterProfiler is unavailable.", width = 8.5, height = 5.8)
    return("skipped: msigdbr or clusterProfiler unavailable")
  }
  ranks <- deg$logFC
  names(ranks) <- deg$gene
  ok <- is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))
  ranks <- sort(tapply(ranks[ok], names(ranks)[ok], mean), decreasing = TRUE)
  hallmark <- tryCatch(msigdbr::msigdbr(species = "Homo sapiens", category = "H"), error = function(e) NULL)
  if (is.null(hallmark) || !nrow(hallmark)) {
    save_table(data.frame(status = "skipped", reason = "Hallmark gene sets unavailable"), file.path(cfg$output_dir, "tables", "Hallmark_GSEA_real.csv"))
    save_placeholder_figure(file.path(outdir, "Hallmark_GSEA_dotplot"), "Hallmark GSEA", "Skipped because Hallmark gene sets are unavailable.", width = 8.5, height = 5.8)
    return("skipped: Hallmark gene sets unavailable")
  }
  term2gene <- unique(hallmark[, c("gs_name", "gene_symbol")])
  gsea <- tryCatch(clusterProfiler::GSEA(ranks, TERM2GENE = term2gene, pAdjustMethod = "BH", verbose = FALSE), error = function(e) NULL)
  gsea_df <- tryCatch(as.data.frame(gsea), error = function(e) data.frame())
  save_table(gsea_df, file.path(cfg$output_dir, "tables", "Hallmark_GSEA_real.csv"))
  if (!nrow(gsea_df)) {
    save_placeholder_figure(file.path(outdir, "Hallmark_GSEA_dotplot"), "Hallmark GSEA", "No Hallmark terms passed filtering.", width = 8.5, height = 5.8)
    return("created empty table; no Hallmark terms")
  }
  plot_df <- head(gsea_df[order(gsea_df$p.adjust), , drop = FALSE], 20)
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(stats::reorder(Description, NES), NES, color = p.adjust, size = abs(NES))) +
    ggplot2::geom_point(alpha = 0.9) + ggplot2::coord_flip() +
    ggplot2::labs(title = "Hallmark GSEA", x = NULL, y = "Normalized enrichment score", color = "FDR", size = "|NES|") +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
  save_ggplot_pdf_png(p, file.path(outdir, "Hallmark_GSEA_dotplot"), width = 8.5, height = 5.8)
  "created"
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
  logfc_cutoff <- if (!is.null(cut$logfc_cutoff)) cut$logfc_cutoff else 0.5
  fdr_cutoff <- if (!is.null(cut$fdr_cutoff)) cut$fdr_cutoff else 0.05
  genes <- split_deg_genes(deg, logfc_cutoff = logfc_cutoff, fdr_cutoff = fdr_cutoff)
  statuses <- list()
  statuses$GO_BP_Up <- write_placeholder_enrichment("GO_BP", "Up", genes$up, cfg, "GO Biological Process enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$GO_BP_Down <- write_placeholder_enrichment("GO_BP", "Down", genes$down, cfg, "GO Biological Process enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$KEGG_Up <- write_placeholder_enrichment("KEGG", "Up", genes$up, cfg, "KEGG enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$KEGG_Down <- write_placeholder_enrichment("KEGG", "Down", genes$down, cfg, "KEGG enrichment placeholder. Install clusterProfiler/org.Hs.eg.db for full enrichment.")
  statuses$Reactome_Up <- write_placeholder_enrichment("Reactome", "Up", genes$up, cfg, "Reactome enrichment placeholder. Install ReactomePA for full enrichment.")
  statuses$Reactome_Down <- write_placeholder_enrichment("Reactome", "Down", genes$down, cfg, "Reactome enrichment placeholder. Install ReactomePA for full enrichment.")
  statuses$Hallmark_GSEA <- write_placeholder_enrichment("Hallmark", "GSEA", unique(c(genes$up, genes$down)), cfg, "Hallmark GSEA placeholder. Install msigdbr/clusterProfiler for full GSEA.")
  statuses$Real_enrichment_requested_outputs <- paste(unlist(run_real_enrichment_outputs(deg, cfg)), collapse = "; ")
  statuses$Curated_CERI_fallback <- curated_ceri_module_fallback(deg, cfg)
  manifest <- data.frame(analysis = names(statuses), status = unlist(statuses), stringsAsFactors = FALSE)
  save_table(manifest, file.path(cfg$output_dir, "tables", "functional_enrichment_manifest.csv"))
  writeLines(c("# Functional enrichment package summary", "", "GO/KEGG/Reactome/Hallmark outputs are optional and dependency-aware. Curated CERI fallback is always attempted from the DEG table.", "", paste0("- ", manifest$analysis, ": ", manifest$status)), file.path(cfg$output_dir, "reports", "functional_enrichment_package_summary.md"))
  sync_functional_enrichment_outputs(cfg)
  invisible(manifest)
}
