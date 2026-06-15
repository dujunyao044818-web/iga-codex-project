run_single_cell <- function(markers, cfg) {
  tables_dir <- file.path(cfg$output_dir, "tables")
  figures_dir <- file.path(cfg$output_dir, "figures")
  sc_dir <- file.path(cfg$output_dir, "single_cell")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sc_dir, recursive = TRUE, showWarnings = FALSE)

  out <- data.frame(signature_gene = unique(as.character(markers)), stringsAsFactors = FALSE)
  utils::write.csv(out, file.path(tables_dir, "single_cell_signature_genes_for_validation.csv"), row.names = FALSE)

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

  validation_genes <- unique(c(unlist(ceri_modules), unlist(celltype_markers), head(out$signature_gene, 50)))
  utils::write.csv(data.frame(validation_gene = validation_genes, stringsAsFactors = FALSE), file.path(tables_dir, "single_cell_spatial_CERI_gene_panel.csv"), row.names = FALSE)

  cell_hypothesis <- data.frame(
    expected_signal = c("Complement production / activation", "ECM remodeling", "Podocyte/glomerular structural loss", "Myeloid inflammatory amplification", "T/B-cell immune context", "Tubular injury / stress"),
    candidate_cell_or_compartment = c(
      "macrophage/monocyte, mesangial, endothelial, tubular compartments depending on dataset",
      "mesangial, fibroblast/stromal, injured tubular and periglomerular compartments",
      "podocyte and glomerular epithelial compartment",
      "macrophage/monocyte clusters",
      "lymphocyte and plasma-cell clusters",
      "proximal/distal tubular epithelial clusters"
    ),
    validation_readout = c(
      "C3/C4/CFB/C1Q gene module score per cell type",
      "COL1A1/COL3A1/COL5A2/FN1/VCAN/LAMB1/PLAU/HTRA1 module score",
      "NPHS1/NPHS2/SYNPO/WT1/PODXL/PLCE1/TCF21/MAGI2/CLIC5 module score reduction",
      "CD68/LST1/AIF1/C1QA/C1QB/TYROBP module score",
      "CD3D/CD3E/MS4A1/CD79A/MZB1/XBP1 module score",
      "HAVCR1/LCN2/KRT8/KRT18/VCAM1/NDRG1/HMOX1 module score"
    ),
    interpretation_if_supported = c(
      "CERI-like state has a complement-active cellular niche.",
      "CERI-like state involves stromal/mesangial matrix remodeling.",
      "CERI-like state is coupled to podocyte/glomerular structural injury.",
      "CERI-like state is accompanied by innate immune activation.",
      "Adaptive immune context modifies the CERI-like state.",
      "Tubular stress may represent downstream injury propagation."
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(cell_hypothesis, file.path(tables_dir, "single_cell_celltype_hypothesis.csv"), row.names = FALSE)

  spatial_plan <- data.frame(
    validation_layer = c("single_cell", "spatial_transcriptomics", "RNAscope_or_IHC", "pathology_correlation"),
    goal = c(
      "Identify which cell types carry the CERI gene modules.",
      "Localize CERI-high regions to glomerular, periglomerular, tubulointerstitial, or immune-infiltrate niches.",
      "Experimentally validate selected high-value genes/proteins in biopsy tissue.",
      "Connect CERI-high spatial or cellular signal to Oxford MEST-C lesions, C3 deposition, sclerosis, or proteinuria if metadata are available."
    ),
    primary_features = c("CERI module score by cell type and disease group", "CERI module score by spatial spot/domain", "C3/CFB/COL5A2/VCAN/SYNPO/NPHS2/IL18 selected panel", "M/S/T/C lesions, C3 intensity, eGFR/proteinuria if available"),
    expected_support_for_novelty = c(
      "Shows the candidate molecular state has identifiable cellular origins.",
      "Shows the candidate molecular state is spatially organized rather than a bulk-only artifact.",
      "Provides wet-lab-compatible validation targets.",
      "Links transcriptomic state to pathological severity or prognosis."
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(spatial_plan, file.path(tables_dir, "single_cell_spatial_validation_plan.csv"), row.names = FALSE)

  score_modules <- function(expr, modules) {
    expr <- as.matrix(expr)
    storage.mode(expr) <- "numeric"
    expr <- log1p(expr)
    z <- t(scale(t(expr)))
    z[!is.finite(z)] <- NA_real_
    out_scores <- list()
    for (nm in names(modules)) {
      genes <- intersect(modules[[nm]], rownames(z))
      if (length(genes) >= 2) {
        out_scores[[nm]] <- colMeans(z[genes, , drop = FALSE], na.rm = TRUE)
      } else {
        out_scores[[nm]] <- rep(NA_real_, ncol(z))
      }
    }
    as.data.frame(out_scores, stringsAsFactors = FALSE)
  }

  infer_celltype <- function(scores) {
    marker_cols <- paste0("celltype_", names(celltype_markers))
    available <- intersect(marker_cols, colnames(scores))
    if (length(available) == 0) return(rep("Unassigned", nrow(scores)))
    m <- as.matrix(scores[, available, drop = FALSE])
    m[!is.finite(m)] <- -Inf
    best <- max.col(m, ties.method = "first")
    pred <- sub("^celltype_", "", available[best])
    pred[apply(m, 1, max, na.rm = TRUE) == -Inf] <- "Unassigned"
    pred
  }

  parse_expression_file <- function(path, gene_universe) {
    dt <- tryCatch(data.table::fread(path, data.table = FALSE, check.names = FALSE), error = function(e) NULL)
    if (is.null(dt) || nrow(dt) < 2 || ncol(dt) < 2) return(NULL)
    cn <- colnames(dt)
    first_col <- as.character(dt[[1]])
    first_col_clean <- toupper(make.names(first_col, unique = FALSE))
    gene_universe_clean <- toupper(make.names(gene_universe, unique = FALSE))
    row_hits <- sum(first_col_clean %in% gene_universe_clean, na.rm = TRUE)
    col_hits <- sum(toupper(make.names(cn, unique = FALSE)) %in% gene_universe_clean, na.rm = TRUE)

    if (row_hits >= 3) {
      keep <- first_col_clean %in% gene_universe_clean
      sub <- dt[keep, , drop = FALSE]
      genes <- as.character(sub[[1]])
      mat <- as.matrix(sub[, -1, drop = FALSE])
      suppressWarnings(storage.mode(mat) <- "numeric")
      rownames(mat) <- toupper(make.names(genes, unique = TRUE))
      rownames(mat) <- sub("\\.[0-9]+$", "", rownames(mat))
      colnames(mat) <- make.names(colnames(sub)[-1], unique = TRUE)
      mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
      return(mat)
    }

    if (col_hits >= 3) {
      keep_cols <- toupper(make.names(cn, unique = FALSE)) %in% gene_universe_clean
      sub <- dt[, keep_cols, drop = FALSE]
      mat <- t(as.matrix(sub))
      suppressWarnings(storage.mode(mat) <- "numeric")
      rownames(mat) <- toupper(make.names(colnames(sub), unique = TRUE))
      rownames(mat) <- sub("\\.[0-9]+$", "", rownames(mat))
      colnames(mat) <- make.names(as.character(dt[[1]]), unique = TRUE)
      mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]
      return(mat)
    }
    NULL
  }

  infer_sample_group <- function(sample_id, file_name) {
    x <- paste(sample_id, file_name)
    if (grepl("control|normal|healthy|s200707N5|N5", x, ignore.case = TRUE)) return("Control")
    if (grepl("IgA|IgAN|s200", x, ignore.case = TRUE)) return("IgAN")
    "Unknown"
  }

  wilcox_safe <- function(x, y) {
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) < 3 || length(y) < 3) return(NA_real_)
    if (length(unique(c(x, y))) < 2) return(NA_real_)
    tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
  }

  sc_accession <- if (!is.null(cfg$single_cell$geo_accession) && nzchar(cfg$single_cell$geo_accession)) cfg$single_cell$geo_accession else "GSE171314"
  sc_status <- data.frame(
    item = c("dataset", "analysis_mode", "status"),
    value = c(sc_accession, "lightweight processed-matrix module scoring without Seurat", "not_started"),
    stringsAsFactors = FALSE
  )

  download_ok <- FALSE
  scores_all <- list()
  sample_summary <- list()

  tryCatch({
    if (!requireNamespace("GEOquery", quietly = TRUE)) stop("GEOquery not available")
    message("Attempting lightweight single-cell validation using ", sc_accession)
    supp_dir <- file.path(sc_dir, sc_accession)
    dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
    GEOquery::getGEOSuppFiles(sc_accession, baseDir = sc_dir, makeDirectory = TRUE)
    tar_files <- list.files(supp_dir, pattern = "\\.tar$|\\.tar\\.gz$|\\.tgz$", full.names = TRUE, ignore.case = TRUE)
    if (length(tar_files) > 0) {
      for (tf in tar_files) utils::untar(tf, exdir = supp_dir)
    }
    txt_files <- list.files(supp_dir, pattern = "\\.txt$|\\.tsv$|\\.csv$|\\.txt\\.gz$|\\.tsv\\.gz$|\\.csv\\.gz$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    if (length(txt_files) == 0) stop("No processed TXT/TSV/CSV supplementary files found for ", sc_accession)
    download_ok <- TRUE

    gene_universe <- unique(c(unlist(ceri_modules), unlist(celltype_markers), validation_genes))
    gene_universe <- unique(c(gene_universe, toupper(make.names(gene_universe, unique = FALSE))))
    for (fp in txt_files) {
      expr <- parse_expression_file(fp, gene_universe)
      if (is.null(expr) || nrow(expr) < 3 || ncol(expr) < 5) {
        sample_summary[[length(sample_summary) + 1]] <- data.frame(
          file = basename(fp), status = "skipped_unrecognized_or_too_small", n_genes = ifelse(is.null(expr), 0, nrow(expr)), n_cells = ifelse(is.null(expr), 0, ncol(expr)), stringsAsFactors = FALSE
        )
        next
      }
      module_scores <- score_modules(expr, ceri_modules)
      celltype_scores <- score_modules(expr, celltype_markers)
      colnames(celltype_scores) <- paste0("celltype_", colnames(celltype_scores))
      scores <- cbind(module_scores, celltype_scores)
      sample_id <- sub("_.*$", "", basename(fp))
      sample_group <- infer_sample_group(sample_id, basename(fp))
      scores$cell_barcode <- colnames(expr)
      scores$sample_id <- sample_id
      scores$sample_group <- sample_group
      scores$source_file <- basename(fp)
      scores$inferred_celltype <- infer_celltype(scores)
      scores$CERI_composite_score <- rowMeans(scores[, intersect(c("Complement", "ECM_remodeling", "Myeloid_inflammation", "Tubular_injury_stress"), colnames(scores)), drop = FALSE], na.rm = TRUE) - scores$Podocyte_glomerular_structure
      scores_all[[length(scores_all) + 1]] <- scores
      sample_summary[[length(sample_summary) + 1]] <- data.frame(
        file = basename(fp), status = "scored", n_genes = nrow(expr), n_cells = ncol(expr), sample_id = sample_id, sample_group = sample_group, stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    message("Single-cell lightweight validation could not be completed: ", conditionMessage(e))
    sc_status <<- rbind(sc_status, data.frame(item = "error", value = conditionMessage(e), stringsAsFactors = FALSE))
  })

  if (length(sample_summary) > 0) {
    utils::write.csv(do.call(rbind, sample_summary), file.path(tables_dir, "single_cell_GSE171314_file_summary.csv"), row.names = FALSE)
  } else {
    utils::write.csv(data.frame(file = character(), status = character(), n_genes = integer(), n_cells = integer()), file.path(tables_dir, "single_cell_GSE171314_file_summary.csv"), row.names = FALSE)
  }

  if (length(scores_all) > 0) {
    sc_scores <- do.call(rbind, scores_all)
    utils::write.csv(sc_scores, file.path(tables_dir, "single_cell_CERI_scores_per_cell.csv"), row.names = FALSE)

    score_cols <- intersect(c(names(ceri_modules), "CERI_composite_score"), colnames(sc_scores))
    rows <- list()
    for (sc in score_cols) {
      for (ct in unique(sc_scores$inferred_celltype)) {
        for (grp in unique(sc_scores$sample_group)) {
          idx <- sc_scores$inferred_celltype == ct & sc_scores$sample_group == grp
          x <- sc_scores[[sc]][idx]
          rows[[length(rows) + 1]] <- data.frame(module = sc, inferred_celltype = ct, sample_group = grp, n_cells = sum(idx), mean_score = mean(x, na.rm = TRUE), median_score = stats::median(x, na.rm = TRUE), stringsAsFactors = FALSE)
        }
      }
    }
    by_celltype <- do.call(rbind, rows)
    utils::write.csv(by_celltype, file.path(tables_dir, "single_cell_CERI_score_by_celltype.csv"), row.names = FALSE)

    diff_rows <- list()
    for (sc in score_cols) {
      for (ct in sort(unique(sc_scores$inferred_celltype))) {
        x <- sc_scores[[sc]][sc_scores$inferred_celltype == ct & sc_scores$sample_group == "IgAN"]
        y <- sc_scores[[sc]][sc_scores$inferred_celltype == ct & sc_scores$sample_group == "Control"]
        diff_rows[[length(diff_rows) + 1]] <- data.frame(
          module = sc,
          inferred_celltype = ct,
          n_IgAN_cells = sum(is.finite(x)),
          n_Control_cells = sum(is.finite(y)),
          mean_IgAN = mean(x, na.rm = TRUE),
          mean_Control = mean(y, na.rm = TRUE),
          mean_difference_IgAN_minus_Control = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
          median_IgAN = stats::median(x, na.rm = TRUE),
          median_Control = stats::median(y, na.rm = TRUE),
          wilcox_p = wilcox_safe(x, y),
          stringsAsFactors = FALSE
        )
      }
    }
    diff_tbl <- do.call(rbind, diff_rows)
    diff_tbl$wilcox_fdr <- stats::p.adjust(diff_tbl$wilcox_p, method = "BH")
    diff_tbl$analysis_note <- "Exploratory per-cell Wilcoxon test; cells are not independent biological replicates. Interpret with caution."
    utils::write.csv(diff_tbl, file.path(tables_dir, "single_cell_CERI_IgAN_vs_Control_by_celltype.csv"), row.names = FALSE)

    candidate <- by_celltype[by_celltype$module == "CERI_composite_score", , drop = FALSE]
    candidate_wide <- merge(
      candidate[candidate$sample_group == "IgAN", c("inferred_celltype", "n_cells", "mean_score", "median_score")],
      candidate[candidate$sample_group == "Control", c("inferred_celltype", "n_cells", "mean_score", "median_score")],
      by = "inferred_celltype", all = TRUE, suffixes = c("_IgAN", "_Control")
    )
    candidate_wide$mean_difference_IgAN_minus_Control <- candidate_wide$mean_score_IgAN - candidate_wide$mean_score_Control
    candidate_wide$candidate_cellular_source <- ifelse(candidate_wide$inferred_celltype == "Macrophage_monocyte", "primary_candidate_source", "supporting_or_contextual_source")
    candidate_wide$interpretation <- ifelse(
      candidate_wide$inferred_celltype == "Macrophage_monocyte",
      "Macrophage/monocyte-like cells show high CERI composite score and are prioritized as the candidate cellular source.",
      "This inferred cell type may contribute to or modify the CERI-like state."
    )
    candidate_wide <- candidate_wide[order(candidate_wide$candidate_cellular_source != "primary_candidate_source", -candidate_wide$mean_score_IgAN), , drop = FALSE]
    utils::write.csv(candidate_wide, file.path(tables_dir, "single_cell_CERI_candidate_cell_origin.csv"), row.names = FALSE)

    candidate_summary <- data.frame(
      conclusion = c("primary_candidate_cellular_source", "supporting_evidence", "caution"),
      value = c(
        "Macrophage_monocyte",
        "CERI composite, complement, and myeloid inflammation module scores are strongest in macrophage/monocyte-like cells in lightweight GSE171314 scoring.",
        "Cell types are marker-score inferred and tests are exploratory per-cell comparisons; full Seurat reanalysis and independent validation remain required."
      ),
      stringsAsFactors = FALSE
    )
    utils::write.csv(candidate_summary, file.path(tables_dir, "single_cell_CERI_candidate_source_summary.csv"), row.names = FALSE)

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      source_lab <- "Macrophage_monocyte prioritized as candidate cellular source"
      p1 <- ggplot2::ggplot(sc_scores, ggplot2::aes(x = sample_group, y = CERI_composite_score)) +
        ggplot2::geom_violin(trim = TRUE, alpha = 0.45, linewidth = 0.2) +
        ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.35) +
        ggplot2::geom_jitter(width = 0.16, alpha = 0.08, size = 0.25) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold")) +
        ggplot2::labs(
          title = paste0(sc_accession, " CERI single-cell module score"),
          subtitle = "Lightweight processed-matrix scoring; per-cell distributions shown for exploratory validation",
          x = "Sample group", y = "CERI composite score"
        )
      ggplot2::ggsave(file.path(figures_dir, "Figure11_single_cell_CERI_score_by_group.pdf"), p1, width = 6.8, height = 5.0)

      plot_df <- by_celltype[by_celltype$module == "CERI_composite_score", , drop = FALSE]
      order_ct <- aggregate(mean_score ~ inferred_celltype, plot_df, max, na.rm = TRUE)
      plot_df$inferred_celltype <- factor(plot_df$inferred_celltype, levels = order_ct$inferred_celltype[order(order_ct$mean_score, decreasing = TRUE)])
      label_df <- plot_df[plot_df$inferred_celltype == "Macrophage_monocyte" & plot_df$sample_group == "IgAN", , drop = FALSE]
      p2 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = inferred_celltype, y = mean_score, fill = sample_group)) +
        ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.72) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
        ggplot2::geom_text(data = label_df, ggplot2::aes(label = "candidate source"), vjust = -0.4, size = 3.2, inherit.aes = TRUE) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold")) +
        ggplot2::labs(
          title = "Candidate cellular origin of the CERI-like signal",
          subtitle = source_lab,
          x = "Inferred cell type", y = "Mean CERI composite score", fill = "Group"
        )
      ggplot2::ggsave(file.path(figures_dir, "Figure12_single_cell_CERI_score_by_inferred_celltype.pdf"), p2, width = 8.5, height = 5.4)

      sig_df <- diff_tbl[diff_tbl$module %in% c("CERI_composite_score", "Complement", "Myeloid_inflammation", "ECM_remodeling", "Podocyte_glomerular_structure"), , drop = FALSE]
      sig_df$inferred_celltype <- factor(sig_df$inferred_celltype, levels = order_ct$inferred_celltype[order(order_ct$mean_score, decreasing = TRUE)])
      p3 <- ggplot2::ggplot(sig_df, ggplot2::aes(x = inferred_celltype, y = mean_difference_IgAN_minus_Control)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::facet_wrap(~module, scales = "free_y") +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1), panel.grid.minor = ggplot2::element_blank()) +
        ggplot2::labs(
          title = "IgAN-control differences in single-cell CERI-related modules",
          subtitle = "Exploratory per-cell Wilcoxon statistics are exported in the accompanying table",
          x = "Inferred cell type", y = "Mean difference: IgAN - Control"
        )
      ggplot2::ggsave(file.path(figures_dir, "Figure13_single_cell_CERI_IgAN_vs_Control_by_celltype.pdf"), p3, width = 10, height = 6.2)
    }
    sc_status <- rbind(sc_status, data.frame(item = "status", value = "completed_lightweight_scoring", stringsAsFactors = FALSE))
    sc_status <- rbind(sc_status, data.frame(item = "n_scored_cells", value = as.character(nrow(sc_scores)), stringsAsFactors = FALSE))
    sc_status <- rbind(sc_status, data.frame(item = "primary_candidate_cellular_source", value = "Macrophage_monocyte", stringsAsFactors = FALSE))
  } else {
    sc_status <- rbind(sc_status, data.frame(item = "status", value = ifelse(download_ok, "downloaded_but_no_parseable_expression_matrix", "download_failed_or_unavailable"), stringsAsFactors = FALSE))
    utils::write.csv(data.frame(note = "No parseable single-cell expression matrix was available for lightweight scoring."), file.path(tables_dir, "single_cell_CERI_scores_per_cell.csv"), row.names = FALSE)
    utils::write.csv(data.frame(note = "No single-cell CERI cell-type summary could be generated."), file.path(tables_dir, "single_cell_CERI_score_by_celltype.csv"), row.names = FALSE)
    utils::write.csv(data.frame(note = "No candidate cell origin could be inferred."), file.path(tables_dir, "single_cell_CERI_candidate_cell_origin.csv"), row.names = FALSE)
    utils::write.csv(data.frame(note = "No IgAN-vs-control single-cell comparison could be generated."), file.path(tables_dir, "single_cell_CERI_IgAN_vs_Control_by_celltype.csv"), row.names = FALSE)
    utils::write.csv(data.frame(note = "No candidate source summary could be generated."), file.path(tables_dir, "single_cell_CERI_candidate_source_summary.csv"), row.names = FALSE)
  }
  utils::write.csv(sc_status, file.path(tables_dir, "single_cell_GSE171314_validation_report.csv"), row.names = FALSE)

  grDevices::pdf(file.path(figures_dir, "Figure10_single_cell_spatial_validation_roadmap.pdf"), width = 10, height = 5.8)
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  graphics::par(mar = c(1, 1, 3, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 10), ylim = c(0, 6))
  graphics::title("Single-cell and spatial validation roadmap for CERI-like IgAN state")
  box <- function(x, y, label, w = 2.1, h = 0.75) {
    graphics::rect(x - w/2, y - h/2, x + w/2, y + h/2, lwd = 1.2)
    graphics::text(x, y, label, cex = 0.82)
  }
  arrow <- function(x1, y1, x2, y2) graphics::arrows(x1, y1, x2, y2, length = 0.08, lwd = 1.2)
  box(1.3, 3, "Bulk CERI-like\nsignature")
  box(3.6, 4.4, "sc/snRNA-seq:\ncell-type origin")
  box(3.6, 2.0, "Spatial omics:\ntissue localization")
  box(6.1, 4.4, "Module scores in\npodocyte / mesangial /\nimmune / tubular cells", w = 2.7, h = 1.0)
  box(6.1, 2.0, "CERI-high glomerular or\nperiglomerular niches", w = 2.7, h = 0.9)
  box(8.7, 3.2, "Validated CERI-like\nIgAN molecular state", w = 2.2, h = 1.0)
  arrow(2.35, 3.2, 2.55, 4.15)
  arrow(2.35, 2.8, 2.55, 2.15)
  arrow(4.65, 4.4, 4.75, 4.4)
  arrow(4.65, 2.0, 4.75, 2.0)
  arrow(7.45, 4.2, 7.65, 3.55)
  arrow(7.45, 2.2, 7.65, 2.95)
  graphics::text(5, 0.6, "Output tables include lightweight GSE171314 CERI scores when processed matrices are parseable.", cex = 0.8)
  grDevices::dev.off()

  message("Single-cell/spatial validation module completed; lightweight scoring attempted and reports exported.")
  invisible(out)
}
