find_single_cell_artifact_dir <- function(cfg) {
  candidates <- unique(c(
    Sys.getenv("GSE171314_REANALYSIS_DIR", unset = NA_character_),
    Sys.getenv("SC_REANALYSIS_DIR", unset = NA_character_),
    if (!is.null(cfg$single_cell$reanalysis_artifact_dir)) cfg$single_cell$reanalysis_artifact_dir else NA_character_,
    "input/GSE171314_reanalysis",
    "data/external/GSE171314_reanalysis",
    "output/GSE171314_reanalysis",
    "output/single_cell_reanalysis"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hits <- candidates[file.exists(candidates) & file.info(candidates)$isdir]
  if (length(hits)) hits[1] else NA_character_
}

copy_sc_file <- function(src_dir, filename, dest_dirs) {
  src <- file.path(src_dir, filename)
  if (!file.exists(src)) return(FALSE)
  invisible(vapply(dest_dirs, function(dest) {
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    file.copy(src, file.path(dest, filename), overwrite = TRUE)
  }, logical(1)))
  TRUE
}

run_single_cell_reanalysis_integration <- function(cfg) {
  message("Integrating independent GSE171314 single-cell reanalysis artifact if available.")
  report_dir <- file.path(cfg$output_dir, "reports")
  main_dir <- file.path(cfg$output_dir, "main_figures")
  supp_sc_dir <- file.path(cfg$output_dir, "supplementary_figures", "single_cell")
  table_dir <- file.path(cfg$output_dir, "tables")
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(supp_sc_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  expected_figures <- c(
    "Figure11_sc_UMAP_celltype.pdf",
    "Figure11b_sc_UMAP_group.pdf",
    "Figure12_sc_CERI_score_umap.pdf",
    "Figure13_sc_CERI_by_celltype_group.pdf",
    "Figure14_sc_CERI_IgAN_vs_Control_by_celltype.pdf"
  )
  main_supporting <- c(
    "Figure12_sc_CERI_score_umap.pdf",
    "Figure13_sc_CERI_by_celltype_group.pdf",
    "Figure14_sc_CERI_IgAN_vs_Control_by_celltype.pdf"
  )
  expected_tables <- c(
    "sc_qc_summary.csv",
    "sc_cell_metadata_with_CERI_scores.csv",
    "sc_celltype_annotation_summary.csv",
    "sc_CERI_candidate_cell_origin.csv",
    "sc_CERI_candidate_source_summary.csv"
  )
  src_dir <- find_single_cell_artifact_dir(cfg)
  status_rows <- list()
  if (is.na(src_dir) || !nzchar(src_dir)) {
    for (f in expected_figures) {
      save_placeholder_figure(file.path(supp_sc_dir, tools::file_path_sans_ext(f)), f, "Independent GSE171314 artifact was not available to this run.", width = 8, height = 5)
      status_rows[[length(status_rows) + 1]] <- data.frame(file = f, type = "figure", status = "placeholder", stringsAsFactors = FALSE)
    }
    save_table(data.frame(
      candidate_source = "macrophage/monocyte-like cells",
      evidence_status = "artifact unavailable in this run",
      caution = "Cell-level tests are exploratory and not independent biological replicates.",
      stringsAsFactors = FALSE
    ), file.path(table_dir, "sc_CERI_candidate_source_summary.csv"))
  } else {
    for (f in expected_figures) {
      copied_supp <- copy_sc_file(src_dir, f, supp_sc_dir)
      if (f %in% main_supporting) copy_sc_file(src_dir, f, main_dir)
      status_rows[[length(status_rows) + 1]] <- data.frame(file = f, type = "figure", status = ifelse(copied_supp, "copied", "missing"), stringsAsFactors = FALSE)
    }
    for (f in expected_tables) {
      copied <- copy_sc_file(src_dir, f, table_dir)
      status_rows[[length(status_rows) + 1]] <- data.frame(file = f, type = "table", status = ifelse(copied, "copied", "missing"), stringsAsFactors = FALSE)
    }
  }
  status <- do.call(rbind, status_rows)
  save_table(status, file.path(table_dir, "single_cell_reanalysis_import_status.csv"))
  report <- c(
    "# GSE171314 single-cell reanalysis integration report",
    "",
    paste0("Artifact source directory: ", ifelse(is.na(src_dir) || !nzchar(src_dir), "not available", src_dir)),
    "",
    "The independent GSE171314 Seurat reanalysis from workflow run 27565463060 is treated as a formal external single-cell validation layer when its artifact files are provided to this pipeline.",
    "",
    "Interpretation: macrophage/monocyte-like cells are the primary candidate cellular source of the CERI-like signal when the imported candidate-source tables and figures support that pattern.",
    "",
    "Caution: cell-level tests are exploratory because individual cells are not independent biological replicates; subject-level or pseudobulk validation remains required before strong causal or diagnostic claims.",
    "",
    "## Import status",
    paste0("- ", status$file, " (", status$type, "): ", status$status)
  )
  writeLines(report, file.path(report_dir, "single_cell_reanalysis_report.md"))
  invisible(status)
}
