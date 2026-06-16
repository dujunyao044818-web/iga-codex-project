copy_if_exists <- function(from, to_dir) {
  if (!file.exists(from)) return(FALSE)
  dir.create(to_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE)
}

run_figure_organization_and_reports <- function(cfg) {
  message("Organizing main and supplementary figures and writing reports.")
  main_dir <- file.path(cfg$output_dir, "main_figures")
  supp_dir <- file.path(cfg$output_dir, "supplementary_figures")
  report_dir <- file.path(cfg$output_dir, "reports")
  dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

  main_map <- c(
    Figure1_study_design_overview = "Figure10_single_cell_spatial_validation_roadmap.pdf",
    Figure2_CERI_subtype_heatmap = "Figure2_subtype_heatmap_publication.pdf",
    Figure3_CERI_DEG_mechanism = "Figure3_volcano_publication.pdf",
    Figure4_external_bulk_validation = "Figure8_external_validation_consensus.pdf",
    Figure5_external_bulk_ML_gene_panel = "FigureS7_external_bulk_ML_model_comparison.pdf",
    Figure6_single_cell_candidate_source = "Figure10_single_cell_spatial_validation_roadmap.pdf"
  )

  for (nm in names(main_map)) {
    src <- file.path(cfg$output_dir, "figures", main_map[[nm]])
    if (!file.exists(src)) src <- file.path(cfg$output_dir, "supplementary", "ml_external_bulk", main_map[[nm]])
    if (!copy_if_exists(src, main_dir)) {
      save_placeholder_figure(file.path(main_dir, nm), nm, "Source figure was unavailable; see pipeline logs and supplementary reports.", width = 8, height = 5)
    }
  }

  fig_files <- list.files(file.path(cfg$output_dir, "figures"), pattern = "\\.(pdf|png)$", full.names = TRUE, recursive = TRUE)
  if (length(fig_files) > 0) invisible(vapply(fig_files, copy_if_exists, logical(1), to_dir = supp_dir))

  main_report <- c(
    "# Main results summary",
    "",
    "Recommended main figures are stored in `output/main_figures/`.",
    "",
    "- Figure 1: study design and dataset overview",
    "- Figure 2: IgAN-only CERI-like subtype heatmap and representative markers",
    "- Figure 3: CERI-like DEG and pathway mechanism",
    "- Figure 4: external bulk validation of CERI-like signature",
    "- Figure 5: external bulk machine-learning gene panel, only if robust enough",
    "- Figure 6: single-cell validation and candidate cellular source",
    "",
    "Subtype2 should be treated as an exploratory small-cluster subtype when it falls below the n=5 threshold."
  )
  writeLines(main_report, file.path(report_dir, "main_results_summary.md"))

  supp_report <- c(
    "# Supplementary results inventory",
    "",
    "Supplementary figures are stored in `output/supplementary_figures/`, `output/supplementary_figures/single_cell/`, `output/supplementary/qc/`, `output/supplementary/ml_external_bulk/`, and `output/supplementary/functional_enrichment/` when generated.",
    "",
    paste0("Number of copied supplementary figure files: ", length(fig_files))
  )
  writeLines(supp_report, file.path(report_dir, "supplementary_results_inventory.md"))

  fig_report <- c(
    "# Figure allocation report",
    "",
    "QC-heavy figures, PCA QC plots, consensus CDF/delta plots, density/boxplots, WGCNA/PPI/ML diagnostics, and small-cluster diagnostics are supplementary, not main text figures.",
    "",
    "Main figure directory: `output/main_figures/`",
    "Supplementary figure directory: `output/supplementary_figures/`"
  )
  writeLines(fig_report, file.path(report_dir, "figure_allocation_report.md"))

  qc_corr <- c(
    "# QC corrections report",
    "",
    "Plotting functions use finite-data checks and placeholder figures to avoid errors such as `need finite 'ylim' values`.",
    "",
    "If a dataset lacks usable metadata or expression values, the pipeline records a skipped status rather than failing."
  )
  writeLines(qc_corr, file.path(report_dir, "qc_corrections_report.md"))

  gaps <- c(
    "# Remaining analysis gaps report",
    "",
    "The main text remains focused on the CERI-like IgAN mechanism: complement activation, ECM remodeling, immune/inflammatory amplification, and glomerular structural injury.",
    "",
    "QC, ML, WGCNA, PPI, and full enrichment outputs are allocated to supplementary results unless external evidence is exceptionally strong.",
    "",
    "The current pipeline avoids claiming a robust diagnostic classifier. External ML analyses are exploratory marker-prioritization and stratification tools.",
    "",
    "Remaining gaps:",
    "- Harmonized external expression matrices are required for full external QC and fully independent model fitting.",
    "- Single-cell GSE171314 conclusions depend on availability of the independent Seurat reanalysis artifact during the integrated pipeline run.",
    "- WGCNA should be restricted to larger external cohorts or clearly labelled exploratory."
  )
  writeLines(gaps, file.path(report_dir, "remaining_analysis_gaps_report.md"))
  invisible(TRUE)
}
