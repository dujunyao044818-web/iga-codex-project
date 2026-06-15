run_single_cell <- function(markers, cfg) {
  tables_dir <- file.path(cfg$output_dir, "tables")
  figures_dir <- file.path(cfg$output_dir, "figures")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  # Lightweight, reproducible planning module: records bulk-derived subtype genes
  # and defines the next single-cell/spatial validation steps without requiring
  # heavy Seurat or spatial dependencies during the main GitHub Actions run.
  out <- data.frame(signature_gene = unique(as.character(markers)), stringsAsFactors = FALSE)
  utils::write.csv(out, file.path(tables_dir, "single_cell_signature_genes_for_validation.csv"), row.names = FALSE)

  ceri_core <- c("C3", "C4A", "C4B", "C1QA", "C1QB", "C1QC", "CFB", "COL1A1", "COL1A2", "COL3A1", "COL5A2", "FN1", "VCAN", "LAMB1", "PLAU", "HTRA1", "SYNPO", "NPHS1", "NPHS2", "WT1", "PODXL", "PLCE1", "TCF21", "MAGI2", "CLIC5", "IL18", "ALOX5AP", "FYN")
  ceri_core <- unique(c(ceri_core, head(out$signature_gene, 30)))
  utils::write.csv(data.frame(validation_gene = ceri_core, stringsAsFactors = FALSE), file.path(tables_dir, "single_cell_spatial_CERI_gene_panel.csv"), row.names = FALSE)

  cell_hypothesis <- data.frame(
    expected_signal = c(
      "Complement production / activation",
      "ECM remodeling",
      "Podocyte/glomerular structural loss",
      "Myeloid inflammatory amplification",
      "T/B-cell immune context",
      "Tubular injury / stress"
    ),
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
    primary_features = c(
      "CERI module score by cell type and disease group",
      "CERI module score by spatial spot/domain",
      "C3/CFB/COL5A2/VCAN/SYNPO/NPHS2/IL18 selected panel",
      "M/S/T/C lesions, C3 intensity, eGFR/proteinuria if available"
    ),
    expected_support_for_novelty = c(
      "Shows the candidate molecular state has identifiable cellular origins.",
      "Shows the candidate molecular state is spatially organized rather than a bulk-only artifact.",
      "Provides wet-lab-compatible validation targets.",
      "Links transcriptomic state to pathological severity or prognosis."
    ),
    stringsAsFactors = FALSE
  )
  utils::write.csv(spatial_plan, file.path(tables_dir, "single_cell_spatial_validation_plan.csv"), row.names = FALSE)

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
  graphics::text(5, 0.6, "Output tables: single_cell_spatial_CERI_gene_panel.csv, single_cell_celltype_hypothesis.csv, single_cell_spatial_validation_plan.csv", cex = 0.8)
  grDevices::dev.off()

  message("Single-cell/spatial validation planning completed; CERI gene panel and validation roadmap exported.")
  invisible(out)
}
