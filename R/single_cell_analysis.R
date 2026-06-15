run_single_cell <- function(markers, cfg) {
  # Lightweight placeholder: records the bulk-derived subtype signature for later scRNA-seq validation.
  # Full Seurat analysis can be added after the bulk pipeline is stable.
  out <- data.frame(signature_gene = markers, stringsAsFactors = FALSE)
  utils::write.csv(out, file.path(cfg$output_dir, "tables", "single_cell_signature_genes_for_validation.csv"), row.names = FALSE)
  message("Single-cell validation placeholder completed; markers exported for future Seurat analysis.")
  invisible(out)
}
