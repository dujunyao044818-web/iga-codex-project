load_sc_geo <- function(cfg) {
  ensure_packages(c("GEOquery", "Seurat", "data.table"))
  supp_dir <- file.path("data/raw", cfg$single_cell$geo_accession)
  dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
  GEOquery::getGEOSuppFiles(cfg$single_cell$geo_accession, makeDirectory = FALSE, baseDir = supp_dir)
  files <- list.files(supp_dir, full.names = TRUE)
  rds <- files[grepl("\\.rds$|\\.RDS$", files)]
  if (length(rds)) return(readRDS(rds[1]))
  mtx <- files[grepl("matrix\\.mtx(\\.gz)?$", basename(files))]
  if (length(mtx)) return(Seurat::CreateSeuratObject(Seurat::Read10X(dirname(mtx[1])), project = cfg$single_cell$geo_accession))
  tab <- files[grepl("\\.txt(\\.gz)?$|\\.tsv(\\.gz)?$|\\.csv(\\.gz)?$", files)]
  if (length(tab)) {
    counts <- data.table::fread(tab[1], data.table = FALSE)
    rownames(counts) <- counts[[1]]
    counts[[1]] <- NULL
    return(Seurat::CreateSeuratObject(as.matrix(counts), project = cfg$single_cell$geo_accession))
  }
  stop("Downloaded GSE171314 supplemental files, but no RDS, 10x matrix, or tabular count matrix was found. Provide a compatible object in data/raw/GSE171314.", call. = FALSE)
}

annotate_cells <- function(seu) {
  markers <- list(Podocyte=c("NPHS1","NPHS2","PODXL"), Mesangial=c("PDGFRB","ITGA8","RGS5"), Endothelial=c("PECAM1","KDR","VWF"), Proximal_tubule=c("LRP2","SLC34A1","ALDOB"), Loop_of_Henle=c("UMOD","SLC12A1"), Distal_tubule=c("SLC12A3","PVALB"), Collecting_duct=c("AQP2","KRT8","KRT18"), T_NK=c("CD3D","NKG7","GNLY"), B_Plasma=c("MS4A1","CD79A","JCHAIN"), Myeloid=c("LYZ","C1QA","CD68"), Fibroblast=c("COL1A1","DCN","LUM"))
  avg <- Seurat::AverageExpression(seu, features = unique(unlist(markers)), assays = DefaultAssay(seu), verbose = FALSE)[[1]]
  cl <- levels(Idents(seu)); lab <- setNames(rep("Unassigned", length(cl)), cl)
  for (k in cl) {
    scores <- vapply(markers, function(g) mean(avg[intersect(g, rownames(avg)), k], na.rm = TRUE), numeric(1))
    lab[k] <- names(which.max(scores))
  }
  seu$cell_type <- lab[as.character(Idents(seu))]
  seu
}

run_single_cell <- function(signatures, cfg) {
  ensure_packages(c("Seurat", "ggplot2"))
  seu <- load_sc_geo(cfg)
  if (!inherits(seu, "Seurat")) stop("Single-cell object must be a Seurat object.")
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu <- subset(seu, subset = nFeature_RNA >= cfg$single_cell$min_features & nFeature_RNA <= cfg$single_cell$max_features & percent.mt <= cfg$single_cell$max_percent_mt)
  if (ncol(seu) > cfg$single_cell$max_cells) seu <- subset(seu, cells = sample(colnames(seu), cfg$single_cell$max_cells))
  seu <- SCTransform(seu, vst.flavor = cfg$single_cell$sct_method, verbose = FALSE)
  seu <- RunPCA(seu, verbose = FALSE) |> FindNeighbors(dims = 1:cfg$single_cell$dims) |> FindClusters(resolution = cfg$single_cell$resolution) |> RunUMAP(dims = 1:cfg$single_cell$dims)
  seu <- annotate_cells(seu)
  ggsave(file.path(cfg$output_dir, "figures", "scrna_umap_clusters.pdf"), DimPlot(seu, reduction = "umap", label = TRUE), width = 7, height = 6)
  ggsave(file.path(cfg$output_dir, "figures", "scrna_cell_type_annotation.pdf"), DimPlot(seu, reduction = "umap", group.by = "cell_type", label = TRUE), width = 8, height = 6)
  for (nm in names(signatures)) seu <- AddModuleScore(seu, features = list(intersect(signatures[[nm]], rownames(seu))), name = paste0(nm, "_Score"))
  score_cols <- grep("Subtype.*_Score1$", colnames(seu@meta.data), value = TRUE)
  p <- VlnPlot(seu, features = score_cols, group.by = "cell_type", pt.size = 0) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(cfg$output_dir, "figures", "scrna_subtype_signature_enrichment.pdf"), p, width = 11, height = 6)
  saveRDS(seu, file.path(cfg$output_dir, "models", "seurat_gse171314_processed.rds"))
  save_table(seu@meta.data, file.path(cfg$output_dir, "tables", "scrna_cell_metadata_with_scores.csv"))
  seu
}
