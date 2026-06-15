# Single-cell and spatial validation plan for the CERI-like IgAN state

## Current status

The bulk transcriptomic discovery and external bulk validation are now the main completed parts of the project.

Completed:

- IgAN-only bulk discovery analysis in GSE104948
- Small-cluster-like CERI signature definition
- External bulk validation across independent IgAN transcriptomic cohorts
- External pathway/immune comparison showing Complement up, ECM remodeling up, and podocyte/glomerular structure down

Not yet fully performed:

- Full single-cell or single-nucleus RNA-seq reanalysis
- Full spatial transcriptomics analysis
- RNAscope, IHC, IF, or other wet-lab tissue validation

The current single-cell and spatial module is therefore a validation-preparation module, not a definitive single-cell or spatial result.

## Why single-cell validation is needed

Bulk transcriptomic data cannot determine which cell type carries the CERI-like signal. Single-cell or single-nucleus RNA-seq validation is needed to identify whether the signal is enriched in:

- podocytes
- mesangial cells
- endothelial cells
- macrophage/monocyte populations
- T cells, B cells, or plasma cells
- tubular epithelial cells
- stromal/fibroblast-like compartments

The key question is whether the CERI-like state reflects one dominant cellular source or a coordinated multicellular tissue niche.

## Why spatial validation is needed

Single-cell data lose tissue position. Spatial validation is needed to determine whether CERI-high signals localize to:

- glomerular regions
- periglomerular regions
- tubulointerstitial regions
- vascular regions
- immune-infiltrate regions
- fibrotic or sclerotic regions

This is important because the CERI-like state is defined by coordinated Complement up, ECM remodeling up, and Podocyte/glomerular structure down. These features may represent a spatially organized injury niche rather than a single-cell-intrinsic state.

## Candidate dataset

The current config keeps GSE171314 as a candidate IgAN sc/snRNA-seq validation dataset.

At this stage, the main pipeline should not download or process the full single-cell dataset because full Seurat processing may be heavy and unstable in GitHub Actions. Instead, the repository exports a validation gene panel and a future-ready Seurat scoring template.

## Suggested CERI validation modules

### Complement module

C3, C4A, C4B, C1QA, C1QB, C1QC, CFB

### ECM remodeling module

COL1A1, COL1A2, COL3A1, COL5A2, FN1, VCAN, LAMB1, PLAU, HTRA1

### Podocyte/glomerular structure module

NPHS1, NPHS2, SYNPO, WT1, PODXL, PLCE1, TCF21, MAGI2, CLIC5

### Inflammatory module

IL18, ALOX5AP, FYN, CD68, LST1, AIF1

### Tubular injury/stress module

HAVCR1, LCN2, KRT8, KRT18, VCAM1, NDRG1, HMOX1

## Future single-cell analysis steps

1. Obtain a processed or processable IgAN sc/snRNA-seq object.
2. Annotate major kidney and immune cell types.
3. Compute CERI module scores by cell.
4. Summarize module scores by cell type and disease group.
5. Identify the main cellular sources of Complement, ECM remodeling, podocyte/glomerular structural loss, and inflammatory signals.
6. Test whether CERI-like signals are concentrated in one cell type or distributed across multiple compartments.

## Future spatial analysis steps

1. Obtain IgAN spatial transcriptomic data if available.
2. Annotate spatial regions or spots as glomerular, tubulointerstitial, vascular, or immune-rich.
3. Compute CERI module score per spot.
4. Map CERI-high regions to tissue compartments.
5. Compare spatial CERI-high regions with histology features such as sclerosis, interstitial fibrosis, cellular crescents, or immune infiltrates if metadata are available.

## Wet-lab-compatible validation alternatives

If public spatial transcriptomics data are not available, the CERI-like state can be validated with RNAscope, IHC, or IF using a targeted panel:

- Complement: C3, CFB, C1QA/C1QB/C1QC
- ECM remodeling: COL5A2, VCAN, FN1, LAMB1, PLAU
- Podocyte/glomerular structure: SYNPO, NPHS1, NPHS2, WT1, PODXL
- Inflammation: IL18, ALOX5AP, CD68

## Manuscript-safe wording

Use:

- We propose a single-cell and spatial validation framework.
- Future cell-type and tissue-localization validation will be required.
- The current data support a candidate CERI-like molecular state but do not yet establish its cellular origin or spatial localization.

Avoid:

- Single-cell analysis confirmed the cellular origin of the CERI-like state.
- Spatial transcriptomics confirmed the localization of the CERI-like state.
- The CERI-like state is definitively localized to podocytes/macrophages/mesangial cells.

## Current project outputs related to this plan

- output/tables/single_cell_signature_genes_for_validation.csv
- output/tables/single_cell_spatial_CERI_gene_panel.csv
- output/tables/single_cell_celltype_hypothesis.csv
- output/tables/single_cell_spatial_validation_plan.csv
- output/figures/Figure10_single_cell_spatial_validation_roadmap.pdf

## Final interpretation

The CERI-like IgAN state is currently supported by bulk discovery and external bulk validation. Single-cell and spatial analyses are the next required steps to determine whether this state reflects a specific cell type, a multicellular tissue niche, or a spatially organized injury microenvironment.
