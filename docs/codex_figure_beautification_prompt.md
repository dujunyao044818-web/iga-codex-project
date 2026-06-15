# Codex instruction: figure-first manuscript polishing and Nature-style figure beautification

Use this prompt in Codex after the single-cell reanalysis workflow has completed and the latest artifacts are available.

## Goal

Convert the current IgAN CERI project into a figure-first manuscript package. Results that can be shown visually should be shown as figures, not as main-text tables. Tables should be used only for supplementary material, QC reports, and exact numerical values.

## Required skills

Use the available figure/publication skills as much as possible, especially:

- nature-skill, if available
- scipilot-figure-skill, if available
- any journal-style plotting or figure-polishing skill available in the Codex environment

If a skill is not installed, do not fail. Continue with publication-quality ggplot2/base R defaults.

## Figure-first principle

Use figures for:

- sample grouping and PCA
- subtype heatmap
- volcano plot
- external validation consensus
- pathway/signature direction summaries
- mechanism model
- single-cell UMAPs
- single-cell CERI module score maps
- IgAN vs Control differences by cell type
- candidate macrophage/monocyte source summary

Use tables only for:

- exact statistics
- supplementary numerical results
- QC reports
- full gene lists
- full differential expression tables
- file parsing reports

## Main figures to produce or polish

### Figure 1: Study design and bulk discovery overview

Create or polish a multi-panel figure showing:

A. Workflow schematic: GSE104948 discovery -> IgAN-only clustering -> external bulk validation -> GSE171314 single-cell reanalysis -> CERI model.
B. PCA of all GSE104948 samples by diagnosis.
C. IgAN-only subtype/sample distribution.

Expected output:

- output/figures/Figure1_study_design_and_bulk_overview.pdf
- output/figures/Figure1_study_design_and_bulk_overview.png

### Figure 2: IgAN CERI-like molecular state discovery

Use the existing subtype heatmap and DEG results to make a clean multi-panel figure:

A. Consensus clustering/subtype heatmap.
B. Volcano plot of small-cluster-like state vs dominant IgAN group.
C. Top marker genes or marker-direction summary.
D. A small annotation box noting that the small subgroup is exploratory.

Expected output:

- output/figures/Figure2_CERI_discovery.pdf
- output/figures/Figure2_CERI_discovery.png

### Figure 3: External bulk validation consensus

Use external_validation_consensus_summary.csv and external_validation_high_low_pathway_comparison.csv.

Show:

A. Directional consensus across datasets for Complement, ECM remodeling, Podocyte/glomerular structure, Macrophage/monocyte, B cell, T cell, Plasma cell, Tubular injury/stress, Cell cycle/proliferation.
B. Effect size dot plot or lollipop plot.
C. Heatmap of CERI-high vs CERI-low module effects.

Expected output:

- output/figures/Figure3_external_bulk_validation.pdf
- output/figures/Figure3_external_bulk_validation.png

### Figure 4: Independent single-cell reanalysis of GSE171314

Use output_single_cell results from the independent workflow.

Show:

A. UMAP by inferred cell type.
B. UMAP by IgAN vs Control.
C. UMAP FeaturePlot of CERI composite score.
D. Bar/violin plot of CERI composite score by inferred cell type and group.
E. Mark macrophage/monocyte-like cells as the candidate cellular source.

Expected output:

- output_single_cell/figures/Figure4_single_cell_CERI_reanalysis.pdf
- output_single_cell/figures/Figure4_single_cell_CERI_reanalysis.png

### Figure 5: Cell-type module differences and candidate cellular source

Use sc_CERI_IgAN_vs_Control_by_celltype.csv and sc_CERI_candidate_cell_origin.csv.

Show:

A. IgAN-Control mean differences for CERI composite, Complement, Myeloid inflammation, ECM remodeling, Podocyte/glomerular structure.
B. Highlight macrophage/monocyte-like cells.
C. Highlight reduced podocyte/glomerular structural program in IgAN podocyte-like cells.
D. Include a concise caution label: exploratory per-cell statistics; marker-score-inferred cell types.

Expected output:

- output_single_cell/figures/Figure5_single_cell_module_differences.pdf
- output_single_cell/figures/Figure5_single_cell_module_differences.png

### Figure 6: Mechanistic model

Convert the current mechanism diagram into a cleaner graphical abstract-style model:

CERI-like IgAN state -> macrophage/monocyte-associated complement and inflammatory activation -> ECM remodeling -> podocyte/glomerular structural loss -> IgAN injury progression.

Expected output:

- output/figures/Figure6_CERI_mechanism_model.pdf
- output/figures/Figure6_CERI_mechanism_model.png

## Style requirements

Use a clean Nature-like style:

- White background
- Minimal grid lines
- Consistent typography
- Clear panel labels: A, B, C, D
- Colorblind-friendly palette
- Avoid rainbow colors
- Use consistent colors across figures:
  - IgAN: red/orange family
  - Control: blue/gray family
  - CERI-high: red/orange
  - CERI-low: blue/gray
  - Complement: purple/red
  - ECM remodeling: brown/orange
  - Podocyte/glomerular structure: blue/teal
  - Macrophage/monocyte: green/dark teal
- Keep font sizes readable for journal figures
- Export both PDF and PNG when possible
- Use 300 dpi for PNG files

## Results writing instruction

Update docs/results_draft_CERI_IgAN.md so that the Results text is figure-led.

Do not describe long tables in the main Results. Instead, write:

- Figure 1 shows...
- Figure 2 demonstrates...
- Figure 3 validates...
- Figure 4 maps...
- Figure 5 supports...
- Figure 6 summarizes...

Keep tables as supplementary references only.

## Manuscript-safe wording

Use cautious wording:

- candidate CERI-like molecular state
- small-cluster-like IgAN molecular state
- supportive external validation
- independent single-cell reanalysis
- macrophage/monocyte-like cells as a candidate cellular source
- marker-score-inferred cell type labels
- exploratory per-cell statistics

Avoid:

- definitive molecular subtype
- proven cellular origin
- spatial confirmation
- clinical validation
- independent cells as patient-level replicates

## Final deliverables

After running, provide:

1. List of newly generated figures.
2. List of updated Results sections.
3. Notes on which tables remain supplementary.
4. Any figure-generation failures and why they failed.
5. A short manuscript-ready summary of the figure set.
