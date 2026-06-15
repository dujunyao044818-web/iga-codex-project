# Results draft: CERI-like IgAN molecular state

## 1. Study design and bulk transcriptomic sample selection

We first constructed a bulk transcriptomic analysis workflow to explore molecular heterogeneity in IgA nephropathy (IgAN). The discovery analysis was based on GSE104948, a human kidney transcriptomic cohort containing multiple renal disease groups. Diagnosis information was extracted from the `diagnosis:ch1` metadata field. Principal component analysis was first performed across all available samples as a quality-control and global structure assessment step. IgAN-specific downstream analyses were then restricted to the IgAN subset, yielding 27 IgAN samples for molecular subtype discovery.

Because the available IgAN discovery cohort was modest in size, the objective of the discovery step was not to define definitive clinical subtypes, but to identify candidate transcriptomic states that could be evaluated through independent datasets and cell-type-level analyses.

## 2. Unsupervised clustering identified a small CERI-like IgAN molecular state

Consensus clustering of the 27 IgAN samples supported a dominant molecular group and a small distinct subgroup. The small subgroup contained four IgAN samples, whereas the remaining samples formed the dominant group. Given the small size of this subgroup, we conservatively refer to it as a small-cluster-like molecular state rather than a robust disease subtype.

Differential expression analysis between the small subgroup and the dominant IgAN group identified a set of marker genes, including NDNF, ABCC3, SRGAP2, CDC42EP2, RRAD, PDLIM2, MAGI2, PLCE1, VCAN, COL5A2, ALOX5AP, FYN, CLIC5, and IL18. These markers suggested coordinated changes in immune activation, complement-related signaling, extracellular matrix remodeling, and glomerular structural programs. Based on this pattern, we termed the candidate state the Complement–ECM Remodeling–Injury-like IgAN molecular state, abbreviated as the CERI-like state.

## 3. External bulk validation supported coordinated complement activation, ECM remodeling, and podocyte structural loss

To test whether the CERI-like signature reflected a reproducible transcriptomic pattern rather than a discovery-cohort artifact, we evaluated the signature across independent IgAN transcriptomic datasets. External validation was performed by scoring the CERI-like signature and comparing pathway-level module activity between CERI-high and CERI-low IgAN samples.

Across the external validation datasets, CERI-high samples consistently showed increased complement and extracellular matrix remodeling signatures. Complement activity was increased in all three evaluable external datasets, with a mean effect size of approximately 1.23 and a best false-discovery-rate-adjusted P value of 7.0 × 10^-6. ECM/stromal remodeling was also increased in all three datasets, with a mean effect size of approximately 0.90 and a best false-discovery-rate-adjusted P value of 4.0 × 10^-5.

In contrast, the podocyte/glomerular structural signature showed the opposite pattern. Podocyte/glomerular structural programs were reduced in all three external datasets, with a mean effect size of approximately -1.12 and a best false-discovery-rate-adjusted P value of 4.0 × 10^-5. This reciprocal pattern suggested that the CERI-like state is characterized by complement and matrix-remodeling activation coupled to loss of glomerular structural programs.

Additional immune and injury-related signals were also observed. Macrophage/monocyte, B-cell, plasma-cell, T-cell, cell-cycle/proliferation, and tubular injury/stress signatures tended to be higher in CERI-high samples, suggesting that the CERI-like state may represent a broader inflammatory and tissue-injury program rather than a single isolated pathway.

## 4. The CERI-like state represents a candidate multicompartment injury program

The external validation results indicated that the CERI-like signature was not simply a marker-gene list derived from a small discovery subgroup. Instead, the signature defined a reproducible pattern across independent IgAN transcriptomic cohorts: complement activation and ECM remodeling were increased, whereas podocyte/glomerular structural programs were decreased.

This pattern is biologically consistent with a multicompartment injury program in IgAN. Complement activation may reflect inflammatory amplification, ECM remodeling may reflect mesangial or stromal injury responses, and reduced podocyte/glomerular structural signatures may reflect glomerular epithelial injury or loss of differentiated podocyte programs. We therefore interpreted the CERI-like state as a candidate inflammatory-remodeling-injury state requiring further cell-type and tissue-localization validation.

## 5. Lightweight single-cell module scoring prioritized macrophage/monocyte-like cells as a candidate source

We next performed a lightweight single-cell module-score validation using GSE171314, an IgAN single-cell transcriptomic dataset containing four IgAN kidney biopsy samples and one control sample. Processed expression matrices were parsed and scored for CERI-related modules, including complement, ECM remodeling, podocyte/glomerular structure, myeloid inflammation, adaptive immune context, and tubular injury/stress.

This lightweight analysis identified macrophage/monocyte-like cells as the leading candidate cellular source of the CERI-like signal. In IgAN samples, macrophage/monocyte-like cells showed the highest mean CERI composite score among inferred cell types. These cells also showed prominent complement and myeloid inflammatory module activity, supporting the hypothesis that the inflammatory and complement-active component of the CERI-like state is linked to macrophage/monocyte-like cells.

In contrast, podocyte-like cells showed a low CERI composite score and reduced podocyte/glomerular structural module activity in IgAN compared with control. This direction was consistent with the external bulk validation pattern, in which podocyte/glomerular structural signatures were consistently reduced in CERI-high IgAN samples.

## 6. Independent Seurat reanalysis of GSE171314 supported a macrophage-associated CERI-like signal

To further evaluate the single-cell findings, we performed an independent Seurat-based reanalysis of GSE171314 in a separate workflow that was not linked to the main bulk pipeline. Five processed matrices were loaded, including four IgAN samples and one control sample. After quality control, 4,002 cells were retained for downstream normalization, variable feature selection, scaling, principal component analysis, clustering, UMAP embedding, marker-score-based cell-type annotation, and CERI module scoring.

After quality control, the retained cells included 2,709 IgAN-derived cells and 1,293 control-derived cells. Marker-score-based annotation identified tubular, T-cell, B/plasma-cell, endothelial, macrophage/monocyte, mesangial/stromal, and podocyte-like compartments. IgAN samples showed relatively increased tubular, macrophage/monocyte, and mesangial/stromal compartments, while the podocyte-like compartment was relatively reduced compared with control.

Consistent with the lightweight scoring analysis, macrophage/monocyte-like cells showed the highest CERI composite score among IgAN cell compartments in the independent Seurat reanalysis. The mean CERI composite score in IgAN macrophage/monocyte-like cells was approximately 0.32. Tubular, B/plasma-cell, T-cell, and mesangial/stromal compartments showed lower but positive CERI composite scores, whereas podocyte-like cells showed a negative CERI composite score.

Module-level comparisons further supported the proposed mechanism. Complement activity was highest in macrophage/monocyte-like cells and was modestly higher in IgAN macrophage/monocyte-like cells than in control macrophage/monocyte-like cells. ECM remodeling signatures were increased across several compartments, including tubular, T-cell, and endothelial-like cells, suggesting that matrix remodeling may represent a multicompartment tissue response rather than a single-cell-type phenomenon. Most notably, podocyte/glomerular structural programs were markedly reduced in IgAN podocyte-like cells compared with control podocyte-like cells, consistent with the loss of glomerular structural programs observed in external bulk validation.

These findings support a model in which the CERI-like IgAN state reflects a macrophage-associated inflammatory and complement-active injury program coupled to multicompartment ECM remodeling and podocyte/glomerular structural loss.

## 7. Proposed working model

Together, the bulk discovery, external bulk validation, and independent single-cell reanalysis support a coherent working model. A small-cluster-like IgAN molecular state identified in the discovery cohort is characterized by a CERI-like transcriptomic program. This program is externally reproducible as a coordinated pattern of complement activation, ECM remodeling, and reduced podocyte/glomerular structural signatures. Single-cell reanalysis suggests that macrophage/monocyte-like cells are the primary candidate cellular source of the inflammatory and complement-active component of this state, while podocyte-like cells exhibit reduced glomerular structural programs.

We therefore propose that the CERI-like state represents a candidate multicompartment IgAN injury state involving macrophage/monocyte-associated inflammatory activation, complement activity, ECM remodeling, and podocyte/glomerular structural loss. Because the discovery subgroup was small and the single-cell cell-type labels were inferred using marker scores, this state should be interpreted as a hypothesis-generating molecular state requiring further validation using independent single-cell datasets, spatial transcriptomics, and tissue-based assays such as RNAscope, immunohistochemistry, or immunofluorescence.

## Suggested Results subsection titles

1. Bulk transcriptomic profiling identifies candidate molecular heterogeneity in IgAN
2. Consensus clustering reveals a small CERI-like IgAN molecular state
3. External validation supports complement activation, ECM remodeling, and podocyte structural loss
4. CERI-high IgAN samples show inflammatory and injury-associated pathway activation
5. Single-cell module scoring prioritizes macrophage/monocyte-like cells as candidate CERI signal carriers
6. Independent Seurat reanalysis supports a macrophage-associated inflammatory-remodeling-injury program
7. A proposed multicompartment CERI-like injury model in IgAN

## Manuscript-safe key claim

The safest central claim is:

"We identified a candidate CERI-like IgAN molecular state characterized by coordinated complement activation, ECM remodeling, and reduced podocyte/glomerular structural programs. External bulk validation and independent single-cell reanalysis support macrophage/monocyte-like cells as a candidate source of the inflammatory and complement-active component of this state."

## Claims to avoid

Avoid claiming that:

- the study discovered definitive IgAN molecular subtypes;
- the CERI-like state is fully validated clinically;
- macrophage/monocyte-like cells are the proven cellular origin of CERI;
- spatial localization has been confirmed;
- the per-cell Wilcoxon tests represent independent patient-level statistical evidence.
