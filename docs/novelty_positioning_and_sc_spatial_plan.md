# Novelty positioning and single-cell/spatial validation plan

## Working biological concept

The current working concept is a candidate **CERI-like IgAN molecular state**:

**CERI = Complement–ECM Remodeling–Injury-like**

This should be described cautiously as a reproducible transcriptomic state or candidate subgroup, not as a definitive new clinical subtype.

## What is not new

The following biological components are already established in IgAN and should not be claimed as first discoveries:

1. Complement activation is involved in IgAN pathogenesis and C3 deposition is a common pathological feature.
2. Extracellular matrix remodeling, matrix expansion, and fibrosis-related processes are common features of progressive IgAN and other chronic kidney diseases.
3. Podocyte and glomerular structural injury are established contributors to proteinuria and glomerular disease progression.

## What may be new in this project

The novelty is the integrated patient-level transcriptomic pattern:

1. Discovery of a small-cluster-like IgAN molecular state in GSE104948.
2. Construction of a marker signature from this small cluster.
3. External detection of high-score samples in independent IgAN transcriptomic datasets.
4. Cross-cohort pathway comparison showing a coordinated pattern:
   - Complement signature increased
   - ECM/stromal remodeling increased
   - Podocyte/glomerular structural signature decreased
5. Framing this coordinated signal as a candidate CERI-like state.

## Recommended novelty statement

Although complement activation, extracellular matrix remodeling, and podocyte injury have each been implicated in IgAN pathogenesis, whether these processes co-define a reproducible transcriptomic patient subgroup remains unclear. Here, we identified a small-cluster-like IgAN molecular state characterized by coordinated complement activation, ECM remodeling, and reduced podocyte/glomerular structural signatures, and validated this pattern across independent IgAN transcriptomic cohorts.

## Claims to avoid

Avoid:

- We discovered complement activation in IgAN.
- We discovered ECM remodeling in IgAN.
- We discovered podocyte injury in IgAN.
- We identified two robust IgAN subtypes.
- This is a definitive clinical subtype.

Use instead:

- candidate molecular state
- exploratory small-cluster-like subgroup
- reproducible transcriptomic signature
- CERI-like IgAN state

## Single-cell validation plan

Goal: determine which cell types carry the CERI-like signature.

Primary analyses:

1. Compute CERI module score in each cell.
2. Compare scores by cell type and disease group.
3. Test whether complement, ECM, podocyte-structure, and immune modules localize to distinct cell populations.
4. Prioritize cell types for wet-lab validation.

Candidate cell-type hypotheses:

- Complement: macrophage/monocyte, mesangial, endothelial, tubular compartments depending on dataset.
- ECM remodeling: mesangial, stromal/fibroblast-like, injured tubular, and periglomerular compartments.
- Podocyte/glomerular structural loss: podocyte and glomerular epithelial compartment.
- Immune activation: macrophage/monocyte, T-cell, B-cell, plasma-cell compartments.

## Spatial validation plan

Goal: determine where the CERI-like state is located in kidney tissue.

Primary analyses:

1. Compute CERI module score per spatial spot/domain.
2. Overlay CERI score with glomerular, tubulointerstitial, vascular, and immune-infiltrate domains.
3. Test whether CERI-high regions are glomerular, periglomerular, or tubulointerstitial.
4. Validate selected genes using RNAscope/IHC if spatial omics is unavailable.

Suggested validation panel:

- Complement: C3, C4A/C4B, CFB, C1QA/C1QB/C1QC
- ECM remodeling: COL1A1, COL3A1, COL5A2, FN1, VCAN, LAMB1, PLAU, HTRA1
- Podocyte/glomerular structure: NPHS1, NPHS2, SYNPO, WT1, PODXL, PLCE1, TCF21, MAGI2, CLIC5
- Inflammation: IL18, ALOX5AP, FYN, CD68, LST1, AIF1

## Figure strategy

Recommended final structure:

- Figure 1: GSE104948 PCA and disease composition
- Figure 2: IgAN-only subtype heatmap
- Figure 3: subtype DEG volcano
- Figure 4: subtype QC and small-cluster characterization
- Figure 5: exploratory LASSO / marker model
- Figure 6: external small-cluster signature validation
- Figure 7: external high-vs-low pathway/immune heatmap
- Figure 8: consensus external validation summary
- Figure 9: CERI-like mechanism model
- Figure 10: single-cell/spatial validation roadmap

## Manuscript-safe conclusion

These findings support a candidate CERI-like IgAN molecular state characterized by complement activation, ECM remodeling, and podocyte/glomerular structural injury. Because the discovery cluster contains only four samples, the result should be framed as an exploratory but externally supported transcriptomic state requiring additional single-cell, spatial, and clinical-pathological validation.
