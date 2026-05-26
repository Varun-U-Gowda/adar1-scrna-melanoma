# Known Limitations

## 1. Pseudoreplication in Condition DE

Differential expression between ADAR KO and Control was performed using
Wilcoxon rank-sum test on individual cells. With only n=2 biological
replicates per condition, this treats cells from 2 mice as independent
observations — a statistical issue known as pseudoreplication.

The gold standard approach is pseudobulk DE: aggregate counts per
sample per cell type (producing a 2 vs 2 matrix), then apply DESeq2
or edgeR. With n=2 this is underpowered, but would produce
conservative, biologically robust results.

All condition DE results in this project should be interpreted as
exploratory. The ISG findings are consistent with the original paper
(Ishizuka et al. 2019) and reproduced across multiple independent
cell type clusters, which partially mitigates the pseudoreplication
concern.

## 2. No Doublet Detection

Scrublet or DoubletFinder was not run. A maximum gene threshold
(max_genes < 6,000) is used as a crude proxy for doublet removal.
At the dataset scale (7,403 cells) the impact is likely minor but
cannot be quantified without explicit doublet scoring.

Recommended addition: run Scrublet before filtering and remove
predicted doublets (score > 0.25 typically).

## 3. No Ambient RNA Correction

SoupX or similar ambient RNA decontamination was not applied.
GEO-sourced 10x data can contain ambient RNA contamination
from lysed cells, which may inflate expression of highly expressed
genes (e.g. haemoglobin genes, ribosomal genes) in all cells.

## 4. Clustering Resolution

Leiden resolution 0.5 produces 10 immune clusters. The original
paper used resolution 0.8 and identified 15 populations. Two
MoDC subtypes (clusters 3 and 9) and potentially M1/M2 macrophage
subtypes are merged at this resolution. Increasing resolution would
provide finer annotation but requires more careful validation.

## 5. TCGA Survival Is Correlational

The ISG survival association reflects patients whose tumours have
high IFN activity. The causal driver — specifically ADAR1 loss —
is not confirmed in TCGA because ADAR1 mutation/expression status
was not linked to the survival analysis. A Cox model adjusting for
ADAR1 expression, tumour stage, age, and sex would provide stronger
evidence of independent prognostic value.