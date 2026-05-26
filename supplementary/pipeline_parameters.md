# Pipeline Parameters

## QC Filtering

| Parameter | Value | Rationale |
|---|---|---|
| min_genes per cell | 200 | Remove empty droplets and dead cells |
| min_cells per gene | 3 | Remove lowly detected genes |
| max_mito_pct | 15% | Remove dying cells (paper used 10%; slightly relaxed) |
| max_genes per cell | 6,000 | Crude doublet proxy |

## Normalisation and Feature Selection

| Step | Parameter | Value |
|---|---|---|
| Normalisation | target_sum | 10,000 (CPM equivalent) |
| Log transform | log1p | applied after normalisation |
| HVG selection | n_top_genes | 2,000 |
| HVG flavour | flavor | seurat |
| Raw storage | adata.raw + adata.layers['lognorm'] | stored before scaling |

## Dimensionality Reduction

| Step | Parameter | Value |
|---|---|---|
| PCA | n_comps | 50 |
| Neighbours | n_neighbors | 15 |
| Neighbours | n_pcs | 30 |
| UMAP | default parameters | — |
| Leiden clustering | resolution | 0.5 |

Paper used resolution 0.8 producing 15 clusters.
Our resolution 0.5 produces 10 immune clusters — some subtypes merged.

## Differential Expression

| Analysis | Method | padj threshold | logFC threshold | Layer |
|---|---|---|---|---|
| Cluster marker DE | Wilcoxon | < 0.05 | > 0.5 | lognorm |
| Condition DE | Wilcoxon | < 0.05 | > 1.0 | lognorm |

All DE uses layer='lognorm', use_raw=False (Scanpy 1.10+ requirement).
Using scaled adata.X inflates fold changes — this was a key fix
identified during development.

Condition DE excludes non-immune clusters (5, 11) before testing.
Minimum 30 cells per condition per cluster required to run DE.

## TCGA Survival Analysis

| Parameter | Value |
|---|---|
| Expression data | RNA-seq V2 RSEM (cBioPortal skcm_tcga) |
| ISG score | mean z-score across 6 genes |
| Z-score method | per-gene across all patients (ddof=1) |
| Group split | median ISG score |
| Survival test | log-rank (lifelines) |
| ISG15 alias | stored as G1P2 in TCGA matrix |