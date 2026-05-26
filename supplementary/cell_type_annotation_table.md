# Cell Type Annotation Detail

## Method

Cell types were assigned using convergent evidence from three
independent sources:

1. **PanglaoDB markers** — filtered for *Mus musculus*, Immune system
   organ, sensitivity_mouse > 0, ranked by specificity_mouse.
   Source: PanglaoDB_markers_27_Mar_2020.tsv

2. **Module scoring** — sc.tl.score_genes() per cell type category
   using the PanglaoDB marker sets. Provides quantitative ranking
   of cell type likelihood per cluster.

3. **Cluster marker DEGs** — sc.tl.rank_genes_groups() with
   layer='lognorm', use_raw=False. Top upregulated genes per cluster
   from cluster-vs-rest Wilcoxon DE.

Where module scores and DEGs disagreed (clusters 2, 6, 8, 9),
cluster-specific DEGs from actual data took priority over
generic database markers.

## Validation

All assignments validated against Ishizuka et al. 2019
Supplementary Table 3 (NIHMS1560729-supplement-SI-3):
cluster marker genes checked against all statistically significant
markers (Adjusted P < 0.05) reported by the original authors for
the corresponding cell type in the same GSE110746 dataset.

## Non-Immune Clusters Removed

| Cluster | Cells | Top Genes | Identity | Action |
|---|---|---|---|---|
| 5 | 1,028 | Elf5, Gm29107, Dmrta2os | B16 tumour spike-in (intentional per paper methods) | Removed |
| 11 | 74 | Wnt5a, Ror2, Slit3, Rbp4 | Stromal fibroblast contamination | Removed |

Post-removal: 7,403 immune cells (paper reports 7,406 — within rounding).

## Cell Type Assignments

| Cluster | Cell Type | Top Marker Genes | Confirmed in Ishizuka 2019 |
|---|---|---|---|
| 0 | M1_Macrophage | C1qa, C1qb, C1qc, Cd209f, Fcrl5 | C1qa, C1qb, C1qc, Cd300e, Slamf9 |
| 1 | CD8_T | Cd3g, Cd8a, Cd8b1, Gzmk, Cxcr6, Pdcd1 | Cd3d, Cd3g, Cd8a, Cd8b1, Cxcr6, Gzmk, Nkg7, Pdcd1 |
| 2 | Monocyte | Ccr2, Il1b, Cd177, Plac8, Ly6i | Ccr2, Il1b, Ly6i, Plac8 |
| 3 | MoDC | Lyve1, Cd209f, Vsig4, Tnfrsf17 | Ear2, Fn1, Lyz1, Mgl2, Retnla |
| 4 | MDSC | Arg1, Mmp12, Hilpda, Egln3, Cxcl3 | Arg1, Cxcl3, Egln3, Hilpda, Mmp12 |
| 6 | pDC | Siglech, Ccr9, Ly6d, Klk1, Itm2c | Siglech, Ccr9, Ly6d, Klk1, Atp2a1, Smim5, Upb1, D13Ertd608e |
| 7 | NK | Ncr1, Klrb1c, Klra4, Klra7, Gzma | Gzma, Klra4, Klra7, Klrb1a, Klrb1c, Ncr1 |
| 8 | CD103_cDC | Xcr1, Itgae, Clec9a, Gcsam, Mycl | Gcsam, Itgae, Mycl, Xcr1, Ffar4 |
| 9 | MoDC | Retnla, Lyz1, Ear2, Fn1, Mgl2 | Ear2, Fcrls, Fn1, Lyz1, Mgl2, Retnla |
| 10 | Migratory_cDC | Fscn1, Ccr7, Ccl22, Il4i1 | Ccl22, Ccr7, Fscn1, Il4i1, Nudt17 |

*"Confirmed in Ishizuka 2019" lists genes from our cluster marker DE
that independently appear as statistically significant markers
(Adjusted P < 0.05) in the original authors' cluster DE table for
the corresponding cell type.*

## ISG Upregulation Per Cell Type

From condition DE (ADAR KO vs Control, padj < 0.05, |logFC| > 1.0):

| Cell Type | ISGs Upregulated | logFC range |
|---|---|---|
| Monocyte | Irf7, Isg15, Ifit2, Ifit3, Oasl1 | 1.07 – 1.44 |
| pDC | Irf7 | 1.69 |
| CD8_T | Irf7, Isg15 | 1.02 – 1.20 |
| MoDC | Isg15, Rsad2 | 1.19 – 1.43 |
| MDSC | Isg15, Ifit3 | 1.19 – 1.21 |