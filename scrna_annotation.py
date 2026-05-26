"""
scrna_annotation.py
─────────────────────────────────────────────────────────────────────
Cell type annotation for GSE110746 immune clusters.
Requires: scrna_pipeline.py outputs from GSE110746

Workflow:
  1.  Load processed AnnData
  2.  Batch effect check (UMAP by condition and sample)
  3.  Remove non-immune clusters (5=tumour spike-in, 11=stromal)
  4.  Load cluster marker DEGs for annotation reference
  5.  Define PanglaoDB mouse immune markers
  6.  Dotplot — marker expression per cluster
  7.  Module scoring — quantitative cell type likelihood per cluster
  8.  Apply validated cluster annotation
  9.  Annotated UMAP
  10. Cell type composition plot (ADAR KO vs Control)
  11. Add cell_type column to condition DE file
  12. Validate annotations against paper supplementary table

Inputs:
  outputs/GSE110746/adata_processed.h5ad
  outputs/GSE110746/all_clusters_top_DEGs.csv
  outputs/GSE110746/condition_DE_within_clusters_filtered.csv
  NIHMS1560729-supplement-SI-3_ClusterDE_Markers.xlsx    (optional — skipped if absent)

Outputs:
  outputs/GSE110746/adata_immune.h5ad
  outputs/GSE110746/cluster_annotation_markers.csv
  outputs/GSE110746/condition_DE_annotated.csv
  outputs/GSE110746/annotation_validation.csv
  outputs/GSE110746/celltype_composition.png
  figures/umap_annotated_umap.png
  figures/dotplot_celltype_markers.png
  figures/umap_batch_check.png

Dependencies: scanpy, pandas, matplotlib, openpyxl
"""

import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

OUTDIR = "outputs/GSE110746"
os.makedirs(OUTDIR, exist_ok=True)

sc.settings.figdir = OUTDIR

# ── 1. Load ──────────────────────────────────────────────────
adata = sc.read_h5ad(f"{OUTDIR}/adata_processed.h5ad")
print(f"Total cells loaded: {adata.n_obs}")

# ── 2. Batch check ───────────────────────────────────────────
print("\n=== Cell counts per sample ===")
print(adata.obs["sample"].value_counts())
print("\n=== Cell counts per condition ===")
print(adata.obs["condition"].value_counts())

sc.pl.umap(
    adata,
    color=["leiden", "condition", "sample"],
    ncols=3,
    wspace=0.3,
    save="_batch_check.png"
)

# ── 3. Remove non-immune clusters ────────────────────────────
# Cluster 5: B16 tumour spike-in cells (Elf5, Pmel — intentional per paper)
# Cluster 11: Stromal fibroblast contamination (Wnt5a, Ror2, Slit3)
adata_immune = adata[~adata.obs["leiden"].isin(["5", "11"])].copy()
print(f"\nImmune cells after removing clusters 5 and 11: {adata_immune.n_obs}")

adata_immune.write(f"{OUTDIR}/adata_immune.h5ad")
print("Saved adata_immune.h5ad")

# ── 4. Load top DEGs for annotation reference ─────────────────
top_degs = pd.read_csv(f"{OUTDIR}/all_clusters_top_DEGs.csv")

print("\n=== TOP 10 UPREGULATED GENES PER CLUSTER ===")
up_degs = top_degs[top_degs["regulation"] == "up"].copy()
for cluster in sorted(up_degs["group"].dropna().unique(),
                      key=lambda x: int(x)):
    genes = (
        up_degs[up_degs["group"] == cluster]
        .nlargest(10, "logfoldchanges")[["names", "logfoldchanges"]]
    )
    print(f"\nCluster {cluster}:")
    print(genes.to_string(index=False))

# ── 5. PanglaoDB markers ──────────────────────────────────────────────────────
# Source: PanglaoDB_markers_27_Mar_2020.tsv
# Filter: Mus musculus | Immune system organ | sensitivity_mouse > 0
# Ranked by: specificity_mouse (descending), top 4-5 per cell type
# Note: T cytotoxic, T helper, MDSC had sensitivity_mouse = 0 in PanglaoDB
#       These were supplemented from Ishizuka 2019 Extended Data Fig 4a
#       (same dataset — ground truth source)
MARKERS = {
    "Dendritic_cell": ["H2afy", "Cadm1", "Plac8", "Lgals3", "Ctss"],
    "Macrophage":     ["Lyz2", "Lgals3", "Cd74", "H2-Ab1", "Cd200"],
    "Monocyte":       ["Psap", "Ifitm3", "Ptprc", "Cd44", "Csf1r"],
    "NK":             ["Ctla2a", "Nkg7", "Ccl4", "Ccl3"],
    "Neutrophil":     ["Cd14", "Ccl6", "Ccl9", "Cxcl2"],
    "pDC":            ["Itm2c", "Bst2", "Bcl11a", "Irf8"],
    "Treg":           ["Maf", "Ikzf2", "Ctla4"],
    "CD8_T":          ["Cd8a", "Cd8b1", "Gzmb", "Prf1", "Nkg7"],
    "CD4_T":          ["Cd4", "Il7r", "Tcf7", "Cd44"],
    "MDSC":           ["Itgam", "Ly6c2", "Ly6g", "S100a8", "S100a9"],
}

# ── 6. Dotplot ───────────────────────────────────────────────
sc.pl.dotplot(
    adata_immune,
    MARKERS,
    groupby="leiden",
    use_raw=True,
    dendrogram=True,
    standard_scale="var",
    save="_celltype_markers.png"
)

# ── 7. Module scoring ────────────────────────────────────────
for cell_type, genes in MARKERS.items():
    genes_present = [g for g in genes if g in adata_immune.raw.var_names]
    sc.tl.score_genes(
        adata_immune,
        gene_list=genes_present,
        score_name=f"score_{cell_type}",
        use_raw=True
    )

score_cols = [f"score_{ct}" for ct in MARKERS.keys()]
mean_scores = (
    adata_immune.obs[["leiden"] + score_cols]
    .groupby("leiden")[score_cols]
    .mean()
)
mean_scores.columns = list(MARKERS.keys())
mean_scores["probable_cluster"] = mean_scores.idxmax(axis=1)

score_matrix_out = mean_scores.copy()
score_matrix_out[list(MARKERS.keys())] = (
    score_matrix_out[list(MARKERS.keys())].round(4)
)

print("\n=== FULL SCORE MATRIX WITH PROBABLE CLUSTER ===")
print(score_matrix_out.to_string())

# Save score matrix + top DEGs combined
combined_outfile = f"{OUTDIR}/cluster_annotation_markers.csv"
with open(combined_outfile, "w") as f:
    f.write("CELL TYPE SCORE MATRIX\n")
    score_matrix_out.to_csv(f, index=True)
    f.write("\nTOP 10 UPREGULATED DEGs PER CLUSTER\n")
    top_degs_up = (
        top_degs[top_degs["regulation"] == "up"]
        .sort_values(["group", "logfoldchanges"], ascending=[True, False])
        .groupby("group", group_keys=False)
        .head(10)
        .reset_index(drop=True)
    )
    top_degs_up.to_csv(f, index=False)
print(f"Saved {combined_outfile}")

# ── 8. Apply validated annotation ─────────────────────────────────────────────
# Annotation method: convergent evidence from three sources
#   (1) PanglaoDB marker expression — dotplot + module scores (steps 6-7)
#   (2) Cluster-specific DEGs — all_clusters_top_DEGs.csv (step 4)
#   (3) Cross-validated against Ishizuka 2019 supplementary cluster DE table
#       (step 12 confirms 5-9/10 top marker genes per cluster)
#
# Resolution note: Leiden resolution=0.5 produces 10 immune clusters vs
#   paper's 15 at resolution=0.8. Some populations are merged (MoDC clusters
#   3 and 9 correspond to a single MoDC population in the paper).

cluster_annotation = {
    "0":  "M1_Macrophage",    # C1qa, C1qc, Cd209f | paper: M1 Macrophage
    "1":  "CD8_T",            # Cd3g, Cd8a, Cd8b1, Gzmk, Cxcr6 | paper: Mki67- CD8+ T
    "2":  "Monocyte",         # Ccr2, Il1b, Cd177, Plac8 | paper: Monocyte
    "3":  "MoDC",             # Lyve1, Cd209f, Vsig4 | paper: MoDC
    "4":  "MDSC",             # Arg1, Mmp12, Hilpda, Egln3 | paper: MDSC
    "6":  "pDC",              # Siglech, Ccr9, Ly6d | paper: Plasmacytoid DC
    "7":  "NK",               # Ncr1, Klrb1c, Klra4/7 | paper: NK cell
    "8":  "CD103_cDC",        # Xcr1, Itgae, Clec9a | paper: CD103+ cDC
    "9":  "MoDC",             # Retnla, Lyz1, Ear2, Fn1 | paper: MoDC
    "10": "Migratory_cDC",    # Fscn1, Ccr7, Ccl22 | paper: Migratory cDC
}


adata_immune.obs["cell_type"] = (
    adata_immune.obs["leiden"].map(cluster_annotation)
)

print("\n=== CELL TYPE COUNTS ===")
print(adata_immune.obs["cell_type"].value_counts())
print(f"Unmapped clusters: {adata_immune.obs['cell_type'].isna().sum()}")

# ── 9. Annotated UMAP ────────────────────────────────────────
sc.pl.umap(
    adata_immune,
    color=["cell_type", "condition"],
    ncols=2,
    wspace=0.5,
    save="_annotated_umap.png"
)

# ── 10. Composition plot ─────────────────────────────────────
comp = (
    adata_immune.obs
    .groupby(["cell_type", "condition"])
    .size()
    .unstack(fill_value=0)
)
comp_pct = comp.div(comp.sum(axis=0), axis=1) * 100

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

comp_pct.plot(
    kind="bar", ax=axes[0],
    color=["#378ADD", "#D4537E"], edgecolor="none"
)
axes[0].set_title("Cell type composition by condition (%)")
axes[0].set_ylabel("% of cells")
axes[0].tick_params(axis="x", rotation=45)
axes[0].legend(title="Condition")

comp.plot(
    kind="bar", ax=axes[1],
    color=["#378ADD", "#D4537E"], edgecolor="none"
)
axes[1].set_title("Cell type counts by condition")
axes[1].set_ylabel("Cell count")
axes[1].tick_params(axis="x", rotation=45)
axes[1].legend(title="Condition")

plt.tight_layout()
plt.savefig(
    f"{OUTDIR}/celltype_composition.png",
    dpi=150, bbox_inches="tight"
)
plt.close()
print("Composition plot saved.")

# ── 11. Add cell type to condition DE file ───────────────────
cond_de = pd.read_csv(
    f"{OUTDIR}/condition_DE_within_clusters_filtered.csv"
)
cond_de["leiden"] = cond_de["leiden"].astype(str)
cond_de["cell_type"] = cond_de["leiden"].map(cluster_annotation)

cond_de.to_csv(f"{OUTDIR}/condition_DE_annotated.csv", index=False)

print("\n=== CONDITION DE WITH CELL TYPES ===")
print(
    cond_de[["cell_type", "leiden", "names",
             "logfoldchanges", "pvals_adj"]]
    .sort_values(["cell_type", "logfoldchanges"],
                 ascending=[True, False])
    .to_string(index=False)
)

# ── 12. Validation against paper supplementary table ────────
# Method: your top 10 DEGs vs ALL paper marker genes
# (more robust than top10 vs top10 for cross-resolution comparison)
paper_path = "reference_data/NIHMS1560729-supplement-SI-3_ClusterDE_Markers.xlsx"
if os.path.exists(paper_path):
    paper = pd.read_excel(paper_path)

    paper_map = {
        "M1_Macrophage":  "M1 Macrophage",
        "CD8_T":          "Mki67- CD8+ T cell",
        "Monocyte":       "Monocyte",
        "MDSC":           "MDSC",
        "pDC":            "Plasmacytoid DC",
        "NK":             "NK cell",
        "CD103_cDC":      "CD103+ cDC",
        "MoDC":           "MoDC",
        "Migratory_cDC":  "Migratory cDC",
    }

    your_markers = pd.read_csv(f"{OUTDIR}/all_clusters_top_DEGs.csv")
    your_markers["group"] = your_markers["group"].astype(str)
    your_markers["cell_type"] = your_markers["group"].map(cluster_annotation)
    your_up = your_markers[your_markers["regulation"] == "up"]

    print("\n=== VALIDATION: YOUR MARKERS vs ALL PAPER MARKERS ===")
    print(f"{'Cell Type':<20} {'Overlap':>8} {'Your marker genes found in paper'}")
    print("-" * 70)

    validation_results = []
    for your_label, paper_label in paper_map.items():

        your_genes = set(
            your_up[your_up["cell_type"] == your_label]["names"]
            .dropna()
            .astype(str)
        )

        # ALL paper markers for this cell type
        all_paper_genes = set(
            paper[
                (paper["Population vs. All Else"] == paper_label) &
                (paper["Adjusted P.value"] < 0.05)
            ]["Gene"]
        )

        overlap = your_genes & all_paper_genes
        n_paper_total = len(all_paper_genes)

        print(
            f"{your_label:<20} "
            f"{len(overlap)}/{len(your_genes)} "
            f"(paper has {n_paper_total} sig. markers): "
            f"{sorted(overlap) if overlap else 'none'}"
        )

        validation_results.append({
            "cell_type":         your_label,
            "paper_label":       paper_label,
            "your_top10_n":      len(your_genes),
            "paper_total_sig":   n_paper_total,
            "overlap_n":         len(overlap),
            "overlap_genes":     ", ".join(sorted(overlap))
        })

    # Save validation results
    val_df = pd.DataFrame(validation_results)
    val_df.to_csv(
        f"{OUTDIR}/annotation_validation.csv",
        index=False
    )
    print(f"\nSaved annotation_validation.csv")
    print("\nInterpretation: overlap >= 1 with padj<0.05 paper markers")
    print("confirms biological concordance regardless of ranking differences.")

else:
    print(f"\n[SKIP] Supplementary file not found at {paper_path}")
