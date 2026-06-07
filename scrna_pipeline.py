"""
scRNA-seq Pipeline: Single-Cell RNA-seq Analysis (Scanpy)
─────────────────────────────────────────────────────────
Author: Varun
Dataset: GSE110746 — ADAR1-null vs Control B16 melanoma TME (Ishizuka et al. 2019, Nature)

Purpose:
  End-to-end scRNA-seq analysis workflow:
  QC → normalisation → HVGs → PCA/UMAP → Leiden clustering →
  cluster marker DE → condition DE (ADAR KO vs Control)

Input formats supported:
  1) 10x H5     : filtered_feature_bc_matrix.h5
  2) 10x MTX v3 : folder with matrix.mtx(.gz), barcodes.tsv(.gz), features.tsv(.gz)
  3) 10x MTX v2 : folder with matrix.mtx(.gz), barcodes.tsv(.gz), genes.tsv(.gz) — GEO format
  4) h5ad       : AnnData file

Outputs (written to --outdir):
  scRNA_Report.pdf                           — QC + UMAP + HVG plots
  cell_qc_metrics.csv                        — per-cell QC metrics
  all_clusters_top_DEGs.csv                  — top 10 up/down DEGs per cluster (cluster vs rest)
  condition_DE_within_clusters_all.csv       — all genes tested (ADAR KO vs Control per cluster)
  condition_DE_within_clusters_filtered.csv  — padj<0.05, |logFC|>1.0 filtered
  adata_processed.h5ad                       — processed AnnData with raw + lognorm layer

Key implementation notes:
  - adata.raw and adata.layers['lognorm'] stored before sc.pp.scale()
    ensures DE uses log-normalised values not scaled values (Scanpy 1.10+)
  - layer='lognorm', use_raw=False passed to all rank_genes_groups calls
  - Non-immune clusters excluded before condition DE
  - Condition DE thresholds: padj<0.05, |logFC|>1.0 (standard 2-fold change)
  - Cluster marker DE thresholds: padj<0.05, |logFC|>0.5 (relaxed for immune subtype DE)

Command:
  python scrna_pipeline.py \
    --input /Users/varunugowda/Documents/My_Documents/Github_uploads/Single_cell/data/GSE110746 \
    --input_type 10x_mtx_v2 \
    --outdir outputs/GSE110746

Dependencies: scanpy, pandas, scipy, matplotlib
"""

import os
import argparse
import scanpy as sc
import pandas as pd
from scipy.io import mmread
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from typing import Optional


def load_data(input_path: str, input_type: str) -> sc.AnnData:
    """
    Load scRNA-seq data from:
      - 10x_h5       : Cell Ranger HDF5
      - 10x_mtx      : 10x v3 MTX (features.tsv with 3 columns)
      - 10x_mtx_v2   : 10x v2 MTX (genes.tsv with 2 columns, common in GEO)
      - h5ad         : AnnData
    """
    if input_type == "10x_h5":
        adata = sc.read_10x_h5(input_path)

    elif input_type == "10x_mtx":
        # 10x v3 format
        adata = sc.read_10x_mtx(input_path, var_names="gene_symbols", cache=True)

    elif input_type == "10x_mtx_v2":
        # Manual loader for 10x v2 (genes.tsv + barcodes.tsv), common in GEO
        X = mmread(f"{input_path}/matrix.mtx.gz").tocsr()

        # GEO often stores genes x cells; AnnData expects cells x genes
        if X.shape[0] > X.shape[1]:
            X = X.T.tocsr()

        barcodes = pd.read_csv(
            f"{input_path}/barcodes.tsv.gz", header=None, sep="\t"
        )[0].astype(str).values

        genes = pd.read_csv(
            f"{input_path}/genes.tsv.gz", header=None, sep="\t"
        )

        gene_ids = genes[0].astype(str).values
        gene_symbols = genes[1].astype(str).values

        adata = sc.AnnData(X)
        adata.obs_names = barcodes
        adata.var_names = gene_symbols
        adata.var["gene_ids"] = gene_ids

    elif input_type == "h5ad":
        adata = sc.read_h5ad(input_path)

    else:
        raise ValueError("input_type must be one of: 10x_h5, 10x_mtx, 10x_mtx_v2, h5ad")

    adata.var_names_make_unique()
    return adata


def add_condition_from_barcodes(adata: sc.AnnData) -> sc.AnnData:
    """
    Parse barcodes like: AAACCTGAGCGCCTCA-ADAR_S1
    Creates:
      - cell_barcode (left part)
      - sample       (right part, e.g., ADAR_S1)
      - condition    (prefix of sample, e.g., ADAR or Control)
    Robust across pandas versions.
    """
    s = pd.Series(adata.obs_names.astype(str), index=adata.obs_names)

    # rsplit into 2 columns (left = barcode, right = sample tag)
    parts = s.str.rsplit("-", n=1, expand=True)

    # If no "-" exists, parts[1] will be NaN
    adata.obs["cell_barcode"] = parts.iloc[:, 0].to_numpy()
    adata.obs["sample"] = parts.iloc[:, 1].to_numpy()

    # condition = prefix before "_"
    cond = pd.Series(adata.obs["sample"]).astype(str).str.split("_", n=1, expand=True).iloc[:, 0]
    # convert "nan" strings back to NA
    cond = cond.replace("nan", pd.NA)

    adata.obs["condition"] = cond.to_numpy()

    return adata


def add_mito_qc(adata: sc.AnnData, mito_prefix: str = "MT-") -> None:
    """Add mitochondrial gene percentage QC."""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mito_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


def basic_filtering(
    adata: sc.AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_mito: float = 15.0,
    max_genes: Optional[int] = 6000,
) -> sc.AnnData:
    """
    Basic scRNA filtering:
      - remove low-quality cells and lowly detected genes
      - remove high-mito cells
      - optionally remove extremely high-gene cells (potential doublets)
    """
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    if "pct_counts_mt" in adata.obs.columns:
        adata = adata[adata.obs["pct_counts_mt"] < max_mito].copy()

    if max_genes is not None and "n_genes_by_counts" in adata.obs.columns:
        adata = adata[adata.obs["n_genes_by_counts"] < max_genes].copy()

    return adata


def condition_de_within_clusters(
    adata: sc.AnnData,
    outdir: str,
    condition_col: str = "condition",
    case: str = "ADAR",
    control: str = "Control",
    min_cells_per_group: int = 30,
    padj_cut: float = 0.05,
    lfc_cut: float = 1.0,
) -> None:
    """
    Condition DE within each Leiden cluster:
      For each cluster: case vs control (e.g., ADAR vs Control)

    Filters results by:
      - pvals_adj < padj_cut
      - |logfoldchanges| > lfc_cut

    Uses lognorm layer (log-normalized pre-scaling values) for DE.

    Writes:
      - condition_DE_within_clusters_all.csv
      - condition_DE_within_clusters_filtered.csv
    """

    if "leiden" not in adata.obs.columns:
        print("[WARN] Leiden clusters not found. Skipping condition DE.")
        return

    if condition_col not in adata.obs.columns:
        print(f"[WARN] {condition_col} not found in adata.obs. Skipping condition DE.")
        return

    if "lognorm" not in adata.layers:
        print("[WARN] lognorm layer not found. Skipping condition DE.")
        return

    results = []

    for clust in sorted(adata.obs["leiden"].unique(), key=lambda x: int(x)):
        ad = adata[adata.obs["leiden"] == clust].copy()

        vc = ad.obs[condition_col].value_counts()
        if case not in vc.index or control not in vc.index:
            continue
        if vc[case] < min_cells_per_group or vc[control] < min_cells_per_group:
            continue

        sc.tl.rank_genes_groups(
            ad,
            groupby=condition_col,
            reference=control,
            method="wilcoxon",
            layer="lognorm",      # FIX — use log-normalized values not scaled
            use_raw=False
        )

        df = sc.get.rank_genes_groups_df(ad, group=case)
        df["leiden"] = clust
        df["n_case"] = int(vc[case])
        df["n_control"] = int(vc[control])
        results.append(df)
        print(f"[INFO] Cluster {clust}: ADAR={vc[case]}, Control={vc[control]}")

    if not results:
        print("[WARN] No clusters met thresholds for condition DE.")
        return

    # Combine all clusters
    cond_de = pd.concat(results, ignore_index=True)

    # Save unfiltered results
    cond_de.to_csv(
        os.path.join(outdir, "condition_DE_within_clusters_all.csv"),
        index=False
    )

    # Apply filters
    cond_de_filt = cond_de[
        (cond_de["pvals_adj"] < padj_cut) &
        (cond_de["logfoldchanges"].abs() > lfc_cut)
    ].copy()

    # Round for readability
    for col in ["logfoldchanges", "pvals", "pvals_adj", "scores"]:
        if col in cond_de_filt.columns:
            cond_de_filt[col] = pd.to_numeric(
                cond_de_filt[col], errors="coerce"
            ).round(3)

    cond_de_filt.to_csv(
        os.path.join(outdir, "condition_DE_within_clusters_filtered.csv"),
        index=False
    )

    print(
        f"[INFO] Condition DE complete | padj < {padj_cut}, |logFC| > {lfc_cut} "
        f"→ {cond_de_filt.shape[0]} DEGs"
    )


def run_pipeline(
    adata: sc.AnnData,
    outdir: str,
    n_top_genes: int = 2000,
    n_pcs: int = 50,
    neighbors_k: int = 15,
    neighbors_pcs: int = 30,
    leiden_res: float = 0.5,
    regress_out: bool = False,
) -> sc.AnnData:
    """
    Run standard scRNA-seq pipeline and write all outputs.

    Steps:
      1. Save QC metrics CSV
      2. Generate PDF report (QC violins, scatter plots, HVG plot,
         PCA variance ratio, UMAP by leiden, UMAP by QC metrics)
      3. Normalise (CPM) → log1p → select HVGs
      4. Store log-normalised values in adata.raw and adata.layers['lognorm']
         before scaling — required for correct DE fold changes in Scanpy 1.10+
      5. Optional regression of total_counts and pct_counts_mt
      6. Scale to unit variance (max_value=10)
      7. PCA (n_comps) → neighbours → UMAP → Leiden clustering
      8. Cluster marker DE: each cluster vs rest (Wilcoxon, lognorm layer)
         Saves all_clusters_top_DEGs.csv
      9. Condition DE: ADAR KO vs Control within each immune cluster
         Excludes non-immune clusters before DE
         Saves condition_DE_within_clusters_all.csv and _filtered.csv
      10. Save adata_processed.h5ad

    Parameters
    ----------
    adata : sc.AnnData
        Filtered AnnData object with condition and sample metadata in obs.
    outdir : str
        Output directory. Created if it does not exist.
    n_top_genes : int
        Number of highly variable genes for PCA. Default 2000.
    n_pcs : int
        Number of PCA components to compute. Default 50.
    neighbors_k : int
        Number of nearest neighbours for graph construction. Default 15.
    neighbors_pcs : int
        Number of PCs to use for neighbour graph. Default 30.
    leiden_res : float
        Leiden clustering resolution. Higher = more clusters. Default 0.5.
    regress_out : bool
        If True, regress out total_counts and pct_counts_mt before scaling.
        Recommended only for datasets with strong technical confounding.

    Returns
    -------
    sc.AnnData
        Processed AnnData with UMAP, leiden, raw, and lognorm layer.

    """
    os.makedirs(outdir, exist_ok=True)

    # Save QC table
    adata.obs.to_csv(os.path.join(outdir, "cell_qc_metrics.csv"))

    report_path = os.path.join(outdir, "scRNA_Report.pdf")
    with PdfPages(report_path) as pdf:
        # QC plots
        qc_keys = ["n_genes_by_counts", "total_counts"]
        if "pct_counts_mt" in adata.obs.columns:
            qc_keys.append("pct_counts_mt")

        sc.pl.violin(adata, qc_keys, jitter=0.4, multi_panel=True, show=False)
        pdf.savefig(bbox_inches="tight")
        plt.close()

        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
        pdf.savefig(bbox_inches="tight")
        plt.close()

        if "pct_counts_mt" in adata.obs.columns:
            sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
            pdf.savefig(bbox_inches="tight")
            plt.close()

        # Normalize + log
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # HVGs
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_top_genes)
        sc.pl.highly_variable_genes(adata, show=False)
        pdf.savefig(bbox_inches="tight")
        plt.close()

        # FIX — store log-normalized values before scaling
        adata.raw = adata
        adata.layers["lognorm"] = adata.raw.X.copy()

        # Optional regression
        if regress_out:
            covars = ["total_counts"]
            if "pct_counts_mt" in adata.obs.columns:
                covars.append("pct_counts_mt")
            sc.pp.regress_out(adata, covars)

        # Scale
        sc.pp.scale(adata, max_value=10)

        # PCA
        sc.pp.pca(adata, n_comps=n_pcs)
        sc.pl.pca_variance_ratio(adata, log=True, show=False)
        pdf.savefig(bbox_inches="tight")
        plt.close()

        # Neighbors + UMAP
        sc.pp.neighbors(adata, n_neighbors=neighbors_k, n_pcs=neighbors_pcs)
        sc.tl.umap(adata)

        # Clustering
        sc.tl.leiden(adata, resolution=leiden_res)
        sc.pl.umap(adata, color=["leiden"], show=False)
        pdf.savefig(bbox_inches="tight")
        plt.close()

        qc_umap = ["total_counts", "n_genes_by_counts"]
        if "pct_counts_mt" in adata.obs.columns:
            qc_umap.append("pct_counts_mt")
        sc.pl.umap(adata, color=qc_umap, show=False)
        pdf.savefig(bbox_inches="tight")
        plt.close()

    # ----------------------------
    # A) Cluster marker DE (cluster vs rest)
    # ----------------------------
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        layer="lognorm",
        use_raw=False
    )

    deg_df = sc.get.rank_genes_groups_df(adata, group=None)

    # Consistent thresholds
    LFC_CUT = 0.5
    PADJ_CUT = 0.05
    LFC_MAX = 10.0
    TOP_N = 10

    # Filter once — apply to all clusters
    deg_filtered = deg_df[
        (deg_df["pvals_adj"] < PADJ_CUT) &
        (deg_df["logfoldchanges"].abs() > LFC_CUT) &
        (deg_df["logfoldchanges"].abs() < LFC_MAX)
        ].copy()

    # Get top N up and down per cluster in one operation
    up = (
        deg_filtered[deg_filtered["logfoldchanges"] > LFC_CUT]
        .groupby("group", group_keys=False)
        .apply(lambda x: x.nlargest(TOP_N, "logfoldchanges"))
        .copy()
    )
    up["regulation"] = "up"

    down = (
        deg_filtered[deg_filtered["logfoldchanges"] < -LFC_CUT]
        .groupby("group", group_keys=False)
        .apply(lambda x: x.nsmallest(TOP_N, "logfoldchanges"))
        .copy()
    )
    down["regulation"] = "down"

    # Combine and sort by cluster then logFC
    final_top = (
        pd.concat([up, down], axis=0)
        .sort_values(["group", "logfoldchanges"], ascending=[True, False])
        .reset_index(drop=True)
    )

    # Round numeric columns
    for col in ["logfoldchanges", "pvals", "pvals_adj", "scores"]:
        if col in final_top.columns:
            final_top[col] = pd.to_numeric(
                final_top[col], errors="coerce"
            ).round(3)

    final_top.to_csv(
        os.path.join(outdir, "all_clusters_top_DEGs.csv"),
        index=False
    )

    print(f"[INFO] Cluster marker DE complete → {len(final_top)} marker genes")
    print(final_top.groupby(["group", "regulation"])["names"].count().unstack(fill_value=0))

    # ----------------------------
    # B) Condition DE within clusters
    # FIX — remove non-immune clusters before condition DE
    # ----------------------------
    non_immune_clusters = ["5", "11"]
    adata_immune_pipe = adata[
        ~adata.obs["leiden"].isin(non_immune_clusters)
    ].copy()
    print(f"[INFO] Immune cells for condition DE: {adata_immune_pipe.n_obs}")

    condition_de_within_clusters(
        adata_immune_pipe,
        outdir=outdir,
        case="ADAR",
        control="Control",
        min_cells_per_group=30
    )

    # Save processed AnnData
    adata.write(os.path.join(outdir, "adata_processed.h5ad"))

    return adata


def main():
    parser = argparse.ArgumentParser(description="scRNA-seq analysis pipeline (Scanpy).")
    parser.add_argument("--input", required=True, help="Path to input file or folder.")
    parser.add_argument(
        "--input_type",
        required=True,
        choices=["10x_h5", "10x_mtx", "10x_mtx_v2", "h5ad"],
        help="Input format: 10x_h5, 10x_mtx (v3 folder), 10x_mtx_v2 (v2 GEO folder), or h5ad."
    )
    parser.add_argument("--outdir", default="scRNA_outputs", help="Output directory.")
    parser.add_argument("--mito_prefix", default="MT-", help="Mito gene prefix (human usually MT-).")
    parser.add_argument("--min_genes", type=int, default=200)
    parser.add_argument("--min_cells", type=int, default=3)
    parser.add_argument("--max_mito", type=float, default=15.0)
    parser.add_argument("--max_genes", type=int, default=6000)
    parser.add_argument("--regress_out", action="store_true", help="Regress out total_counts and pct_counts_mt.")
    parser.add_argument("--leiden_res", type=float, default=0.5)
    args = parser.parse_args()

    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=100, facecolor="white")

    adata = load_data(args.input, args.input_type)

    # Add condition/sample from barcodes (works for GSE110746-style barcodes)
    adata = add_condition_from_barcodes(adata)

    # Quick sanity check
    if "condition" in adata.obs.columns:
        print("\n=== CONDITION/SAMPLE CHECK ===")
        print(adata.obs["condition"].value_counts(dropna=False))
        if "sample" in adata.obs.columns:
            print(adata.obs["sample"].value_counts(dropna=False).head(10))

    # QC metrics
    add_mito_qc(adata, mito_prefix=args.mito_prefix)

    # Filtering
    adata = basic_filtering(
        adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mito=args.max_mito,
        max_genes=args.max_genes,
    )

    # Run pipeline
    run_pipeline(
        adata,
        outdir=args.outdir,
        leiden_res=args.leiden_res,
        regress_out=args.regress_out,
    )

    print(f"Done. Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()
