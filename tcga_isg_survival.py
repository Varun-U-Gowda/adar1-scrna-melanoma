"""
tcga_isg_survival.py
────────────────────────────────────────────────────────────────────
Translational validation of the ADAR1-KO interferon-stimulated gene
(ISG) signature in TCGA Skin Cutaneous Melanoma (SKCM) patient data.

Purpose
-------
Scores TCGA SKCM patients on a 6-gene ISG signature derived from
ADAR1-KO scRNA-seq condition DE (see scrna_annotation.py outputs).
Tests whether ISG-high patients have significantly different overall
survival using Kaplan-Meier analysis and log-rank test.

ISG Genes
---------
IRF7, ISG15, IFIT2, IFIT3, OASL, RSAD2
Note: ISG15 is stored under the alias G1P2 in the TCGA RSEM matrix.
Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC8834048/

Workflow
--------
1.  Load TCGA SKCM RNA-seq expression matrix (RSEM)
2.  Extract 6 ISG genes with alias handling (ISG15 stored as G1P2)
3.  Load clinical data
4.  Harmonise patient IDs between expression and clinical data
5.  Compute ISG score: log2(RSEM+1) → per-gene z-score → mean
6.  Median split into ISG High vs ISG Low groups
7.  Parse overall survival (months and event)
8.  Log-rank test
9.  Generate 3-panel diagnostic figure (KM + score distribution + boxplot)
10. Generate KM-only figure for README
11. Save patient-level scores CSV

Inputs
------
tcga_skcm/data_mrna_seq_v2_rsem.txt
    RNA-seq RSEM expression matrix from cBioPortal (skcm_tcga)

skcm_tcga_clinical_data.tsv
    Clinical data including overall survival from cBioPortal (skcm_tcga)

Outputs
-------
outputs/GSE110746_run3/tcga_isg_survival_diagnostics.png/.pdf
    Multi-panel survival diagnostic figure

outputs/GSE110746_run3/tcga_isg_km_survival.png/.pdf
    Kaplan-Meier survival figure

outputs/GSE110746/tcga_isg_patient_scores.csv
    Per-patient ISG score, group assignment, and survival data

Dependencies
------------
pandas, numpy, matplotlib, scipy, lifelines
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy import stats
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings("ignore")

# -- 0. Paths ------------------------------------------------------------------
EXPR_FILE = "/Users/varunugowda/Documents/My_Documents/Github_uploads/Single_cell/tcga_skcm/data_mrna_seq_v2_rsem.txt"
CLINICAL_FILE = "/Users/varunugowda/Documents/My_Documents/Github_uploads/Single_cell/tcga_skcm/skcm_tcga_clinical_data.tsv"
OUT_DIR = "outputs1/GSE110746"

# TCGA stores ISG15 as G1P2 in this RSEM matrix
ISG_GENES_TCGA = {
    "IRF7":  "IRF7",
    "ISG15": "G1P2",
    "IFIT2": "IFIT2",
    "IFIT3": "IFIT3",
    "OASL":  "OASL",
    "RSAD2": "RSAD2",
}

COLORS = {"High": "#d62728", "Low": "#1f77b4"}

# -- 1. Load & clean expression ------------------------------------------------
print("Loading expression matrix ...")
expr_raw = pd.read_csv(EXPR_FILE, sep="\t", index_col=0)

if expr_raw.columns[0] == "Entrez_Gene_Id":
    expr_raw = expr_raw.drop(columns=["Entrez_Gene_Id"])
    print("  Dropped Entrez_Gene_Id column")

n_before = len(expr_raw)
expr_raw = expr_raw[~expr_raw.index.duplicated(keep="first")]
print(f"  Dropped {n_before - len(expr_raw)} duplicate gene rows")
print(f"  Expression: {expr_raw.shape[0]} genes x {expr_raw.shape[1]} samples")

# -- 2. Extract ISG genes with TCGA alias handling -----------------------------
found_map = {
    display_name: tcga_name
    for display_name, tcga_name in ISG_GENES_TCGA.items()
    if tcga_name in expr_raw.index
}

missing_map = {
    display_name: tcga_name
    for display_name, tcga_name in ISG_GENES_TCGA.items()
    if tcga_name not in expr_raw.index
}

print(f"  ISG genes found  : {list(found_map.keys())}")
if missing_map:
    print(f"  Genes NOT found  : {missing_map}")

assert len(found_map) >= 3, "Too few ISG genes found — check gene symbol format"

# Extract using TCGA names, then rename to standard gene symbols
isg_expr = expr_raw.loc[list(found_map.values())].copy()
isg_expr.index = list(found_map.keys())
found = list(found_map.keys())

# -- 3. Load clinical ----------------------------------------------------------
print("\nLoading clinical data ...")
clin = pd.read_csv(CLINICAL_FILE, sep="\t")
clin.columns = clin.columns.str.strip()

OS_MONTH_COL = "Overall Survival (Months)"
OS_STATUS_COL = "Overall Survival Status"
PT_ID_COL = "Patient ID"

print(f"  Clinical: {clin.shape[0]} rows x {clin.shape[1]} cols")

# -- 4. Harmonise sample IDs ---------------------------------------------------
expr_cols_map = {col: col[:12] for col in isg_expr.columns}
isg_expr = isg_expr.rename(columns=expr_cols_map)
isg_expr = isg_expr.loc[:, ~isg_expr.columns.duplicated(keep="first")]

clin[PT_ID_COL] = clin[PT_ID_COL].astype(str).str.strip().str[:12]
clin = clin.drop_duplicates(subset=PT_ID_COL)

common_pts = sorted(set(isg_expr.columns) & set(clin[PT_ID_COL]))
print(f"  Overlapping patients: {len(common_pts)}")

assert len(common_pts) >= 50, \
    f"Only {len(common_pts)} overlapping patients — check ID format"

isg_sub = isg_expr[common_pts]
clin_sub = clin.set_index(PT_ID_COL).loc[common_pts].copy()

# -- 5. ISG signature score ----------------------------------------------------
# log2(RSEM+1) -> z-score per gene across patients -> mean across genes
isg_log = np.log2(isg_sub + 1)

isg_z_arr = np.apply_along_axis(
    lambda row: stats.zscore(row, ddof=1),
    axis=1,
    arr=isg_log.values
)

isg_z = pd.DataFrame(
    isg_z_arr,
    index=isg_log.index,
    columns=isg_log.columns
)

isg_score = isg_z.mean(axis=0)
clin_sub["ISG_score"] = isg_score

print(
    f"\n  ISG score  min={isg_score.min():.3f}  "
    f"median={isg_score.median():.3f}  max={isg_score.max():.3f}"
)

# -- 6. Median split -----------------------------------------------------------
median_score = isg_score.median()
clin_sub["ISG_group"] = np.where(
    clin_sub["ISG_score"] >= median_score,
    "High",
    "Low"
)

print(
    f"  ISG High: {(clin_sub['ISG_group'] == 'High').sum()}   "
    f"ISG Low: {(clin_sub['ISG_group'] == 'Low').sum()}"
)

# -- 7. Parse OS ---------------------------------------------------------------
clin_sub["OS_months"] = pd.to_numeric(clin_sub[OS_MONTH_COL], errors="coerce")

def parse_status(s):
    s = str(s).upper().strip()
    return 1 if ("DECEASED" in s or s in ("1", "1:DECEASED")) else 0

clin_sub["OS_event"] = clin_sub[OS_STATUS_COL].apply(parse_status)
clin_sub = clin_sub.dropna(subset=["OS_months", "OS_event"])

print(
    f"  Patients after dropna: {len(clin_sub)}   "
    f"Events (deaths): {int(clin_sub['OS_event'].sum())}"
)

# -- 8. Log-rank test ----------------------------------------------------------
hi = clin_sub[clin_sub["ISG_group"] == "High"]
lo = clin_sub[clin_sub["ISG_group"] == "Low"]

lr = logrank_test(
    hi["OS_months"],
    lo["OS_months"],
    event_observed_A=hi["OS_event"],
    event_observed_B=lo["OS_event"]
)

pval = lr.p_value
p_str = f"p = {pval:.2e}" if pval < 0.001 else f"p = {pval:.4f}"

print(f"\n  Log-rank p = {pval:.4e}")

# -- 9. Full 3-panel figure ----------------------------------------------------
fig = plt.figure(figsize=(14, 10))
gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.38)

# Panel A: Kaplan-Meier
ax_km = fig.add_subplot(gs[0, :])
kmf = KaplanMeierFitter()

for grp, df_g in [("High", hi), ("Low", lo)]:
    kmf.fit(
        df_g["OS_months"],
        df_g["OS_event"],
        label=f"ISG {grp} (n={len(df_g)})"
    )
    kmf.plot_survival_function(
        ax=ax_km,
        ci_show=True,
        color=COLORS[grp],
        linewidth=2.5
    )

ax_km.set_title(
    f"TCGA SKCM — 6-Gene ISG Signature Kaplan-Meier Overall Survival\n"
    f"Genes: {', '.join(found)} | Log-rank {p_str}",
    fontsize=12,
    fontweight="bold",
    pad=10
)

ax_km.set_xlabel("Time (months)", fontsize=12)
ax_km.set_ylabel("Survival Probability", fontsize=12)
ax_km.set_xlim(left=0)
ax_km.set_ylim(0, 1.05)
ax_km.legend(fontsize=11, loc="upper right")
ax_km.grid(axis="y", linestyle="--", alpha=0.35)

# Panel B: ISG score distribution
ax_dist = fig.add_subplot(gs[1, 0])

for grp, df_g in [("High", hi), ("Low", lo)]:
    ax_dist.hist(
        df_g["ISG_score"],
        bins=25,
        alpha=0.65,
        color=COLORS[grp],
        label=f"ISG {grp}"
    )

ax_dist.axvline(
    median_score,
    color="k",
    linestyle="--",
    linewidth=1.5,
    label=f"Median = {median_score:.2f}"
)

ax_dist.set_xlabel("ISG Score (mean z-score)", fontsize=11)
ax_dist.set_ylabel("Patients", fontsize=11)
ax_dist.set_title("ISG Score Distribution", fontsize=11, fontweight="bold")
ax_dist.legend(fontsize=9)
ax_dist.grid(axis="y", linestyle="--", alpha=0.35)

# Panel C: per-gene boxplot
ax_box = fig.add_subplot(gs[1, 1])

pos_hi = np.arange(len(found)) * 2.5
pos_lo = pos_hi + 0.9

bp_hi = ax_box.boxplot(
    [isg_log.loc[g, hi.index].values for g in found],
    positions=pos_hi,
    widths=0.75,
    patch_artist=True,
    medianprops=dict(color="white", linewidth=2),
    flierprops=dict(marker="o", markersize=2, alpha=0.4)
)

bp_lo = ax_box.boxplot(
    [isg_log.loc[g, lo.index].values for g in found],
    positions=pos_lo,
    widths=0.75,
    patch_artist=True,
    medianprops=dict(color="white", linewidth=2),
    flierprops=dict(marker="o", markersize=2, alpha=0.4)
)

for patch in bp_hi["boxes"]:
    patch.set_facecolor(COLORS["High"])
    patch.set_alpha(0.75)

for patch in bp_lo["boxes"]:
    patch.set_facecolor(COLORS["Low"])
    patch.set_alpha(0.75)

ax_box.set_xticks(pos_hi + 0.45)
ax_box.set_xticklabels(found, fontsize=9, rotation=30, ha="right")
ax_box.set_ylabel("log2(RSEM+1)", fontsize=11)
ax_box.set_title("Per-Gene Expression by ISG Group", fontsize=11, fontweight="bold")
ax_box.legend(
    handles=[
        Patch(color=COLORS["High"], label="ISG High"),
        Patch(color=COLORS["Low"], label="ISG Low")
    ],
    fontsize=9
)
ax_box.grid(axis="y", linestyle="--", alpha=0.35)

for fmt in ("pdf", "png"):
    plt.savefig(
        f"{OUT_DIR}/tcga_isg_survival_diagnostics.{fmt}",
        bbox_inches="tight",
        dpi=150
    )

print("  Saved: tcga_isg_survival_diagnostics.pdf / .png")
plt.close()

# -- 9b. README-friendly Kaplan-Meier-only figure ------------------------------
fig_km, ax_km = plt.subplots(figsize=(10, 5.5))

kmf = KaplanMeierFitter()

for grp, df_g in [("High", hi), ("Low", lo)]:
    kmf.fit(
        df_g["OS_months"],
        df_g["OS_event"],
        label=f"ISG {grp} (n={len(df_g)})"
    )
    kmf.plot_survival_function(
        ax=ax_km,
        ci_show=True,
        color=COLORS[grp],
        linewidth=2.8
    )

ax_km.set_title(
    "TCGA SKCM — ISG Signature and Overall Survival\n"
    f"Median OS: ISG High = {hi['OS_months'].median():.1f} months | "
    f"ISG Low = {lo['OS_months'].median():.1f} months | "
    f"Log-rank {p_str}",
    fontsize=12,
    fontweight="bold",
    pad=10
)

ax_km.set_xlabel("Time (months)", fontsize=12)
ax_km.set_ylabel("Survival Probability", fontsize=12)
ax_km.set_xlim(left=0)
ax_km.set_ylim(0, 1.05)
ax_km.legend(fontsize=11, loc="upper right")
ax_km.grid(axis="y", linestyle="--", alpha=0.35)

plt.tight_layout()

for fmt in ("pdf", "png"):
    plt.savefig(
        f"{OUT_DIR}/tcga_isg_km_survival.{fmt}",
        bbox_inches="tight",
        dpi=200
    )

print("  Saved: tcga_isg_km_survival.pdf / .png")
plt.close()

# -- 10. Save patient results --------------------------------------------------
out_cols = [
    "ISG_score",
    "ISG_group",
    "OS_months",
    "OS_event",
    "Diagnosis Age",
    "Sex",
    "Neoplasm Disease Stage American Joint Committee on Cancer Code"
]

out_cols = [c for c in out_cols if c in clin_sub.columns]

clin_sub[out_cols].to_csv(
    f"{OUT_DIR}/tcga_isg_patient_scores.csv"
)

print("  Saved: tcga_isg_patient_scores.csv")

# -- 11. Summary ---------------------------------------------------------------
print(f"""
== TCGA ISG Survival Summary ================================================
  ISG genes used       : {found}
  TCGA lookup mapping  : {found_map}
  Patients analysed    : {len(clin_sub)}
  ISG High n           : {len(hi)}   median OS = {hi['OS_months'].median():.1f} months
  ISG Low  n           : {len(lo)}   median OS = {lo['OS_months'].median():.1f} months
  Log-rank p-value     : {pval:.4e}
  Outputs              : tcga_isg_survival_diagnostics.pdf / .png
                         tcga_isg_km_survival.pdf / .png
                         tcga_isg_patient_scores.csv
============================================================================
""")
