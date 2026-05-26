"""
panglodb_markers.py
────────────────────────────────────────────────────────────────────
Extract top mouse immune cell type markers from PanglaoDB for use
in scRNA-seq cell type annotation (scrna_annotation.py).

Purpose
-------
Filters PanglaoDB marker database for mouse immune cell types and
returns the top markers by specificity score. Output informs the
MARKERS dict used in scrna_annotation.py for dotplot and module
scoring.

Filters Applied
---------------
- Species: Mus musculus (contains "Mm")
- Organ: Immune system
- sensitivity_mouse > 0 (gene must be detectable in mouse)
- Top 5 per cell type ranked by specificity_mouse (descending)

Note on Missing Cell Types
--------------------------
T cytotoxic cells, T helper cells, and Myeloid-derived suppressor
cells returned zero results after sensitivity_mouse > 0 filter.
Markers for these populations were supplemented from the original
paper's Extended Data Fig 4a (Ishizuka et al. 2019, same dataset).

Note on Gene Symbol Capitalisation
-----------------------------------
PanglaoDB returns human uppercase symbols (CD8A).
Mouse convention uses first-letter capitalisation only (Cd8a).
Hyphens are preserved correctly: H2-Ab1 not H2-ab1.

Input
-----
PanglaoDB_markers_27_Mar_2020.tsv
    Downloaded from https://panglaodb.se
    Tab-separated file with columns:
    species, official gene symbol, cell type, organ,
    specificity_mouse, sensitivity_mouse, canonical marker

Output
------
Prints filtered marker table to console.
Copy the mouse_symbol column values into the MARKERS dict
in scrna_annotation.py.

Usage
-----
python panglodb_markers.py

Dependencies
------------
pandas
"""

import pandas as pd

df = pd.read_csv("reference_data/PanglaoDB_markers_27_Mar_2020.tsv", sep="\t")

target_cell_types = [
    "Macrophages",
    "Monocytes",
    "Plasmacytoid dendritic cells",
    "Dendritic cells",
    "NK cells",
    "T cytotoxic cells",
    "T helper cells",
    "T regulatory cells",
    "Myeloid-derived suppressor cells",
    "Neutrophils",
]

mouse_immune = df[
    (df["species"].str.contains("Mm", na=False)) &
    (df["cell type"].isin(target_cell_types)) &
    (df["organ"].str.contains("Immune system", na=False)) &
    (df["specificity_mouse"].notna()) &
    (df["sensitivity_mouse"] > 0.0)
].copy()

mouse_immune = mouse_immune.sort_values(
    ["cell type", "specificity_mouse"],
    ascending=[True, False]
)

top_markers = (
    mouse_immune
    .groupby("cell type")
    .head(5)[["cell type", "official gene symbol",
              "specificity_mouse", "sensitivity_mouse"]]
    .reset_index(drop=True)
)

# Proper mouse capitalisation
# Rule: first letter of each word capital, rest lowercase
# But preserve hyphens — H2-Ab1 not H2-ab1
def to_mouse_symbol(gene):
    parts = gene.split("-")
    result = []
    for i, part in enumerate(parts):
        if i == 0:
            result.append(part.capitalize())
        else:
            # After hyphen — capitalise only if it starts with a letter
            if part and part[0].isalpha():
                result.append(part.capitalize())
            else:
                result.append(part)
    return "-".join(result)

top_markers["mouse_symbol"] = top_markers["official gene symbol"].apply(to_mouse_symbol)

print(top_markers[["cell type", "mouse_symbol",
                    "specificity_mouse", "sensitivity_mouse"]].to_string())

print(f"\nUnique cell types found: {top_markers['cell type'].unique()}")
print(f"Total markers: {len(top_markers)}")
