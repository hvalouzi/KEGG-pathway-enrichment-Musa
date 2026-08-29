# ============================================================
# KEGG PATHWAY ENRICHMENT ANALYSIS
# ============================================================
#
# Purpose:
#   Perform KEGG pathway enrichment analysis using:
#       1. Fisher's exact test
#       2. Benjamini-Hochberg false discovery rate (BH-FDR)
#
# Statistical universe:
#   Only KEGG genes with at least one KEGG pathway annotation
#   are included in the enrichment universe.
#
# Final statistical comparison:
#   Target KEGG genes with pathway annotation
#       versus
#   Background KEGG genes with pathway annotation
#
# Statistical test:
#   One-sided Fisher's exact test (alternative = "greater")
#
# Multiple-testing correction:
#   Benjamini-Hochberg false discovery rate (BH-FDR)
#
# Requirements:
#   Python 3.x
#   pandas
#   numpy
#   scipy
#   openpyxl
#
# Repository structure:
#   KEGG-analysis/
#   ├── input/
#   │   ├── KEGG_background_final_enrichment_genes.xlsx
#   │   ├── KEGG_target_final_enrichment_genes.xlsx
#   │   └── KEGG_mus_gene_pathway_mapping.xlsx
#   ├── results/
#   └── scripts/
#       └── 02_kegg_enrichment.py
#
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

# ============================================================
# PROJECT DIRECTORIES
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent

INPUT_DIR = PROJECT_DIR / "input"
RESULTS_DIR = PROJECT_DIR / "results"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# INPUT FILES
# ============================================================

BACKGROUND_FILE = INPUT_DIR / "KEGG_background_final_enrichment_genes.xlsx"
TARGET_FILE = INPUT_DIR / "KEGG_target_final_enrichment_genes.xlsx"
PATHWAY_FILE = INPUT_DIR / "KEGG_mus_gene_pathway_mapping.xlsx"

# ============================================================
# OUTPUT FILE
# ============================================================

OUTPUT_FILE = RESULTS_DIR / "KEGG_enrichment_results.xlsx"

# ============================================================
# CHECK INPUT FILES
# ============================================================

required_files = [
    BACKGROUND_FILE,
    TARGET_FILE,
    PATHWAY_FILE
]

for file_path in required_files:
    if not file_path.exists():
        raise FileNotFoundError(
            f"Required input file was not found:\n{file_path}"
        )

# ============================================================
# LOAD DATA
# ============================================================

print("=" * 80)
print("KEGG PATHWAY ENRICHMENT ANALYSIS")
print("=" * 80)

print("\nProject directory:")
print(PROJECT_DIR)

print("\nInput directory:")
print(INPUT_DIR)

print("\nResults directory:")
print(RESULTS_DIR)

background_df = pd.read_excel(BACKGROUND_FILE)
target_df = pd.read_excel(TARGET_FILE)
pathway_df = pd.read_excel(PATHWAY_FILE)

# ============================================================
# VALIDATE REQUIRED COLUMNS
# ============================================================

required_background_columns = {"kegg_gene"}
required_target_columns = {"kegg_gene"}
required_pathway_columns = {"gene", "pathway"}

if not required_background_columns.issubset(background_df.columns):
    raise ValueError("Background file is missing the required 'kegg_gene' column.")

if not required_target_columns.issubset(target_df.columns):
    raise ValueError("Target file is missing the required 'kegg_gene' column.")

if not required_pathway_columns.issubset(pathway_df.columns):
    raise ValueError("Pathway file must contain both 'gene' and 'pathway' columns.")

# ============================================================
# CLEAN KEGG GENE IDENTIFIERS
# ============================================================

background_genes = set(
    background_df["kegg_gene"]
    .dropna()
    .astype(str)
    .str.strip()
)

target_genes = set(
    target_df["kegg_gene"]
    .dropna()
    .astype(str)
    .str.strip()
)

# Remove empty and literal "nan" entries
background_genes = {g for g in background_genes if g and g.lower() != "nan"}
target_genes = {g for g in target_genes if g and g.lower() != "nan"}

# ============================================================
# CLEAN PATHWAY DATA
# ============================================================

pathway_df["gene"] = pathway_df["gene"].astype(str).str.strip()
pathway_df["pathway"] = pathway_df["pathway"].astype(str).str.strip()

# Remove empty and literal "nan" entries
pathway_df = pathway_df[
    (pathway_df["gene"] != "") &
    (pathway_df["gene"].str.lower() != "nan") &
    (pathway_df["pathway"] != "") &
    (pathway_df["pathway"].str.lower() != "nan")
].copy()

# Remove duplicate gene-pathway records
pathway_df = pathway_df.drop_duplicates(subset=["gene", "pathway"]).reset_index(drop=True)

# ============================================================
# INITIAL SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("INPUT GENE SETS")
print("=" * 80)

print("\nUnique background KEGG genes:", len(background_genes))
print("Unique target KEGG genes:", len(target_genes))
print("Unique KEGG genes with pathway records:", pathway_df["gene"].nunique())
print("Unique KEGG pathways:", pathway_df["pathway"].nunique())

# ============================================================
# CHECK TARGET ⊂ BACKGROUND
# ============================================================

missing_from_background = target_genes - background_genes

if missing_from_background:
    raise ValueError(
        "The following target KEGG genes are missing from the background gene set:\n"
        + "\n".join(sorted(missing_from_background))
    )

# ============================================================
# FINAL STATISTICAL UNIVERSE
# ============================================================

pathway_genes = set(pathway_df["gene"])
background_universe = background_genes & pathway_genes
target_universe = target_genes & pathway_genes

print("\n" + "=" * 80)
print("FINAL STATISTICAL UNIVERSE")
print("=" * 80)

print("\nTarget KEGG genes used for enrichment:", len(target_universe))
print("Background KEGG genes used for enrichment:", len(background_universe))

# ============================================================
# VALIDATE STATISTICAL UNIVERSE
# ============================================================

if len(target_universe) == 0:
    raise ValueError("No target KEGG genes with pathway annotations were found.")

if len(background_universe) == 0:
    raise ValueError("No background KEGG genes with pathway annotations were found.")

if not target_universe.issubset(background_universe):
    raise ValueError("Target statistical universe is not a subset of the background statistical universe.")

# ============================================================
# CREATE PATHWAY → GENE SET DICTIONARY
# ============================================================

pathway_groups = pathway_df.groupby("pathway")["gene"].apply(set).to_dict()

# ============================================================
# FISHER'S EXACT TEST
# ============================================================

results = []

for pathway, genes in pathway_groups.items():
    pathway_background = genes & background_universe
    pathway_target = genes & target_universe

    target_in = len(pathway_target)
    target_out = len(target_universe) - target_in

    background_in = len(pathway_background)
    background_out = len(background_universe) - background_in

    a = target_in
    b = target_out
    c = background_in - target_in
    d = background_out - target_out

    contingency_table = [
        [a, b],
        [c, d]
    ]

    odds_ratio, p_value = fisher_exact(contingency_table, alternative="greater")

    results.append({
        "pathway": pathway,
        "target_genes": target_in,
        "target_total": len(target_universe),
        "background_genes": background_in,
        "background_total": len(background_universe),
        "target_percent": (target_in / len(target_universe)) * 100,
        "background_percent": (background_in / len(background_universe)) * 100,
        "odds_ratio": odds_ratio,
        "p_value": p_value,
        "target_gene_list": "; ".join(sorted(pathway_target))
    })

results_df = pd.DataFrame(results)

# ============================================================
# BENJAMINI-HOCHBERG FDR
# ============================================================

def benjamini_hochberg(p_values):
    """
    Calculate Benjamini-Hochberg adjusted P-values.
    """
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)

    if n == 0:
        return np.array([])

    order = np.argsort(p_values)
    ranked_p = p_values[order]

    adjusted = ranked_p * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)

    q_values = np.empty(n)
    q_values[order] = adjusted

    return q_values

results_df["FDR"] = benjamini_hochberg(results_df["p_value"])
results_df["significant_FDR_0.05"] = results_df["FDR"] < 0.05

# ============================================================
# SORT RESULTS
# ============================================================

results_df = results_df.sort_values(["FDR", "p_value"]).reset_index(drop=True)

fdr_significant = results_df[results_df["FDR"] < 0.05].copy()
nominal_significant = results_df[results_df["p_value"] < 0.05].copy()

# ============================================================
# SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("ENRICHMENT SUMMARY")
print("=" * 80)

print("\nTotal pathways tested:", len(results_df))
print("Nominally significant pathways (P < 0.05):", len(nominal_significant))
print("FDR-significant pathways (FDR < 0.05):", len(fdr_significant))

display_columns = [
    "pathway",
    "target_genes",
    "background_genes",
    "target_percent",
    "background_percent",
    "odds_ratio",
    "p_value",
    "FDR",
    "target_gene_list"
]

print("\nTop pathways:")
print(results_df[display_columns].head(20).to_string(index=False))

# ============================================================
# SAVE RESULTS
# ============================================================

with pd.ExcelWriter(OUTPUT_FILE, engine="openpyxl") as writer:
    results_df[display_columns + ["significant_FDR_0.05"]].to_excel(
        writer, sheet_name="All_pathways", index=False
    )
    fdr_significant[display_columns].to_excel(
        writer, sheet_name="FDR_significant", index=False
    )
    nominal_significant[display_columns].to_excel(
        writer, sheet_name="Pvalue_significant", index=False
    )
    target_df.to_excel(writer, sheet_name="Target_genes", index=False)
    background_df.to_excel(writer, sheet_name="Background_genes", index=False)

    summary_df = pd.DataFrame({
        "Metric": [
            "Target KEGG genes",
            "Target KEGG genes with pathways",
            "Background KEGG genes",
            "Background KEGG genes with pathways",
            "KEGG pathways tested",
            "Nominal P < 0.05",
            "FDR < 0.05"
        ],
        "Value": [
            len(target_genes),
            len(target_universe),
            len(background_genes),
            len(background_universe),
            len(results_df),
            len(nominal_significant),
            len(fdr_significant)
        ]
    })

    summary_df.to_excel(writer, sheet_name="Summary", index=False)

# ============================================================
# FINAL REPORT
# ============================================================

print("\n" + "=" * 80)
print("KEGG ENRICHMENT ANALYSIS FINISHED")
print("=" * 80)

print("\nFinal statistical comparison:")
print("Target KEGG genes with pathway annotation:", len(target_universe))
print("Background KEGG genes with pathway annotation:", len(background_universe))

print("\nNominally significant pathways:", len(nominal_significant))
print("FDR-significant pathways:", len(fdr_significant))

print("\nOutput file:")
print(OUTPUT_FILE)
