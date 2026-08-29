# ============================================================
# KEGG ENRICHMENT DOT PLOT
# ============================================================
#
# Purpose:
#   Generate a publication-quality dot plot from the final
#   KEGG pathway enrichment results.
#
# Plot:
#   - Top 10 pathways ranked by nominal P-value
#   - X-axis: -log10(nominal P-value)
#   - Dot size: number of target genes
#   - Dashed vertical line: nominal P = 0.05
#
# Statistical note:
#   The figure displays nominal P-values for visualization.
#   FDR-adjusted P-values are retained in the enrichment
#   results and are not represented as significant in this
#   figure.
#
# Requirements:
#   Python 3.x
#   pandas
#   numpy
#   matplotlib
#   openpyxl
#
# Repository structure:
#
#   KEGG-analysis/
#   ├── input/
#   ├── results/
#   │   └── KEGG_enrichment_results.xlsx
#   └── scripts/
#       └── 03_kegg_dotplot.py
#
# Run from repository root:
#
#   python scripts/03_kegg_dotplot.py
#
# ============================================================


from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# PROJECT DIRECTORIES
# ============================================================

PROJECT_DIR = Path.cwd()

RESULTS_DIR = PROJECT_DIR / "results"

RESULTS_DIR.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# INPUT FILE
# ============================================================

INPUT_FILE = (
    RESULTS_DIR /
    "KEGG_enrichment_results.xlsx"
)


# ============================================================
# OUTPUT FILE
# ============================================================

FIGURE_FILE = (
    RESULTS_DIR /
    "KEGG_enrichment_dotplot.png"
)


# ============================================================
# LOAD ENRICHMENT RESULTS
# ============================================================

print("=" * 80)
print("KEGG ENRICHMENT DOT PLOT")
print("=" * 80)

print("\nReading enrichment results:")
print(INPUT_FILE)


if not INPUT_FILE.exists():

    raise FileNotFoundError(
        f"Input file not found:\n{INPUT_FILE}\n\n"
        "Please run 02_kegg_enrichment.py first."
    )


df = pd.read_excel(
    INPUT_FILE,
    sheet_name="All_pathways"
)


# ============================================================
# VALIDATE REQUIRED COLUMNS
# ============================================================

required_columns = {
    "pathway",
    "target_genes",
    "p_value",
    "FDR"
}


missing_columns = (
    required_columns -
    set(df.columns)
)


if missing_columns:

    raise ValueError(
        "The enrichment results are missing "
        "the following required columns:\n"
        + "\n".join(
            sorted(missing_columns)
        )
    )


# ============================================================
# CLEAN DATA
# ============================================================

df["target_genes"] = pd.to_numeric(
    df["target_genes"],
    errors="coerce"
)

df["p_value"] = pd.to_numeric(
    df["p_value"],
    errors="coerce"
)

df["FDR"] = pd.to_numeric(
    df["FDR"],
    errors="coerce"
)


df["pathway"] = (
    df["pathway"]
    .astype(str)
    .str.strip()
)


# ============================================================
# KEEP PATHWAYS CONTAINING TARGET GENES
# ============================================================

plot_df = df[
    df["target_genes"] > 0
].copy()


# Remove missing P-values

plot_df = plot_df[
    plot_df["p_value"].notna()
].copy()


# Remove zero or negative P-values

plot_df = plot_df[
    plot_df["p_value"] > 0
].copy()


# ============================================================
# RANK PATHWAYS BY NOMINAL P-VALUE
# ============================================================

plot_df = (
    plot_df
    .sort_values(
        "p_value",
        ascending=True
    )
    .head(10)
    .copy()
)


if plot_df.empty:

    raise ValueError(
        "No pathways containing target genes "
        "are available for plotting."
    )


# ============================================================
# CALCULATE -LOG10(NOMINAL P-VALUE)
# ============================================================

plot_df["minus_log10_p"] = (
    -np.log10(
        plot_df["p_value"]
    )
)


# ============================================================
# SORT FOR HORIZONTAL PLOT
# ============================================================
#
# The first category appears at the bottom.
# Therefore, ascending -log10(P) places the pathway with
# the smallest P-value at the top of the figure.
#
# ============================================================

plot_df = (
    plot_df
    .sort_values(
        "minus_log10_p",
        ascending=True
    )
    .copy()
)


# ============================================================
# REPORT PATHWAYS INCLUDED IN FIGURE
# ============================================================

print("\n" + "=" * 80)
print("PATHWAYS INCLUDED IN FIGURE")
print("=" * 80)

print(
    plot_df[
        [
            "pathway",
            "target_genes",
            "p_value",
            "FDR"
        ]
    ].to_string(
        index=False
    )
)


# ============================================================
# MARKER SIZE FUNCTION
# ============================================================
#
# Dot size represents the number of target genes.
#
# Square-root scaling prevents pathways containing several
# target genes from becoming disproportionately large.
#
# ============================================================

def calculate_marker_size(
    gene_count
):
    """
    Convert target-gene count to marker area.
    """

    return (
        70 +
        np.sqrt(gene_count) *
        65
    )


sizes = (
    plot_df["target_genes"]
    .apply(
        calculate_marker_size
    )
)


# ============================================================
# CREATE FIGURE
# ============================================================

fig, ax = plt.subplots(
    figsize=(10, 6.5)
)


# ============================================================
# PLOT DOTS
# ============================================================

ax.scatter(

    plot_df["minus_log10_p"],

    plot_df["pathway"],

    s=sizes,

    alpha=0.85,

    linewidth=0.8,

    edgecolors="black"
)


# ============================================================
# NOMINAL P-VALUE THRESHOLD
# ============================================================
#
# P = 0.05 corresponds to:
#
#   -log10(0.05) = 1.301
#
# This is a nominal threshold only.
# It does NOT represent FDR significance.
#
# ============================================================

threshold = (
    -np.log10(0.05)
)


ax.axvline(

    threshold,

    linestyle="--",

    linewidth=1.2
)


# ============================================================
# THRESHOLD LABEL
# ============================================================

ax.text(

    threshold + 0.03,

    0.98,

    "P = 0.05",

    transform=ax.get_xaxis_transform(),

    fontsize=9,

    va="top"
)


# ============================================================
# AXIS LABELS
# ============================================================

ax.set_xlabel(

    r"$-\log_{10}$(nominal P-value)",

    fontsize=12
)


ax.set_ylabel(

    "KEGG pathway",

    fontsize=12
)


# ============================================================
# FIGURE TITLE
# ============================================================

ax.set_title(

    "KEGG Pathway Enrichment",

    fontsize=14,

    pad=12
)


# ============================================================
# TICK SETTINGS
# ============================================================

ax.tick_params(

    axis="both",

    labelsize=10
)


# ============================================================
# REMOVE TOP AND RIGHT SPINES
# ============================================================

ax.spines[
    "top"
].set_visible(False)


ax.spines[
    "right"
].set_visible(False)


# ============================================================
# X-AXIS GRID
# ============================================================

ax.grid(

    axis="x",

    linestyle=":",

    linewidth=0.7,

    alpha=0.5
)


ax.set_axisbelow(
    True
)


# ============================================================
# TARGET-GENE SIZE LEGEND
# ============================================================

unique_counts = sorted(

    plot_df[
        "target_genes"
    ]
    .dropna()
    .astype(int)
    .unique()
)


legend_handles = []


for count in unique_counts:

    marker_size = (
        calculate_marker_size(
            count
        )
    )

    handle = ax.scatter(

        [],

        [],

        s=marker_size,

        edgecolors="black",

        linewidth=0.8,

        alpha=0.85,

        label=str(count)
    )

    legend_handles.append(
        handle
    )


if legend_handles:

    ax.legend(

        handles=legend_handles,

        title="Target genes",

        loc="lower right",

        frameon=False,

        fontsize=9,

        title_fontsize=10,

        scatterpoints=1
    )


# ============================================================
# LAYOUT
# ============================================================

plt.tight_layout()


# ============================================================
# SAVE HIGH-RESOLUTION FIGURE
# ============================================================

plt.savefig(

    FIGURE_FILE,

    dpi=600,

    bbox_inches="tight"
)


plt.show()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("DOT PLOT FINISHED")
print("=" * 80)

print(
    "\nPathways plotted:",
    len(plot_df)
)


print(
    "Nominally significant pathways among plotted:",
    int(
        (
            plot_df["p_value"] < 0.05
        ).sum()
    )
)


print(
    "FDR-significant pathways among plotted:",
    int(
        (
            plot_df["FDR"] < 0.05
        ).sum()
    )
)


print(
    "\nImportant:"
)


print(
    "The figure displays nominal P-values. "
    "FDR significance is assessed separately."
)


print(
    "\nFigure saved:"
)


print(
    FIGURE_FILE
)
