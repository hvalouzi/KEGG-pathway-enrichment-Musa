# KEGG Pathway Enrichment Analysis of lncRNA-Associated Target Genes

## Overview

This repository contains Python scripts used for KEGG pathway annotation and enrichment analysis of genes located upstream of lncRNAs identified in *Musa*.

The analysis was performed to characterize the biological pathways associated with the predicted lncRNA-associated target genes.

The workflow includes:

1. Mapping UniProt protein identifiers to KEGG gene identifiers.
2. Retrieving KEGG pathway annotations.
3. Defining the final statistical background and target gene sets.
4. Performing KEGG pathway enrichment using Fisher's exact test.
5. Correcting for multiple testing using the Benjamini-Hochberg false discovery rate (BH-FDR).
6. Generating a publication-quality KEGG enrichment dot plot.

---

## Repository structure

```text
KEGG-analysis/
│
├── README.md
│
├── scripts/
│   ├── 01_kegg_mapping.py
│   ├── 02_kegg_enrichment.py
│   └── 03_kegg_dotplot.py
│
├── input/
│   ├── KEGG_background_final_enrichment_genes.xlsx
│   ├── KEGG_target_final_enrichment_genes.xlsx
│   └── KEGG_mus_gene_pathway_mapping.xlsx
│
└── results/
    ├── KEGG_enrichment_results.xlsx
    └── KEGG_enrichment_dotplot.png
```

---

## Analysis workflow

### Step 1 — KEGG mapping and pathway annotation

`01_kegg_mapping.py`

This script:

* reads the target and background UniProt protein identifiers;
* queries the UniProt REST API;
* retrieves corresponding KEGG gene identifiers;
* retrieves KEGG gene-to-pathway relationships using the KEGG REST API;
* retrieves KEGG pathway names;
* generates the mapping tables required for downstream enrichment analysis.

The script produces UniProt-to-KEGG mapping files and a KEGG gene-to-pathway annotation table.

---

### Step 2 — KEGG pathway enrichment

`02_kegg_enrichment.py`

KEGG pathway enrichment is performed using a one-sided Fisher's exact test.

For each pathway, a 2 × 2 contingency table is constructed:

|                  | In pathway | Not in pathway |
| ---------------- | ---------: | -------------: |
| Target genes     |          a |              b |
| Background genes |          c |              d |

The alternative hypothesis is enrichment of target genes in the pathway.

The resulting P-values are corrected for multiple testing using the Benjamini-Hochberg procedure.

A pathway is considered FDR-significant when:

```text
FDR < 0.05
```

Nominal enrichment is reported separately using:

```text
P < 0.05
```

---

## Statistical universe

Genes without KEGG pathway annotations were excluded from the final enrichment universe.

Therefore, the statistical analysis was performed using only KEGG genes with pathway annotations.

The final statistical comparison was:

```text
Target KEGG genes with pathway annotation:       6
Background KEGG genes with pathway annotation:   2,014
```

Thus, the enrichment analysis was performed using:

**6 target KEGG genes versus 2,014 background KEGG genes.**

The larger numbers of mapped KEGG genes before pathway filtering were not used directly as the statistical universe.

---

## Final KEGG mapping summary

The mapping stage produced the following numbers:

| Category                                      | Number |
| --------------------------------------------- | -----: |
| Background UniProt proteins                   | 34,967 |
| Background unique KEGG genes                  |  8,016 |
| Background KEGG genes with pathway annotation |  2,014 |
| Target UniProt proteins                       |     76 |
| Target unique KEGG genes                      |     17 |
| Target KEGG genes with pathway annotation     |      6 |

The remaining target KEGG genes did not have a corresponding pathway annotation and were therefore excluded from the enrichment universe.

---

## Enrichment results

The final enrichment analysis tested the KEGG pathways represented in the background statistical universe.

The analysis identified:

* **2 pathways with nominal P < 0.05**
* **0 pathways with FDR < 0.05**

The two pathways with nominal P < 0.05 were:

| KEGG pathway            | P-value |    FDR |
| ----------------------- | ------: | -----: |
| Fatty acid biosynthesis |  0.0353 | 1.0000 |
| Propanoate metabolism   |  0.0439 | 1.0000 |

Therefore, these pathways were **nominally significant but did not remain significant after multiple-testing correction**.

No pathway met the predefined FDR threshold of 0.05.

---

## Dot plot

`03_kegg_dotplot.py`

The dot plot displays the ten pathways with the lowest nominal P-values among pathways containing at least one target gene.

Plot characteristics:

* X-axis: `−log10(nominal P-value)`
* Y-axis: KEGG pathway
* Dot size: number of target genes
* Dashed vertical line: nominal `P = 0.05`

The figure is intended to visualize the nominal enrichment results.

Importantly, the figure does **not** indicate FDR significance.

FDR-adjusted results should be interpreted from the enrichment results table.

---

## Requirements

The analysis requires Python 3.x and the following packages:

```text
pandas
numpy
scipy
matplotlib
requests
openpyxl
```

The packages can be installed using:

```bash
pip install pandas numpy scipy matplotlib requests openpyxl
```

---

## Running the analysis

Run the scripts from the repository root.

### 1. KEGG mapping

```bash
python scripts/01_kegg_mapping.py
```

This retrieves UniProt-to-KEGG mappings and KEGG pathway annotations.

### 2. KEGG enrichment

```bash
python scripts/02_kegg_enrichment.py
```

This performs Fisher's exact tests and BH-FDR correction.

### 3. Generate the dot plot

```bash
python scripts/03_kegg_dotplot.py
```

The resulting figure is saved in the `results/` directory.

---

## Reproducibility

The scripts use publicly accessible UniProt and KEGG REST APIs for identifier mapping and pathway annotation.

The enrichment analysis itself is performed locally using Python and does not depend on a graphical enrichment-analysis platform.

The statistical procedure consists of:

1. defining the target and background KEGG gene sets;
2. restricting both sets to genes with KEGG pathway annotations;
3. constructing pathway-specific contingency tables;
4. applying a one-sided Fisher's exact test;
5. applying Benjamini-Hochberg FDR correction;
6. reporting both nominal P-values and FDR-adjusted P-values.

---

## Interpretation

Nominal P-values are provided to describe the strongest pathway-level signals observed in the analysis.

However, because no pathway remained significant after BH-FDR correction, the results should not be interpreted as statistically significant KEGG pathway enrichment at an FDR threshold of 0.05.

The nominally lowest P-values may be considered exploratory observations and should be interpreted cautiously, particularly given the small number of target KEGG genes in the final statistical universe.

---

## Software and databases

The workflow uses:

* Python 3.x
* pandas
* NumPy
* SciPy
* Matplotlib
* Requests
* UniProt REST API
* KEGG REST API

The scripts are provided to facilitate transparency and reproducibility of the computational analysis.

---

## Citation

If this repository or the associated scripts are used, please cite the accompanying research article.

The relevant database and software resources should also be cited according to their respective citation requirements.

---

## License

This repository is intended to provide the computational workflow associated with the analysis. A license can be added according to the authors' preferred terms.
