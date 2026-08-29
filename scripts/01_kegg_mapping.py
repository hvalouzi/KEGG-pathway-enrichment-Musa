# ============================================================
# KEGG MAPPING AND PATHWAY ANNOTATION
# ============================================================
#
# Purpose:
# 1. Map UniProt protein identifiers to KEGG gene identifiers
# 2. Retrieve KEGG pathway annotations
# 3. Retrieve KEGG pathway names
# 4. Generate the mapping tables used for KEGG enrichment
#
# Organism:
# Musa-related KEGG entries (KEGG organism code: mus)
#
# Requirements:
# Python 3.x
# pandas
# requests
# openpyxl
#
# Repository structure:
#
# KEGG-analysis/
# ├── input/
# ├── results/
# └── scripts/
#
# ============================================================

import time
from pathlib import Path

import pandas as pd
import requests

# ============================================================
# PROJECT DIRECTORIES
# ============================================================

# This script is expected to be located in:
# KEGG-analysis/scripts/
# Therefore, the repository root is one level above "scripts".

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

# ============================================================
# OUTPUT FILES
# ============================================================

TARGET_MAPPING_FILE = RESULTS_DIR / "KEGG_target_uniprot_kegg_mapping.xlsx"
BACKGROUND_MAPPING_FILE = RESULTS_DIR / "KEGG_background_uniprot_kegg_mapping.xlsx"
PATHWAY_MAPPING_FILE = RESULTS_DIR / "KEGG_mus_gene_pathway_mapping.xlsx"

# ============================================================
# KEGG / UNIPROT API SETTINGS
# ============================================================

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/search"
KEGG_LINK_URL = "https://rest.kegg.jp/link/pathway/mus"
KEGG_PATHWAY_LIST_URL = "https://rest.kegg.jp/list/pathway/mus"

# ============================================================
# HELPER FUNCTION: CLEAN UNIPROT IDENTIFIER
# ============================================================

def clean_uniprot_id(value):
    """
    Return a cleaned UniProt identifier.
    """
    if pd.isna(value):
        return ""
    return str(value).strip()

# ============================================================
# HELPER FUNCTION: UNIPROT → KEGG MAPPING
# ============================================================

def get_uniprot_kegg_mapping(uniprot_ids, batch_size=100, pause=0.2):
    """
    Map UniProt protein identifiers to KEGG identifiers
    using the UniProt REST API.

    Parameters
    ----------
    uniprot_ids : iterable
        UniProt identifiers.
    batch_size : int
        Number of UniProt identifiers queried per request.
    pause : float
        Pause between API requests in seconds.

    Returns
    -------
    pandas.DataFrame
        UniProt-to-KEGG mapping table.
    """
    uniprot_ids = [clean_uniprot_id(x) for x in uniprot_ids]
    uniprot_ids = sorted(set(x for x in uniprot_ids if x))

    print(f"Unique UniProt IDs to map: {len(uniprot_ids)}")

    mappings = []
    total_batches = (len(uniprot_ids) + batch_size - 1) // batch_size

    for batch_number, start in enumerate(range(0, len(uniprot_ids), batch_size), start=1):
        batch = uniprot_ids[start:start + batch_size]
        query = " OR ".join(f"accession:{x}" for x in batch)

        params = {
            "query": query,
            "format": "tsv",
            "fields": "accession,id,protein_name,database(Kegg)",
            "size": batch_size
        }

        print(f"Processing batch {batch_number}/{total_batches}...")

        response = requests.get(UNIPROT_URL, params=params, timeout=120)
        response.raise_for_status()

        text = response.text.strip()
        if not text:
            time.sleep(pause)
            continue

        lines = text.splitlines()
        if len(lines) < 2:
            time.sleep(pause)
            continue

        header = lines[0].split("\t")

        for line in lines[1:]:
            fields = line.split("\t")
            row = dict(zip(header, fields))

            accession = row.get("Entry", "").strip()
            entry_name = row.get("Entry Name", "").strip()
            protein_name = row.get("Protein names", "").strip()
            kegg_value = row.get("KEGG", "").strip()

            mappings.append({
                "uniprot_id": accession,
                "entry_name": entry_name,
                "protein_name": protein_name,
                "kegg_ids": kegg_value
            })

        time.sleep(pause)

    mapping_df = pd.DataFrame(mappings)

    if mapping_df.empty:
        mapping_df = pd.DataFrame(
            columns=["uniprot_id", "entry_name", "protein_name", "kegg_ids"]
        )

    return mapping_df

# ============================================================
# READ INPUT DATA
# ============================================================

print("=" * 80)
print("KEGG MAPPING AND PATHWAY ANNOTATION")
print("=" * 80)

print("\nInput directory:")
print(INPUT_DIR)

print("\nResults directory:")
print(RESULTS_DIR)

background_df = pd.read_excel(BACKGROUND_FILE)
target_df = pd.read_excel(TARGET_FILE)

# ============================================================
# IDENTIFY UNIPROT COLUMNS
# ============================================================

possible_columns = [
    "uniprot_id",
    "UniProt_ID",
    "uniprot",
    "UniProt"
]

background_uniprot_column = None
target_uniprot_column = None

for column in possible_columns:
    if column in background_df.columns:
        background_uniprot_column = column
        break

for column in possible_columns:
    if column in target_df.columns:
        target_uniprot_column = column
        break

if background_uniprot_column is None:
    raise ValueError(
        "Could not identify the UniProt ID column in the background input file."
    )

if target_uniprot_column is None:
    raise ValueError(
        "Could not identify the UniProt ID column in the target input file."
    )

print("\nBackground UniProt column:", background_uniprot_column)
print("Target UniProt column:", target_uniprot_column)

# ============================================================
# CLEAN UNIPROT IDENTIFIERS
# ============================================================

background_uniprot = (
    background_df[background_uniprot_column]
    .dropna()
    .astype(str)
    .str.strip()
)

target_uniprot = (
    target_df[target_uniprot_column]
    .dropna()
    .astype(str)
    .str.strip()
)

background_uniprot = sorted(set(x for x in background_uniprot if x))
target_uniprot = sorted(set(x for x in target_uniprot if x))

print("\nUnique background UniProt IDs:", len(background_uniprot))
print("Unique target UniProt IDs:", len(target_uniprot))

# ============================================================
# MAP BACKGROUND PROTEINS
# ============================================================

print("\n" + "=" * 80)
print("BACKGROUND UNIPROT → KEGG MAPPING")
print("=" * 80)

background_mapping = get_uniprot_kegg_mapping(background_uniprot)

# ============================================================
# MAP TARGET PROTEINS
# ============================================================

print("\n" + "=" * 80)
print("TARGET UNIPROT → KEGG MAPPING")
print("=" * 80)

target_mapping = get_uniprot_kegg_mapping(target_uniprot)

# ============================================================
# SAVE UNIPROT → KEGG MAPPINGS
# ============================================================

background_mapping.to_excel(BACKGROUND_MAPPING_FILE, index=False)
target_mapping.to_excel(TARGET_MAPPING_FILE, index=False)

print("\nBackground mapping saved:")
print(BACKGROUND_MAPPING_FILE)

print("\nTarget mapping saved:")
print(TARGET_MAPPING_FILE)

# ============================================================
# HELPER FUNCTION: EXTRACT UNIQUE KEGG GENE IDS
# ============================================================

def extract_kegg_ids(mapping_df):
    """
    Extract individual KEGG gene identifiers from the
    semicolon-separated 'kegg_ids' column.
    """
    kegg_ids = []

    for _, row in mapping_df.iterrows():
        value = row.get("kegg_ids", "")
        if pd.isna(value):
            continue

        value = str(value).strip()
        if not value:
            continue

        for kegg_id in value.split(";"):
            kegg_id = kegg_id.strip()
            if not kegg_id:
                continue
            kegg_ids.append(kegg_id)

    return sorted(set(kegg_ids))

background_kegg_ids = extract_kegg_ids(background_mapping)
target_kegg_ids = extract_kegg_ids(target_mapping)

print("\nUnique background KEGG genes:", len(background_kegg_ids))
print("Unique target KEGG genes:", len(target_kegg_ids))

# ============================================================
# RETRIEVE KEGG GENE → PATHWAY LINKS
# ============================================================

print("\n" + "=" * 80)
print("RETRIEVING KEGG PATHWAY ANNOTATIONS")
print("=" * 80)

response = requests.get(KEGG_LINK_URL, timeout=120)
response.raise_for_status()

pathway_rows = []

for line in response.text.splitlines():
    if not line.strip():
        continue

    parts = line.split("\t")
    if len(parts) != 2:
        continue

    gene_id = parts[0].strip()
    pathway_id = parts[1].strip()

    pathway_rows.append({
        "gene": gene_id,
        "pathway": pathway_id
    })

pathway_df = pd.DataFrame(pathway_rows)

# ============================================================
# CLEAN PATHWAY MAPPING
# ============================================================

pathway_df["gene"] = pathway_df["gene"].astype(str).str.strip()
pathway_df["pathway"] = pathway_df["pathway"].astype(str).str.strip()
pathway_df = pathway_df.drop_duplicates().reset_index(drop=True)

# ============================================================
# RETRIEVE KEGG PATHWAY NAMES
# ============================================================

print("\nRetrieving KEGG pathway names...")

response = requests.get(KEGG_PATHWAY_LIST_URL, timeout=120)
response.raise_for_status()

pathway_names = {}

for line in response.text.splitlines():
    if not line.strip():
        continue

    parts = line.split("\t", maxsplit=1)
    if len(parts) != 2:
        continue

    pathway_id = parts[0].strip()
    pathway_name = parts[1].strip()

    pathway_names[pathway_id] = pathway_name

# ============================================================
# ADD PATHWAY NAMES
# ============================================================

pathway_df["pathway_name"] = pathway_df["pathway"].map(pathway_names)

# ============================================================
# SAVE KEGG PATHWAY MAPPING
# ============================================================

pathway_df.to_excel(PATHWAY_MAPPING_FILE, index=False)

# ============================================================
# CALCULATE FINAL ENRICHMENT UNIVERSE
# ============================================================

background_with_pathway = set(background_kegg_ids) & set(pathway_df["gene"])
target_with_pathway = set(target_kegg_ids) & set(pathway_df["gene"])

# ============================================================
# FINAL MAPPING SUMMARY
# ============================================================

print("\n" + "=" * 80)
print("FINAL MAPPING SUMMARY")
print("=" * 80)

print("\nBackground UniProt proteins:", len(background_uniprot))
print("Background unique KEGG genes:", len(background_kegg_ids))
print("Background KEGG genes with pathway annotation:", len(background_with_pathway))

print("\nTarget UniProt proteins:", len(target_uniprot))
print("Target unique KEGG genes:", len(target_kegg_ids))
print("Target KEGG genes with pathway annotation:", len(target_with_pathway))

print("\nTotal KEGG gene → pathway records:", len(pathway_df))
print("Unique KEGG pathways:", pathway_df["pathway"].nunique())

# ============================================================
# FINAL ENRICHMENT UNIVERSE
# ============================================================

print("\n" + "=" * 80)
print("FINAL STATISTICAL UNIVERSE")
print("=" * 80)

print("Target KEGG genes used for enrichment:", len(target_with_pathway))
print("Background KEGG genes used for enrichment:", len(background_with_pathway))

print("\nExpected final statistical comparison:")
print(f"{len(target_with_pathway)} target KEGG genes vs {len(background_with_pathway)} background KEGG genes")

# ============================================================
# OUTPUT FILES
# ============================================================

print("\n" + "=" * 80)
print("MAPPING STAGE COMPLETED")
print("=" * 80)

print("\nOutput files:")
print(TARGET_MAPPING_FILE)
print(BACKGROUND_MAPPING_FILE)
print(PATHWAY_MAPPING_FILE)
