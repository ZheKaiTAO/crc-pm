from functools import cache
from pathlib import Path

import pandas as pd


MAPPING_PATH = Path.home() / "project/crc-pm/tables/GTF"

BIOTYPE_PRIORITY = {
    "protein_coding": 0,
    "lncRNA": 1,
    "miRNA": 2,
    "snRNA": 2,
    "snoRNA": 2,
    "scaRNA": 2,
    "scRNA": 2,
    "ribozyme": 2,
    "vault_RNA": 2,
    "IG_V_gene": 3,
    "IG_D_gene": 3,
    "IG_J_gene": 3,
    "IG_C_gene": 3,
    "TR_V_gene": 3,
    "TR_D_gene": 3,
    "TR_J_gene": 3,
    "TR_C_gene": 3,
    "rRNA": 4,
    "Mt_rRNA": 4,
    "Mt_tRNA": 4,
    "translated_processed_pseudogene": 5,
    "transcribed_unprocessed_pseudogene": 6,
    "transcribed_processed_pseudogene": 6,
    "transcribed_unitary_pseudogene": 6,
    "unprocessed_pseudogene": 7,
    "processed_pseudogene": 7,
    "unitary_pseudogene": 7,
    "rRNA_pseudogene": 7,
    "IG_V_pseudogene": 7,
    "IG_C_pseudogene": 7,
    "IG_J_pseudogene": 7,
    "TR_V_pseudogene": 7,
    "TR_J_pseudogene": 7,
    "TEC": 8,
    "misc_RNA": 8,
}


def get_versions():
    return [p.name for p in MAPPING_PATH.iterdir() if p.is_dir()]

@cache
def get_gene_mappings(version):
    mapping = pd.read_csv(
        MAPPING_PATH / version / "gene_id_gene_name_gene_biotype.csv",
        dtype=str,
    )
    mapping["biotype_priority"] = (
        mapping["gene_biotype"].map(BIOTYPE_PRIORITY).fillna(9)
    )
    mapping = mapping.sort_values("biotype_priority")
    mapping = mapping.drop_duplicates("gene_name", keep="first")

    gene2id = dict(zip(mapping["gene_name"], mapping["gene_id"]))
    gene2biotype = dict(
        zip(mapping["gene_name"], mapping["gene_biotype"])
    )
    return gene2id, gene2biotype
