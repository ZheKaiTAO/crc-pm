from functools import cache
from pathlib import Path
import re
import gzip as gz
import warnings
import pandas as pd

# settings
GTF_VERSION = "GRCh38.111"
HGNC_PATH = Path(__file__).resolve().parents[2] / "tables/HGNC/hgnc_complete_set.txt"
MAPPING_PATH = Path(__file__).resolve().parents[2]/ "tables/GTF"

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

# def _load_gtf(version:str):
#     # read the preprocessed GTF
#     mapping = pd.read_csv(
#         MAPPING_PATH / version / "gene_id_gene_name_gene_biotype.csv",
#         dtype=str,
#     )
#     # map priority
#     mapping["biotype_priority"] = (
#         mapping["gene_biotype"].map(BIOTYPE_PRIORITY).fillna(9)
#     )
#     return mapping.sort_values("biotype_priority", kind="stable")

def _load_gtf(version:str = GTF_VERSION):
    GTF_PATH = MAPPING_PATH / version / f"Homo_sapiens.{version}.gtf.gz"
    df_dict = dict(
        chromosome= [],
        start = [],
        end = [],
        strand = [],
        gene_id = [],
        gene_name = [],
        gene_biotype = [],
    )
    with gz.open(GTF_PATH, mode='rt', encoding='utf-8') as gtf_file:
        for text in gtf_file:
            if text.startswith("#") or not text.strip():
                continue
            seg = text.rstrip('\n').split('\t')
            if seg[2] != 'gene':
                continue
            df_dict['chromosome'].append(seg[0])
            df_dict['start'].append(int(seg[3]))
            df_dict['end'].append(int(seg[4]))
            df_dict['strand'].append(seg[6])
            attr_dict  = dict(re.findall(r'(\w+)\s+"([^"]*)"', seg[8]))
            df_dict['gene_id'].append(attr_dict['gene_id'])
            df_dict['gene_name'].append(attr_dict.get('gene_name'))
            df_dict['gene_biotype'].append(attr_dict.get('gene_biotype'))

    GTF_DF = pd.DataFrame(df_dict)
    # map priority
    GTF_DF["biotype_priority"] = (
        GTF_DF["gene_biotype"].map(BIOTYPE_PRIORITY).fillna(9)
    )
    return GTF_DF.sort_values("biotype_priority", kind="stable")

@cache
def get_history_mapping(version):
    """
    HGNC current/previous/alias names → ENSG, using GTF biotype priority.
    """
    gtf = _load_gtf(version).drop_duplicates("gene_id")
    id2priority = dict(zip(gtf["gene_id"], gtf["biotype_priority"]))
    # the structure of HGNC file is described here: https://www.genenames.org/help/statistics-and-files/
    hgnc = pd.read_csv(HGNC_PATH, sep="\t", dtype=str, keep_default_na=False)
    # items approved and with valid ENSG IDs
    hgnc = hgnc[
        (hgnc["status"] == "Approved")
        &
        hgnc["ensembl_gene_id"].str.startswith("ENSG")]
    rows = [] # store tuples of (name, gene_id, priority)
    for idx, gene in hgnc.iterrows(): # take the whole row as a tuple (idx, serie)
        gene_id = gene["ensembl_gene_id"]
        for col in ("symbol", "prev_symbol", "alias_symbol"):
            for name in gene[col].split("|"): # some with multiple names
                if name == '':
                    continue
                rows.append((name.strip(), gene_id, id2priority.get(gene_id, 9)))
    history = pd.DataFrame(rows, columns=["name", "gene_id", "priority"])
    history = history.sort_values("priority", kind="stable").drop_duplicates("name")
    return dict(zip(history["name"], history["gene_id"]))

@cache
def _get_gene_mappings(version):
    """Build consistent name/ID lookups using the GTF biotype priority."""
    mapping = _load_gtf(version)
    by_name = mapping.dropna(subset=["gene_name"]).drop_duplicates("gene_name")
    by_id = mapping.drop_duplicates("gene_id")
    return (
        dict(zip(by_name["gene_name"], by_name["gene_id"])),
        {gene_id: name if pd.notna(name) else None
         for gene_id, name in zip(by_id["gene_id"], by_id["gene_name"])},
        {gene_id: biotype if pd.notna(biotype) else None
         for gene_id, biotype in zip(by_id["gene_id"], by_id["gene_biotype"])},
    )


def name2id(genes, Warning=True, version=GTF_VERSION):
    """Return IDs in input order, including duplicates and None for missing names.

    Resolve GTF names first, then HGNC current names, previous names and aliases.
    """
    normal_mapping, _, _ = _get_gene_mappings(version)
    result = []
    missing = []
    for name in genes:
        gene_id = normal_mapping.get(name)
        if gene_id is None and pd.notna(name):
            gene_id = get_history_mapping(version).get(name)
        result.append(gene_id)
        if gene_id is None:
            missing.append(name)
    if Warning and missing:
        print(f"Warning: {len(missing)} names not found: {', '.join(map(str, missing))}")
    return result


def id2name(ids, version=GTF_VERSION):
    """Return names in input order, including duplicates and None for missing IDs."""
    _, mapping, _ = _get_gene_mappings(version)
    return [mapping.get(gene_id) for gene_id in ids]


def id2biotype(ids, version=GTF_VERSION):
    """Return biotypes in input order; IDs absent from this GTF return None."""
    _, _, mapping = _get_gene_mappings(version)
    return [mapping.get(gene_id) for gene_id in ids]


def name2biotype(names, version=GTF_VERSION):
    """Resolve names/aliases to IDs, then return their biotypes in input order.

    Missing names or IDs absent from the selected GTF return None.
    """
    return id2biotype(name2id(names, version=version), version=version)


@cache
def _get_gene_positions(version):
    mapping = _load_gtf(version).drop_duplicates("gene_id", keep="first")
    return {
        row.gene_id: (
            row.chromosome, int(row.start), int(row.end), row.strand
        )
        for row in mapping.itertuples(index=False)
    }


def id2pos(ids: list, version=GTF_VERSION, flavour=None):
    """Return (chr, start, end, strand) tuples in input order.

    Coordinates retain the GTF 1-based, inclusive convention. Unknown IDs
    return None and are reported together in one warning.
    With flavour=None, chromosome names retain their GTF spelling. With
    flavour="chr", names have a single chr prefix and mitochondria use chrM.
    """
    if flavour not in (None, "chr"):
        raise ValueError("flavour must be None or 'chr'")
    positions = _get_gene_positions(version)
    result = []
    missing = {}
    for gene_id in ids:
        position = positions.get(gene_id)
        if position is not None and flavour == "chr":
            chromosome, start, end, strand = position
            chromosome = chromosome.removeprefix("chr")
            chromosome = "chrM" if chromosome in {"MT", "M"} else f"chr{chromosome}"
            position = (chromosome, start, end, strand)
        result.append(position)
        if gene_id not in positions:
            missing[gene_id] = None
    if missing:
        warnings.warn(
            f"{len(missing)} gene IDs not found in {version}: "
            + ", ".join(str(gene_id) for gene_id in missing),
            UserWarning,
            stacklevel=2,
        )
    return result
