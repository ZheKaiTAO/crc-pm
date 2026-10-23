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

# using the previous dicts to get ENSG ID
def name2id(genes, fill_NA = False, Warning = True):
    """
    resolve a list of single gene names eg. ['TP53','FABP']
    first we check if it could converted using GTF mapping
    if we encouter some aliases or previous names, check hgnc
    if that doesn't help, we tell it to fuck off and warns about it 
    """
    GTF_DF = _load_gtf(version).drop_duplicates("gene_name", keep="first")
    normal_mapping = dict(zip(GTF_DF["gene_name"], GTF_DF["gene_id"]))
    history_mapping = get_history_mapping(GTF_VERSION)
    rows = []
    notfound_lst = []
    for item in genes:
        gene_id = normal_mapping.get(item, 'GTF not found')
        if gene_id == 'GTF not found':
            gene_id = history_mapping.get(item, 'GTF and HGNC not found')
            if gene_id == 'GTF and HGNC not found':
                notfound_lst.append(item)
                if fill_NA:
                    rows.append([item, 'NA'])
                continue
        rows.append([item, gene_id])
    if Warning and len(notfound_lst) > 0:
        print(f"Warning: {len(notfound_lst)} names not found: {', '.join(notfound_lst)}")
    return pd.DataFrame(rows, columns=['gene_name', 'gene_id'])

def id2name(ids:list, version = GTF_VERSION):
    """
    Transform IDs to gene names, preserving unknown IDs and input order.
    """
    mapping = _load_gtf(version)
    id2name_dict = dict(zip(mapping["gene_id"],mapping["gene_name"]))
    return [id2name_dict.get(gene_id, gene_id) for gene_id in ids]

@cache
def _get_gene_positions(version):
    mapping = _load_gtf(version).drop_duplicates("gene_id", keep="first")
    return {
        row.gene_id: (
            row.chromosome, int(row.start), int(row.end), row.strand
        )
        for row in mapping.itertuples(index=False)
    }


def id2pos(ids: list, version=GTF_VERSION):
    """Return (chr, start, end, strand) tuples in input order.

    Coordinates retain the GTF 1-based, inclusive convention. Unknown IDs
    return None and are reported together in one warning.
    """
    positions = _get_gene_positions(version)
    result = []
    missing = {}
    for gene_id in ids:
        result.append(positions.get(gene_id))
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
