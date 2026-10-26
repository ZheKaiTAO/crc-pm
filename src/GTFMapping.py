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


def resolve(src, dest, *, version=GTF_VERSION, Warning=True, flavour=None):
    """
    src是一个可迭代的对象, 里面装着ENSG ID(注意不检查是否去除版本号)或gene name
    dest 可选 'id', 'name', 'biotype', 'pos'
    flavour 'chr'返回infercnv能用的格式
    """
    return _resolve(src, dest, version, Warning, flavour, source="auto")


def _resolve(src, dest, version, Warning, flavour=None, *, source):
    if dest not in ("id", "name", "biotype", "pos"):
        raise ValueError("dest must be 'id', 'name', 'biotype', or 'pos'")
    if flavour not in (None, "chr"):
        raise ValueError("flavour must be None or 'chr'")

    normal_mapping = id_names = id_biotypes = None
    if source != "id" or dest != "pos":
        normal_mapping, id_names, id_biotypes = _get_gene_mappings(version)
    positions = _get_gene_positions(version) if dest == "pos" else None
    target = {"name": id_names, "biotype": id_biotypes, "pos": positions}.get(dest)

    result = []
    missing_names = []
    missing_ids = {}
    for value in src:
        is_id = source == "id" or (
            source == "auto" and isinstance(value, str) and value.startswith("ENSG")
        )
        if is_id:
            gene_id = value
        else:
            gene_id = normal_mapping.get(value)
            if gene_id is None and pd.notna(value):
                gene_id = get_history_mapping(version).get(value)
            if gene_id is None:
                missing_names.append(value)

        if dest == "id":
            resolved = gene_id if not is_id or gene_id in id_names else None
        else:
            resolved = target.get(gene_id)
        if dest == "pos" and resolved is not None and flavour == "chr":
            chromosome, start, end, strand = resolved
            chromosome = chromosome.removeprefix("chr")
            chromosome = "chrM" if chromosome in {"MT", "M"} else f"chr{chromosome}"
            resolved = (chromosome, start, end, strand)
        result.append(resolved)

        if source == "auto":
            if gene_id is not None and gene_id not in id_names:
                missing_ids[gene_id] = None
        elif source == "id" and dest == "pos" and gene_id not in positions:
            missing_ids[gene_id] = None

    if Warning and missing_names:
        print(f"Warning: {len(missing_names)} names not found: {', '.join(map(str, missing_names))}")
    if Warning and missing_ids:
        warnings.warn(
            f"{len(missing_ids)} gene IDs not found in {version}: "
            + ", ".join(str(gene_id) for gene_id in missing_ids),
            UserWarning,
            stacklevel=3,
        )
    return result


def name2id(genes, Warning=True, version=GTF_VERSION):
    return _resolve(genes, "id", version, Warning, source="name")


def id2name(ids, version=GTF_VERSION):
    return _resolve(ids, "name", version, False, source="id")


def id2biotype(ids, version=GTF_VERSION):
    return _resolve(ids, "biotype", version, False, source="id")


def name2biotype(names, version=GTF_VERSION):
    return _resolve(names, "biotype", version, True, source="name")


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
    """
    Return (chr, start, end, strand) tuples

    flavour="chr", names have a single chr prefix and mitochondria use chrM.
    """
    return _resolve(ids, "pos", version, True, flavour, source="id")
