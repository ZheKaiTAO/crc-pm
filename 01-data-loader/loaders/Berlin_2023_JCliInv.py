###
#  Import modules
###

from pathlib import Path
import re
import anndata as ad
import scanpy as sc
import scipy as sp
import src

DATADIR = Path.home() / "project/crc-pm/data/Berlin_2023_JCliInv/GSE234804"
STUDY = "Berlin_2023_JCliInv"
GTF_VERSION = "GRCh38.98"


def load():
    adata_list = []

    for path in sorted(DATADIR.glob("*.h5ad")):
        sample = path.stem.split("_")[1]
        adata = sc.read_h5ad(path)

        # only keep raw counts
        for layer_name in list(adata.layers.keys()):
            if layer_name in {"scale.data", "data"}:
                del adata.layers[layer_name]
                continue
            elif layer_name == "counts":
                adata.layers[layer_name] = sp.sparse.csr_matrix(
                    src.ConvDtype.to_int32_counts(adata.layers[layer_name])
                )

        cell_ids = adata.obs_names.astype(str)
        prefix, number = re.findall(r"^(\w+)(\d+)$", sample)[0]
        tissue = {"PC": "pm", "CRC": "crc"}[prefix]

        adata.obs["Study"] = STUDY
        adata.obs["Sample"] = sample
        adata.obs["CellID"] = cell_ids
        adata.obs["Batch"] = sample
        adata.obs["Assay"] = "scRNA"
        adata.obs["Protocol"] = "3prime"
        adata.obs["Tissue"] = tissue
        adata.obs["Patient"] = "patient-" + number
        adata.obs["Glob_CellID"] = STUDY + "-" + sample + "_" + cell_ids
        adata.obs["Glob_Sample"] = STUDY + "-" + sample
        adata.obs["Glob_Patient"] = STUDY + "-" + adata.obs["Patient"]
        adata.obs = adata.obs.set_index("Glob_CellID", drop=False)

        adata_list.append(adata)

    if not adata_list:
        raise FileNotFoundError(f"No h5ad files found in {DATADIR}")

    adata = ad.concat(
        adata_list,
        join="inner",
        merge="first",
        uns_merge=None,
    )
    adata.obsm = {}
    adata.obsp = {}

    _, gene2biotype = src.GTFMapping.get_gene_mappings(GTF_VERSION)
    adata.var["gene_name"] = adata.var_names.astype(str)
    adata.var = adata.var.rename(columns={"gene_ids": "gene_id"})
    adata.var = adata.var.set_index("gene_id", drop=False)
    adata.var.index.name = None
    adata.var["gene_biotype"] = adata.var["gene_name"].map(gene2biotype)
    return adata
