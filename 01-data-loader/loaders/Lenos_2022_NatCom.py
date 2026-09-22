from pathlib import Path

import pandas as pd
import scanpy as sc

import src


DATADIR = Path.home() / "project/crc-pm/data/Lenos_2022_NatCom/GSE183916/GSE183916.h5ad"
STUDY = "Lenos_2022_NatCom"
GTF_VERSION = "GRCh38.98"


def load():
    adata = sc.read_h5ad(DATADIR)
    del adata.layers["data"]
    adata.layers["counts"] = src.ConvDtype.to_int32_counts(
        adata.layers["counts"]
    )

    obs = adata.obs[["loc", "pat_id"]].copy()
    obs = obs.rename(columns={"pat_id": "Patient", "loc": "Tissue"})
    obs["Tissue"] = obs["Tissue"].map(
        {"metastasis": "pm", "primary": "crc"}
    )
    obs["Study"] = STUDY
    obs["Sample"] = obs["Patient"] + "_" + obs["Tissue"]
    obs["Batch"] = obs["Sample"]
    obs["Assay"] = "snRNA"
    obs["Protocol"] = "5prime"
    obs["CellID"] = obs.index.astype(str)
    obs["Glob_CellID"] = STUDY + "-" + obs["CellID"]
    obs["Glob_Sample"] = STUDY + "-" + obs["Sample"]
    obs["Glob_Patient"] = STUDY + "-" + obs["Patient"]
    obs = obs.set_index("Glob_CellID", drop=False)
    adata.obs = obs

    gene_names = adata.var_names.astype(str)
    adata.var = pd.DataFrame({"gene_name": gene_names}, index=gene_names)

    gene2id, gene2biotype = src.GTFMapping.get_gene_mappings(GTF_VERSION)
    adata.var["gene_id"] = adata.var["gene_name"].map(gene2id)
    adata.var["gene_biotype"] = adata.var["gene_name"].map(gene2biotype)
    adata = adata[:, adata.var["gene_id"].notna()].copy()
    adata.var = adata.var.set_index("gene_id", drop=False)
    adata.var.index.name = None
    return adata
