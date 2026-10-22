from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp

import src


DATADIR = Path.home() / "project/crc-pm/data/Pelka_2021_Cell/pelka-2021.h5ad"
STUDY = "Pelka_2021_Cell"
GTF_VERSION = "GRCh38.84"


def load():
    adata = sc.read_h5ad(DATADIR)
    adata.layers["counts"] = src.ConvDtype.to_int32_counts(adata.X.copy())
    del adata.X

    obs_columns = [
        "batchID",
        "clTopLevel",
        "clMidwayPr",
        "cl295v11SubFull",
        "SPECIMEN_TYPE",
        "PROCESSING_TYPE",
        "SINGLECELL_TYPE",
        "PID",
    ]
    obs = adata.obs[obs_columns].copy()
    # transform categorical to str
    for column in obs_columns:
        obs[column] = obs[column].astype(str)

    obs["CellID"] = obs.index.astype(str)
    rename_dict = {
            "PID": "Patient",
            "batchID": "Sample",
            "SPECIMEN_TYPE": "Tissue",
            "clMidwayPr": "OrigCellType",
        }
    obs = obs.rename(columns=rename_dict)
    obs["Tissue"] = obs["Tissue"].map({"N": "normal", "T": "crc"})
    obs["Study"] = STUDY
    obs["Batch"] = obs["Sample"]
    obs["Assay"] = "scRNA"
    obs["Protocol"] = "3prime"
    obs["Glob_Sample"] = STUDY + "-" + obs["Sample"]
    obs["Glob_Patient"] = STUDY + "-" + obs["Patient"]
    obs["Glob_CellID"] = STUDY + "-" + obs["CellID"]
    obs["OrigCellTypeFull"] = (
        obs["clTopLevel"] + "_" + obs["OrigCellType"] + "_"+ obs["cl295v11SubFull"]
        )
    obs = obs.set_index("Glob_CellID", drop=False)
    adata.obs = obs

    adata = adata[adata.obs["PROCESSING_TYPE"] == "unsorted"].copy()
    adata.obs = adata.obs.drop(
        columns=["clTopLevel", "cl295v11SubFull", "PROCESSING_TYPE","SINGLECELL_TYPE"]
    )

    adata.var["gene_name"] = adata.var_names
    adata.var = adata.var.rename(columns={"gene_ids": "gene_id"})
    adata = adata[:, adata.var["gene_id"].str.startswith("ENSG")].copy()
    adata.var["gene_id"] = adata.var["gene_id"].str.split(".").str[0]
    #####
    # since ENSG ids are duplicated (lucky enough only 28 of them))
    # we need to add them up. The following is written by ChatGPT
    #####
    # map duplicated gene_id -> one column
    codes, unique_gene_ids = pd.factorize(adata.var["gene_id"], sort=False)
    mapping = sp.sparse.csr_matrix(
        (
            np.ones(adata.n_vars, dtype=np.int8),
            (np.arange(adata.n_vars), codes),
        ),
        shape=(adata.n_vars, len(unique_gene_ids)),
    )
    # sum counts of duplicated genes
    counts = adata.layers["counts"] @ mapping
    # construct new var
    var = (
        adata.var.groupby("gene_id", sort=False)
        .first()
        .loc[unique_gene_ids]
    )
    var["gene_id"] = var.index
    var.index.name = None

    adata = ad.AnnData(
        X=None,
        obs=adata.obs.copy(),
        var=var,
    )
    adata.layers["counts"] = counts.astype(np.int32)
    adata.var = adata.var.set_index("gene_id", drop=False)
    adata.var.index.name = None
    adata.var = adata.var.drop(columns=["genome", "feature_types"])

    _, gene2biotype = src.GTFMapping.get_gene_mappings(GTF_VERSION)
    adata.var["gene_biotype"] = adata.var["gene_name"].map(gene2biotype)
    return adata
