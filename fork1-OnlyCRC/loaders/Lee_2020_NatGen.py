from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

import src


DATADIR = Path.home() / "project/crc-pm/data/Lee_2020_NatGen"
GTF_VERSION = "GRCh38.84"


def _prepare(adata, study):
    adata = adata[adata.obs["Class"] != "Border"].copy()
    adata.layers["counts"] = src.ConvDtype.to_int32_counts(adata.X.copy())
    del adata.X

    obs = adata.obs.copy()
    obs = obs.rename(columns={"Cell_type": "OrigCellType", "Class": "Tissue"})
    for column in obs.columns:
        obs[column] = obs[column].astype(str)

    obs["Study"] = study
    obs["CellID"] = obs.index.astype(str)
    obs["Glob_Sample"] = study + "-" + obs["Sample"]
    obs["Glob_Patient"] = study + "-" + obs["Patient"]
    obs["Glob_CellID"] = study + "-" + obs["CellID"]
    obs["OrigCellTypeFull"] = obs["OrigCellType"] + "_" + obs["Cell_subtype"]
    obs["Batch"] = obs["Sample"]
    obs["Assay"] = "scRNA"
    obs["Protocol"] = "3prime"
    obs["Tissue"] = obs["Tissue"].map({"Tumor": "crc", "Normal": "normal"})
    obs = obs.drop(columns=["study", "Cell_subtype"], errors="ignore")
    obs = obs.set_index("Glob_CellID", drop=False)
    adata.obs = obs
    return adata


def load():
    kul3 = sc.read_h5ad(DATADIR / "lee-2020-KUL3.h5ad")
    smc = sc.read_h5ad(DATADIR / "lee-2020-SMC.h5ad")

    kul3 = _prepare(kul3, "Lee_2020_NatGen_KUL3")
    smc = _prepare(smc, "Lee_2020_NatGen_SMC")

    adata = ad.concat(
        [kul3, smc],
        join="outer",
        merge="first",
        uns_merge=None,
    )
    gene_names = adata.var_names.astype(str)
    adata.var = pd.DataFrame({"gene_name": gene_names}, index=gene_names)

    gene2id, gene2biotype = src.GTFMapping.get_gene_mappings(GTF_VERSION)
    adata.var["gene_id"] = adata.var["gene_name"].map(gene2id)
    adata.var["gene_biotype"] = adata.var["gene_name"].map(gene2biotype)
    adata = adata[:, adata.var["gene_id"].notna()].copy()
    adata.var = adata.var.set_index("gene_id", drop=False)
    adata.var.index.name = None
    return adata
