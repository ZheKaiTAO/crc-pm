"""Data integration: crc01 + crc03 + crc06 -> merged.h5ad."""
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import anndata as ad
import sys
import re

sys.path.insert(0, "..")
import utils as ut

DataRootDir = Path("../../data")
outdir = Path("../../result/curated")
outdir.mkdir(parents=True, exist_ok=True)

# ── crc01 ──────────────────────────────────────────────────────────────
print("Loading crc01...")
adata_crc01 = sc.read_h5ad(DataRootDir / "crc01/HRA003293.h5ad")
adata_crc01.X = ut.to_float32(adata_crc01.layers["data"])
del adata_crc01.layers["data"]
adata_crc01.layers["counts"] = ut.to_int32_counts(adata_crc01.layers["counts"])

crc01_obs = adata_crc01.obs.copy()
crc01_obs["Study"] = "crc01"
crc01_obs = crc01_obs.rename(
    columns={"patient": "Patient", "tissue": "Tissue", "sample": "Sample"}
)
crc01_obs["Glob_CellID"] = crc01_obs["Study"] + "-" + crc01_obs["cells"]
crc01_obs = crc01_obs.drop(
    columns=["orig.ident", "tissuegroup", "tumorgroup", "site", "cells"]
)
crc01_obs["Glob_Sample"] = crc01_obs["Study"] + "-" + crc01_obs["Sample"]
crc01_obs["Glob_Patient"] = crc01_obs["Study"] + "-" + crc01_obs["Patient"]
crc01_obs["Tech"] = "scRNA"
crc01_obs = crc01_obs.set_index("Glob_CellID", drop=False)
adata_crc01.obs = crc01_obs.copy()
crc01_obs_keep = ["Study_cells_lm"]  # (unused var, deliberately not referenced)
del crc01_obs

# Remove liver metastasis
adata_crc01 = adata_crc01[adata_crc01.obs["Tissue"] != "lm"].copy()
adata_crc01.obsm = {}
adata_crc01.obsp = {}
adata_crc01.var = pd.DataFrame(index=adata_crc01.var.index)
print(f"  crc01: {adata_crc01.shape}")

# ── crc03 ──────────────────────────────────────────────────────────────
print("Loading crc03...")
CurrentDir = DataRootDir / "crc03/GSE234804"
FileList = sorted(f for f in os.listdir(CurrentDir) if f.endswith(".h5ad"))
adata_list = []

for file in FileList:
    name = Path(file).stem.split("_")[1]
    adata = sc.read_h5ad(CurrentDir / file)

    for layer_name in list(adata.layers.keys()):
        if layer_name == "scale.data":
            del adata.layers[layer_name]
        elif layer_name == "counts":
            adata.layers[layer_name] = ut.to_int32_counts(adata.layers[layer_name])
        elif layer_name == "data":
            adata.layers[layer_name] = sp.sparse.csr_matrix(
                ut.to_float32(adata.layers[layer_name])
            )

    adata.X = adata.layers["data"]
    del adata.layers["data"]

    adata.obs["Study"] = "crc03"
    adata.obs["Sample"] = name
    adata.obs["Tech"] = "scRNA"
    adata.obs["Glob_CellID"] = "crc03-" + name + "_" + adata.obs.index
    adata.obs["Glob_Sample"] = "crc03-" + name

    prefix, num = re.findall(r"(\w+)(\d+)", name)[0]
    tissue = "pm" if prefix == "PC" else "crc"
    patient = "patient" + num
    adata.obs["Tissue"] = tissue
    adata.obs["Patient"] = patient
    adata.obs["Glob_Patient"] = "crc03-" + adata.obs["Patient"]
    adata.obs = adata.obs.set_index("Glob_CellID", drop=False)

    adata_list.append(adata)

adata_crc03 = ad.concat(adata_list, join="outer", merge="first", uns_merge=None)
adata_crc03.obsm = {}
adata_crc03.obsp = {}
print(f"  crc03: {adata_crc03.shape}")

# ── crc06 ──────────────────────────────────────────────────────────────
print("Loading crc06...")
adata_crc06 = sc.read_h5ad(DataRootDir / "crc06/GSE183916/GSE183916.h5ad")
adata_crc06.X = sp.sparse.csr_matrix(ut.to_float32(adata_crc06.layers["data"]))
del adata_crc06.layers["data"]
adata_crc06.layers["counts"] = ut.to_int32_counts(adata_crc06.layers["counts"])

crc06_obs = adata_crc06.obs.copy()
crc06_obs["Study"] = "crc06"
crc06_obs = crc06_obs.rename(columns={"pat_id": "Patient", "loc": "Tissue"})
crc06_obs["Tissue"] = crc06_obs["Tissue"].map({"metastasis": "pm", "primary": "crc"})
crc06_obs["Sample"] = crc06_obs["Patient"] + "_" + crc06_obs["Tissue"]
crc06_obs["Tech"] = "snRNA"
crc06_obs["Glob_CellID"] = "crc06-" + crc06_obs["Sample"] + "_" + crc06_obs.index
crc06_obs["Glob_Sample"] = "crc06-" + crc06_obs["Sample"]
crc06_obs["Glob_Patient"] = "crc06-" + crc06_obs["Patient"]
crc06_obs = crc06_obs.drop(columns=["orig.ident", "nCount_RNA", "nFeature_RNA"])
crc06_obs = crc06_obs.set_index("Glob_CellID", drop=False)
adata_crc06.obs = crc06_obs.copy()
del crc06_obs

adata_crc06.obsm = {}
adata_crc06.obsp = {}
adata_crc06.var = pd.DataFrame(index=adata_crc06.var.index)
print(f"  crc06: {adata_crc06.shape}")

# ── Merge ──────────────────────────────────────────────────────────────
print("\nMerging...")
adata_crc01.var["gene_name"] = adata_crc01.var_names
adata_crc03.var["gene_name"] = adata_crc03.var_names
adata_crc06.var["gene_name"] = adata_crc06.var_names
adata_crc01.var["gene_id"] = ""
adata_crc03.var = adata_crc03.var.rename(columns={"gene_ids": "old_gene_id"})

adata_merged = ad.concat(
    [adata_crc03, adata_crc06, adata_crc01],
    join="outer",
    merge="first",
    uns_merge=None,
)

# Fix NaN in layers introduced by ad.concat outer join
for layer_name in adata_merged.layers:
    layer = adata_merged.layers[layer_name]
    if sp.sparse.issparse(layer):
        nan_mask = np.isnan(layer.data)
        if nan_mask.any():
            layer.data[nan_mask] = 0
            print(f"  fixed {nan_mask.sum():,} NaN in layers['{layer_name}']")

# ── Gene ID mapping (GRCh38.111) ──────────────────────────────────────
print("Mapping gene IDs...")
mapping_path = Path.home() / "t_rna/database/GRCh38.111/gene_id_gene_name.csv"
gene_mapping = pd.read_csv(mapping_path, dtype=str)
gene_mapping["gene_id"] = gene_mapping["gene_id"].str.replace(r"\.\d+$", "", regex=True)

symbol_to_gene_ids = (
    gene_mapping.dropna(subset=["gene_name", "gene_id"])
    .groupby("gene_name")["gene_id"]
    .agg(lambda ids: ",".join(sorted(ids.unique())))
)

adata_merged.var["gene_id"] = adata_merged.var["gene_name"].map(symbol_to_gene_ids)
mapped = adata_merged.var["gene_id"].notna()
print(f"  mapped to GRCh38.111: {mapped.sum():,} / {adata_merged.n_vars:,} ({mapped.mean():.2%})")

# Gene biotype
biotype_path = Path.home() / "t_rna/database/GRCh38.111/gene_id_gene_name_gene_biotype.csv"
biotype_mapping = pd.read_csv(biotype_path, dtype=str)
biotype_mapping["gene_id"] = biotype_mapping["gene_id"].str.replace(r"\.\d+$", "", regex=True)

symbol_to_biotype = (
    biotype_mapping.groupby("gene_name")["gene_biotype"]
    .agg(lambda b: ",".join(sorted(b.unique())))
)
adata_merged.var["Class"] = adata_merged.var["gene_name"].map(symbol_to_biotype)

# ── Preview ────────────────────────────────────────────────────────────
ut.preview_adata(adata_merged)
print(f"\nX dtype: {adata_merged.X.dtype}")
print(f"counts dtype: {adata_merged.layers['counts'].dtype}")
print(f"\nStudy counts:\n{adata_merged.obs['Study'].value_counts()}")
print(f"\nTissue counts:\n{adata_merged.obs['Tissue'].value_counts()}")

# ── Save ───────────────────────────────────────────────────────────────
export_path = outdir / "merged.h5ad"
adata_merged.write_h5ad(export_path)
print(f"\nSaved to {export_path}")
print("DONE")
