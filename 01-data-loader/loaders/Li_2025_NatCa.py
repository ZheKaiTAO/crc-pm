###
#  Import modules
###
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp
import anndata as ad
from anndata import AnnData
import sys
import re
import src

STUDY = 'Li_2025_NatCa'
GTF_VERSION = 'GRCh38.111'
DATADIR = Path.home()/"project/crc-pm/data/Li_2025_NatCa/HRA003293.h5ad"

def load():
    # remove log data and tranfer to int32
    adata = sc.read_h5ad(DATADIR)
    del adata.layers["data"] 
    del adata.X
    adata.layers["counts"] = src.ConvDtype.to_int32_counts(adata.layers["counts"])

    # modify obs
    obs = adata.obs.copy()
    obs_list = ['sample','cells','tissue','patient','method','celltype']
    obs = obs[obs_list]
    obs['Study'] = STUDY
    rename_dict = {"patient": "Patient","tissue":"Tissue","sample":"Sample","cells":"CellID", "method":"Protocol","celltype":"OrigCellType"}
    obs = obs.rename(columns=rename_dict)
    obs['Batch'] = obs['Sample']
    obs['Assay'] = 'scRNA'
    obs['Protocol'] = obs['Protocol'].map({'end3': '3prime', 'end5': '5prime'})
    obs['Glob_CellID'] = obs['Study'] + '-' + obs['CellID']
    obs['Glob_Sample'] = obs['Study'] + '-' + obs['Sample']
    obs['Glob_Patient'] = obs['Study'] + '-' + obs['Patient']
    obs = obs.set_index("Glob_CellID", drop=False)
    adata.obs = obs.copy()
    del obs
    # remove lm
    lm_mask = adata.obs['Tissue'] == 'lm'
    adata = adata[~lm_mask,:].copy()
    adata.obsm = {}
    adata.obsp = {}
    gene_names = adata.var_names.astype(str)
    adata.var = pd.DataFrame({'gene_name': gene_names}, index=gene_names)

    gene2id, gene2biotype = src.GTFMapping.get_gene_mappings(GTF_VERSION)
    adata.var['gene_id'] = adata.var['gene_name'].map(gene2id)
    adata.var['gene_biotype'] = adata.var['gene_name'].map(gene2biotype)
    adata = adata[:, adata.var['gene_id'].notna()].copy()
    adata.var = adata.var.set_index('gene_id', drop=False)
    adata.var.index.name = None
    return adata
