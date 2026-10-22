"""从 Stromal_compartment.ipynb 提取的最终分群与人工注释。"""

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


CLUSTER_ANNOTATIONS = {
    '0': 'Endothelium venous',
    '1': 'Fibroblast S3',
    '2': 'Pericyte',
    '3': 'Fibroblast S3',
    '4': 'Fibroblast S3',
    '5': 'Endothelium arterial',
    '6': 'Fibroblast S3',
    '7': 'Endothelium lymphatic',
    '8': 'Schwann Cell',
    '9': 'Fibroblast S2',
    '10': 'Fibroblast S1',
}


def annotate_stromal(adata):
    """按原参数重算 Stromal 分群，返回以细胞 ID 索引的群号和注释表。

    不修改输入，不复制完整表达矩阵，不计算 UMAP，也不读写文件。
    群数或群号与原注释字典不符时停止，避免输出缺失注释。
    """
    if not adata.obs_names.is_unique:
        raise ValueError("细胞 ID 必须唯一。")
    mask = adata.obs["Compartment"].eq("Stromal").to_numpy()
    if not mask.any():
        raise ValueError("没有 Stromal 细胞，无法应用现有分群注释。")
    latent = np.asarray(adata.obsm["X_scVI"][mask]).copy()
    if latent.ndim != 2 or not np.isfinite(latent).all():
        raise ValueError("Stromal 的 X_scVI 必须是有限值二维矩阵。")
    # 只复制细胞索引和 latent，避免复制 counts/log1p 等完整矩阵。
    work = AnnData(obs=pd.DataFrame(index=adata.obs_names[mask].copy()))
    work.obsm["X_scVI"] = latent
    sc.pp.neighbors(
        work,
        use_rep="X_scVI",
        n_neighbors=20,
        metric="cosine",
        random_state=42,
    )
    sc.tl.leiden(
        work,
        key_added="leiden_stromal_res0.25",
        resolution=0.25,
        flavor="igraph",
        n_iterations=2,
        random_state=0,
    )
    clusters = work.obs["leiden_stromal_res0.25"].astype(str)
    observed = set(clusters)
    expected = set(CLUSTER_ANNOTATIONS)
    if observed != expected:
        missing = sorted(expected - observed)
        unmapped = sorted(observed - expected)
        raise ValueError(
            f"Stromal 群号与原注释不一致：预期 {len(expected)} 群，"
            f"实际 {len(observed)} 群；缺少群号 {missing}；未映射群号 {unmapped}。"
        )
    result = work.obs[["leiden_stromal_res0.25"]].copy()
    result["man_anno"] = pd.Categorical(clusters.map(CLUSTER_ANNOTATIONS))
    return result
