"""从 Immune_compartment.ipynb 提取的最终分群与人工注释。"""

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


CLUSTER_ANNOTATIONS = {
    '0': 'Monocyte',
    '1': 'Macrophage',
    '2': 'Macrophage',
    '3': 'Macrophage',
    '4': 'Monocyte',
    '5': 'Monocyte',
    '6': 'Neutrophil',
    '7': 'Dendritic Cell',
    '8': 'Mast Cell',
    '9': 'Macrophage',
    '10': 'Treg',
    '11': 'CD4+T',
    '12': 'CD8+T',
    '13': 'Ambigous',
    '14': 'Naïve T',
    '15': 'CD4+T',
    '16': 'CD4+T',
    '17': 'CD4+T',
    '18': 'CD8+T',
    '19': 'CD8+T',
    '20': 'CD8+T',
    '21': 'CD4+T',
    '22': 'Ambigous',
    '23': 'Gamma-Delta T',
    '24': 'NK-ILC',
    '25': 'Memory B',
    '26': 'Naïve B',
    '27': 'Memory B',
    '28': 'GC B',
    '29': 'Plasma',
    '30': 'Ambigous',
}


def annotate_immune(adata):
    """按原参数重算 Immune 分群，返回以细胞 ID 索引的群号和注释表。

    不修改输入，不复制完整表达矩阵，不计算 UMAP，也不读写文件。
    群数或群号与原注释字典不符时停止，避免输出缺失注释。
    """
    if not adata.obs_names.is_unique:
        raise ValueError("细胞 ID 必须唯一。")
    mask = adata.obs["Compartment"].eq("Immune").to_numpy()
    if not mask.any():
        raise ValueError("没有 Immune 细胞，无法应用现有分群注释。")
    latent = np.asarray(adata.obsm["X_scVI"][mask]).copy()
    if latent.ndim != 2 or not np.isfinite(latent).all():
        raise ValueError("Immune 的 X_scVI 必须是有限值二维矩阵。")
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
        key_added="leiden_immune_res1.5",
        resolution=1.5,
        flavor="igraph",
        n_iterations=2,
        random_state=0,
    )
    clusters = work.obs["leiden_immune_res1.5"].astype(str)
    observed = set(clusters)
    expected = set(CLUSTER_ANNOTATIONS)
    if observed != expected:
        missing = sorted(expected - observed)
        unmapped = sorted(observed - expected)
        raise ValueError(
            f"Immune 群号与原注释不一致：预期 {len(expected)} 群，"
            f"实际 {len(observed)} 群；缺少群号 {missing}；未映射群号 {unmapped}。"
        )
    result = work.obs[["leiden_immune_res1.5"]].copy()
    result["man_anno"] = pd.Categorical(clusters.map(CLUSTER_ANNOTATIONS))
    return result
