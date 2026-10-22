"""
辅助分群的helper模块; 纯AI生成我还没看

# 找marker 
markers = ha.find_markers(adata, "leiden")

# 从marker中提取 top marker
top = ha.top_markers(markers)
genes = ha.top_markers(markers, as_dict=True)

# 点图看maker
plot = ha.dotplot(adata, genes, "leiden", standard_scale="var")
plot.show()

# 显示cluster的组成
counts, fractions = ha.cluster_composition(adata, "leiden", "sample")
"""

from collections.abc import Mapping
from copy import deepcopy
from numbers import Integral
import re
import warnings

import numpy as np
import pandas as pd

from . import GTFMapping

__all__ = ["find_markers", "top_markers", "dotplot", "cluster_composition"]

def _display_names(index):
    # 将一列index中所有ENSG者转换为gene name
    names = list(index)
    positions = [i for i, name in enumerate(names)
                 if re.fullmatch(r"ENSG\d{11}", name)]
    if positions:
        mapped = GTFMapping.id2name([names[i] for i in positions])
        for i, name in zip(positions, mapped):
            names[i] = name
    return names


def find_markers(adata, groupby, *, layer="log1p", groups=None,
                 reference="rest", method="wilcoxon",
                 corr_method="benjamini-hochberg", key_added=None, **kwargs):
    """
    每群与参照群比较，返回所有基因的长表，并保存 Scanpy 原始结果到 uns。

    默认使用 log1p 层、Wilcoxon、BH 校正；layer=None 使用 X，始终不使用 raw。
    groups=None 表示全部群，也可传单个标签或标签列表；reference='rest'
    表示其余所有细胞（包括未选入 groups 的群），也可指定一个参照群。
    默认存入uns的键为'helpanno_markers_<groupby>'
    method 支持 wilcoxon、t-test、t-test_overestim_var；kwargs 透传
    tie_correct 等检验参数

    不允许截断或屏蔽基因

    返回列：group、names（原始索引）、gene_name（展示名）、scores、
    logfoldchanges（Scanpy 近似 log2FC）、pvals、pvals_adj、pct_nz_group、
    pct_nz_reference。表达比例按所用矩阵的非零值计算，范围为 [0, 1]。

    示例：find_markers(adata, 'leiden', groups=['1', '2'], reference='0')。
    """
    import scanpy as sc
    from anndata import AnnData

    matrix = adata.X if layer is None else adata.layers[layer]
    labels = adata.obs[groupby]
    categories = (labels.cat.remove_unused_categories().cat.categories.tolist()
                  if isinstance(labels.dtype, pd.CategoricalDtype)
                  else list(pd.unique(labels)))
    category_names = list(map(str, categories))
    selected = (category_names if groups is None or
                (isinstance(groups, str) and groups == "all") else
                [str(groups)] if isinstance(groups, (str, Integral)) else
                list(map(str, groups)))
    reference = str(reference)
    if reference != "rest":
        selected = [group for group in selected if group != reference]
    # Scanpy 会将 obs 字符串转为分类变量，使用独立元数据避免修改调用方。
    obs = pd.DataFrame({groupby: pd.Categorical(labels.astype(str),
                                               categories=category_names)},
                       index=adata.obs_names)
    work = AnnData(X=matrix, obs=obs, var=adata.var.copy())
    if layer is not None:
        work.layers[layer] = matrix
    if "log1p" in adata.uns:
        work.uns["log1p"] = deepcopy(adata.uns["log1p"])
    key = key_added if key_added is not None else f"helpanno_markers_{groupby}"
    sc.tl.rank_genes_groups(work, groupby, groups=selected, reference=reference,
                            layer=layer, use_raw=False, method=method,
                            corr_method=corr_method, n_genes=adata.n_vars,
                            pts=True, key_added=key, **kwargs)
    result = work.uns[key]
    display = pd.Series(_display_names(adata.var_names), index=adata.var_names)
    frames = []
    for group in result["names"].dtype.names:
        frame = pd.DataFrame({field: result[field][group] for field in
                              ("names", "scores", "logfoldchanges", "pvals", "pvals_adj")})
        frame.insert(0, "group", group)
        frame.insert(2, "gene_name", frame["names"].map(display))
        frame["pct_nz_group"] = frame["names"].map(result["pts"][group])
        # Scanpy 在指定参照群时不生成 pts_rest，使用参照群的 pts。
        ref_fraction = (result["pts_rest"][group] if reference == "rest"
                        else result["pts"][reference])
        frame["pct_nz_reference"] = frame["names"].map(ref_fraction)
        frames.append(frame)
    table = pd.concat(frames, ignore_index=True)
    adata.uns[key] = result
    return table


def top_markers(markers, n=10, *, pvals_adj=0.05, logfoldchanges=0.5,
                min_pct=0.25, as_dict=False, gene_column="names"):
    """
    按阈值筛选上调基因，再按 scores 降序取每群前 n 个（并列保留原顺序）。

    使用严格的 pvals_adj < 阈值、logfoldchanges > 阈值，以及
    pct_nz_group >= min_pct。默认返回包含全部字段的 DataFrame；
    as_dict=True 返回 {群: [基因]}，gene_column 可选 names 或 gene_name。
    默认字典使用原始索引，以避免展示名歧义；无合格基因的群不进入字典。

    示例：top_markers(markers, n=5, as_dict=True)。
    """
    filtered = markers.loc[(markers["pvals_adj"] < pvals_adj)
                           & (markers["logfoldchanges"] > logfoldchanges)
                           & (markers["pct_nz_group"] >= min_pct)]
    # 分组后排序，保持群首次出现的顺序。
    parts = [group.sort_values("scores", ascending=False, kind="stable").head(n)
             for _, group in filtered.groupby("group", sort=False, observed=True)]
    top = (pd.concat(parts) if parts else filtered).reset_index(drop=True).copy()
    if as_dict:
        return {group: frame[gene_column].tolist()
                for group, frame in top.groupby("group", sort=False, observed=True)}
    return top


def dotplot(adata, markers, groupby, *, layer="log1p", standard_scale=None,
            categories_order=None, show=False, **kwargs):
    """
    绘制并返回一个 marker 的 dotplot

    默认使用 log1p 层，颜色不标准化; 始终不使用 raw; standard_scale='var' 逐基因标准化。
    categories_order 指定全部有细胞的群的顺序；忽略未使用的分类类别。
    其余 kwargs 透传 Scanpy 绘图参数， 如 cmap、figsize、swap_axes、dendrogram。
    
    原始 var_names、obs 和表达矩阵均不修改；图中标签使用展示名。

    示例：dotplot(adata, {'T': ['CD3D', 'CD3E']}, 'leiden').show()。
    """
    import scanpy as sc
    from anndata import AnnData

    matrix = adata.X if layer is None else adata.layers[layer]
    names = _display_names(adata.var_names)
    originals = set(adata.var_names)
    aliases = {}
    for gene_id, name in zip(adata.var_names, names):
        aliases.setdefault(name, []).append(gene_id)

    def resolve(genes):
        if isinstance(genes, str):
            genes = [genes]
        resolved = []
        for gene in genes:
            if gene in originals:
                resolved.append(gene)
            elif gene in aliases:
                if len(aliases[gene]) > 1:
                    raise ValueError(f"展示名 {gene!r} 对应多个索引 {aliases[gene]}，请使用原 ID。")
                resolved.append(aliases[gene][0])
            else:
                warnings.warn(f"缺失 marker {gene!r}，已跳过。", UserWarning, stacklevel=3)
        return list(dict.fromkeys(resolved))

    if isinstance(markers, Mapping):
        resolved = {}
        for label, genes in markers.items():
            hits = resolve(genes)
            if hits:
                resolved[label] = hits
            else:
                warnings.warn(f"marker 分组 {label!r} 为空，已跳过。", UserWarning, stacklevel=2)
        flat = [gene for genes in resolved.values() for gene in genes]
    else:
        resolved = resolve(markers)
        flat = resolved
    if not flat:
        raise ValueError("没有可绘制的 marker。")
    unique = list(dict.fromkeys(flat))
    positions = adata.var_names.get_indexer(unique)
    work = AnnData(X=matrix[:, positions].copy(), obs=adata.obs.copy(),
                   var=adata.var.iloc[positions].copy())
    if isinstance(work.obs[groupby].dtype, pd.CategoricalDtype):
        work.obs[groupby] = work.obs[groupby].cat.remove_unused_categories()
    if layer is not None:
        work.layers[layer] = work.X
    plot = sc.pl.dotplot(work, resolved, groupby, layer=layer, use_raw=False,
                         standard_scale=standard_scale, categories_order=categories_order,
                         return_fig=True, show=False, **kwargs)
    # 在解析完原始 ID 后替换图标签，允许多个原 ID 展示同一个名称。
    display = dict(zip(adata.var_names, names))
    plot.var_names = [display[gene] for gene in plot.var_names]
    plot.dot_color_df.columns = [display[gene] for gene in plot.dot_color_df.columns]
    plot.dot_size_df.columns = [display[gene] for gene in plot.dot_size_df.columns]
    if show:
        plot.show()
    return plot


def cluster_composition(adata, groupby, by, normalize="cluster"):
    """
    返回 (counts, fractions)，行是细胞群 groupby，列是来源 by；不绘图。

    normalize='cluster'：每行除以该群总细胞数，即群内来源组成。
    normalize='source'：每列除以该来源总细胞数，即来源内部的群比例。
    保留分类变量的类别顺序与未使用类别；零分母对应比例为 0。

    示例：counts, fractions = cluster_composition(adata, 'leiden', 'sample')。
    """
    counts = pd.crosstab(adata.obs[groupby], adata.obs[by], dropna=False)
    if normalize == "cluster":
        fractions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
    else:
        fractions = counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1)
    return counts, fractions.fillna(0.0)
