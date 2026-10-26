"""
辅助分群的helper模块
"""
import re
import numpy as np
import pandas as pd
from . import GTFMapping
import scanpy as sc
from anndata import AnnData
from pathlib import Path

def get_markers(DIR, exist_name = None):
    """
    需求注释marker文件格式: csv
    cell_identity,marker
    cell1,marker1
    exist_name是一个list或serie, 往往是adata.var的列, 用于检查分群marker是否存在于数据中
    DIR是一个Path对象
    """
    marker_df = pd.read_csv(DIR)
    marker_dict = {}
    for idx, row in marker_df.iterrows():
        cell_id = row.iloc[0]
        marker = row.iloc[1]
        if exist_name is not None and marker not in list(exist_name):
            continue
        if marker_dict.get(cell_id, None) is None:
            marker_dict[cell_id] = [marker]
        else:
            marker_dict[cell_id].append(marker)
    return marker_dict

def dotplot(adata, marker_dict, groupby, gene_symbols = 'gene_name',layer='log1p',cmap = "Reds"):
    dp = sc.pl.dotplot(
        adata,
        var_names= marker_dict,
        gene_symbols = gene_symbols,
        groupby = groupby,
        layer = layer,
        var_group_rotation=90,  # 顶部文字旋转
        cmap=cmap,
        figsize=(30, 6),
        return_fig=True,
        categories_order=None
    )
    dp.add_totals()
    dp.show()
    return dp

def dotplotplus():
    pass
    # TODO 自动对亚群进行排序, 创建一个阶梯状的图像; 思考是改categories_order还是var_names字典呢?
    # TODO 对于不同的marker组(竖列), 比如都是Fibroblast S1, 统一使用一种颜色,就用GenCmap中提供的各种颜色, 同时删除那个颜色图注, 以免每个颜色都有图注
    # TODO 对于右侧的totals, 能改颜色吗? 探索一下

def find_top_markers(adata, groupby:str, use_cluster:list[str],return_figure = False):
    pass
    # TODO 对于给定的数据, 探索给定列(grouby)给定亚群的(use_cluster(比如['1','4','7']))的高表达Marker
    # 使用sc.tl.rank_genes_groups函数作为核心
    # return_figure时, 返回figure
    # 除此之外, 返回一个dict, {'cluster1':['marker1','maker2']}这样
    # 但是总是要fig.show()


