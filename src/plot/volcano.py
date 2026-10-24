"""绘制差异表达火山图：点面积按绝对 log2FoldChange 的幂映射, 显著性越高颜色越深 """

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

def _color_levels(y):
    """
    第 10–90 分位数色阶；分位数范围为零时使用中间色。"""
    lower, upper = np.percentile(y, [10, 90])
    if lower < upper:
        return Normalize(vmin=lower, vmax=upper, clip=True)(y)
    return np.full(y.size, .5)
def _decorate_axes(ax, conditions, lfc_threshold, cutoff, figarea, title):
    """
    添加阈值线
    添加<->箭头装饰
    添加标题
    添加图例
    """
    # 竖直阈值线
    for threshold in (-lfc_threshold, lfc_threshold):
        ax.axvline(threshold, color='dimgrey', linestyle='--', linewidth=1)
    # 横向阈值线
    ax.axhline(cutoff, color='dimgrey', linestyle='--', linewidth=1)
    # 坐标轴
    ax.set(xlabel='log2 Fold Change', ylabel='−log10(adjusted p-value)')
    if figarea is not None:
        ax.set(xlim=figarea[:2], ylim=figarea[2:])

    # <->箭头装饰
    if conditions is not None:
        ax.annotate('', xy=(0.98, 1.04), xytext=(0.02, 1.04),
                xycoords='axes fraction', annotation_clip=False, clip_on=False,
                arrowprops={'arrowstyle': '<->', 'color': 'black', 'lw': 1})
        for xpos, condition, alignment in zip((0, 1), conditions, ('left', 'right')):
            ax.text(xpos, 1.07, f'Up in {condition}', transform=ax.transAxes,
                    ha=alignment, va='bottom', fontsize=9)
    if title is not None:
        ax.set_title(title, y=1.16, fontsize=11)
    # 使用固定大小的示意圆点, handles是一个列表装了两个空的线对象
    handles = [
        Line2D(
            [], [], linestyle='None', # mock一条线出来
            marker='o', markersize=6, markerfacecolor=color, markeredgewidth=0, label=status)
        for status, color in [('up', '#E6554F'), ('down', '#4292C6')]
    ]
    ax.legend(handles=handles, frameon=False, loc='upper left', fontsize=9)


def _annotate_genes(ax, x, y, genes, sizes, labels,):
    """
    使用textalloc标注基因
    """
    from textalloc import allocate
    # AIGC, 我不懂就硬用呗
    # textalloc 的避让尺寸是像素半径, scatter 的 s 则为平方磅 
    radii = np.sqrt(sizes) * ax.figure.dpi / 72 / 2
    allocate(
        ax, x[labels], y[labels], [str(gene) for gene in genes[labels]],
        x_scatter=x, y_scatter=y, scatter_sizes=radii,
        text_scatter_sizes=radii[labels],
        textsize=8, textcolor='black', linecolor='gray', linewidth=.5,
        draw_lines=True, draw_all=True, nbr_candidates=500,
        margin=.005, min_distance=.015, max_distance=.3,
        priority_strategy=None, avoid_label_lines_overlap=True,
        avoid_crossing_label_lines=True,
        clip_on=True, zorder=4, plot_kwargs={'clip_on': True, 'zorder': 3},
    )
    
def plot_volcano(
    df: pd.DataFrame,
    conditions: list[str],
    lfc_threshold=0.5,
    p_threshold=0.05,
    anno=None,
    figarea=None,
    save_dir=None,
    title=None,
    figsize=(4, 6),
    size_factor=10,
    size_power=1.5,
    ):
    """
    绘制火山图, 返回 (fig, ax)

    df 的索引为基因名, 包含数值列 log2FoldChange 和 padj,
    需要满足:两列无缺失且有限, padj 在 [0, 1] 内
    
    conditions 为 [左侧条件, 右侧条件]
    
    lfc_threshold, p_threshold: 两个阈值

    仅对显著点按 −log10(padj) 上色
    点面积为 abs(log2FoldChange) * size_factor ** size_power
    
    anno 为 None 不标注 'all' 标注显著基因 或自己传入基因名列表 
    figarea 可传入为 (xmin, xmax, ymin, ymax)

    save_dir 指定时创建目录并保存 volcano.png
    """
    x = df['log2FoldChange'].to_numpy()
    p = df['padj'].to_numpy()
    # padj=0 映射到有限的纵坐标
    y = -np.log10(np.clip(p, np.finfo(float).tiny, None))
    # 显著点
    significant = (p < p_threshold) & (np.abs(x) > lfc_threshold)
    # 点大小
    sizes = size_factor * np.abs(x) ** size_power
    # 显著值的cutoff
    cutoff = -np.log10(p_threshold)
    # 分组
    groups = {
        'normal': ~significant,
        'down': significant & (x < -lfc_threshold),
        'up': significant & (x > lfc_threshold),
    }
    # 颜色映射
    cmaps = {
        'down': LinearSegmentedColormap.from_list('deg_down', ['#A6CBE3', '#4292C6', '#2166AC']),
        'up': LinearSegmentedColormap.from_list('deg_up', ['#F3AAA4', '#E6554F', '#C62828']),
    }

    fig, ax = plt.subplots(figsize=figsize)
    for status, mask in groups.items():
        colors = (
            '#9E9E9E' if status == 'normal'
            else cmaps[status](_color_levels(y[mask]))
        )
        # 绘制边缘
        ax.scatter(
            x[mask], y[mask],
            s=sizes[mask], linewidths=1, c='black',
            )
        ax.scatter(
            x[mask], y[mask],
            s=sizes[mask], linewidths=0, c='white',
            )
        # 正常
        ax.scatter(
            x[mask], y[mask],
            s=sizes[mask], linewidths=0, c=colors, label=status,
            )

    _decorate_axes(ax, conditions, lfc_threshold, cutoff, figarea, title)
    fig.tight_layout()
    match anno:
        case None:
            labels = None
        case 'all':
            labels = significant
        case [*genes]:
            labels = df.index.isin(genes)

    if labels is not None:
        _annotate_genes(ax, x, y, df.index, sizes, labels,)

    if save_dir is not None:
        directory = Path(save_dir)
        directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(directory / 'volcano.png', dpi=300, bbox_inches='tight')
    return fig, ax




