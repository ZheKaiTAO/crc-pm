import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from pathlib import Path
from scipy import sparse


def preview_adata(
    adata: AnnData,
    max_levels: int = 5,
    n_show: int = 5,
    obs_keys: list[str] | None = None,
) -> None:
    print(f"shape: {adata.shape}")
    if adata.X is None:
        print("X: None")
    else:
        print(
            f"X: {type(adata.X).__name__}, sparse={sparse.issparse(adata.X)}, "
            f"dtype={adata.X.dtype}"
        )
    print(
        "unique index: "
        f"obs={adata.obs_names.is_unique}, var={adata.var_names.is_unique}"
    )
    cols = obs_keys if obs_keys is not None else adata.obs.columns
    print("obs:")
    for col in cols:
        s = adata.obs[col]
        n_nan = int(s.isna().sum())
        if pd.api.types.is_numeric_dtype(s):
            detail = f"min={s.min():.3g}, median={s.median():.3g}, max={s.max():.3g}"
        else:
            vc = s.value_counts(dropna=False)
            if len(vc) <= max_levels:
                detail = str(dict(vc))
            else:
                detail = f"top={vc.index[:max_levels].tolist()}"
        if n_nan:
            detail += f", nan={n_nan}"
        print(f"  {col}: {s.dtype} n_unique={s.nunique()} {detail}")
    print(f"var: {adata.var.columns.tolist()}")
    print(f"var_names[:{n_show}] = {adata.var_names[:n_show].tolist()}")
    print(f"obs_names[:{n_show}] = {adata.obs_names[:n_show].tolist()}")
    print(f"layers: {[f'{k} ({adata.layers[k].dtype})' for k in adata.layers]}")
    print(f"raw: {adata.raw is not None}")
    print(f"obsm: {[f'{k} {adata.obsm[k].shape}' for k in adata.obsm]}")
    print(f"obsp: {list(adata.obsp)}")
    print(f"uns: {list(adata.uns)}")


def to_float32(matrix: np.ndarray | sparse.spmatrix) -> np.ndarray | sparse.spmatrix:
    if not isinstance(matrix, np.ndarray | sparse.spmatrix):
        raise TypeError("matrix must be a NumPy array or SciPy sparse matrix")
    if not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(matrix.dtype, np.bool_):
        raise ValueError("matrix dtype must be numeric and non-boolean")
    if np.issubdtype(matrix.dtype, np.complexfloating):
        raise ValueError("matrix dtype must be real")

    values = matrix.data if sparse.issparse(matrix) else matrix
    if not np.all(np.isfinite(values)):
        raise ValueError("matrix contains non-finite values")
    float32_info = np.finfo(np.float32)
    if np.any(values < float32_info.min) or np.any(values > float32_info.max):
        raise ValueError("matrix contains values outside the float32 range")

    return matrix.astype(np.float32, copy=False)


def to_int32_counts(matrix: np.ndarray | sparse.spmatrix) -> np.ndarray | sparse.spmatrix:
    if not isinstance(matrix, np.ndarray | sparse.spmatrix):
        raise TypeError("matrix must be a NumPy array or SciPy sparse matrix")
    if not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(matrix.dtype, np.bool_):
        raise ValueError("counts dtype must be numeric and non-boolean")
    if np.issubdtype(matrix.dtype, np.complexfloating):
        raise ValueError("counts dtype must be real")

    values = matrix.data if sparse.issparse(matrix) else matrix
    if not np.all(np.isfinite(values)):
        raise ValueError("counts contain non-finite values")
    if np.any(values < 0):
        raise ValueError("counts contain negative values")
    if np.any(values != np.floor(values)):
        raise ValueError("counts contain non-integer values")
    if np.any(values > np.iinfo(np.int32).max):
        raise ValueError("counts exceed the int32 range")

    return matrix.astype(np.int32, copy=False)


def preview_cell_quality(
    adata: AnnData,
    metric: str,
    value: object = None,
    enumeration: bool = False,
    calc_all: bool = False,
    outdir: str | Path = Path("result/QC-results"),
) -> None:
    """Preview QC metrics (total_counts, pct_counts_mt, counts vs genes) by group.

    enumeration: when True, go through all values in metric and store fig
    """
    if calc_all:
        p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
        p2 = sc.pl.violin(adata, "pct_counts_mt")
        p3 = sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    elif not enumeration:
        data = adata[adata.obs[metric] == value]
        p1 = sns.displot(data.obs["total_counts"], bins=50, kde=True)
        p1.fig.suptitle(f"{metric}: {value} - total counts")
        p1.fig.tight_layout()

        p2 = sc.pl.violin(data, "pct_counts_mt", show=False)
        p2.set_title(f"{metric}: {value} - mitochondrial percentage")

        p3 = sc.pl.scatter(
            data, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False
        )
        p3.set_title(f"{metric}: {value} - counts vs genes")

        plt.show()
    else:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        for i, name in enumerate(adata.obs[metric].unique(), start=1):
            data = adata[adata.obs[metric] == name]

            group_dir = outdir / f"{i:02d}_{name}"
            group_dir.mkdir(parents=True, exist_ok=True)

            p1 = sns.displot(data.obs["total_counts"], bins=50, kde=True)
            p1.fig.suptitle(f"{metric}: {name} - total counts")
            p1.fig.tight_layout()
            p1.fig.savefig(group_dir / "01_total_counts.png", dpi=300, bbox_inches="tight")
            plt.close(p1.fig)

            p2 = sc.pl.violin(data, "pct_counts_mt", show=False)
            p2.set_title(f"{metric}: {name} - mitochondrial percentage")
            p2.figure.savefig(group_dir / "02_pct_counts_mt.png", dpi=300, bbox_inches="tight")
            plt.close(p2.figure)

            p3 = sc.pl.scatter(
                data, "total_counts", "n_genes_by_counts", color="pct_counts_mt", show=False
            )
            p3.set_title(f"{metric}: {name} - counts vs genes")
            p3.figure.savefig(group_dir / "03_counts_vs_genes.png", dpi=300, bbox_inches="tight")
            plt.close(p3.figure)
