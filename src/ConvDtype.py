import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData
from pathlib import Path
from scipy import sparse


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