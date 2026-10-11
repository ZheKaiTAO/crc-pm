import numpy as np
import pandas as pd
from formulaic import model_matrix
from formulaic.errors import FormulaicError


def IsFullRank(metadata: pd.DataFrame, design: str | None = None) -> bool:
    """
    检查 DEG 设计矩阵是否列满秩，返回 bool。
    metadata 行是样本，列是协变量；design 是 Formulaic 单侧公式。
    字符串列由 Formulaic 编码为分类变量，不自动转换为连续数值。
    """
    data = metadata.copy(deep=True)

    # 默认使用 metadata 中所有的列，并包含截距。
    if design is None:
        # 内部列名避免原始列名中的特殊字符影响公式解析。
        data.columns = [f"__deg_covariate_{i}" for i in range(data.shape[1])]
        design = "~ " + " + ".join(data.columns)

    try:
        matrix = model_matrix(
            design, data, context={},
            na_action="raise", output="pandas",
        )
    except (FormulaicError, ValueError, TypeError) as exc:
        raise ValueError(f"Cannot construct design matrix: {exc}") from exc
    values = matrix.to_numpy(dtype=float)
    return bool(np.linalg.matrix_rank(values) == values.shape[1])
