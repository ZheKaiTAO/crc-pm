# 基因集名称解析与评分

基因集评分现在独立于数据 loader 的 GTF 映射。没有修改 loader、整合数据或表达矩阵。

## 数据与规则

- `tables/geneset/2026_Nat_Canellas_TME-HR.csv`：12 条 `|||` 记录拆分，2,431 行变为 2,449 行；按名称保序去重后为 2,447 行。原始 CSV 保存在 `tables/geneset_original/`，不会被 `get_names()` 枚举。
- `tables/HGNC/hgnc_complete_set.txt`：本地 HGNC complete set 快照。`metadata.json` 记录来源、下载时间和 SHA-256；解析器首次加载时核对校验值，不在运行评分时联网。
- 先识别 HGNC 正式名称，再联合查询历史名和别名。正式名称优先于其作为其他基因别名的用途.
- HGNC 候选与本地 GRCh38.111、98、84 的名称→ENSG 对应关系核对；HGNC 没有对应时允许唯一 GTF 对应回退。任何多个基因或多个 ENSG 候选都标记 `ambiguous`，不会按生物类型或矩阵中是否存在来猜测。这可能使此前精确匹配成功的名称也进入冲突报告。
- ENSG 输入去掉版本后缀；symbol 与 ENSG 可以混合。有效结果按 ENSG 保序去重。
- `38231`、`36951`、`37316` 在 2022 原始补充表的 TME-HR 列中本来就是数字字符串（Excel 行 623、1485、1785），暂不推断原基因，标记 `invalid`。
- 评分仅使用矩阵中存在的已解析 ENSG。缺失比例分母为唯一已解析 ENSG 数加唯一未解析名称数；分子为缺失 ENSG 数加未解析名称数。超过 25% 报错，恰好 25% 允许；空集或完全没有可用基因报错。

## 使用

Notebook 已加载旧模块时先重启 kernel，然后继续使用原来的调用：

```python
src.Eval.sEval(adata, layer='log1p', score_name='2026_Nat_Canellas_TME-HR')
```

查看映射和实际数据覆盖：

```python
ids, report = src.Eval.read_gene_set('2026_Nat_Canellas_TME-HR', return_report=True)
available, report, summary = src.GeneSetMapping.assess_coverage(report, adata.var_names)
report.loc[report.status.ne('mapped') | report.in_dataset.eq(False)]
```

`read_gene_set()` 默认仍返回 ENSG 列表；传入 `return_report=True` 返回 `(ids, report)`。报告列为原始条目、拆分名称、HGNC 标准名称、ENSG、来源、状态、候选 ENSG；覆盖检查另加 `in_dataset`。`unresolved`、`invalid`、`ambiguous` 均计入未解析数量。

重新生成全部基因集报告（只读取 H5AD 的 var 索引）：

```sh
python code/07-gene_set_eval/audit_gene_sets.py
```

输出位于 `result/gene-set-mapping/`；`summary.csv` 为汇总，每个基因集另有逐项 CSV。可用 `--adata` 和 `--output` 指定输入与输出路径。

验证命令（已在项目 pm 环境的 Scanpy 1.12.1 下通过，包括真实小矩阵评分）：

```sh
MPLCONFIGDIR=/tmp/crc-pm-mpl NUMBA_CACHE_DIR=/tmp/crc-pm-numba /home/zhekai/miniconda3/envs/pm/bin/python -m unittest discover -s code/tests -p 'test_gene_set_mapping.py' -v
```

## 参考

- HGNC multi-symbol checker：https://www.genenames.org/tools/multi-symbol-checker/
- HGNC 字段说明：https://www.genenames.org/help/custom-downloads/
- HGNChelper：https://waldronlab.io/HGNChelper/reference/checkGeneSymbols.html

更新 HGNC 快照时同时更新 metadata 的下载时间和校验值，重新运行审计并检查映射变化；运行中的 kernel 需重启以清除缓存。
