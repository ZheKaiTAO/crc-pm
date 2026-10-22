# PM signature 分析记录

更新：2026-10-19。

## 已确定的分析方向

- U-CAN 原发灶 bulk：114 位 IV 期患者，PM Yes 24、No 90。
- 沿用当前 notebook 的模型 `~ site + treatment + msi + histo + pm`，对比 PM Yes / No。
- 单向探索集：未 shrink 的 DE 结果，`padj < 0.1` 且 `log2FoldChange > 0.5`，40 个基因；其中 22 个通过 `padj < 0.05`。
- 阈值固定，不以基因数量、预期重合度或生存显著性继续调整。
- 将 PM-up 拆分成互斥的 epi、immune、stromal；每个可归属基因只出现一次，同时记录来源证据强弱，不能将唯一分类解释为仅在该 compartment 表达。
- 生存分析终点尚未确定，不用于当前基因筛选。
- 单独在原发灶恶性上皮中识别 CLDN18 相关细胞状态，从完整转录组独立构建其 signature，再比较与 PM_epi 的重合和富集。验证“PM 患者原发灶中更多”需要患者级 PM 标签；PM 组织表达不能替代该证据。

## 下载与预处理

数据目录：`~/project/crc-pm/bulk-data/GSE39395/`。

- `GSE39395_series_matrix.txt.gz`：54,675 probes × 24 samples，已处理的 RMA log2 信号。
- `GSE39395_family.soft.gz`：样本和平台信息。
- `GPL570.annot.gz`：GPL570 平台注释。
- `download_manifest.json`：下载 URL、字节数、SHA256；下载后已读取完整 gzip 验证。

8 位患者，每位患者都有 epi、immune、stromal 各一个样本。直接用处理好的 log2 信号；不重复取 log、不当作 counts 输入 PyDESeq2。

`prepare_gse39395.py` 整理样本、表达和 coverage。多基因映射探针排除；同一 gene symbol 的多个可用探针，预先按所有 24 个样本的平均表达最高选择，不按 compartment 效应或显著性挑探针。以 PM 基因的 HGNC approved symbol 优先匹配，未匹配时使用 GRCh38.98 symbol。来源分类后续仍需探针一致性检查。

派生文件在数据目录的 `processed/`：

- `samples.tsv`：sample_id、patient、compartment。
- `expression_gene_log2.tsv.gz`：20,848 个唯一 symbol × 24 个样本。
- `selected_probes.tsv`：保留探针及注释。
- `PM_signature_coverage.tsv`：40 基因 DE 信息、覆盖标记、选定探针、三个 compartment 均值、严格 FDR 标记。
- `coverage_summary.json`：覆盖汇总。
- `CLDN18_probes_by_sample.tsv`：CLDN18 全部 4 个探针的逐样本表达。

## Coverage 结果

40 个 PM-up 中 36 个有可用的唯一基因探针，可靠覆盖率 90%。

- ERFL、NPIPB13：GPL570 下载注释中没有对应 symbol。
- PRSS1：探针同时映射 PRSS2 / PRSS1，未作为可靠唯一覆盖。
- RGPD4：探针同时映射多个 RGPD 家族基因，未作为可靠唯一覆盖。

“未可靠覆盖”不表示基因不表达。未覆盖基因保留在原 PM signature，后续可用单细胞补充来源，不能凭芯片强制分配。

仅按选定探针的 compartment 平均 log2 信号最大值，36 个基因暂时为 epi 23、stromal 10、immune 3。这是描述性分配，尚未进行 limma 比较，不是最终细胞来源结论。

## CLDN18：需要特别核对

GPL570 中有 4 个 CLDN18 探针：214135_at、221132_at、221133_s_at、232578_at。

按预先规则选中的 214135_at，epi / immune / stromal 平均 log2 信号分别为 5.313 / 4.566 / 5.615。多数样本信号接近，patient V 的上皮信号显著升高，patient IV、V 的基质也较高。其他探针不能假定方向相同。因此不能仅凭单探针均值将 CLDN18 断言为基质来源，也不能为符合上皮假设手工更换探针或分类。

后续保留所有探针核对，结合配对检验和原发灶单细胞表达判断 CLDN18 是否标记部分患者的上皮状态；需要进一步证据解释基质信号。

## 下一步的三分法

1. 先查看 PM 基因的配对表达热图及 CLDN18 四探针逐患者图。
2. 对全部可用基因用 limma 拟合 `~ 0 + compartment + patient`，进行 epi-immune、epi-stromal、immune-stromal 比较。先拟合完整表达矩阵，再提取 PM 子集。
3. 对覆盖基因按最高 compartment 分配唯一标签，同时报告与另两组的差值、统计支持和患者一致性；接近或探针不一致的结果标记低置信度，不声称特异性。
4. 用原发灶单细胞数据复核来源；PM 和正常组织分别补充验证。Unknown 标签不强制归入三类。
5. 独立构建 CLDN18 状态 signature；比较重合时使用两种数据共同可检测且可进入筛选的基因作背景，并检查去掉 CLDN18 后的一致性。

## 可复用命令

```bash
python 06-bulk-process/prepare_gse39395.py \
  --data-dir ~/project/crc-pm/bulk-data/GSE39395 \
  --signature /path/to/PM_signature.tsv \
  --output-dir ~/project/crc-pm/bulk-data/GSE39395/processed
```

输入 signature 表至少包含 `gene_id`、`gene_name`、`padj`；当前已生成的 coverage 表可作为该输入。当前 notebook 磁盘版本仍是 0.05 阈值；本次未修改 notebook，在独立进程重建同一模型并按 0.1 提取 40 基因。

资料：[GSE39395](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39395)、[Cañellas-Socias 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC7616986/)、[limma 用户指南](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)。
