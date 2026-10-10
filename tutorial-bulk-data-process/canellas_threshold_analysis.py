"""Repeat the PM analysis with nominal P<0.05 and positive PM effect.

The threshold follows Canellas-Socias 2022; the phenotype remains PM status,
not relapse. All variant outputs are isolated from the prior strict analysis.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/crc-pm-bulk-mpl")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import hashlib
import json
import shutil
import textwrap
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from pipeline import Pipeline, COHORTS, HERE, write_json


def nominal_pm_up(result):
    """Raw-P and positive-effect rule; adjusted P and effect magnitude are not filters."""
    return result.pvalue.lt(.05) & result.log2FoldChange.gt(0)


class ThresholdPipeline(Pipeline):
    def __init__(self):
        super().__init__()
        self.source_out = self.out
        self.out = self.source_out / "canellas_threshold"
        for cohort in COHORTS + ["signature"]:
            (self.out / cohort).mkdir(parents=True, exist_ok=True)

    def matrix(self, cohort, name):
        path = self.out / cohort / f"{name}.parquet"
        if not path.exists():
            path = self.source_out / cohort / f"{name}.parquet"
        return pd.read_parquet(path)

    def clinical(self, cohort, name="clinical"):
        path = self.out / cohort / f"{name}.tsv"
        if not path.exists():
            path = self.source_out / cohort / f"{name}.tsv"
        return pd.read_csv(path, sep="\t", index_col=0)

    def volcano(self, results, cohort, effect, label):
        super().volcano(results, cohort, effect, label, significance="pvalue", title_note="Nominal P < 0.05; no fold-change minimum")

    def heatmap(self, expression, meta, genes, cohort, plot_name="signature_heatmap"):
        genes = [g for g in genes if g in expression.index]
        super().heatmap(expression, meta, genes, cohort, plot_name=plot_name,
            signature_title=f"PM-up set, nominal P < 0.05 ({len(genes)} genes displayed)")

    def select(self):
        self.annotation()
        main = self.clinical("U-CAN", "DE_main")
        mask = nominal_pm_up(main)
        signature = main.loc[mask].sort_values(["pvalue", "log2FoldChange"], ascending=[True, False]).copy()
        sensitivity = self.clinical("U-CAN", "DE_with_sensitivity")
        extra = [c for c in sensitivity if c.startswith(("untreated_", "colon_", "expanded_"))]
        signature = signature.join(sensitivity[extra])
        signature["annotation_status"] = np.where(signature.symbol.notna(), "unique_standard_symbol", "Ensembl_only_no_unique_standard_symbol")
        self.table(signature, "signature", "PM_up_signature")
        (self.out / "signature/PM_up_ensembl.txt").write_text("".join(f"{g}\n" for g in signature.index))
        (self.out / "signature/PM_up_symbols.txt").write_text("".join(f"{g}\n" for g in signature.symbol.dropna()))
        self.table(main.loc[main.pvalue.lt(.05) & main.log2FoldChange.lt(0)], "signature", "PM_down_genes")
        main["nominal_signature_member"] = mask
        self.table(main, "U-CAN", "DE_main")
        shutil.copyfile(self.source_out / "signature/ranking_background.tsv", self.out / "signature/ranking_background.tsv")
        definition = dict(selection="raw P < 0.05 AND log2FoldChange > 0", multiple_testing_filter=False,
            minimum_positive_fold_change=False, model="PyDESeq2, site + pre-treatment + MSI + PM", genes=len(signature),
            standard_symbols=int(signature.symbol.notna().sum()), phenotype="stage IV PM versus other metastatic sites, primary tumor",
            reference="https://pmc.ncbi.nlm.nih.gov/articles/PMC7616986/", adaptation="Positive log2FC in PM comparison substitutes for HR>1 in relapse Cox screening; no survival-based selection",
            scoring="Same sample-wise percentile-rank mean and fixed shared background as prior analysis",
            analysis_status="post_hoc_threshold_extension; AC-ICAM and TCGA were already inspected",
            source_DE_sha256=hashlib.sha256((self.source_out / "U-CAN/DE_main.tsv").read_bytes()).hexdigest())
        write_json(self.out / "signature/definition.json", definition)
        self.volcano(main, "U-CAN", "log2FoldChange", "log2 fold change (PM / non-PM)")
        u = np.log2(self.matrix("U-CAN", "tpm")+1)
        m = self.clinical("U-CAN")
        self.heatmap(u, m, signature.index[:50], "U-CAN")
        self.heatmap(u, m, signature.index, "U-CAN", plot_name="signature_heatmap_full")
        return signature

    def overrepresentation(self, signature):
        result = self.clinical("U-CAN", "DE_main")
        background = set(result.loc[result.pvalue.notna(), "symbol"].dropna())
        selected = set(signature.symbol.dropna()) & background
        N, K = len(background), len(selected)
        manifests = []
        for label in self.config["gene_sets"]:
            path = self.source_out / "signature" / f"{label}.gmt"
            rows = []
            for line in path.read_text().splitlines():
                parts = line.split("\t")
                members = set(parts[2:]) & background
                if not 15 <= len(members) <= 500:
                    continue
                overlap = sorted(members & selected)
                n, k = len(members), len(overlap)
                rows.append(dict(term=parts[0], background_size=N, signature_size=K, set_size=n,
                    overlap_size=k, fold_enrichment=(k/K)/(n/N), pvalue=float(hypergeom.sf(k-1,N,n,K)), overlap_genes=";".join(overlap)))
            table = pd.DataFrame(rows)
            table["padj"] = multipletests(table.pvalue, method="fdr_bh")[1]
            table = table.sort_values(["padj", "pvalue"])
            self.table(table, "signature", f"ORA_{label}", index=False)
            plot = table.loc[table.padj.lt(.05)].head(15).iloc[::-1]
            fig, ax = plt.subplots(figsize=(8, max(4,.45*len(plot))))
            if len(plot):
                labels = [textwrap.fill(t.replace("HALLMARK_", "").replace("GOBP_", "").replace("_", " "),40) for t in plot.term]
                ax.barh(range(len(plot)), -np.log10(plot.padj.clip(lower=1e-300)), color="#c84343")
                ax.set_yticks(range(len(plot)), labels, fontsize=8)
            else:
                ax.text(.5,.5,"No terms at BH-FDR < 0.05",ha="center",transform=ax.transAxes)
            ax.set(xlabel="−log10(BH-FDR)", title=f"PM-up nominal set: {label} overrepresentation")
            self.figure(fig, "signature", f"ORA_{label}")
            manifests.append(dict(library=label, source=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                tested_terms=len(table), significant_terms=int(table.padj.lt(.05).sum()), background_size=N, signature_size=K,
                test="hypergeometric, BH across all tested terms per library, including zero overlaps"))
        write_json(self.out / "signature/ORA_definition.json", manifests)
        # A threshold change cannot alter whole-transcriptome preranked GSEA.
        for label in self.config["gene_sets"]:
            for extension in ["tsv", "pdf", "png"]:
                shutil.copyfile(self.source_out / "U-CAN" / f"GSEA_{label}.{extension}", self.out / "U-CAN" / f"GSEA_{label}.{extension}")
        write_json(self.out / "U-CAN/GSEA_reuse.json", dict(reason="Same genes and Wald-statistic ranking; signature-selection threshold is not an input to preranked GSEA", source=str(self.source_out / "U-CAN")))

    def variant_report(self, signature):
        validation = json.loads((self.out / "AC-ICAM/signature_validation.json").read_text())
        survival = json.loads((self.out / "TCGA-COAD/OS_summary.json").read_text())
        scoring = pd.read_csv(self.out / "signature/scoring_genes.tsv", sep="\t")
        adjusted = survival["adjusted"]
        lines = ["# Cañellas阈值的PM gene set分析", "", "## 定义", "",
            "采用原始P<0.05且PM效应为正（log2FC>0），不要求FDR<0.05或log2FC≥0.5。保持U-CAN主模型及其部位、预治疗、MSI调整；只重新筛选已有模型结果，没有必要重拟合相同模型。",
            "[Cañellas-Socias 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC7616986/)用复发Cox模型的HR>1、原始P<0.05。本项目将这一阈值思路用于PM差异表达，以正log2FC对应正效应；统计模型与原论文不同。",
            "本分析为查看过AC-ICAM与TCGA结果后的阈值扩展，跨队列结果应视为追加验证证据，不能称为新的盲法独立验证。", "",
            f"得到{len(signature)}个唯一Ensembl基因，其中{signature.symbol.notna().sum()}个有唯一标准symbol。未能唯一对应标准symbol的条目仍按Ensembl ID保留在完整清单，并明确标注。",
            f"其中{signature.padj.lt(.05).sum()}个达到BH-FDR<0.05；其余{signature.padj.ge(.05).sum()}个只通过名义P门槛。原64基因严格signature全部保留于宽gene set。",
            f"评分使用共有背景内{len(scoring)}个基因，覆盖率{len(scoring)/len(signature):.1%}；评分公式沿用样本内百分位秩等权平均。完整清单保持{len(signature)}个，不按验证表现增删。", "",
            "## U-CAN敏感性分析", ""]
        for name in ["untreated", "colon", "expanded"]:
            c = f"{name}_log2FoldChange"
            available = signature[c].notna()
            lines.append(f"{name}：{available.sum()}个可评估基因中，{signature.loc[available,c].gt(0).sum()}个保持正效应。")
        lines += ["", "## AC-ICAM追加验证", "",
            f"7 PM对53其他转移。评分均值差异={validation['mean_difference']:.5f}，患者bootstrap 95%CI={validation['bootstrap_CI']}，双侧置换P={validation['permutation_two_sided_p']:.5g}。",
            f"逐基因可评估{validation['genes_assessable']}个，{validation['genes_same_direction']}个在PM组表达更高。", "",
            "## TCGA-COAD OS", "",
            f"OS n={survival['n']}、死亡{survival['events']}，高/低评分={survival['high_n']}/{survival['low_n']}，中位数分组log-rank P={survival['logrank_p']:.5g}。",
            f"调整年龄、性别和分期的Cox：n={adjusted['n']}，HR/评分SD={adjusted['HR']:.3f}，95%CI {adjusted['CI_low']:.3f}–{adjusted['CI_high']:.3f}，P={adjusted['p']:.5g}；PH标记={adjusted['PH_violations']}。",
            "维持OS，不搜索最优截点，也不按TCGA生存结果修改基因集。", "", "## 富集和交付", "",
            "新增Hallmark和GO BP超几何富集（ORA），背景为U-CAN有有效原始P且可注释的被检验基因；每个库对全部检验通路作BH校正。完整统计和重叠基因均输出。",
            "全转录组GSEA的排序没有变化，复用此前相同排序的结果，并记录复用依据。热图同时交付完整gene set版本和前50基因的可读版本。", ""]
        for label in self.config["gene_sets"]:
            t = pd.read_csv(self.out / "signature" / f"ORA_{label}.tsv", sep="\t")
            significant = t.loc[t.padj.lt(.05)]
            lines += [f"{label}：{len(significant)}项BH-FDR<0.05；前列为"+"；".join(significant.term.head(8))+"。", ""]
        lines += ["## 与严格版本比较", "", "| 版本 | gene set基因数 | 评分基因数 | AC置换P | TCGA KM P | TCGA调整Cox P |", "|---|---:|---:|---:|---:|---:|"]
        old_v=json.loads((self.source_out / "AC-ICAM/signature_validation.json").read_text())
        old_s=json.loads((self.source_out / "TCGA-COAD/OS_summary.json").read_text())
        old_score=pd.read_csv(self.source_out / "signature/scoring_genes.tsv",sep="\t")
        old_sig=pd.read_csv(self.source_out / "signature/PM_up_signature.tsv",sep="\t")
        lines += [f"| 严格FDR+FC | {len(old_sig)} | {len(old_score)} | {old_v['permutation_two_sided_p']:.4g} | {old_s['logrank_p']:.4g} | {old_s['adjusted']['p']:.4g} |",
            f"| Cañellas名义P阈值 | {len(signature)} | {len(scoring)} | {validation['permutation_two_sided_p']:.4g} | {survival['logrank_p']:.4g} | {adjusted['p']:.4g} |", "",
            "结果隔离存放于result/bulk/canellas_threshold；原严格版本未覆盖。宽gene set适合下游来源探索，但名义P不等于1158个基因均已被多重检验确认。", ""]
        (self.out / "analysis_report_zh.md").write_text("\n".join(lines))


if __name__ == "__main__":
    pipeline = ThresholdPipeline()
    signature = pipeline.select()
    print(f"Selected {len(signature)} nominal-P PM-up genes", flush=True)
    pipeline.validate()
    pipeline.heatmap(pipeline.matrix("AC-ICAM","expression"), pipeline.clinical("AC-ICAM"), signature.index,
        "AC-ICAM", plot_name="signature_heatmap_full")
    pipeline.survival()
    pipeline.overrepresentation(signature)
    pipeline.variant_report(signature)
    hashes = {str(p.relative_to(pipeline.source_out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
        [pipeline.source_out/"U-CAN/DE_main.tsv",pipeline.source_out/"signature/PM_up_signature.tsv",
         pipeline.source_out/"signature/ranking_background.tsv",pipeline.source_out/"TCGA-COAD/clinical.tsv"]}
    write_json(pipeline.out / "run_manifest.json", dict(status="complete", completed_utc=datetime.now(timezone.utc).isoformat(),
        source_inputs_sha256=hashes, base_pipeline_config=pipeline.config,
        selection_override=dict(raw_p_less_than=.05, positive_log2fc=True, FDR_filter=False, minimum_log2fc_filter=False),
        script_sha256=hashlib.sha256((HERE/"canellas_threshold_analysis.py").read_bytes()).hexdigest()))
    print((pipeline.out / "analysis_report_zh.md").read_text(), flush=True)
