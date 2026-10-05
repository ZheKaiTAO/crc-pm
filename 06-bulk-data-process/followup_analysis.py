"""Requested single-gene OS tests and a descriptive signature-threshold audit.

These analyses never change the frozen PM signature or choose cutoffs by outcome.
"""
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/crc-pm-bulk-mpl")
import argparse
import json
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, proportional_hazard_test
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from pipeline import Pipeline, write_json


def single_gene_os(pipeline, genes):
    pipeline.annotation()
    expression = pipeline.matrix("TCGA-COAD", "tpm")
    clinical = pipeline.clinical("TCGA-COAD")
    folder = pipeline.out / "TCGA-COAD/single_gene"
    folder.mkdir(parents=True, exist_ok=True)
    summaries = []
    for gene in genes:
        gene_id, status = pipeline.resolve(gene)
        if status != "mapped" or gene_id not in expression.index:
            raise ValueError(f"Gene {gene} not uniquely present in TCGA")
        d = clinical.join(expression.loc[gene_id].rename("TPM"))
        d = d.loc[d.OS_exclusion.fillna("").eq("") & d.TPM.notna()].copy()
        d["log2_TPM_plus_1"] = np.log2(d.TPM + 1)
        mean, sd = d.log2_TPM_plus_1.mean(), d.log2_TPM_plus_1.std(ddof=1)
        if sd <= 0 or not np.isfinite(sd):
            raise ValueError(f"Constant or invalid expression for {gene}")
        d["expression_z"] = (d.log2_TPM_plus_1 - mean) / sd
        cutoff = float(d.TPM.median())
        d["expression_group"] = np.where(d.TPM.gt(cutoff), "High", "Low")
        if d.expression_group.nunique() != 2:
            raise ValueError(f"Median split produces one group for {gene}")
        pipeline.table(d, "TCGA-COAD/single_gene", f"{gene}_OS_patients")
        high, low = d.loc[d.expression_group.eq("High")], d.loc[d.expression_group.eq("Low")]
        lr = logrank_test(high.OS_days, low.OS_days, high.OS_event, low.OS_event)
        fig, ax = plt.subplots(figsize=(7, 6))
        fitters = []
        for name, group, color in [("High", high, "#c84343"), ("Low", low, "#307ca8")]:
            fit = KaplanMeierFitter(label=f"{name} (n={len(group)})").fit(group.OS_days/365.25, group.OS_event)
            fit.plot_survival_function(ax=ax, color=color, ci_show=True)
            fitters.append(fit)
            pipeline.table(fit.survival_function_, "TCGA-COAD/single_gene", f"{gene}_KM_{name}")
        add_at_risk_counts(*fitters, ax=ax)
        ax.set(xlabel="Years from diagnosis", ylabel="Overall survival probability", ylim=(0, 1.02),
            title=f"TCGA-COAD {gene}: log-rank P={lr.p_value:.3g}\nMedian TPM={cutoff:.4g}; high > median")
        pipeline.figure(fig, "TCGA-COAD/single_gene", f"{gene}_OS_KM")
        summary = dict(gene=gene, ensembl_id=gene_id, endpoint="OS", n=len(d), events=int(d.OS_event.sum()),
            high_n=len(high), low_n=len(low), high_events=int(high.OS_event.sum()), low_events=int(low.OS_event.sum()),
            median_TPM=cutoff, logrank_p=float(lr.p_value), expression_transform="log2(TPM+1)",
            expression_mean=float(mean), expression_SD=float(sd), HR_units="per OS-population SD of log2(TPM+1)")
        for label, cols, formula in [("unadjusted", ["OS_days", "OS_event", "expression_z"], "expression_z"),
            ("adjusted", ["OS_days", "OS_event", "expression_z", "age", "sex", "stage"], "expression_z + age + C(sex) + C(stage)")]:
            data = d[cols].dropna()
            model = CoxPHFitter().fit(data, duration_col="OS_days", event_col="OS_event", formula=formula)
            pipeline.table(model.summary, "TCGA-COAD/single_gene", f"{gene}_Cox_{label}")
            pipeline.table(data, "TCGA-COAD/single_gene", f"{gene}_Cox_{label}_patients")
            ph = proportional_hazard_test(model, data, time_transform="rank").summary
            pipeline.table(ph, "TCGA-COAD/single_gene", f"{gene}_PH_{label}")
            row = model.summary.loc["expression_z"]
            summary.update({f"{label}_n": len(data), f"{label}_events": int(data.OS_event.sum()),
                f"{label}_HR": float(row["exp(coef)"]), f"{label}_CI_low": float(row["exp(coef) lower 95%"]),
                f"{label}_CI_high": float(row["exp(coef) upper 95%"]), f"{label}_p": float(row["p"]),
                f"{label}_PH_expression_p": float(ph.loc["expression_z", "p"]),
                f"{label}_PH_violations": ph.index[ph.p.lt(.05)].tolist()})
        summaries.append(summary)
    result = pd.DataFrame(summaries)
    # Each test/model family contains the two requested genes; keep raw and BH P.
    for family in ["logrank", "unadjusted", "adjusted"]:
        result[f"{family}_BH_q"] = multipletests(result[f"{family}_p"], method="fdr_bh")[1]
    pipeline.table(result, "TCGA-COAD/single_gene", "OS_summary", index=False)
    write_json(folder / "OS_summary.json", result.to_dict("records"))
    return result


def threshold_audit(pipeline):
    result = pipeline.clinical("U-CAN", "DE_main")
    rules = [
        ("nominal_P_0.05_up", "P < 0.05, log2FC > 0", result.pvalue.lt(.05) & result.log2FoldChange.gt(0)),
        ("nominal_P_0.05_up_FC", "P < 0.05, log2FC >= 0.5", result.pvalue.lt(.05) & result.log2FoldChange.ge(.5)),
        ("FDR_0.05_up", "FDR < 0.05, log2FC > 0", result.padj.lt(.05) & result.log2FoldChange.gt(0)),
        ("fixed_signature", "FDR < 0.05, log2FC >= 0.5", result.padj.lt(.05) & result.log2FoldChange.ge(.5)),
        ("FDR_0.10_up", "FDR < 0.10, log2FC > 0", result.padj.lt(.1) & result.log2FoldChange.gt(0)),
    ]
    rows = []
    for name, description, mask in rules:
        annotated = mask & result.symbol.notna() & result.symbol.ne("")
        rows.append(dict(rule=name, threshold=description, total_genes=int(mask.sum()),
            uniquely_annotated_genes=int(annotated.sum()), status="frozen_primary" if name == "fixed_signature" else "descriptive_exploration_only"))
    table = pd.DataFrame(rows)
    pipeline.table(table, "signature", "threshold_audit", index=False)
    broad = result.loc[rules[0][2]].sort_values("pvalue").copy()
    broad["status"] = "exploratory_nominal_candidates_not_FDR_confirmed_signature"
    pipeline.table(broad, "signature", "PM_up_nominal_candidates_exploratory")
    return table


def followup_report(pipeline, genes, audit):
    lines = ["# CLDN18/CD55单基因OS及signature数量说明", "", "## 单基因OS", "",
        "沿用原TCGA患者去重与OS构建：原发灶、全分期、可靠死亡或随访时间。中位数分组（高表达>中位数）；连续Cox使用log2(TPM+1)，HR按OS人群的一个标准差报告；调整年龄、性别和分期。未搜索最佳截点。", "",
        "本次只检验OS，没有新增PFI/DFI终点，因此不能据此判断原signature的OS阴性是否由终点选择导致。", ""]
    lines += ["[TCGA-CDR原始研究](https://pmc.ncbi.nlm.nih.gov/articles/PMC6066282/)支持COAD使用OS、PFI、DFI和DSS，因此不能把OS阴性直接归因于OS不适用。若重点是疾病进展/术后复发，建议预先指定补充PFI或适用患者中的DFI及其人群；这些终点仍不是PM特异性结局。", ""]
    for r in genes.to_dict("records"):
        lines += [f"### {r['gene']}", "", f"OS n={r['n']}，死亡{r['events']}；高/低表达各{r['high_n']}/{r['low_n']}例；中位TPM={r['median_TPM']:.4g}。",
            f"KM log-rank P={r['logrank_p']:.4g}，两基因BH q={r['logrank_BH_q']:.4g}。",
            f"未调整连续Cox HR={r['unadjusted_HR']:.3f}（95%CI {r['unadjusted_CI_low']:.3f}–{r['unadjusted_CI_high']:.3f}），P={r['unadjusted_p']:.4g}。",
            f"调整Cox n={r['adjusted_n']}，HR={r['adjusted_HR']:.3f}（95%CI {r['adjusted_CI_low']:.3f}–{r['adjusted_CI_high']:.3f}），P={r['adjusted_p']:.4g}，两基因BH q={r['adjusted_BH_q']:.4g}。",
            f"PH检验：表达项未调整P={r['unadjusted_PH_expression_p']:.4g}、调整P={r['adjusted_PH_expression_p']:.4g}；调整模型所有标记项为{r['adjusted_PH_violations']}。",
            "PH检验若提示违反，常数HR仅是全随访期间的概括，不应解释为恒定风险比。", ""]
    lines += ["## 基因数量差异", "", "完整PM signature是64个基因；跨队列评分用57个共同可用基因；热图只展示按FDR排序的前50个。50是展示限制，不是筛选限制。", "",
        "Cañellas-Socias 2022整合了7个队列的1830个原发CRC样本，以复发相关Cox模型的HR>1、原始P<0.05构建2530基因All-HR；之后才筛选出99基因的上皮EpiHR。本项目发现集只有24 PM/90其他转移，比较的是IV期转移部位，并使用BH-FDR<0.05及log2FC≥0.5。人群、问题、统计模型及阈值均不相同。",
        "来源：[论文正文与方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC7616986/)。其原始P阈值不能等同于多重检验校正的FDR阈值。", "",
        "| U-CAN规则 | 基因数 |", "|---|---:|"]
    lines += [f"| {r.threshold} | {r.total_genes} |" for r in audit.itertuples()]
    nominal_count = int(audit.loc[audit.rule.eq("nominal_P_0.05_up"), "total_genes"].iloc[0])
    lines += ["", f"推荐保留64基因作为已固定的严格signature；用于下游细胞来源探索时，同时利用完整Wald排序/GSEA和单独标注的{nominal_count}个名义P候选。后者不是{nominal_count}个已确认PM基因，也不替换现有signature。",
        "不建议按期待数量手挑基因，或因为TCGA OS不显著而调整阈值。若决定采用更宽筛选规则，应明确标为新的探索性分析，并评估患者重采样稳定性及外部验证；不能把已经查看过的验证结果再当作未经选择的独立证据。",
        "本次仅生成阈值审计和候选表，未改动原64基因清单或原signature评分、生存结果。", ""]
    (pipeline.out / "followup_report_zh.md").write_text("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes", nargs="+", default=["CLDN18", "CD55"])
    args = parser.parse_args()
    if len(set(args.genes)) != len(args.genes):
        parser.error("Gene names must be unique")
    pipeline = Pipeline()
    genes = single_gene_os(pipeline, args.genes)
    audit = threshold_audit(pipeline)
    followup_report(pipeline, genes, audit)
    print(genes[["gene", "n", "logrank_p", "adjusted_HR", "adjusted_p", "adjusted_PH_expression_p"]].to_string(index=False))
    print(audit.to_string(index=False))
