"""PM bulk RNA-seq discovery, independent validation and prognosis analysis."""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import re
import tarfile
import urllib.request
import warnings
import textwrap
import numpy as np
import pandas as pd
import scipy.stats as st
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
COHORTS = ["U-CAN", "AC-ICAM", "TCGA-COAD"]


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=str) + "\n")


def numeric(series):
    return pd.to_numeric(series, errors="coerce")


def unique_value(series):
    values = series.dropna().unique()
    return values[0] if len(values) == 1 else np.nan


def stage_number(value):
    m = re.match(r"(?:Stage\s*)?(IV|III|II|I)(?:[ABC]|\b)", str(value), flags=re.I)
    return {"I": 1, "II": 2, "III": 3, "IV": 4}.get(m[1].upper(), np.nan) if m else np.nan


def choose_tcga(sheet):
    """Deterministic patient-level selection, preserving reasons for every file."""
    s = sheet.copy()
    primary = s["Tumor Descriptor"].eq("Primary") & s["Sample ID"].str[13:15].eq("01")
    s["selection_reason"] = "not_primary_tumor"
    p = s.loc[primary].copy()
    p["priority"] = (~p["Sample ID"].str.endswith("01A")).astype(int)
    p = p.sort_values(["Case ID", "priority", "Sample ID", "File ID"])
    chosen = p.drop_duplicates("Case ID").index
    s.loc[primary, "selection_reason"] = "other_primary_sample_or_file"
    s.loc[chosen, "selection_reason"] = "selected"
    return s, s.loc[chosen].sort_values("Case ID")


def build_os(clinical, follow):
    """One patient, one OS record; no guessed death times or endpoint imputation."""
    rows = []
    for patient, d in clinical.groupby("cases.submitter_id", sort=True):
        primary = d.loc[d["diagnoses.diagnosis_is_primary_disease"].astype(str).str.lower().eq("true") &
                        numeric(d["diagnoses.days_to_diagnosis"]).eq(0)].drop_duplicates()
        reason = ""
        statuses = set(d["demographic.vital_status"].dropna()) - {"Not Reported", "'--"}
        deaths = numeric(d["demographic.days_to_death"]).dropna().unique()
        deaths = deaths[deaths >= 0]
        vital = next(iter(statuses)) if len(statuses) == 1 else None
        f = follow.loc[follow["cases.submitter_id"].eq(patient)]
        times = pd.concat([numeric(primary["diagnoses.days_to_last_follow_up"]),
                           numeric(f["follow_ups.days_to_follow_up"])])
        times = times[times >= 0]
        if primary.empty:
            reason = "no_unambiguous_primary_diagnosis_at_day_zero"
        elif vital not in {"Alive", "Dead"}:
            reason = "unknown_or_conflicting_vital_status"
        elif len(deaths) > 1:
            reason = "conflicting_death_times"
        elif vital == "Alive" and len(deaths):
            reason = "alive_with_death_time"
        elif vital == "Dead" and len(deaths) != 1:
            reason = "death_without_time"
        event = int(vital == "Dead")
        time = float(deaths[0]) if event and len(deaths) == 1 else (float(times.max()) if len(times) else np.nan)
        if event and len(deaths) == 1 and len(times) and times.max() > time:
            reason = reason or "follow_up_after_death"
        if not np.isfinite(time) or time <= 0:
            reason = reason or "nonpositive_or_missing_time"
        age = unique_value(numeric(primary["diagnoses.age_at_diagnosis"])) / 365.25
        stage = unique_value(primary["diagnoses.ajcc_pathologic_stage"].map(stage_number))
        sex = unique_value(d["demographic.sex_at_birth"].replace({"'--": np.nan, "Not Reported": np.nan}))
        rows.append(dict(patient_id=patient, OS_days=time, OS_event=event, age=age, stage=stage,
                         sex=sex, OS_exclusion=reason, vital_status=vital,
                         primary_diagnosis_rows=len(primary), follow_up_rows=len(f)))
    return pd.DataFrame(rows).set_index("patient_id")


def rank_score(expression, signature, background, minimum=0.8):
    """Genes x samples. Scores use a fixed shared background and observed genes only."""
    if not signature:
        raise ValueError("Empty signature")
    background = pd.Index(background)
    if not background.isin(expression.index).all():
        raise ValueError("Fixed ranking background missing from expression")
    available = [g for g in signature if g in expression.index and g in background]
    coverage = len(available) / len(signature)
    if coverage < minimum:
        raise ValueError(f"Signature coverage {coverage:.3f} < {minimum}")
    ranks = expression.loc[background].rank(axis=0, method="average", pct=True)
    return ranks.loc[available].mean(axis=0), coverage


class Pipeline:
    def __init__(self):
        self.config = json.loads((HERE / "config.json").read_text())
        self.raw = PROJECT / "bulk-data"
        self.out = PROJECT / "result" / "bulk"
        for name in COHORTS + ["signature"]:
            (self.out / name).mkdir(parents=True, exist_ok=True)
        self.mapping = None
        sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)

    def table(self, d, cohort, name, index=True):
        d.to_csv(self.out / cohort / f"{name}.tsv", sep="\t", index=index)

    def save_matrix(self, d, cohort, name):
        assert d.index.is_unique and d.columns.is_unique
        assert np.isfinite(d.to_numpy()).all() and (d.to_numpy() >= 0).all()
        d.index.name = "ensembl_id"
        d.to_parquet(self.out / cohort / f"{name}.parquet")

    def matrix(self, cohort, name):
        return pd.read_parquet(self.out / cohort / f"{name}.parquet")

    def clinical(self, cohort, name="clinical"):
        return pd.read_csv(self.out / cohort / f"{name}.tsv", sep="\t", index_col=0)

    def annotation(self):
        h = pd.read_csv(PROJECT / "tables/HGNC/hgnc_complete_set.txt", sep="\t", dtype=str).fillna("")
        h = h.loc[h.status.eq("Approved")]
        official, alternate, gtf = {}, {}, {}
        def add(target, key, value):
            if key and value and value.startswith("ENSG"):
                target.setdefault(key, set()).add(value.split(".")[0])
        for r in h.itertuples():
            add(official, r.symbol, r.ensembl_gene_id)
            for field in [r.prev_symbol, r.alias_symbol]:
                for alias in field.split("|"):
                    add(alternate, alias, r.ensembl_gene_id)
        for version in ["111", "98", "84"]:
            d = pd.read_csv(PROJECT / f"tables/GTF/GRCh38.{version}/gene_id_gene_name_gene_biotype.csv")
            for r in d.itertuples():
                add(gtf, str(r.gene_name), str(r.gene_id))
        ids = h.loc[h.ensembl_gene_id.str.startswith("ENSG"), ["ensembl_gene_id", "symbol"]].copy()
        ids = ids.loc[~ids.ensembl_gene_id.duplicated(False) & ~ids.symbol.duplicated(False)]
        self.id_symbols = ids.set_index("ensembl_gene_id").symbol.to_dict()
        self.maps = official, alternate, gtf

    def resolve(self, symbol):
        official, alternate, gtf = self.maps
        candidates = official.get(symbol, alternate.get(symbol, gtf.get(symbol, set())))
        # Official HGNC wins; historical names must not contradict GTF.
        if symbol not in official:
            candidates = candidates | gtf.get(symbol, set())
        return (next(iter(candidates)), "mapped") if len(candidates) == 1 else (None, "ambiguous" if candidates else "unresolved")

    def figure(self, fig, cohort, name):
        fig.savefig(self.out / cohort / f"{name}.pdf", bbox_inches="tight")
        fig.savefig(self.out / cohort / f"{name}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    def qc(self, expression, meta, cohort):
        from sklearn.decomposition import PCA
        x = np.log2(expression + 1) if cohort != "AC-ICAM" else expression
        variable = x.var(axis=1).nlargest(min(2000, len(x))).index
        pca = PCA(n_components=2, random_state=self.config["seed"])
        coords = pca.fit_transform(x.loc[variable].T)
        d = pd.DataFrame(coords, index=x.columns, columns=["PC1", "PC2"]).join(meta)
        self.table(d, cohort, "PCA_samples")
        fig, ax = plt.subplots(figsize=(6, 5))
        sns.scatterplot(data=d, x="PC1", y="PC2", hue="PM" if "PM" in d else "stage", ax=ax)
        ax.set(xlabel=f"PC1 ({pca.explained_variance_ratio_[0]:.1%})", ylabel=f"PC2 ({pca.explained_variance_ratio_[1]:.1%})", title=f"{cohort} primary tumors")
        self.figure(fig, cohort, "PCA")
        self.table(pd.DataFrame({"sum": expression.sum(), "detected_genes": expression.gt(0).sum()}), cohort, "sample_QC")

    def prepare(self):
        self.annotation()
        udir = self.raw / "U-CAN"
        u = pd.read_excel(udir / "Supplementary_Table_01.xlsx", header=2)
        self.table(u, "U-CAN", "clinical_source", index=False)
        meta = pd.DataFrame(dict(sample_id=u["RNA Tumor Sample Barcode"], patient_id=u["Sample ID"],
            PM=u["Peritoneum metastases"].map({"Yes": 1, "No": 0}), site=u["Primary Site Disease"],
            treatment=u["Pre-Treated"], MSI=u["MSI Status"], age=u["Age at diagnosis"], sex=u.Sex,
            specimen=u["Tumour Specimen Type"], stage=u["Tumour Stage"].map(stage_number),
            CMS=u["CMS Tumour"], purity=u["Tumour Cell Content Pathology"])) .set_index("sample_id")
        meta["exclusion"] = np.where(meta.PM.isna(), "PM_not_applicable", "")
        self.table(meta, "U-CAN", "clinical_all")
        eligible = meta.loc[meta.PM.notna() & meta.stage.eq(4)].copy()
        eligible.PM = eligible.PM.astype(int)
        assert eligible.PM.value_counts().to_dict() == {0: 90, 1: 24}
        assert eligible.patient_id.is_unique
        self.table(eligible, "U-CAN", "clinical")
        for kind in ["count", "TPM"]:
            d = pd.read_csv(udir / f"CRC.SW.mRNA.{kind}.txt.gz", sep="\t", index_col=0)
            d.index = d.index.str.split(".").str[0]
            d = d.loc[:, eligible.index]
            if kind == "count":
                assert np.equal(d.to_numpy(), np.floor(d.to_numpy())).all()
            self.save_matrix(d, "U-CAN", "counts" if kind == "count" else "tpm")
        um = self.matrix("U-CAN", "tpm")
        self.table(pd.DataFrame({"symbol": [self.id_symbols.get(g, "") for g in um.index]}, index=um.index), "U-CAN", "gene_mapping")
        self.qc(um, eligible, "U-CAN")

        adir = self.raw / "AC-ICAM"
        patients = pd.read_csv(adir / "data_clinical_patient.txt", sep="\t", comment="#")
        samples = pd.read_csv(adir / "data_clinical_sample.txt", sep="\t", comment="#")
        assert patients.PATIENT_ID.is_unique and samples.SAMPLE_ID.is_unique
        a = samples.merge(patients, on="PATIENT_ID", validate="one_to_one").set_index("SAMPLE_ID")
        site_recorded = a.LOCATION_METASTASIS_AT_DX.notna()
        a["PM"] = a.LOCATION_METASTASIS_AT_DX.fillna("").str.contains(r"\bPeritoneum\b", case=False).astype(int)
        a["exclusion"] = np.where(~a.AJCC_PATH_STAGE.eq(4), "not_stage_IV", np.where(~site_recorded, "metastatic_site_missing", ""))
        self.table(a, "AC-ICAM", "clinical_all")
        am = pd.DataFrame(dict(patient_id=a.PATIENT_ID, PM=a.PM, age=a.AGE_AT_DX, sex=a.SEX,
                              stage=a.AJCC_PATH_STAGE, CMS=a.CMS, site=a.TUMOR_ANATOMIC_LOCATION,
                              treatment=a.NEOADJUVANT_TREATED)).loc[a.exclusion.eq("")]
        assert am.PM.value_counts().to_dict() == {0: 53, 1: 7}
        self.table(am, "AC-ICAM", "clinical")
        expression = pd.read_csv(adir / "data_mrna_seq_expression.txt", sep="\t", index_col=0)
        expression = expression.groupby(level=0, sort=False).mean()
        mapping = pd.DataFrame([dict(symbol=s, ensembl_id=self.resolve(s)[0], status=self.resolve(s)[1]) for s in expression.index])
        self.table(mapping, "AC-ICAM", "gene_mapping", index=False)
        good = mapping.loc[mapping.status.eq("mapped")]
        e = expression.loc[good.symbol].copy()
        e.index = good.ensembl_id.values
        e = e.groupby(level=0, sort=False).mean()
        e = e.loc[e.gt(0).any(axis=1), am.index]
        self.save_matrix(e, "AC-ICAM", "expression")
        write_json(self.out / "AC-ICAM/expression_scale.json", dict(scale="public EDASeq-normalized continuous expression; log-like but transform constant unconfirmed", additional_log_transform=False, duplicate_symbols="mean on public scale", minimum=float(e.min().min()), maximum=float(e.max().max())))
        self.qc(e, am, "AC-ICAM")

        tdir = self.raw / "tcga-coad"
        sheet = pd.read_csv(tdir / "gdc_sample_sheet.2026-07-12.tsv", sep="\t")
        audit, chosen = choose_tcga(sheet)
        paths = {r["File ID"]: tdir / "data" / r["File ID"] / r["File Name"] for _, r in audit.iterrows()}
        duplicates = []
        for sample, group in audit.loc[audit["Tumor Descriptor"].eq("Primary")].groupby("Sample ID"):
            if len(group) > 1:
                reference = None
                for _, r in group.iterrows():
                    expr = pd.read_csv(paths[r["File ID"]], sep="\t", comment="#").set_index("gene_id")
                    v = expr[["unstranded", "tpm_unstranded"]]
                    same = reference is None or v.equals(reference)
                    duplicates.append(dict(sample_id=sample, file_id=r["File ID"], identical_expression=same))
                    if reference is None:
                        reference = v
        self.table(pd.DataFrame(duplicates), "TCGA-COAD", "duplicate_file_comparison", index=False)
        self.table(audit, "TCGA-COAD", "sample_selection", index=False)
        values, count_values, gene_info = {}, {}, None
        for _, r in chosen.iterrows():
            d = pd.read_csv(paths[r["File ID"]], sep="\t", comment="#")
            d = d.loc[d.gene_id.str.startswith("ENSG")].copy()
            # PAR_Y is a distinct GENCODE row, not an Ensembl version suffix.
            d.gene_id = d.gene_id.str.replace(r"\.\d+(?=_|$)", "", regex=True)
            d = d.set_index("gene_id")
            assert d.index.is_unique
            if gene_info is None:
                gene_info = d[["gene_name", "gene_type"]]
            values[r["Case ID"]] = d.tpm_unstranded
            count_values[r["Case ID"]] = d.unstranded
        tm = pd.DataFrame(values)
        self.save_matrix(tm, "TCGA-COAD", "tpm")
        self.save_matrix(pd.DataFrame(count_values), "TCGA-COAD", "counts")
        gene_info["symbol"] = [self.id_symbols.get(g, "") for g in gene_info.index]
        self.table(gene_info, "TCGA-COAD", "gene_mapping")
        with tarfile.open(tdir / "clinical.cart.2026-07-12.tar.gz") as t:
            clinical = pd.read_csv(t.extractfile("clinical.tsv"), sep="\t", low_memory=False)
            follow = pd.read_csv(t.extractfile("follow_up.tsv"), sep="\t", low_memory=False)
        os_table = build_os(clinical, follow).reindex(tm.columns)
        self.table(os_table, "TCGA-COAD", "clinical")
        self.qc(tm, os_table, "TCGA-COAD")
        common = sorted(set(um.index) & set(e.index) & set(tm.index) & set(self.id_symbols))
        pd.DataFrame({"ensembl_id": common, "symbol": [self.id_symbols[g] for g in common]}).to_csv(self.out / "signature/ranking_background.tsv", sep="\t", index=False)
        write_json(self.out / "preparation_summary.json", {"U-CAN": eligible.PM.value_counts().to_dict(), "AC-ICAM": am.PM.value_counts().to_dict(), "TCGA_patients": len(tm.columns), "TCGA_OS_eligible": int(os_table.OS_exclusion.fillna("").eq("").sum()), "common_ranking_genes": len(common)})
        self.provenance()

    def provenance(self):
        packages = ["numpy", "pandas", "scipy", "statsmodels", "pydeseq2", "lifelines", "gseapy", "matplotlib", "seaborn", "pyarrow", "scikit-learn"]
        write_json(self.out / "versions.json", {p: importlib.metadata.version(p) for p in packages})
        write_json(self.out / "analysis_config.json", self.config)
        files = sorted(p for p in self.raw.rglob("*") if p.is_file() and (p.parent.name in ["U-CAN", "AC-ICAM"] or p.suffix in [".tsv", ".json", ".gz"]) and "renv" not in p.parts)
        files += [PROJECT / "tables/HGNC/hgnc_complete_set.txt"] + list((PROJECT / "tables/GTF").glob("*/gene_id_gene_name_gene_biotype.csv"))
        rows = []
        for p in files:
            with p.open("rb") as f:
                sha = hashlib.file_digest(f, "sha256").hexdigest()
            rows.append(dict(path=str(p), bytes=p.stat().st_size, sha256=sha))
        pd.DataFrame(rows).to_csv(self.out / "input_checksums.tsv", sep="\t", index=False)

    def deseq(self, counts, meta, formula, name):
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
        from formulaic import model_matrix
        assert counts.columns.equals(meta.index)
        m = meta.copy()
        m["PM_group"] = pd.Categorical(np.where(m.PM.eq(1), "PM", "non_PM"), categories=["non_PM", "PM"])
        for c in ["site", "treatment", "MSI", "sex", "specimen"]:
            m[c] = m[c].astype("category")
        # Center/scale age for numerical stability; the PM contrast is invariant.
        age_mean = float(m.age.mean())
        m["age"] = (m.age - age_mean) / 10
        design = model_matrix(formula, m)
        if np.linalg.matrix_rank(design.to_numpy()) != design.shape[1]:
            raise ValueError(f"Nonidentifiable design: {name}")
        self.table(design, "U-CAN", f"design_{name}")
        self.table(m, "U-CAN", f"samples_{name}")
        keep = counts.ge(self.config["min_count"]).sum(axis=1).ge(self.config["min_samples"])
        dds = DeseqDataSet(counts=counts.loc[keep].T.astype(int), metadata=m, design=formula,
                           n_cpus=self.config["workers"], refit_cooks=True)
        dds.deseq2()
        ds = DeseqStats(dds, contrast=["PM_group", "PM", "non_PM"], alpha=self.config["fdr"], n_cpus=self.config["workers"])
        ds.summary()
        self.table(dds.var, "U-CAN", f"fit_diagnostics_{name}")
        write_json(self.out / "U-CAN" / f"model_{name}.json", dict(formula=formula, age_center=age_mean,
            age_units="centered decades", patients=len(m), genes=int(keep.sum()),
            convergence={c: {str(k): int(v) for k, v in dds.var[c].value_counts(dropna=False).items()} for c in dds.var if "converged" in c}))
        results = ds.results_df.copy()
        results["CI_low"] = results.log2FoldChange - st.norm.ppf(.975) * results.lfcSE
        results["CI_high"] = results.log2FoldChange + st.norm.ppf(.975) * results.lfcSE
        results["symbol"] = [self.id_symbols.get(g, "") for g in results.index]
        self.table(results, "U-CAN", f"DE_{name}")
        self.save_matrix(pd.DataFrame(dds.layers["normed_counts"].T, index=dds.var_names, columns=dds.obs_names), "U-CAN", f"normalized_counts_{name}")
        return results

    def differential(self):
        self.annotation()
        counts, m = self.matrix("U-CAN", "counts"), self.clinical("U-CAN")
        main = self.deseq(counts, m, "~ site + treatment + MSI + PM_group", "main")
        subsets = {"untreated": (m.treatment.eq("Untreated"), "~ site + MSI + PM_group"),
                   "colon": (m.site.eq("Colon"), "~ treatment + MSI + PM_group"),
                   "expanded": (pd.Series(True, index=m.index), "~ site + treatment + MSI + age + sex + specimen + PM_group")}
        sensitivity = []
        for name, (keep, formula) in subsets.items():
            result = self.deseq(counts.loc[:, m.index[keep]], m.loc[keep], formula, name)
            sensitivity.append(result[["log2FoldChange", "padj"]].rename(columns=lambda c: f"{name}_{c}"))
        joined = main.join(sensitivity)
        signature = joined.loc[joined.padj.lt(self.config["fdr"]) & joined.log2FoldChange.ge(self.config["min_log2fc"]) & joined.symbol.ne("")].sort_values(["padj", "symbol"])
        joined["signature_member"] = joined.index.isin(signature.index)
        self.table(joined, "U-CAN", "DE_with_sensitivity")
        self.table(signature, "signature", "PM_up_signature")
        (self.out / "signature/PM_up_ensembl.txt").write_text("".join(f"{g}\n" for g in signature.index))
        (self.out / "signature/PM_up_symbols.txt").write_text("".join(f"{g}\n" for g in signature.symbol))
        self.table(joined.loc[joined.padj.lt(self.config["fdr"]) & joined.log2FoldChange.le(-self.config["min_log2fc"])], "signature", "PM_down_genes")
        write_json(self.out / "signature/definition.json", dict(source="U-CAN main only", formula="site + treatment + MSI + PM", fdr=self.config["fdr"], min_log2fc=self.config["min_log2fc"], genes=len(signature), annotation_requirement="unique approved HGNC Ensembl-symbol mapping", score="mean sample-wise percentile rank over fixed shared background", threshold_fallback=False))
        self.volcano(main, "U-CAN", "log2FoldChange", "log2 fold change (PM / non-PM)")
        if len(signature):
            self.heatmap(np.log2(self.matrix("U-CAN", "tpm") + 1), m, signature.index[:50], "U-CAN")

    def volcano(self, results, cohort, effect, label, significance="padj", title_note=""):
        d = results.copy()
        d["minus_log10_p"] = -np.log10(d.pvalue.clip(lower=1e-300))
        d["group"] = "other"
        d.loc[d[significance].lt(.05) & d[effect].gt(0), "group"] = "PM higher"
        d.loc[d[significance].lt(.05) & d[effect].lt(0), "group"] = "PM lower"
        fig, ax = plt.subplots(figsize=(7, 5))
        sns.scatterplot(data=d, x=effect, y="minus_log10_p", hue="group", palette={"other": "#b8bec4", "PM higher": "#c84343", "PM lower": "#307ca8"}, s=12, linewidth=0, ax=ax)
        labels = d.dropna(subset=["pvalue"]).nsmallest(8, "pvalue").sort_values("minus_log10_p")
        last_y = -np.inf
        gap = max(.15, float(d.minus_log10_p.max()) * .035)
        for _, r in labels.iterrows():
            text_y = max(r.minus_log10_p + gap/3, last_y + gap)
            ax.annotate(str(r.symbol), (r[effect], r.minus_log10_p), fontsize=7,
                xytext=(r[effect]+.08, text_y), textcoords="data",
                arrowprops={"arrowstyle": "-", "color": "#777777", "lw": .5})
            last_y = text_y
        ax.set(xlabel=label, ylabel="−log10(raw P)", title=f"{cohort}: PM versus other metastases" + (f"\n{title_note}" if title_note else ""))
        if significance == "pvalue":
            ax.axhline(-np.log10(.05), linestyle="--", linewidth=.8, color="#777777")
        self.figure(fig, cohort, "volcano")

    def heatmap(self, expression, meta, genes, cohort, plot_name="signature_heatmap", signature_title="frozen U-CAN PM-up signature"):
        genes = [g for g in genes if g in expression.index]
        if not genes:
            return
        ordered = meta.sort_values(["PM", "patient_id"]).index
        x = expression.loc[genes, ordered]
        z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1).replace(0, 1), axis=0)
        z.index = [self.id_symbols.get(g, g) for g in z.index]
        colors = pd.DataFrame({"PM": meta.PM.map({0: "#307ca8", 1: "#c84343"})}, index=meta.index)
        legends = {"PM": {0: "#307ca8", 1: "#c84343"}}
        for col in ["site", "treatment", "MSI"]:
            if col in meta:
                cats = sorted(meta[col].dropna().astype(str).unique())
                palette = dict(zip(cats, sns.color_palette("Set2", len(cats)).as_hex()))
                colors[col] = meta[col].map(palette).fillna("#dddddd")
                legends[col] = palette
        grid = sns.clustermap(z, cmap="vlag", center=0, vmin=-2.5, vmax=2.5, row_cluster=len(z)>1,
                              col_cluster=False, col_colors=colors.loc[ordered], xticklabels=False,
                              yticklabels=len(z)<=100, dendrogram_ratio=(.16,.02), colors_ratio=(.02,.025),
                              figsize=(12, max(4, .19 * len(z) + 2) if len(z)<=100 else 14), cbar_pos=(1.04,.58,.02,.16),
                              cbar_kws={"label": "gene-wise z score"})
        grid.fig.suptitle(f"{cohort}: {signature_title}", y=1.02)
        from matplotlib.patches import Patch
        handles = []
        for field, palette in legends.items():
            for value, color in palette.items():
                label = {0: "Other metastases", 1: "PM"}.get(value, value) if field == "PM" else value
                handles.append(Patch(facecolor=color, label=f"{field}: {label}"))
        grid.fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.03,.98), frameon=False, fontsize=9)
        self.figure(grid.fig, cohort, plot_name)
        write_json(self.out / cohort / "heatmap_annotation_colors.json", legends)
        self.table(z, cohort, "heatmap_values" if plot_name == "signature_heatmap" else f"{plot_name}_values")

    def figures(self):
        """Regenerate DE figures without refitting statistical models."""
        self.annotation()
        signature = self.clinical("signature", "PM_up_signature")
        self.volcano(self.clinical("U-CAN", "DE_main"), "U-CAN", "log2FoldChange", "log2 fold change (PM / non-PM)")
        self.volcano(self.clinical("AC-ICAM", "DE_validation"), "AC-ICAM", "effect", "Difference on public expression scale (PM − other)")
        if len(signature):
            self.heatmap(np.log2(self.matrix("U-CAN", "tpm")+1), self.clinical("U-CAN"), signature.index[:50], "U-CAN")
            self.heatmap(self.matrix("AC-ICAM", "expression"), self.clinical("AC-ICAM"), signature.index[:50], "AC-ICAM")

    def score_cohort(self, cohort):
        signature = self.clinical("signature", "PM_up_signature").index.tolist()
        if not signature:
            write_json(self.out / cohort / "signature_status.json", dict(status="skipped_empty_discovery_signature"))
            return None
        bg = pd.read_csv(self.out / "signature/ranking_background.tsv", sep="\t").ensembl_id.tolist()
        scoring_genes = [g for g in signature if g in set(bg)]
        pd.DataFrame({"ensembl_id": scoring_genes}).to_csv(self.out / "signature/scoring_genes.tsv", sep="\t", index=False)
        expr = self.matrix(cohort, "expression" if cohort == "AC-ICAM" else "tpm")
        try:
            score, coverage = rank_score(expr, signature, bg, self.config["min_signature_coverage"])
        except ValueError as e:
            write_json(self.out / cohort / "signature_status.json", dict(status="not_scoreable", reason=str(e)))
            return None
        write_json(self.out / cohort / "signature_status.json", dict(status="scored", coverage=coverage, missing=[g for g in signature if g not in bg]))
        m = self.clinical(cohort).join(score.rename("PM_signature_score"))
        self.table(m, cohort, "signature_scores")
        return m

    def validate(self):
        self.annotation()
        expression, m = self.matrix("AC-ICAM", "expression"), self.clinical("AC-ICAM")
        p, n = expression.loc[:, m.index[m.PM.eq(1)]], expression.loc[:, m.index[m.PM.eq(0)]]
        # Exact OLS two-group coefficient/t test, vectorized; public scale is not relabeled as log2FC.
        effect = p.mean(axis=1) - n.mean(axis=1)
        df = p.shape[1] + n.shape[1] - 2
        variance = ((p.shape[1]-1)*p.var(axis=1) + (n.shape[1]-1)*n.var(axis=1)) / df
        se = np.sqrt(variance * (1/p.shape[1] + 1/n.shape[1]))
        statistic = effect.div(se.replace(0, np.nan))
        result = pd.DataFrame(dict(effect=effect, SE=se, stat=statistic, pvalue=2*st.t.sf(abs(statistic), df)))
        result["padj"] = np.nan
        valid = result.pvalue.notna()
        result.loc[valid, "padj"] = multipletests(result.loc[valid, "pvalue"], method="fdr_bh")[1]
        result["CI_low"] = effect - st.t.ppf(.975, df)*se
        result["CI_high"] = effect + st.t.ppf(.975, df)*se
        result["symbol"] = [self.id_symbols.get(g, "") for g in result.index]
        self.table(result, "AC-ICAM", "DE_validation")
        self.volcano(result, "AC-ICAM", "effect", "Difference on public expression scale (PM − other)")
        signature = self.clinical("signature", "PM_up_signature")
        gene_validation = signature.join(result.add_prefix("AC_"))
        gene_validation["same_direction"] = gene_validation.AC_effect.gt(0).where(gene_validation.AC_effect.notna())
        self.table(gene_validation, "signature", "independent_gene_validation")
        self.score_cohort("U-CAN")
        scores = self.score_cohort("AC-ICAM")
        if scores is None:
            return
        v1 = scores.loc[scores.PM.eq(1), "PM_signature_score"].to_numpy()
        v0 = scores.loc[scores.PM.eq(0), "PM_signature_score"].to_numpy()
        observed = float(v1.mean() - v0.mean())
        rng = np.random.default_rng(self.config["seed"])
        boot = [rng.choice(v1,len(v1)).mean()-rng.choice(v0,len(v0)).mean() for _ in range(self.config["bootstrap_iterations"])]
        combined = np.concatenate([v1, v0])
        null = []
        for _ in range(self.config["permutation_iterations"]):
            perm = rng.permutation(combined)
            null.append(perm[:len(v1)].mean()-perm[len(v1):].mean())
        perm_p = (1+np.sum(np.abs(null) >= abs(observed)))/(len(null)+1)
        write_json(self.out / "AC-ICAM/signature_validation.json", dict(PM_n=len(v1), other_n=len(v0), mean_difference=observed,
            bootstrap_CI=np.quantile(boot,[.025,.975]).tolist(), permutation_two_sided_p=perm_p,
            genes_assessable=int(gene_validation.AC_effect.notna().sum()), genes_same_direction=int(gene_validation.AC_effect.gt(0).sum())))
        fig, ax = plt.subplots(figsize=(5, 5))
        sns.boxplot(data=scores, x="PM", y="PM_signature_score", showfliers=False, ax=ax, color="#dce4ea")
        sns.stripplot(data=scores, x="PM", y="PM_signature_score", jitter=.15, ax=ax, color="#3b536a")
        ax.set(xticks=[0,1], xticklabels=["Other metastases", "PM"], xlabel="", ylabel="Frozen PM-up score", title=f"AC-ICAM: permutation P={perm_p:.3g}")
        self.figure(fig, "AC-ICAM", "signature_score_validation")
        self.heatmap(expression, m, signature.index[:50], "AC-ICAM")

    def survival(self):
        from lifelines import KaplanMeierFitter, CoxPHFitter
        from lifelines.statistics import logrank_test, proportional_hazard_test
        from lifelines.plotting import add_at_risk_counts
        scores = self.score_cohort("TCGA-COAD")
        if scores is None:
            return
        eligible = scores.OS_exclusion.fillna("").eq("") & scores.PM_signature_score.notna()
        d = scores.loc[eligible].copy()
        cutoff = float(d.PM_signature_score.median())
        d["score_group"] = np.where(d.PM_signature_score.gt(cutoff), "High", "Low")
        sd = d.PM_signature_score.std(ddof=1)
        if sd == 0 or d.score_group.nunique() != 2:
            raise ValueError("No signature score variation for survival analysis")
        d["score_z"] = (d.PM_signature_score-d.PM_signature_score.mean())/sd
        self.table(d, "TCGA-COAD", "OS_analysis_patients")
        hi, lo = d.loc[d.score_group.eq("High")], d.loc[d.score_group.eq("Low")]
        lr = logrank_test(hi.OS_days, lo.OS_days, hi.OS_event, lo.OS_event)
        fig, ax = plt.subplots(figsize=(7, 6))
        fitters = []
        for label, group, color in [("High", hi, "#c84343"),("Low", lo, "#307ca8")]:
            fit = KaplanMeierFitter(label=f"{label} (n={len(group)})").fit(group.OS_days/365.25, group.OS_event)
            fit.plot_survival_function(ax=ax, color=color, ci_show=True)
            fitters.append(fit)
        add_at_risk_counts(*fitters, ax=ax)
        ax.set(xlabel="Years from diagnosis", ylabel="Overall survival probability", title=f"TCGA-COAD: log-rank P={lr.p_value:.3g}", ylim=(0,1.02))
        self.figure(fig, "TCGA-COAD", "OS_KM")
        summary = dict(endpoint="OS", n=len(d), events=int(d.OS_event.sum()), excluded=int((~eligible).sum()), cutoff=cutoff,
            cutoff_population="OS-eligible patients", high_n=len(hi), low_n=len(lo), logrank_p=float(lr.p_value), score_SD=float(sd))
        for label, cols, formula in [("unadjusted", ["OS_days","OS_event","score_z"], "score_z"),
            ("adjusted", ["OS_days","OS_event","score_z","age","sex","stage"], "score_z + age + C(sex) + C(stage)")]:
            data = d[cols].dropna()
            model = CoxPHFitter().fit(data, duration_col="OS_days", event_col="OS_event", formula=formula)
            self.table(model.summary, "TCGA-COAD", f"OS_Cox_{label}")
            self.table(data, "TCGA-COAD", f"OS_Cox_{label}_patients")
            ph = proportional_hazard_test(model, data, time_transform="rank")
            self.table(ph.summary, "TCGA-COAD", f"OS_PH_{label}")
            summary[label] = dict(n=len(data), events=int(data.OS_event.sum()), HR=float(model.summary.loc["score_z","exp(coef)"]),
                CI_low=float(model.summary.loc["score_z","exp(coef) lower 95%"]), CI_high=float(model.summary.loc["score_z","exp(coef) upper 95%"]),
                p=float(model.summary.loc["score_z","p"]), PH_violations=ph.summary.index[ph.summary.p.lt(.05)].tolist())
            with warnings.catch_warnings(record=True) as captured:
                import contextlib, io
                output = io.StringIO()
                with contextlib.redirect_stdout(output):
                    model.check_assumptions(data, p_value_threshold=.05, show_plots=False)
                (self.out / "TCGA-COAD" / f"OS_PH_{label}.txt").write_text(output.getvalue()+"\n"+"\n".join(str(w.message) for w in captured))
        write_json(self.out / "TCGA-COAD/OS_summary.json", summary)

    def enrichment(self):
        import gseapy as gp
        r = self.clinical("U-CAN", "DE_main").dropna(subset=["stat", "pvalue"])
        r = r.loc[r.symbol.notna() & r.symbol.ne("")].sort_values("stat", ascending=False)
        ranks = r.drop_duplicates("symbol").set_index("symbol").stat
        self.table(ranks.to_frame(), "U-CAN", "GSEA_ranking")
        manifests = []
        for label, url in self.config["gene_sets"].items():
            path = self.out / "signature" / f"{label}.gmt"
            if not path.exists():
                request = urllib.request.Request(url, headers={"User-Agent": "crc-pm-research/1.0"})
                with urllib.request.urlopen(request, timeout=60) as response:
                    path.write_bytes(response.read())
            genesets = {}
            for line in path.read_text().splitlines():
                parts = line.split("\t")
                genesets[parts[0]] = [g.split(",")[0] for g in parts[2:] if g]
            if not genesets or len(genesets) < 10:
                raise ValueError(f"Empty gene set library: {label}")
            manifests.append(dict(library=label, url=url, sha256=hashlib.sha256(path.read_bytes()).hexdigest(), sets=len(genesets)))
            result = gp.prerank(rnk=ranks, gene_sets=genesets, outdir=None, min_size=15, max_size=500,
                permutation_num=self.config["gsea_permutations"], threads=self.config["workers"], seed=self.config["seed"], no_plot=True)
            table = result.res2d.sort_values("FDR q-val")
            self.table(table, "U-CAN", f"GSEA_{label}", index=False)
            plot = pd.concat([table.loc[table.NES.gt(0)].head(8), table.loc[table.NES.lt(0)].head(8)]).drop_duplicates("Term").sort_values("NES")
            fig, ax = plt.subplots(figsize=(8, max(4,.45*len(plot))))
            labels = [textwrap.fill(re.sub(r"^(HALLMARK_|GOBP_)", "", term).replace("_", " "), 40) for term in plot.Term]
            ax.barh(range(len(plot)), plot.NES, color=np.where(plot.NES.gt(0), "#c84343", "#307ca8"))
            ax.set_yticks(range(len(plot)), labels=labels, fontsize=8)
            for i, (_, row) in enumerate(plot.iterrows()):
                ax.text(row.NES, i, f"  q={row['FDR q-val']:.2g}", va="center", ha="left" if row.NES>0 else "right", fontsize=7)
            ax.set(xlabel="Normalized enrichment score", title=f"U-CAN {label}: top ranked pathways (see table for FDR)")
            self.figure(fig, "U-CAN", f"GSEA_{label}")
        write_json(self.out / "signature/gene_set_sources.json", manifests)

    def report(self):
        prep = json.loads((self.out / "preparation_summary.json").read_text())
        r = self.clinical("U-CAN", "DE_with_sensitivity")
        sig = self.clinical("signature", "PM_up_signature")
        lines = ["# 腹膜转移相关 bulk 转录组分析", "", "## 设计与数据", "",
            "主分析比较IV期患者的原发灶：U-CAN发现24 PM对90非PM；AC-ICAM独立验证7 PM对53其他转移部位。未知PM标签、正常组织和非IV期样本不进入发现/验证比较。", "",
            f"TCGA按患者去重后保留{prep['TCGA_patients']}例原发肿瘤，OS数据规则通过者{prep['TCGA_OS_eligible']}例。", "",
            "## 差异表达与 signature", "",
            "U-CAN用PyDESeq2负二项模型，调整结肠/直肠、预治疗和MSI。保留至少18个样本counts≥10的基因；BH-FDR以PyDESeq2独立过滤后的可检验基因为范围。",
            f"可检验结果{r.pvalue.notna().sum()}个；FDR<0.05基因{r.padj.lt(.05).sum()}个；固定PM上调signature为{len(sig)}个基因（FDR<0.05、log2FC≥0.5且唯一标准注释）。", ""]
        duplicates = pd.read_csv(self.out / "TCGA-COAD/duplicate_file_comparison.tsv", sep="\t")
        lines += [f"TCGA有{duplicates.sample_id.nunique()}个样本存在多个表达文件；其中{duplicates.loc[~duplicates.identical_expression, 'sample_id'].nunique()}个与首个文件表达不一致。按预先确定的样本/文件ID规则选择一个文件，未将技术文件当独立患者。具体选择及比较见TCGA的sample_selection和duplicate_file_comparison表。", ""]
        if len(sig):
            lines += ["前列基因："+"、".join(sig.symbol.head(15)), ""]
            score_path = self.out / "signature/scoring_genes.tsv"
            if score_path.exists():
                lines += [f"完整signature为{len(sig)}个基因；跨队列可比评分使用共同背景内的{len(pd.read_csv(score_path, sep=chr(9)))}个基因。完整清单不因验证覆盖而缩减，便于下游解析细胞来源。", ""]
            for name in ["untreated", "colon", "expanded"]:
                assessed = sig[f"{name}_log2FoldChange"].notna()
                lines += [f"{name}敏感性分析：{assessed.sum()}个可评估signature基因中，{sig.loc[assessed,f'{name}_log2FoldChange'].gt(0).sum()}个效应仍为正。"]
        else:
            lines += ["无基因通过预设signature阈值；未放宽阈值。signature评分和KM分析跳过。"]
        vpath = self.out / "AC-ICAM/signature_validation.json"
        if vpath.exists():
            v = json.loads(vpath.read_text())
            lines += ["", "## 独立验证", "", f"AC-ICAM平均评分差异（PM−其他）{v['mean_difference']:.4f}，患者bootstrap 95%区间{v['bootstrap_CI']}，双侧置换P={v['permutation_two_sided_p']:.4g}。",
                f"可评估signature基因{v['genes_assessable']}个，其中{v['genes_same_direction']}个方向一致。验证集只有7例PM，应结合效应区间理解不确定性。",
                "AC-ICAM使用公开连续表达尺度，不重新取对数，也不将效应称为log2FC；逐基因为两组OLS对比的等价向量化计算。"]
        spath = self.out / "TCGA-COAD/OS_summary.json"
        if spath.exists():
            s = json.loads(spath.read_text())
            lines += ["", "## TCGA-COAD OS", "", f"OS分析{ s['n']}例，{s['events']}例死亡。按OS可用患者评分中位数分组，高分{s['high_n']}例、低分{s['low_n']}例，log-rank P={s['logrank_p']:.4g}。"]
            for label in ["unadjusted", "adjusted"]:
                x=s[label]
                lines += [f"{label} Cox（每标准差评分）：n={x['n']}，HR={x['HR']:.3f}，95%CI {x['CI_low']:.3f}–{x['CI_high']:.3f}，P={x['p']:.4g}；PH检查标记：{x['PH_violations']}。"]
            lines += ["调整模型包含年龄、性别和分期，仅使用完整协变量患者；评分标准差固定为主OS人群的标准差。PH违反时，常数HR的解释需谨慎，不能将其解释为全随访期间稳定效应。"]
            if s["adjusted"]["p"] >= .05 and s["logrank_p"] >= .05:
                lines += ["当前TCGA结果未显示该固定PM signature与OS存在显著关联。未为获得显著性而改变基因集、终点或分组阈值。"]
        for label in self.config["gene_sets"]:
            path = self.out / "U-CAN" / f"GSEA_{label}.tsv"
            if path.exists():
                d = pd.read_csv(path, sep="\t")
                positive = d.loc[d["FDR q-val"].lt(.05) & d.NES.gt(0)]
                negative = d.loc[d["FDR q-val"].lt(.05) & d.NES.lt(0)]
                lines += ["", f"## {label}富集", "", f"FDR<0.05：PM正向富集{len(positive)}项，负向{len(negative)}项。", "正向前列："+"；".join(positive.Term.head(8)), "负向前列："+"；".join(negative.Term.head(8))]
                if label == "Hallmark":
                    sig_terms = d.loc[d["FDR q-val"].lt(.05)].set_index("Term").NES
                    high = [t for t in ["HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_COMPLEMENT", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_HYPOXIA"] if sig_terms.get(t, 0)>0]
                    low = [t for t in ["HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_G2M_CHECKPOINT"] if sig_terms.get(t, 0)<0]
                    lines += ["", "生物学解读：PM原发灶中正向富集的程序包括"+"、".join(t.replace("HALLMARK_", "") for t in high)+"；负向包括"+"、".join(t.replace("HALLMARK_", "") for t in low)+"。这提示炎症/补体及相关应激信号与增殖程序存在差异；bulk结果尚不能区分细胞组成变化与细胞内状态变化。"]
                    lines += ["GSEA采用预排序的基因集置换，和患者标签置换检验不同。表中出现的0为有限置换下的软件估计，并非真实概率为0。"]
        lines += ["", "## 解释边界与可重复性", "", "该gene set描述已发生远处转移患者中的PM相关整体组织表达，不能直接视为未来PM预测器或因果基因。TME、免疫和肿瘤来源信号均保留，供下游解析。TCGA缺少可靠PM标签，OS关联不能证明PM特异性。",
            "signature只由U-CAN确定，AC-ICAM和TCGA不参与基因选择。排名背景固定为三队列共有的唯一标准映射基因；平均样本内百分位秩评分，最低覆盖80%。无显著结果也是有效结果。",
            "详见各队列TSV、PDF/PNG、signature定义、版本、配置和输入校验值。预先存在的其他项目改动未触碰。", "",
            "## 方法与基因集来源", "",
            "采用[Cañellas-Socias 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC7616986/)先定义整体关联基因集、再解析细胞来源的顺序；原论文研究复发，本研究以IV期PM标签为发现依据，不是对其复发模型或基因数量的直接复现。",
            "Counts建模：[PyDESeq2](https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_step_by_step.html)；生存建模：[lifelines](https://lifelines.readthedocs.io/en/stable/Quickstart.html)；通路基因集：[Broad MSigDB Human 2025.1](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/)。下载URL和散列保存在signature/gene_set_sources.json。", ""]
        (self.out / "analysis_report_zh.md").write_text("\n".join(lines))
