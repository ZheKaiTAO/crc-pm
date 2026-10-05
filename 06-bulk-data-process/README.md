# PM bulk transcriptome analysis (Python)

Use the project's `pm` Python environment. Run from `code/`:

```bash
/home/zhekai/miniconda3/envs/pm/bin/python 06-bulk-data-process/run_analysis.py --stage all
```

Stages, in dependency order: `prepare`, `differential`, `validate`, `survival`, `enrichment`, `report`.
Individual stages reuse preceding outputs; rerun downstream stages after changing any upstream input or config.
`--stage figures` regenerates volcano plots and heatmaps without refitting models.
The original inputs in `../bulk-data` and annotation in `../tables` are read-only. Generated artifacts are in `../result/bulk`.

Dependencies: numpy, pandas, scipy, statsmodels, pydeseq2, lifelines, gseapy, matplotlib, seaborn,
scikit-learn, formulaic, openpyxl, pyarrow. Exact runtime versions are written to `versions.json`.
The user manages installations. No R, rpy2 or shared-environment modifications are used.

## Fixed statistical design

- U-CAN stage IV, explicit PM labels only (24 PM / 90 non-PM). PyDESeq2 adjusts site, pre-treatment and MSI.
- Sensitivity: untreated only, colon only, and expanded age/sex/specimen adjustment. Nonidentifiable models fail explicitly.
- Retain genes with counts >=10 in at least 18 samples. Signature: main BH-FDR <0.05, log2FC >=0.5,
  uniquely annotated HGNC Ensembl-symbol pairs. Gene set is frozen before either validation cohort is inspected for association.
- PyDESeq2 uses Cook's filtering and independent filtering. Reported FDR is its adjusted `padj`, not a custom recalculation.
- AC-ICAM public expression is treated as continuous, with no extra log. Duplicate symbols/IDs are averaged on that scale.
  Vectorized two-group OLS coefficients and pooled-residual t tests are mathematically equivalent to an intercept+PM statsmodels OLS fit.
- Score: mean within-sample percentile rank of available signature genes over a fixed common ranking background; minimum coverage 80%.
  Ranking background uses available gene IDs, not validation labels or survival outcomes. Ties use average ranks.
- AC score difference: 5,000 patient bootstrap draws and 10,000 two-sided label permutations; gene-level validation never edits the signature.
- TCGA: one primary sample/patient, prioritize 01A then sample/file IDs; all selection decisions and duplicate comparisons exported.
  OS uses reliable death time or maximum available follow-up; conflicts are excluded, never replaced by follow-up for deaths.
- KM uses median among OS-eligible patients (high > median; low <= median). Cox reports score per baseline-population SD,
  unadjusted and age/sex/stage adjusted, with PH diagnostics. No best-cutpoint search or survival-driven gene selection.
- GSEA ranks U-CAN genes by signed Wald statistic, uses 10,000 permutations and fixed downloaded GMT snapshots.
  Libraries are official MSigDB Human 2025.1 Hallmark and GO BP GMT files from Broad Institute.
  First enrichment run requires network access to Broad Institute. Subsequent runs reuse the hashed snapshot.

Configuration is `config.json`; output `analysis_config.json`, input hashes and package versions preserve the exact run.
Technical QC includes finite/nonnegative matrices, unique IDs, clinical alignment, PCA, detected-gene totals and model-rank checks.
PCA outliers alone are not removed. No significant genes or survival findings is a valid outcome.

## Tests

```bash
MPLCONFIGDIR=/tmp/crc-pm-bulk-mpl /home/zhekai/miniconda3/envs/pm/bin/python -m unittest discover -s 06-bulk-data-process -p 'test_pipeline.py' -v
```

Inspect `../result/bulk/analysis_report_zh.md`, cohort tables and PDF/PNG plots.
The notebook `bulk-readin.ipynb` is a read-only output browser, not a second implementation.

## Requested single-gene follow-up

```bash
/home/zhekai/miniconda3/envs/pm/bin/python 06-bulk-data-process/followup_analysis.py --genes CLDN18 CD55
```

Uses the existing deduplicated TCGA primary tumors and OS endpoint. Each gene has a fixed median TPM split for KM;
continuous Cox uses standardized log2(TPM+1), unadjusted and age/sex/stage adjusted. The baseline OS-population SD
is reused for adjusted models. Raw P and BH q across the two genes in each test family are exported; no best-cutpoint search.
Results, patient tables, PH diagnostics and KM plots are in `TCGA-COAD/single_gene/`.

`signature/threshold_audit.tsv` compares U-CAN selection rules. `PM_up_nominal_candidates_exploratory.tsv` contains
nominal-P candidates, explicitly exploratory. These files do not change the frozen 64-gene signature or original scoring.
The follow-up rationale and results are in `../result/bulk/followup_report_zh.md`.

## Canellas nominal-P threshold variant

```bash
/home/zhekai/miniconda3/envs/pm/bin/python 06-bulk-data-process/canellas_threshold_analysis.py
```

Reuses the fitted U-CAN PM model, selects raw P <0.05 and positive log2FC without an FDR or fold-change magnitude
filter, then recalculates cohort scores, AC validation and TCGA OS. All artifacts stay under
`../result/bulk/canellas_threshold/`; the earlier strict 64-gene outputs are not overwritten.
The full 1,158-gene Ensembl list retains one ID without a unique standard symbol, with annotation status recorded.
Full-set heatmaps and readable top-50 heatmaps are both exported. ORA uses the tested annotated U-CAN genes as
background and BH correction per library. Whole-transcriptome GSEA is reused because its ranking is unchanged.
The threshold variant is an extension made after validation outcomes were inspected; its results are exploratory.
It adapts the paper's nominal-P/positive-effect rule to PM differential expression, not to a relapse Cox endpoint.
