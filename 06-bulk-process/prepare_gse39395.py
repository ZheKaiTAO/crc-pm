"""Prepare the processed GSE39395 matrix and audit a supplied PM signature.

Uses the downloaded RMA log2 values as supplied. No DE is performed here.
"""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path

import pandas as pd


def prepare(data_dir, signature_path, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {}
    with gzip.open(data_dir / "GSE39395_series_matrix.txt.gz", "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                row = next(csv.reader([line], delimiter="\t"))
                if row[0] == "!Sample_characteristics_ch1":
                    key = row[1].split(": ", 1)[0]
                    metadata[key] = [value.split(": ", 1)[1] for value in row[1:]]
                else:
                    metadata[row[0]] = row[1:]
            if line.startswith("!series_matrix_table_begin"):
                expression = pd.read_csv(handle, sep="\t", comment="!", index_col=0)
                break
    samples = pd.DataFrame({
        "sample_id": metadata["!Sample_geo_accession"],
        "title": metadata["!Sample_title"],
        "patient": metadata["individual"],
        "compartment": [
            {"epithelial cells": "epi", "leukocytes": "immune", "stromal cells": "stromal"}[value]
            for value in metadata["cell type"]
        ],
    }).set_index("sample_id").loc[expression.columns]
    with gzip.open(data_dir / "GPL570.annot.gz", "rt") as handle:
        for line in handle:
            if line.startswith("ID\t"):
                annotation = pd.read_csv(io.StringIO(line + handle.read()), sep="\t", comment="!", dtype=str)
                break
    # Reject probes annotated to multiple genes, instead of counting them as unique coverage.
    unique = annotation.loc[
        annotation["Gene symbol"].notna()
        & ~annotation["Gene symbol"].str.contains("///", na=False)
        & ~annotation["Gene ID"].str.contains("///", na=False)
    ].copy()
    unique["gene_name"] = unique["Gene symbol"].str.strip()
    unique = unique.loc[unique.gene_name.ne("") & unique.ID.isin(expression.index)]
    unique["mean_log2_expression"] = unique.ID.map(expression.mean(axis=1))
    # Predefined, compartment-independent probe selection.
    selected = unique.sort_values(
        ["mean_log2_expression", "ID"], ascending=[False, True]
    ).drop_duplicates("gene_name")
    genes = expression.loc[selected.ID].copy()
    genes.index = pd.Index(selected.gene_name, name="gene_name")

    signature = pd.read_csv(signature_path, sep="\t", dtype={"gene_id": str})
    hgnc = pd.read_csv(
        Path.home() / "project/crc-pm/tables/HGNC/hgnc_complete_set.txt",
        sep="\t", dtype=str, usecols=["symbol", "ensembl_gene_id"],
    ).dropna().drop_duplicates("ensembl_gene_id", keep=False)
    approved = hgnc.set_index("ensembl_gene_id").symbol
    signature["approved_symbol"] = signature.gene_id.map(approved)
    signature["array_symbol"] = [
        approved_symbol if pd.notna(approved_symbol) and approved_symbol in genes.index
        else original_symbol if pd.notna(original_symbol) and original_symbol in genes.index
        else None
        for original_symbol, approved_symbol in zip(signature.gene_name, signature.approved_symbol)
    ]
    signature["covered"] = signature.array_symbol.notna()
    probe_map = selected.set_index("gene_name").ID
    signature["selected_probe"] = signature.array_symbol.map(probe_map)
    means = genes.T.groupby(samples.compartment).mean().T
    for compartment in ["epi", "immune", "stromal"]:
        signature[f"mean_log2_{compartment}"] = signature.array_symbol.map(means[compartment])
    signature["highest_mean_compartment"] = signature.array_symbol.map(means.idxmax(axis=1))
    signature["strict_fdr005"] = signature.padj.lt(0.05)
    samples.to_csv(output_dir / "samples.tsv", sep="\t")
    genes.to_csv(output_dir / "expression_gene_log2.tsv.gz", sep="\t")
    selected.to_csv(output_dir / "selected_probes.tsv", sep="\t", index=False)
    signature.to_csv(output_dir / "PM_signature_coverage.tsv", sep="\t", index=False)
    expression.loc[unique.loc[unique.gene_name.eq("CLDN18"), "ID"]].T.join(samples).to_csv(
        output_dir / "CLDN18_probes_by_sample.tsv", sep="\t"
    )
    summary = {
        "probe_matrix_shape": list(expression.shape),
        "unique_gene_matrix_shape": list(genes.shape),
        "patients": int(samples.patient.nunique()),
        "samples_per_compartment": samples.compartment.value_counts().to_dict(),
        "signature_genes": len(signature),
        "covered": int(signature.covered.sum()),
        "coverage_fraction": float(signature.covered.mean()),
        "uncovered": signature.loc[~signature.covered, ["gene_id", "gene_name"]].to_dict("records"),
        "note": "Highest mean is descriptive, not a test of compartment specificity. Multi-gene probes excluded.",
    }
    (output_dir / "coverage_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print(signature.loc[signature.gene_name.eq("CLDN18")].to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.data_dir, args.signature, args.output_dir)
