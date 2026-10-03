"""Audit mappings using only the H5AD var index (never the expression matrix)."""
import argparse
import csv
from pathlib import Path
import sys

import h5py
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'code/src'))
from GeneSetMapping import assess_coverage, get_resolver


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--adata', type=Path, default=ROOT / 'result/dr-anno/ready_for_anno.h5ad')
    parser.add_argument('--output', type=Path, default=ROOT / 'result/gene-set-mapping')
    args = parser.parse_args()
    with h5py.File(args.adata, 'r') as handle:
        var = handle['var']
        names = var[var.attrs.get('_index', '_index')].asstr()[:]
    args.output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for path in sorted((ROOT / 'tables/geneset').glob('*.csv')):
        with path.open() as handle:
            genes = [row[0] for row in csv.reader(handle) if row]
        report = get_resolver().resolve_genes(genes)
        _, report, summary = assess_coverage(report, names)
        report.to_csv(args.output / f'{path.stem}.csv', index=False)
        summaries.append({'gene_set': path.stem, 'input_rows': len(genes), **summary})
    summary = pd.DataFrame(summaries)
    summary.to_csv(args.output / 'summary.csv', index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
