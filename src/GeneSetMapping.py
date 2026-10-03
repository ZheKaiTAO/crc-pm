# AIGC

"""Offline, auditable gene-set mapping; independent of loader mappings."""
from collections import defaultdict
from functools import cache
from pathlib import Path
import csv
import hashlib
import json
import re

import pandas as pd

TABLES = Path(__file__).resolve().parents[2] / 'tables'
ENSG = re.compile(r'ENSG\d{11}(?:\.\d+)?\Z')
REPORT_COLUMNS = ['original', 'name', 'standard_symbol', 'gene_id', 'source', 'status', 'candidates']


class GeneSetResolver:
    def __init__(self, hgnc_path, gtf_paths):
        self.approved = {}
        self.aliases = defaultdict(set)
        self.gtf = defaultdict(lambda: defaultdict(set))
        with Path(hgnc_path).open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row['status'] != 'Approved':
                    continue
                symbol = row['symbol']
                self.approved[symbol] = {
                    value.split('.')[0] for value in row['ensembl_gene_id'].split('|')
                    if ENSG.fullmatch(value)
                }
                for field in ('prev_symbol', 'alias_symbol'):
                    for alias in row[field].split('|'):
                        if alias:
                            self.aliases[alias].add((symbol, field))
        for version, path in gtf_paths.items():
            with Path(path).open() as handle:
                for row in csv.DictReader(handle):
                    gene_id = row['gene_id'].split('.')[0]
                    if row['gene_name'] and ENSG.fullmatch(gene_id):
                        self.gtf[row['gene_name']][gene_id].add(version)

    def resolve(self, name):
        if ENSG.fullmatch(name):
            return '', name.split('.')[0], 'input_ENSG', 'mapped', ''
        if not name or name.isdecimal():
            return '', '', '', 'invalid', ''
        # Approved symbols take precedence over their use as another gene's alias.
        matches = {(name, 'symbol')} if name in self.approved else self.aliases.get(name, set())
        symbols = {symbol for symbol, _ in matches}
        candidates, sources = set(), set()
        for symbol, field in matches:
            candidates.update(self.approved[symbol])
            sources.add('HGNC:' + field)
        # Retain every candidate; loader biotype priorities must not hide ambiguity.
        for symbol in symbols | {name}:
            for gene_id, versions in self.gtf.get(symbol, {}).items():
                candidates.add(gene_id)
                sources.update('GTF:' + version for version in versions)
        source, standard = '|'.join(sorted(sources)), '|'.join(sorted(symbols))
        if len(symbols) > 1 or len(candidates) > 1:
            return standard, '', source, 'ambiguous', '|'.join(sorted(candidates))
        if candidates:
            return standard, next(iter(candidates)), source, 'mapped', ''
        return standard, '', source, 'unresolved', ''

    def resolve_genes(self, genes):
        rows = []
        for original in genes:
            original = str(original).strip()
            for name in original.split('|||'):
                name = name.strip()
                standard, gene_id, source, status, candidates = self.resolve(name)
                rows.append([original, name, standard, gene_id, source, status, candidates])
        return pd.DataFrame(rows, columns=REPORT_COLUMNS)


@cache
def get_resolver():
    snapshot = TABLES / 'HGNC/hgnc_complete_set.txt'
    metadata = json.loads((snapshot.parent / 'metadata.json').read_text())
    if hashlib.sha256(snapshot.read_bytes()).hexdigest() != metadata['sha256']:
        raise ValueError('HGNC snapshot checksum does not match metadata.json')
    paths = {
        version: TABLES / 'GTF' / version / 'gene_id_gene_name_gene_biotype.csv'
        for version in ('GRCh38.111', 'GRCh38.98', 'GRCh38.84')
    }
    return GeneSetResolver(snapshot, paths)


def assess_coverage(report, var_names):
    """Count unique resolved IDs plus unique unresolved names, never None entries."""
    report = report.copy()
    present = set(var_names)
    mapped = report['status'].eq('mapped')
    report['in_dataset'] = pd.Series(pd.NA, index=report.index, dtype='boolean')
    report.loc[mapped, 'in_dataset'] = report.loc[mapped, 'gene_id'].isin(present)
    ids = list(dict.fromkeys(report.loc[mapped, 'gene_id']))
    unresolved = set(report.loc[~mapped, 'name'])
    available = [gene for gene in ids if gene in present]
    total = len(ids) + len(unresolved)
    summary = {
        'resolved_unique_ids': len(ids),
        'unresolved_unique_names': len(unresolved),
        'ambiguous_unique_names': int(report.loc[report.status.eq('ambiguous'), 'name'].nunique()),
        'missing_dataset_ids': len(ids) - len(available),
        'available_ids': len(available),
        'total': total,
        'missing_fraction': (total - len(available)) / total if total else None,
    }
    return available, report, summary
