import csv
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import Eval
from src.GeneSetMapping import GeneSetResolver, assess_coverage, get_resolver


class MappingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.ids = [f'ENSG{i:011d}' for i in range(1, 5)]
        hgnc = self.root / 'hgnc.tsv'
        with hgnc.open('w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['status', 'symbol', 'ensembl_gene_id', 'prev_symbol', 'alias_symbol'])
            writer.writerow(['Approved', 'A', self.ids[0], 'OLD_A', 'SHARED'])
            writer.writerow(['Approved', 'B', self.ids[1], '', 'A|SHARED'])
            writer.writerow(['Approved', 'C', self.ids[2], '', ''])
        gtf = self.root / 'gtf.csv'
        with gtf.open('w') as f:
            writer = csv.writer(f)
            writer.writerow(['gene_name', 'gene_id'])
            writer.writerows([['OLD_A', self.ids[0]], ['LEGACY', self.ids[3]], ['CONFLICT', self.ids[0]], ['CONFLICT', self.ids[1]]])
        self.resolver = GeneSetResolver(hgnc, {'test': gtf})

    def test_symbols_aliases_conflicts_and_direct_ids(self):
        self.assertEqual(self.resolver.resolve('A')[1], self.ids[0])
        self.assertEqual(self.resolver.resolve('OLD_A')[1], self.ids[0])
        self.assertEqual(self.resolver.resolve('LEGACY')[1], self.ids[3])
        for name in ('SHARED', 'CONFLICT'):
            self.assertEqual(self.resolver.resolve(name)[3], 'ambiguous')
        self.assertEqual(self.resolver.resolve('38231')[3], 'invalid')
        self.assertEqual(self.resolver.resolve('NO_MATCH')[3], 'unresolved')
        report = self.resolver.resolve_genes([' A ||| OLD_A ', self.ids[0] + '.2', 'B'])
        available, _, summary = assess_coverage(report, self.ids)
        self.assertEqual(available, self.ids[:2])
        self.assertEqual(summary['total'], 2)

    def write_set(self, genes):
        with (self.root / 'fixture.csv').open('w') as f:
            csv.writer(f).writerows([[gene] for gene in genes])

    def test_real_scoring_and_missing_boundaries(self):
        rng = np.random.default_rng(42)
        matrix = rng.random((8, 100))
        var = self.ids[:3] + [f'ENSG{i:011d}' for i in range(100, 197)]
        data = ad.AnnData(matrix, var=pd.DataFrame(index=var))
        data.layers['log1p'] = matrix.copy()
        with patch.object(Eval, 'GENE_SET_ROOT_DIR', self.root), patch.object(Eval, 'get_resolver', return_value=self.resolver):
            self.write_set(['A', 'OLD_A', 'B', 'C', 'UNKNOWN'])
            ids, report = Eval.read_gene_set('fixture', return_report=True)
            self.assertEqual(ids, self.ids[:3])
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                Eval.sEval(data, score_name='fixture')
            self.assertTrue(any('1 unresolved' in str(w.message) for w in caught))
            sc.tl.score_genes(data, gene_list=self.ids[:3], score_name='expected', layer='log1p', use_raw=False, ctrl_as_ref=False)
            np.testing.assert_allclose(data.obs.fixture, data.obs.expected)
            self.write_set(['A', 'B', 'C', 'UNKNOWN', 'ANOTHER'])
            with self.assertRaisesRegex(ValueError, 'more than 25%'):
                Eval.sEval(data, score_name='fixture')
            self.write_set(['UNKNOWN'])
            with self.assertRaisesRegex(ValueError, 'no available'):
                Eval.sEval(data, score_name='fixture')
            self.write_set([])
            with self.assertRaisesRegex(ValueError, 'empty'):
                Eval.sEval(data, score_name='fixture')
            self.write_set(['A', 'B', 'C', 'LEGACY'])
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                Eval.sEval(data, score_name='fixture')
            with self.assertRaisesRegex(ValueError, 'Unknown gene set'):
                Eval.read_gene_set('absent')

    def test_real_old_symbols(self):
        names = ['FAM129B', 'CYR61', 'CTGF', 'LHFP', 'WISP1', 'C10orf10', 'GPR124', 'PTRF', 'FAM198B', 'GUCY1B3', 'PPAP2A', 'C15orf52', 'FAM101B', 'PTPLAD2', 'MRVI1', 'KIAA0430', 'KIRREL', 'GUCY1A3', 'MIR143HG', 'PPAP2B', 'C16orf62', 'ST5', 'SPG20', 'FAM196B', 'PCDHGC5', 'PCDHGA8', 'PCDHGA11', 'PCDHGC3']
        report = get_resolver().resolve_genes(names)
        self.assertTrue(report.status.eq('mapped').all(), report.to_string())
        self.assertEqual(get_resolver().resolve('WISP1-OT1')[3], 'ambiguous')

    def test_expansion_matches_original_and_is_idempotent(self):
        tables = Path(__file__).resolve().parents[2] / 'tables'
        name = '2026_Nat_Canellas_TME-HR.csv'
        with (tables / 'geneset_original' / name).open() as f:
            original = [row[0] for row in csv.reader(f)]
        expanded = list(dict.fromkeys(x.strip() for row in original for x in row.split('|||')))
        with (tables / 'geneset' / name).open() as f:
            actual = [row[0] for row in csv.reader(f)]
        self.assertEqual(actual, expanded)
        self.assertEqual(len(actual), 2447)
        self.assertEqual(actual, list(dict.fromkeys(x.strip() for row in actual for x in row.split('|||'))))
        self.assertTrue(all(actual))


if __name__ == '__main__':
    unittest.main()
