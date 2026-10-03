"""Gene-set scoring with offline historical-symbol resolution."""
from pathlib import Path
import warnings

import pandas as pd

from .GeneSetMapping import assess_coverage, get_resolver

GENE_SET_ROOT_DIR = Path(__file__).resolve().parents[2] / 'tables/geneset'


def get_names():
    return sorted(path.stem for path in GENE_SET_ROOT_DIR.glob('*.csv'))


def read_gene_set(score_name=None, *, return_report=False):
    """Return unique ENSG IDs, optionally with a per-input mapping report."""
    if score_name not in get_names():
        raise ValueError(f'Unknown gene set: {score_name!r}')
    path = GENE_SET_ROOT_DIR / f'{score_name}.csv'
    try:
        genes = pd.read_csv(path, header=None, dtype=str, keep_default_na=False).iloc[:, 0].tolist()
    except pd.errors.EmptyDataError:
        genes = []
    report = get_resolver().resolve_genes(genes)
    ids = list(dict.fromkeys(report.loc[report.status.eq('mapped'), 'gene_id']))
    if return_report:
        return ids, report
    unresolved = report.loc[~report.status.eq('mapped'), 'name'].unique()
    if len(unresolved):
        warnings.warn(
            f'{score_name}: {len(unresolved)} names unresolved or ambiguous: '
            + ', '.join(unresolved) + '. Use return_report=True for details.',
            UserWarning, stacklevel=2,
        )
    return ids


def sEval(adata, flavour='scanpy', layer='log1p', score_name=None):
    import scanpy as sc
    if not score_name:
        raise ValueError('score_name is not provided')
    _, report = read_gene_set(score_name, return_report=True)
    available, report, summary = assess_coverage(report, adata.var_names)
    if not summary['total']:
        raise ValueError(f'Gene set {score_name} is empty')
    if not available:
        raise ValueError(f'Gene set {score_name} has no available genes in the dataset')
    message = (
        f"{score_name}: {summary['unresolved_unique_names']} unresolved names "
        f"({summary['ambiguous_unique_names']} ambiguous), "
        f"{summary['missing_dataset_ids']} resolved ENSG IDs absent from dataset; "
        f"{summary['available_ids']}/{summary['total']} available."
    )
    if summary['missing_fraction'] > 0.25:
        raise ValueError('more than 25% of the genes are unavailable. ' + message)
    if summary['missing_fraction']:
        warnings.warn(message, UserWarning, stacklevel=2)
    ### checked
    
    match flavour:
        case 'scanpy':
            sc.tl.score_genes(
            adata, gene_list=available, score_name=score_name,
            layer=layer, use_raw=False, ctrl_as_ref=False,
            )
        case 'AUCell':
            raise NotImplementedError('AUCell scoring is not implemented yet')
            # TODO: implement AUCell scoring
        case _:
            raise ValueError(f'Unknown method: {flavour!r}')
    
