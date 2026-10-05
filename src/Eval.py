from pathlib import Path
import warnings
import pandas as pd

from .GTFMapping import resolve_gene_name
# . means import relatively

GENE_SET_ROOT_DIR = Path(__file__).resolve().parents[2] / 'tables/geneset'


def get_names():
    return sorted(path.stem for path in GENE_SET_ROOT_DIR.glob('*.csv'))

def parse_multiple(genes:[str], sep = '|'):
    """
    pass a list of gene names, in which there could be multiple name in one str
    """
    single_gene_list = []
    for item in genes:
        if item == '':
            raise ValueError('Empty gene name in the list')
        single_gene_list = single_gene_list + [i.strip() for i in item.split(sep)]
    return single_gene_list

def read_gene_set(score_name=None):
    """
    pass a name in dir
    return (1) a list of ENSG id (2) generated dict
    """
    if score_name not in get_names():
        raise ValueError(f'Unknown gene set: {score_name!r}')
    path = GENE_SET_ROOT_DIR / f'{score_name}.csv'
    genes = parse_multiple(pd.read_csv(
        path,header=None,
        keep_default_na=False,).iloc[:, 0].tolist())
    resolved_ENSG = resolve_gene_name(genes, fill_NA=False) # a pd dataframe
    nona_ids = resolved_ENSG['gene_id'].ne('NA')
    # we don't care about the values actually; list(dict) is just to deduplicated since dict keys will be unique
    ids = list(dict.fromkeys(resolved_ENSG.loc[nona_ids, 'gene_id']))
    return ids, resolved_ENSG

def sEval(adata, flavour='scanpy', layer='log1p', score_name=None):
    if not score_name:
        raise ValueError('score_name is not provided')
    ids, resolved_ENSG = read_gene_set(score_name)
    # need to check if ENSG ids are available in adata
    percentage = 100 * len([notfound for notfound in ids if not notfound in adata.var_names]) / len(ids)
    if percentage > 5:
        warnings.warn(f'Unmatched ENSG ids in adata exceed {percentage}%  ', UserWarning, stacklevel=2)

    match flavour:
        case 'AUCell':
            import decoupler as dc
            # generate net(gene set)
            net = pd.DataFrame({'source': score_name, 'target': ids})
            dc.mt.aucell(adata, layer = layer, net = net, tmin=3)
        case 'scanpy':
            import scanpy as sc
            sc.tl.score_genes(
                adata, gene_list=ids, score_name=score_name,
                layer=layer, use_raw=False, ctrl_as_ref=False,
                )
        case _:
            raise ValueError(f'Unknown method: {flavour!r}')
    
