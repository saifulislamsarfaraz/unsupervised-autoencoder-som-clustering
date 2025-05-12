import scanpy as sc

def load_and_preprocess_pbmc3k():
    adata = sc.datasets.pbmc3k()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    return adata, X
