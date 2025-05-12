import scanpy as sc
from sklearn.metrics import silhouette_score

def visualize_umap(adata, cluster_key):
    sc.pp.neighbors(adata, use_rep='X')
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=[cluster_key])

def evaluate_clustering(adata, encoded_X, cluster_labels):
    adata.obs['SOM_cluster'] = cluster_labels.astype(str)
    silhouette = silhouette_score(encoded_X, cluster_labels)
    print(f"Silhouette Score: {silhouette:.4f}")
    return silhouette
