from data.load_data import load_and_preprocess_pbmc3k
from models.autoencoder import build_autoencoder
from clustering.som_clustering import run_som_clustering
from utils.visualize import visualize_umap, evaluate_clustering

def main():
    adata, X = load_and_preprocess_pbmc3k()
    autoencoder, encoder = build_autoencoder(X.shape[1])

    print("Training Autoencoder...")
    autoencoder.fit(X, X, epochs=50, batch_size=32, shuffle=True, verbose=1)

    print("Encoding data...")
    encoded_X = encoder.predict(X)

    print("Clustering with SOM...")
    cluster_labels = run_som_clustering(encoded_X)

    print("Evaluating clustering...")
    adata.obsm['X'] = encoded_X
    silhouette = evaluate_clustering(adata, encoded_X, cluster_labels)
    print(f"Silhouette Score: {silhouette:.4f}")
    adata.obs['SOM_cluster'] = cluster_labels.astype(str)
    print("Saving results...")
    print("Visualizing...")
    visualize_umap(adata, 'SOM_cluster')

    

if __name__ == "__main__":
    main()
