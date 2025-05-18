# Unsupervised Autoencoder + SOM Clustering on PBMC3k Single-Cell Data

This project demonstrates an end-to-end pipeline for unsupervised clustering of single-cell RNA-seq data using a deep autoencoder and Self-Organizing Maps (SOM). A combination of advanced dimensionality reduction, clustering, and visualization techniques is employed to extract and explore meaningful structures in the data.

## Overview

- **Data Preprocessing:**  
  The PBMC3k dataset from Scanpy is loaded, filtered, normalized, log-transformed, and reduced to its top 2000 highly variable genes. Exploratory analysis is performed by plotting highly expressed genes and PCA variance ratios.

- **Autoencoder for Dimensionality Reduction:**  
  A deep autoencoder is built using the Keras functional API. Its encoder component compresses the high-dimensional gene expression data into a lower-dimensional latent space suitable for downstream tasks.

- **SOM Clustering with Hyperparameter Tuning:**  
  The latent representations are clustered using a Self-Organizing Map (SOM) on a 4x4 grid. Hyperparameters such as sigma and learning rate are tuned by computing the silhouette score to identify optimal settings. The best SOM parameters are then used to generate final cluster assignments and numeric cluster IDs.

- **Additional Clustering Methods and Evaluation:**  
  - **KMeans Clustering:** Applied to the latent data with silhouette scoring.
  - **PCA + Leiden & Louvain Clustering:** PCA is performed, and both Leiden and Louvain clustering are evaluated using the PCA-reduced representation.
  - **SOM Evaluation:** The silhouette score for the SOM clusters is computed using numerical cluster IDs.

- **Visualization with UMAP:**  
  UMAP is used to generate 2D visualizations of the clustering results. Multiple UMAP plots are generated and saved, displaying SOM, KMeans, Leiden, and Louvain clusterings for comparative analysis.

## Installation

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/saifulislamsarfaraz/unsupervised-autoencoder-som-clustering.git
   cd unsupervised-autoencoder-som-clustering
   ```

2. **Install Dependencies:**  
   Use the following command to install required packages:
   ```bash
   pip install minisom scanpy tensorflow matplotlib seaborn igraph leidenalg louvain
   ```

## Notebook Workflow

The main steps in the workflow are detailed in the Jupyter notebook (`single_cell_Autoencoder_and_SOM.ipynb`):

1. **Data Loading & Preprocessing:**  
   - Load the PBMC3k dataset.
   - Perform filtering, normalization, log-transformation, and selection of highly variable genes.
   - Scale the data with StandardScaler.

2. **Autoencoder Construction and Training:**  
   - Build the autoencoder and encoder models.
   - Train the autoencoder with early stopping and track the training/validation losses.
   - Encode the data into a latent space.

3. **SOM Clustering with Hyperparameter Tuning:**  
   - Iterate over different sigma and learning rate values.
   - Evaluate using silhouette scores.
   - Retrain the SOM with the best parameters and record cluster assignments.

4. **Additional Clustering and Evaluation:**  
   - Perform KMeans, PCA + Leiden, and Louvain clustering.
   - Compute silhouette scores for all methods to assess clustering quality.

5. **Visualization:**  
   - Generate UMAP embeddings on the latent space.
   - Visualize the clustering results with separate UMAP plots for SOM, KMeans, Leiden, and Louvain clusters.

   Example UMAP code snippet:
   ```python
   # Generate UMAP plots
   sc.pp.neighbors(adata, use_rep='X_latent', n_neighbors=15)
   sc.tl.umap(adata)
   sc.pl.umap(adata, color=['som_cluster'], title="SOM Clustering on Autoencoder Latent Space", save='_som_umap.png')
   sc.pl.umap(adata, color=['kmeans'], save='_kmeans_umap.png')
   sc.pl.umap(adata, color=['leiden_pca'], save='_pca_leiden_umap.png')
   sc.pl.umap(adata, color=['louvain'], save='_louvain_umap.png')
   ```

## Running the Notebook

Open the notebook in your favorite Jupyter environment or VS Code with the Python extension. Run the cells in order to reproduce the analysis and generate the clustering and UMAP visualizations.


## Project Structure

- **data/load_data.py:**  
  Handles loading and preprocessing of the PBMC3k dataset.

- **models/autoencoder.py:**  
  Contains the definition and compilation of the autoencoder model.

- **clustering/som_clustering.py:**  
  Runs the SOM clustering on the latent features and converts the SOM output into cluster labels.

- **utils/visualize.py:**  
  Provides functions for UMAP visualization and evaluation of clustering performance.

- **main.py:**  
  Integrates all the modules to execute the training, clustering, visualization, and evaluation pipeline.

- **notebooks/:**  
  Contains Jupyter notebooks for exploratory data analysis and prototyping of the workflow.


## Installation

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/saifulislamsarfaraz/unsupervised-autoencoder-som-clustering.git
   cd unsupervised-autoencoder-som-clustering
   python main.py