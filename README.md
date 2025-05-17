# Unsupervised Autoencoder + SOM Clustering on PBMC3k

This project demonstrates an end-to-end pipeline for unsupervised clustering of single-cell RNA-seq data using deep learning techniques. The pipeline incorporates a deep autoencoder for dimensionality reduction and a Self-Organizing Map (SOM) for clustering. Additionally, UMAP visualization is used to assess the quality of the clustering, and the silhouette score is computed for quantitative evaluation.

## Overview

- **Autoencoder for Dimensionality Reduction:**  
  The autoencoder compresses the high-dimensional gene expression data from the PBMC3k dataset into a lower-dimensional latent space. This step reduces noise and preserves the most informative features of the data.

- **Self-Organizing Maps (SOM) for Clustering:**  
  Using the latent representation from the autoencoder, SOM is employed to cluster cells into groups. SOMs effectively highlight patterns and similar structures in the data without supervision.

- **Visualization with UMAP:**  
  The UMAP algorithm is applied to project the latent space to two dimensions, enabling visual interpretation of the clusters.

- **Clustering Evaluation:**  
  The silhouette score is used to quantitatively assess the clustering quality.

## Features

- **Dimensionality Reduction:**  
  Uses a deep autoencoder to compress input features.

- **Unsupervised SOM Clustering:**  
  Applies a Self-Organizing Map to cluster cells based on their encoded features.

- **UMAP Visualization:**  
  Provides a 2D visual representation of clusters for qualitative evaluation.

- **Silhouette Score Evaluation:**  
  Measures clustering performance with silhouette metrics.

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