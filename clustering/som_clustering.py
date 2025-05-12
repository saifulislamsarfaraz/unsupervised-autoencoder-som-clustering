from minisom import MiniSom
import numpy as np

def run_som_clustering(encoded_X, grid_size=4, sigma=1.0, learning_rate=0.5):
    som = MiniSom(grid_size, grid_size, encoded_X.shape[1], sigma=sigma, learning_rate=learning_rate)
    som.random_weights_init(encoded_X)
    som.train_random(encoded_X, 100)

    winner_coordinates = np.array([som.winner(x) for x in encoded_X])
    cluster_labels = np.ravel_multi_index(winner_coordinates.T, (grid_size, grid_size))
    return cluster_labels
