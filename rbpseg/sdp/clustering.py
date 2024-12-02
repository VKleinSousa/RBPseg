from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from joblib import Parallel, delayed
from .umap_processing import sdp_umap
import hdbscan
from scipy.spatial.distance import cdist
import numpy as np
import logging
import sys  
import pandas as pd
from sklearn.cluster import SpectralClustering
from rbpseg.strclust import *
from scipy.ndimage import uniform_filter1d
logger = logging.getLogger(__name__)

def kmeans_clustering(min_k, max_k, umap_result):
    try:
        silhouette_scores = Parallel(n_jobs=-1)(delayed(silhouette_score)(
            umap_result, KMeans(n_clusters=k, random_state=42, n_init='auto').fit_predict(umap_result)
        ) for k in range(min_k, max_k))

        optimal_k = silhouette_scores.index(max(silhouette_scores)) + min_k
        kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init='auto')
        clusters = kmeans.fit_predict(umap_result)
        logger.info(f"KMeans clustering completed with optimal_k={optimal_k}")
        return clusters
    except Exception as e:
        logger.error(f"Error during KMeans clustering: {e}", exc_info=True)
        raise

def hdbscan_clustering(min_cluster_size, umap_result):
    try:
        hdbscan_labels = hdbscan.HDBSCAN(min_samples=None, min_cluster_size=min_cluster_size, leaf_size=40).fit_predict(umap_result)
        valid_clusters = np.unique(hdbscan_labels[hdbscan_labels >= 0])

        if len(valid_clusters) < 2:
            hdbscan_labels[:] = 1
            logger.warning("HDBSCAN found fewer than 2 unique clusters. All points assigned to a single cluster.")
        else:
            for i, label in enumerate(hdbscan_labels):
                if label == -1:
                    nearest_cluster = min(valid_clusters, key=lambda x: cdist([umap_result[i]], umap_result[hdbscan_labels == x]).min())
                    hdbscan_labels[i] = nearest_cluster

        logger.info("HDBSCAN clustering completed")
        return hdbscan_labels
    except Exception as e:
        logger.error(f"Error during HDBSCAN clustering: {e}", exc_info=True)
        raise


def spectral_clustering_sdp(sdp_matrix, min_k=3, max_k=20, smoothing_window=3, gradient_threshold=0.005):
    """
    Perform spectral clustering with silhouette score optimization.
    
    Parameters:
    sdp_matrix (numpy array or pandas DataFrame): The similarity matrix.
    min_k (int): Minimum number of clusters.
    max_k (int): Maximum number of clusters.
    smoothing_window (int): Window size for moving average smoothing.
    gradient_threshold (float): Threshold for silhouette score gradient change.
    
    Returns:
    clusters (array): The predicted cluster labels.
    optimal_k (int): The optimal number of clusters based on silhouette score stabilization.
    """
    
    # If the sdp_matrix is a numpy array, convert it to a DataFrame
    if isinstance(sdp_matrix, np.ndarray):
        sdp_matrix = pd.DataFrame(sdp_matrix)

    silhouette_scores = []
    
    # Loop over different values of k to compute silhouette scores
    for k in range(min_k, max_k + 1):
        clustering_labels = SpectralClustering(n_clusters=k, affinity='precomputed').fit_predict(sdp_matrix)
        score = silhouette_score(sdp_matrix, clustering_labels)
        silhouette_scores.append(score)
        print(f"Silhouette score for k={k}: {score:.4f}")  # Print the silhouette score for each k

    # Convert silhouette scores to a numpy array for further processing
    silhouette_scores = np.array(silhouette_scores)
    
    # Step 2: Apply Moving Average Smoothing
    smoothed_scores = uniform_filter1d(silhouette_scores, size=smoothing_window)
    
    # Step 3: Compute the first derivative (gradient)
    gradients = np.diff(smoothed_scores)  # Rate of change of silhouette scores
    print(f"Gradients: {gradients}")
    
    # Step 4: Find the stabilization point using the gradient
    optimal_k = min_k
    for i in range(1, len(gradients)):
        # If the change in gradient is smaller than the threshold, stop
        if abs(gradients[i]) < gradient_threshold:
            optimal_k = min_k + i  # Because gradients are one element shorter than k range
            print(f"Stabilization detected at k={optimal_k}, gradient={gradients[i]:.4f}")
            break
    
    # If no stabilization found, take the k with max silhouette score
    if optimal_k == min_k:
        optimal_k = np.argmax(smoothed_scores) + min_k
        print(f"No stabilization found. Taking optimal_k based on max silhouette score: {optimal_k}")

    # Step 5: Perform spectral clustering with the optimal number of clusters
    spectral_clustering = SpectralClustering(n_clusters=optimal_k, affinity='precomputed')
    clusters = spectral_clustering.fit_predict(sdp_matrix)
    print(clusters)
    # Log the optimal number of clusters
    logger.info(f"Spectral clustering completed with optimal_k={optimal_k}")

    return clusters

    
def perform_clustering(sdp_matrix, clustering_method, min_k, max_k, min_cluster_size, optimize_umap):
    """
    Performs clustering on the given sdp_matrix using the specified clustering method.

    :param sdp_matrix: The SDP matrix to be clustered.
    :param clustering_method: Clustering method to be used ('kmeans', 'hdbscan', or 'spectral').
    :param min_k: Minimum number of clusters for KMeans.
    :param max_k: Maximum number of clusters for KMeans.
    :param min_cluster_size: Minimum cluster size for HDBSCAN.
    :param optimize_umap: Whether to optimize UMAP parameters or not.
    :return: A tuple containing the clusters and the UMAP result.
    """

        
        
    try:
        optimal_limit = int(len(sdp_matrix)//min_cluster_size)
        if max_k > optimal_limit:
            if optimal_limit < min_k:
                max_k = min_k
                print(f"max_k changed to {max_k}")
            else:
                max_k = optimal_limit + 2
                print(f"max_k changed to {max_k}")
        umap_result = sdp_umap(sdp_matrix, optimize_umap, min_k, max_k)
        if clustering_method.lower() == 'spectral':
            clusters = spectral_clustering_sdp(sdp_matrix, min_k, max_k)
        else:
            if clustering_method.lower() == 'hdbscan':
                clusters = hdbscan_clustering(min_cluster_size, umap_result)
            elif clustering_method.lower() == 'kmeans':
                clusters = kmeans_clustering(min_k, max_k, umap_result)
            else:
                logger.error("Invalid clustering method. Use 'hdbscan', 'kmeans', or 'spectral'.")
                sys.exit(1)
        
        return clusters, umap_result
    except Exception as e:
        logger.error(f"Error during clustering: {e}", exc_info=True)
        raise
