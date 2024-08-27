from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from joblib import Parallel, delayed
from .umap_processing import sdp_umap
import hdbscan
from scipy.spatial.distance import cdist
import numpy as np
import logging


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
def perform_clustering(sdp_matrix, clustering_method, min_k, max_k, min_cluster_size, optimize_umap):
    """
    Performs clustering on the given sdp_matrix using the specified clustering method.

    :param sdp_matrix: The SDP matrix to be clustered.
    :param clustering_method: Clustering method to be used ('kmeans' or 'hdbscan').
    :param min_k: Minimum number of clusters for KMeans.
    :param max_k: Maximum number of clusters for KMeans.
    :param min_cluster_size: Minimum cluster size for HDBSCAN.
    :param optimize_umap: Whether to optimize UMAP parameters or not.
    :return: A tuple containing the clusters and the UMAP result.
    """
    umap_result = sdp_umap(sdp_matrix, optimize_umap, min_k, max_k)
    
    if clustering_method.lower() == 'hdbscan':
        clusters = hdbscan_clustering(min_cluster_size, umap_result)
    elif clustering_method.lower() == 'kmeans':
        clusters = kmeans_clustering(min_k, max_k, umap_result)
    else:
        logger.error("Invalid clustering method. Use 'hdbscan' or 'kmeans'.")
        sys.exit(1)
        
    return clusters, umap_result
