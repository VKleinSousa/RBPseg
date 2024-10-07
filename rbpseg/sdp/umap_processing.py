from umap import UMAP
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from joblib import Parallel, delayed
import numpy as np
import logging

logger = logging.getLogger(__name__)

def umap_optimizer(min_k, max_k, sdp_matrix):
    try:
        parameters = []
        avg_sil = []

        def umap_silhouette(n_neighbors, min_dist, negative_sample_rate, repulsion_strength):
            umap = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, negative_sample_rate=negative_sample_rate,
                        repulsion_strength=repulsion_strength, random_state=42)
            umap_result = umap.fit_transform(sdp_matrix)
            silhouette_scores = [silhouette_score(umap_result, KMeans(n_clusters=k, random_state=42, n_init='auto').fit_predict(umap_result))
                                 for k in range(min_k, max_k)]
            return np.mean(silhouette_scores), [n_neighbors, min_dist, negative_sample_rate, repulsion_strength]

        results = Parallel(n_jobs=-1)(delayed(umap_silhouette)(
            np.random.randint(2, 50), np.random.random(), np.random.randint(1, 50), np.random.randint(1, 50)
        ) for _ in range(50))

        avg_sil, parameters = zip(*results)
        
        if np.mean(avg_sil) > max(avg_sil):
            logger.info('UMAP optimizer could not converge to a set of parameters better than the standard')
            return [0, 0, 0, 0]
        else:
            optimal_params = parameters[np.argmax(avg_sil)]
            logger.info(f'Optimal UMAP parameters: {optimal_params}')
            return optimal_params

    except Exception as e:
        logger.error(f"Error during UMAP optimization: {e}", exc_info=True)
        raise

def sdp_umap(sdp_matrix, optimize_umap, min_k, max_k):
    try:
        if optimize_umap:
            n_neighbors, min_dist, negative_sample_rate, repulsion_strength = umap_optimizer(min_k, max_k, sdp_matrix)
            if n_neighbors == 0:
                umap = UMAP(n_components=2, random_state=42)
            else:
                umap = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, negative_sample_rate=negative_sample_rate,
                            repulsion_strength=repulsion_strength, random_state=42)
        else:
            umap = UMAP(n_components=2, random_state=42)

        umap_result = umap.fit_transform(sdp_matrix)
        logger.info("UMAP projection completed")
        return umap_result
    except Exception as e:
        logger.error(f"Error during UMAP projection: {e}", exc_info=True)
        raise
