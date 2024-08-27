# rbpseg/sdp/__init__.py

import sys
from .sdp_calculation import (
    sdp,
    sdp_t,
    sdp_d,
    sdp_matrix_calculation
)

from .umap_processing import (
    sdp_umap,
    umap_optimizer
)


from .structure_processing import (
    load_structure,
    preprocess_structure,
    extract_chain,
    normalize_bfactor
)

from .clustering import (
    kmeans_clustering,
    hdbscan_clustering,
    perform_clustering
)

from .results_processing import (
    save_structure,
    save_results,
    save_overlap_domains_to_fasta,
    save_pdb_clusters,
    visualize_results,
    three_to_one,
    split_cluster
)

import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),  # Log to console
        logging.FileHandler("protein_analysis.log")  # Log to a file
    ]
)

logger = logging.getLogger(__name__)
