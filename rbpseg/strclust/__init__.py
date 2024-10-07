# rbpseg/strclust/__init__.py

import sys
from .strclust import (
    prepare_diagonal_matrix,
    calculate_cluster_metrics,
    calculate_cluster_means_final,
    calculate_cluster_means,
    calculate_sm,
    exponential_decay,
    exponential_decay_derivative,
    exponential_decay_double_derivative,
    fit_psm,
    calculate_optimal_epsilon,
    calculate_optimal_classes,
    calculate_random_class_mean,
    test_distribution_mean_tm,
    exclude_bad_classes,
    relabel_bad_clusters,
    save_plots,

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
