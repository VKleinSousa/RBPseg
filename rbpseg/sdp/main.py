import argparse
import sys
import logging

from .sdp_calculation import sdp_matrix_calculation
from .structure_processing import load_structure, preprocess_structure
from .clustering import perform_clustering
from .results_processing import save_results, visualize_results

logger = logging.getLogger(__name__)

def main(pdb_file,experimental_structure=False, pair_distance_constant=1, min_k=3, max_k=6, clustering_method='spectral', min_cluster_size=120, symmetry=3, save_overlap_domains=False, save_pdb_domains=True, random_cut=False, optimize_umap=False, only_one_chain=False, normalize_plddt=True):
    logger.info("Starting main process")
    try:
        structure = load_structure(pdb_file, True)
        structure = preprocess_structure(structure, experimental_structure, normalize_plddt)
        sdp_matrix, total_plddt, residues = sdp_matrix_calculation(structure, pair_distance_constant)
        clusters, umap_result = perform_clustering(sdp_matrix, clustering_method, min_k, max_k, min_cluster_size, optimize_umap)
        
        structure_full = load_structure(pdb_file, only_one_chain)
        save_results(structure_full, clusters, residues, save_overlap_domains, save_pdb_domains, pdb_file, symmetry, min_cluster_size)
        visualize_results(sdp_matrix, umap_result, clusters, pdb_file)
        logger.info("Main process completed successfully")
    except Exception as e:
        logger.critical(f"Critical error in main process: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    _defaults = {
        'clustering_method': 'spectral',
        'min_k': 3,
        'max_k': 6,
        'min_domain_size': 120,
        'number_of_chains': 3,
        'pair_distance_constant': 1,
        'save_overlap_domains': False,
        'save_pdb_domains': True,
        'optimize_umap': False,
        'only_single_chain': False,
        'normalize_plddt': True,
        'experimental_structure': False,
    }

    parser = argparse.ArgumentParser(description='Segment PDB files based on the sDp analysis')
    parser.add_argument("-p", "--pdb", required=True, help="PDB file")
    parser.add_argument("-c", "--clustering_method", type=str, default=_defaults['clustering_method'], help=f'Clustering Method. Options: kmeans, hdbscan. Default: {_defaults["clustering_method"]}')
    parser.add_argument("-k", "--max_k", type=int, default=_defaults['max_k'], help=f'Maximum number of possible kmean clusters. Default: {_defaults["max_k"]}')
    parser.add_argument("-mk", "--min_k", type=int, default=_defaults['min_k'], help=f'Minimum number of possible kmean clusters. Default: {_defaults["min_k"]}')
    parser.add_argument("-s", "--min_domain_size", type=int, default=_defaults['min_domain_size'], help=f'Minimal possible domain size. Default: {_defaults["min_domain_size"]}')
    parser.add_argument("-pd", "--pair_distance_constant", type=float, default=_defaults['pair_distance_constant'], help=f'The pair distance constant will influence how the sDp is calculated. Default: {_defaults["pair_distance_constant"]}')
    parser.add_argument("-n", "--number_of_chains", type=int, default=_defaults['number_of_chains'], help=f'Tail fibre symmetry. Default: {_defaults["number_of_chains"]}')
    parser.add_argument("-e", "--experimental_structure", action='store_true', help=f'Add this flag if the input is not an Alphafold/ESMfold prediction with bfactors set as the pLDDT value.')
    parser.add_argument("-so", "--save_overlap_domains", action='store_true', help=f'Add flag to include overlap domains per fasta file.')
    parser.add_argument("-sv", "--save_pdb_domains", action='store_true', help=f'Add flag to save the pdb file of the domains.')
    parser.add_argument("-u", "--optimize_umap", action='store_true', help=f'Add flat to self-optimize UMAP parameters by a Monte Carlo Approach.')
    parser.add_argument("-o", "--only_single_chain", action='store_true', help=f'Add flag of model has a single chain.')
    parser.add_argument("-np", "--normalize_plddt", action='store_true', help=f'Add flag to normalize pLDDT')

    args = parser.parse_args()

    main(
        pdb_file=args.pdb,
        experimental_structure=args.experimental_structure,
        pair_distance_constant=args.pair_distance_constant,
        min_k=args.min_k,
        max_k=args.max_k,
        clustering_method=args.clustering_method,
        min_cluster_size=args.min_domain_size,
        symmetry=args.number_of_chains,
        save_overlap_domains=args.save_overlap_domains,
        save_pdb_domains=args.save_pdb_domains,
        random_cut=False,
        optimize_umap=args.optimize_umap,
        only_one_chain=args.only_single_chain,
        normalize_plddt=args.normalize_plddt
    )
