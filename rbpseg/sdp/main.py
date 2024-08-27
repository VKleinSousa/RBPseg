import argparse
import sys
import logging

from .sdp_calculation import sdp_matrix_calculation
from .structure_processing import load_structure, preprocess_structure
from .clustering import perform_clustering
from .results_processing import save_results, visualize_results

logger = logging.getLogger(__name__)

def main(pdb_file, pair_distance_constant, min_k, max_k, clustering_method, min_cluster_size, symmetry, save_overlap_domains, save_pdb_domains, random_cut, optimize_umap, only_one_chain, normalize_plddt):
    logger.info("Starting main process")
    try:
        structure = load_structure(pdb_file, only_one_chain)
        structure = preprocess_structure(structure, normalize_plddt)
        sdp_matrix, total_plddt, residues = sdp_matrix_calculation(structure, pair_distance_constant)
        clusters, umap_result = perform_clustering(sdp_matrix, clustering_method, min_k, max_k, min_cluster_size, optimize_umap)
        save_results(structure, clusters, residues, save_overlap_domains, save_pdb_domains, pdb_file, symmetry, min_cluster_size)
        visualize_results(sdp_matrix, umap_result, clusters, pdb_file)
        logger.info("Main process completed successfully")
    except Exception as e:
        logger.critical(f"Critical error in main process: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    _defaults = {
        'clustering_method': 'hdbscan',
        'min_k': 3,
        'max_k': 6,
        'min_domain_size': 180,
        'number_of_chains': 3,
        'pair_distance_constant': 1,
        'save_overlap_domains': False,
        'save_pdb_domains': True,
        'optimize_umap': False,
        'only_single_chain': False,
        'normalize_plddt': True
    }

    parser = argparse.ArgumentParser(description='Segment PDB files based on the sDp analysis')
    parser.add_argument("-p", "--pdb", required=True, help="PDB file")
    parser.add_argument("-c", "--clustering_method", type=str, default=_defaults['clustering_method'], help=f'Clustering Method. Options: kmeans, hdbscan. Default: {_defaults["clustering_method"]}')
    parser.add_argument("-k", "--max_k", type=int, default=_defaults['max_k'], help=f'Maximum number of possible kmean clusters. Default: {_defaults["max_k"]}')
    parser.add_argument("-mk", "--min_k", type=int, default=_defaults['min_k'], help=f'Minimum number of possible kmean clusters. Default: {_defaults["min_k"]}')
    parser.add_argument("-s", "--min_domain_size", type=int, default=_defaults['min_domain_size'], help=f'Minimal possible domain size. Default: {_defaults["min_domain_size"]}')
    parser.add_argument("-pd", "--pair_distance_constant", type=float, default=_defaults['pair_distance_constant'], help=f'The pair distance constant will influence how the sDp is calculated. Default: {_defaults["pair_distance_constant"]}')
    parser.add_argument("-n", "--number_of_chains", type=int, default=_defaults['number_of_chains'], help=f'Tail fibre symmetry. Default: {_defaults["number_of_chains"]}')
    parser.add_argument("-so", "--save_overlap_domains", action='store_true', help=f'Include overlap domains per fasta file.')
    parser.add_argument("-sv", "--save_pdb_domains", action='store_true', help=f'Save the pdb file of the domains. Default: {_defaults["save_pdb_domains"]}')
    parser.add_argument("-u", "--optimize_umap", action='store_true', help=f'Self-optimize UMAP parameters.')
    parser.add_argument("-o", "--only_single_chain", action='store_true', help=f'Model has a single chain.')
    parser.add_argument("-np", "--normalize_plddt", action='store_true', help=f'Normalize pLDDT. Default: {_defaults["normalize_plddt"]}')

    args = parser.parse_args()

    main(
        pdb_file=args.pdb,
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

