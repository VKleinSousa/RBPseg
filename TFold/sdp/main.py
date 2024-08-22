import argparse
import sys
import logging

from sdp_calculation import sdp_matrix_calculation
from structure_processing import load_structure, preprocess_structure
from clustering import perform_clustering
from results_processing import save_results, visualize_results

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
    _defaults = {
        'cluster': 'hdbscan',
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
    parser.add_argument("-p", "--pdb", help="PDB file")
    parser.add_argument("-c", "--clustering_method", type=str, default=_defaults['cluster'], help=f'Clustering Method. Options: kmeans, HDBSCAN. Default: {_defaults["cluster"]}')
    parser.add_argument("-k", "--max_k", type=int, default=_defaults['max_k'], help=f'Maximum number of possible kmean clusters. Used only when -c kmeans. Smaller max_k will usually result in larger fasta files. Default: {_defaults["max_k"]}')
    parser.add_argument("-mk", "--min_k", type=int, default=_defaults['min_k'], help=f'Minimum number of possible kmean clusters. Used only when -c kmeans. Smaller max_k will usually result in larger fasta files. Default: {_defaults["min_k"]}')
    parser.add_argument("-s", "--min_domain_size", type=int, default=_defaults['min_domain_size'], help=f'Minimal possible domain size. Default: {_defaults["min_domain_size"]}')
    parser.add_argument("-pd", "--pair_distance_constant", type=float, default=_defaults['pair_distance_constant'], help=f'The pair_distance_constant will influence how the sDp is calculated: >1 values will increase the weight given to the plddt values of distant residues. It should always be >0 . Default: {_defaults["pair_distance_constant"]}')
    parser.add_argument("-n", "--number_of_chains", type=int, default=_defaults['number_of_chains'], help=f'Tail fibre symmetry. Default: {_defaults["number_of_chains"]}' )
    parser.add_argument("-so", "--save_overlap_domains", type=bool, default=_defaults['save_overlap_domains'], help=f'If including only one (False) or two domains per fasta file (True). Default: {_defaults["save_overlap_domains"]}' )
    parser.add_argument("-sv", "--save_pdb_domains", type=bool, default=_defaults['save_pdb_domains'], help=f'Whether to save the pdb file of the domains. Default: {_defaults["save_pdb_domains"]}' )
    parser.add_argument("-u", "--optimize_umap", type=bool, default=_defaults['optimize_umap'], help=f'Whether to self optimize umap parameters. Default: {_defaults["optimize_umap"]}' )
    parser.add_argument("-o", "--only_single_chain", type=bool, default=_defaults['only_single_chain'], help=f'Whether the model has a single chain. Default: {_defaults["only_single_chain"]}' )
    parser.add_argument("-np", "--normalize_plddt", type=bool, default=_defaults['normalize_plddt'], help=f'Whether to normalize plddt. Default: {_defaults["normalize_plddt"]}' )

    args = parser.parse_args()

    main(args.pdb, args.pair_distance_constant, args.min_k, args.max_k, args.clustering_method, args.min_domain_size, args.number_of_chains, args.save_overlap_domains, args.save_pdb_domains, False, args.optimize_umap, args.only_single_chain, args.normalize_plddt)
