import argparse
import sys
import logging
import textwrap
from .sdp_calculation import sdp_matrix_calculation
from .structure_processing import load_structure, preprocess_structure
from .clustering import perform_clustering
from .results_processing import save_results, visualize_results

logger = logging.getLogger(__name__)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
            ___________________________________________
                        RBPseg-sdp v0.1.3
            ___________________________________________
            Segment PDB files based on the sDp analysis. 
            '''), 
    
    epilog=textwrap.dedent('''\
            _______________________________
            Developed by Victor Klein-Sousa.
            
            If you used this script, please consider citing us:

            Towards a complete phage tail fiber structure atlas. Victor Klein-Sousa, Aritz Roa-Eguiara, Claudia Sybille Kielkopf, Nicholas Sofos, Nicholas M. I. Taylor bioRxiv 2024.10.28.620165; doi: https://doi.org/10.1101/2024.10.28.620165
            ______________________________
            
            Please don't hesitate to contact us if you have any problems.
            
            '''))
    
    _defaults = {
        'clustering_method': 'HDBSCAN',
        'min_k': 3,
        'max_k': 20,
        'min_domain_size': 120,
        'min_ov_size': 50,
        'number_of_chains': 3,
        'pair_distance_constant': 1,
        'save_overlap_domains': False,
        'save_pdb_domains': True,
        'optimize_umap': False,
        'only_single_chain': False,
        'normalize_plddt': True,
        'experimental_structure': False,
    }

    # Define arguments
    parser.add_argument("-p", "--pdb", required=True, help="PDB file")
    parser.add_argument("-c", "--clustering_method", type=str, default=_defaults['clustering_method'], help=f'Clustering Method. Options: kmeans, hdbscan. Default: {_defaults["clustering_method"]}')
    parser.add_argument("-k", "--max_k", type=int, default=_defaults['max_k'], help=f'Maximum number of possible kmean clusters. Used only -c is spectral or kmeans. Optimal value will depend on the size of the fiber. Default: {_defaults["max_k"]}')
    parser.add_argument("-mk", "--min_k", type=int, default=_defaults['min_k'], help=f'Minimum number of possible kmean clusters. Default: {_defaults["min_k"]}')
    parser.add_argument("-s", "--min_domain_size", type=int, default=_defaults['min_domain_size'], help=f'Minimal possible domain size. Default: {_defaults["min_domain_size"]}')
    parser.add_argument("-ovs", "--min_ov_size", type=int, default=_defaults['min_ov_size'], help=f'Minimal possible overhang size. Default: {_defaults["min_ov_size"]}')
    parser.add_argument("-pd", "--pair_distance_constant", type=float, default=_defaults['pair_distance_constant'], help=f'The pair distance constant will influence how the sDp is calculated. Default: {_defaults["pair_distance_constant"]}')
    parser.add_argument("-n", "--number_of_chains", type=int, default=_defaults['number_of_chains'], help=f'Tail fibre symmetry. Default: {_defaults["number_of_chains"]}')
    parser.add_argument("-e", "--experimental_structure", action='store_true', help=f'Add this flag if the input is not an Alphafold/ESMfold prediction with bfactors set as the pLDDT value.')
    parser.add_argument("-so", "--save_overlap_domains", action='store_true', help=f'Add flag to include overlap domains per fasta file.')
    parser.add_argument("-sv", "--save_pdb_domains", action='store_true', help=f'Add flag to save the pdb file of the domains.')
    parser.add_argument("-u", "--optimize_umap", action='store_true', help=f'Add flag to self-optimize UMAP parameters.')
    parser.add_argument("-o", "--only_single_chain", action='store_true', help=f'Add flag if the model has a single chain.')
    parser.add_argument("-np", "--normalize_plddt", action='store_true', help=f'Add flag to normalize pLDDT.')

    args = parser.parse_args()

    logger.info("Starting main process")
    try:
        structure = load_structure(args.pdb, args.only_single_chain)
        structure = preprocess_structure(structure, args.experimental_structure, args.normalize_plddt)
        sdp_matrix, total_plddt, residues = sdp_matrix_calculation(structure, args.pair_distance_constant)
        clusters, umap_result = perform_clustering(sdp_matrix, args.clustering_method, args.min_k, args.max_k, args.min_domain_size, args.optimize_umap)
        
        structure_full = load_structure(args.pdb, args.only_single_chain)
        save_results(structure_full, clusters, residues, args.save_overlap_domains, args.save_pdb_domains, args.pdb, args.number_of_chains, args.min_domain_size, args.min_ov_size)
        visualize_results(sdp_matrix, umap_result, clusters, args.pdb)
        logger.info("Main process completed successfully")
    except Exception as e:
        logger.critical(f"Critical error in main process: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    main()
