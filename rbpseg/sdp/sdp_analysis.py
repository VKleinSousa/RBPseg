import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from sklearn.cluster import KMeans
from umap import UMAP
from sklearn.metrics import silhouette_score
from sklearn import preprocessing
import sys
import csv
import argparse
import hdbscan
from random import randrange, random
import heapq
import matplotlib.pyplot as plt
import logging
from joblib import Parallel, delayed
from scipy.spatial.distance import cdist

# Configure logging
logging.basicConfig(
    level=logging.INFO,  # Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),  # Log to console
    ]
)

logger = logging.getLogger(__name__)

def sdp(plddt1, plddt2, d, pair_distance_constant):
    if d == 0:
        return 1
    else:
        return 1 / (1 + np.exp(-1 * pair_distance_constant * (plddt1 + plddt2) / (d * d)))

def sdp_t(plddt1, plddt2, d, pair_distance_constant):
    if d == 0:
        return 1
    else:
        return 1 / (1 + np.exp(-1 * pair_distance_constant * (plddt1 * plddt2) / (d * d)))

def sdp_d(plddt1, plddt2, d, pair_distance_constant):
    if d == 0:
        return 1
    else:
        return 1 / (1 + np.exp(-1 * pair_distance_constant * (plddt1 + plddt2) / d))

def sdp_matrix_calculation(structure, pair_distance_constant, sdp_mode=0):
    try:
        residues = [residue for model in structure for chain in model for residue in chain]
        num_residues = len(residues)
        sdp_matrix = np.zeros((num_residues, num_residues))

        # Extract coordinates and pLDDT values for all residues
        coords = np.array([residue['CA'].get_coord() for residue in residues])
        plddts = np.array([residue['CA'].get_bfactor() for residue in residues])

        distances = cdist(coords, coords)

        if sdp_mode == 0:
            sdp_matrix = 1 / (1 + np.exp(-1 * pair_distance_constant * (plddts[:, None] + plddts[None, :]) / (distances ** 2)))
        elif sdp_mode == 1:
            sdp_matrix = 1 / (1 + np.exp(-1 * pair_distance_constant * (plddts[:, None] * plddts[None, :]) / (distances ** 2)))
        elif sdp_mode == 2:
            sdp_matrix = 1 / (1 + np.exp(-1 * pair_distance_constant * (plddts[:, None] + plddts[None, :]) / distances))

        np.fill_diagonal(sdp_matrix, 1)
        return sdp_matrix, np.mean(plddts), residues

    except KeyError as e:
        logger.error(f"Error accessing CA atom or B-factor for residues: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during SDP matrix calculation: {e}", exc_info=True)
        raise

def split_cluster(cluster_residues, max_residue_jump=10):
    new_cluster_residues = []
    current_subcluster = []
    for cluster in cluster_residues:
        for i in range(len(cluster)):
            current_residue_id = cluster[i].id[1]
            if current_residue_id == cluster[len(cluster)-1].id[1]:
                current_subcluster.append(cluster[i])
                break
            next_residue_id = cluster[i + 1].id[1]

            if next_residue_id - current_residue_id > 1:
                current_subcluster.append(cluster[i])
                if current_subcluster:
                    new_cluster_residues.append(current_subcluster)
                current_subcluster = []
            else:
                current_subcluster.append(cluster[i])

        new_cluster_residues.append(current_subcluster)
        current_subcluster = []
    return new_cluster_residues

def cluster_bfactor(structure, labels):
    try:
        residues = [residue for residue in structure.get_residues()]
        assert len(residues) == len(labels), "Number of labels must match the number of residues in the structure"

        for residue, label in zip(residues, labels):
            for atom in residue:
                atom.set_bfactor(label)

        return structure
    except AssertionError as e:
        logger.error(f"AssertionError: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error in cluster_bfactor: {e}", exc_info=True)
        raise

def save_structure(structure, filename, selector=None):
    try:
        io = PDBIO()
        io.set_structure(structure)
        if selector is None:
            io.save(filename)
        else:
            io.save(filename, selector)
        logger.info(f"Structure saved to {filename}")
    except Exception as e:
        logger.error(f"Error saving structure to {filename}: {e}", exc_info=True)
        raise

def three_to_one(aa):
    aa_mapping = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_mapping.get(aa, aa)

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

def extract_chain(structure, chain_id):
    from Bio import PDB
    try:
        extracted_structure = PDB.Structure.Structure('extracted')
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    extracted_chain = chain.copy()
                    extracted_model = PDB.Model.Model(model.id)
                    extracted_model.add(extracted_chain)
                    extracted_structure.add(extracted_model)

        logger.info(f"Extracted chain {chain_id}")
        return extracted_structure
    except Exception as e:
        logger.error(f"Error extracting chain {chain_id}: {e}", exc_info=True)
        raise

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
            randrange(2, 50, 1), random(), randrange(1, 50, 1), randrange(1, 50, 1)
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

def normalize_bfactor(structure):
    try:
        max_bfactor = max(atom.get_bfactor() for atom in structure.get_atoms())
        min_bfactor = min(atom.get_bfactor() for atom in structure.get_atoms())

        for atom in structure.get_atoms():
            old_bfactor = atom.get_bfactor()
            normalized_bfactor = 100 * (old_bfactor - min_bfactor) / (max_bfactor - min_bfactor)
            atom.set_bfactor(normalized_bfactor)
        logger.info("B-factors normalized")
        return structure
    except Exception as e:
        logger.error("Error normalizing B-factors: {e}", exc_info=True)
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

def load_structure(pdb_file, only_one_chain):
    logger.info(f"Loading structure from {pdb_file}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    if only_one_chain:
        chain_id = 'A'  # Modify as needed
        structure = extract_chain(structure, chain_id)
    return structure

def preprocess_structure(structure, normalize_plddt):
    if normalize_plddt:
        structure = normalize_bfactor(structure)
    return structure

def calculate_sdp(structure, pair_distance_constant):
    sdp_matrix, total_plddt, residues = sdp_matrix_calculation(structure, pair_distance_constant)
    return sdp_matrix, total_plddt, residues

def perform_clustering(sdp_matrix, clustering_method, min_k, max_k, min_cluster_size, optimize_umap):
    umap_result = sdp_umap(sdp_matrix, optimize_umap, min_k, max_k)
    if clustering_method.lower() == 'hdbscan':
        clusters = hdbscan_clustering(min_cluster_size, umap_result)
    elif clustering_method.lower() == 'kmeans':
        clusters = kmeans_clustering(min_k, max_k, umap_result)
    else:
        logger.error("Invalid clustering method. Use 'hdbscan' or 'kmeans'.")
        sys.exit(1)
    return clusters, umap_result

def save_results(structure, clusters, residues, save_overlap_domains, save_pdb_domains, pdb_file, symmetry, min_cluster_size):
    cluster_residues = [[] for _ in range(len(np.unique(clusters)))]
    for i, cluster in enumerate(clusters):
        cluster_residues[cluster].append(residues[i])

    cluster_residues = split_cluster(cluster_residues, max_residue_jump=1)
    sorted_cluster_indices = np.argsort([np.mean([residue.id[1] for residue in cluster]) for cluster in cluster_residues])

    if save_overlap_domains:
        save_overlap_domains_to_fasta(cluster_residues, sorted_cluster_indices, pdb_file, symmetry, min_cluster_size)

    if save_pdb_domains:
        save_pdb_clusters(structure, cluster_residues, pdb_file)

def save_overlap_domains_to_fasta(cluster_residues, sorted_cluster_indices, pdb_file, symmetry, min_cluster_size):
    overlap_lengths = []
    k = 0
    concatenated_sequence = ''
    for i in range(len(sorted_cluster_indices) - 1):
        seq_records = []
        current_cluster_residue_list = cluster_residues[sorted_cluster_indices[i]]
        next_cluster_residue_list = cluster_residues[sorted_cluster_indices[i + 1]]

        current_cluster_sequence = ''.join(three_to_one(residue.get_resname()) for residue in current_cluster_residue_list)
        next_cluster_sequence = ''.join(three_to_one(residue.get_resname()) for residue in next_cluster_residue_list)
        concatenated_sequence = concatenated_sequence + current_cluster_sequence + next_cluster_sequence

        if len(concatenated_sequence) >= min_cluster_size:
            for j in range(symmetry):
                seq_record = SeqRecord(Seq(concatenated_sequence), id=f'domain_{sorted_cluster_indices[i] + 1}-{sorted_cluster_indices[i + 1] + 1}_copy{j + 1}', description='')
                seq_records.append(seq_record)
            overlap_length = len(next_cluster_sequence)
            overlap_lengths.append(overlap_length)
            fasta_filename = f'{pdb_file.split(".")[0]}_seq_{k}.fasta'
            SeqIO.write(seq_records, fasta_filename, "fasta")
            concatenated_sequence = ''
            k += 1

    csv_filename = f'{pdb_file.split(".")[0]}_overlaps_2.csv'
    with open(csv_filename, mode='w', newline='') as csvfile:
        fieldnames = ['Pair', 'Overlap Length']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for i, overlap_length in enumerate(overlap_lengths):
            pair_name = f'domain_{sorted_cluster_indices[i]+1}-{sorted_cluster_indices[i+1]+1}'
            writer.writerow({'Pair': pair_name, 'Overlap Length': overlap_length})

def save_pdb_clusters(structure, cluster_residues, pdb_file):
    for cluster_idx, cluster_residue_list in enumerate(cluster_residues):
        cluster_structure = structure.copy()

        class ResidueSelector(Select):
            def accept_residue(self, residue):
                return residue in cluster_residue_list

        filename = f'{pdb_file.split(".")[0]}_domain_{cluster_idx + 1}.pdb'
        save_structure(cluster_structure, filename, ResidueSelector())

def visualize_results(sdp_matrix, umap_result, clusters, pdb_file):
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.imshow(sdp_matrix, cmap='viridis', origin='upper')
    plt.colorbar(label='SDp')
    plt.xlabel('Residue Index')
    plt.ylabel('Residue Index')
    plt.title('SDp Matrix')

    plt.subplot(1, 2, 2)
    scatter = plt.scatter(umap_result[:, 0], umap_result[:, 1], c=clusters, cmap='viridis', s=5)
    plt.title('Clusters on UMAP Projection')
    plt.xlabel('UMAP Dimension 1')
    plt.ylabel('UMAP Dimension 2')

    combined_filename = f'{pdb_file.split(".")[0]}_combined_plots.png'
    plt.savefig(combined_filename, dpi=300)
    logger.info(f"Visualization saved to {combined_filename}")

def main(pdb_file, pair_distance_constant, min_k, max_k, clustering_method, min_cluster_size, symmetry, save_overlap_domains, save_pdb_domains, random_cut, optimize_umap, only_one_chain, normalize_plddt):
    logger.info("Starting main process")
    try:
        structure = load_structure(pdb_file, only_one_chain)
        structure = preprocess_structure(structure, normalize_plddt)
        sdp_matrix, total_plddt, residues = calculate_sdp(structure, pair_distance_constant)
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
    parser.add_argument("-pd", "--pair_distance_constant", type=float, default=_defaults['pair_distance_constant'], help=f'The pair_distance_constant will influcence on how the sDp is calculated: >1 values will increase the weight given to the plddt values of distant residues. It should always be >0 . Default: {_defaults["pair_distance_constant"]}')
    parser.add_argument("-n", "--number_of_chains", type=int, default=_defaults['number_of_chains'], help=f'Tail fibre symmetry. Default: {_defaults["number_of_chains"]}' )
    parser.add_argument("-so", "--save_overlap_domains", type=bool, default=_defaults['save_overlap_domains'], help=f'If including only one (False) or two domains per fasta file (True). Default: {_defaults["save_overlap_domains"]}' )
    parser.add_argument("-sv", "--save_pdb_domains", type=bool, default=_defaults['save_pdb_domains'], help=f'Whether to save the pdb file of the domains. Default: {_defaults["save_pdb_domains"]}' )
    parser.add_argument("-u", "--optimize_umap", type=bool, default=_defaults['optimize_umap'], help=f'Whether to self optimize umap parameters. Default: {_defaults["optimize_umap"]}' )
    parser.add_argument("-o", "--only_single_chain", type=bool, default=_defaults['only_single_chain'], help=f'Whether the model has a single chain. Default: {_defaults["only_single_chain"]}' )
    parser.add_argument("-np", "--normalize_plddt", type=bool, default=_defaults['normalize_plddt'], help=f'Whether to normalize plddt. Default: {_defaults["normalize_plddt"]}' )

    args = parser.parse_args()

    main(args.pdb, args.pair_distance_constant, args.min_k, args.max_k, args.clustering_method, args.min_domain_size, args.number_of_chains, args.save_overlap_domains, args.save_pdb_domains, False, args.optimize_umap, args.only_single_chain, args.normalize_plddt)
