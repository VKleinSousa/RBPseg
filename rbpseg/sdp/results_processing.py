from Bio.PDB import PDBIO, Select
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import csv
import matplotlib.pyplot as plt
import logging
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

def three_to_one(aa):
    aa_mapping = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    return aa_mapping.get(aa, aa)

def split_cluster(cluster_residues, max_residue_jump=50):
    new_cluster_residues = []
    current_subcluster = []
    for cluster in cluster_residues:
        for i in range(len(cluster)):
            current_residue_id = cluster[i].id[1]
            if current_residue_id == cluster[len(cluster)-1].id[1]:
                current_subcluster.append(cluster[i])
                break
            next_residue_id = cluster[i + 1].id[1]

            if next_residue_id - current_residue_id > max_residue_jump:
                current_subcluster.append(cluster[i])
                if current_subcluster:
                    new_cluster_residues.append(current_subcluster)
                current_subcluster = []
            else:
                current_subcluster.append(cluster[i])

        new_cluster_residues.append(current_subcluster)
        current_subcluster = []
    return new_cluster_residues


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

def save_results(structure, clusters, residues, save_overlap_domains, save_pdb_domains, pdb_file, symmetry, min_cluster_size, min_ov_size):
    cluster_residues = [[] for _ in range(len(np.unique(clusters)))]
    for i, cluster in enumerate(clusters):
        cluster_residues[cluster].append(residues[i])

    cluster_residues = split_cluster(cluster_residues)
    sorted_cluster_indices = np.argsort([np.mean([residue.id[1] for residue in cluster]) for cluster in cluster_residues])
    # Check the last sorted cluster for minimum size
    last_cluster = cluster_residues[sorted_cluster_indices[-1]]
    
    while len(last_cluster) < min_cluster_size:
        # Concatenate with the previous cluster
        if len(sorted_cluster_indices) < 2:
            raise ValueError("Not enough clusters to meet the minimum cluster size requirement.")
        previous_cluster_index = sorted_cluster_indices[-2]
        last_cluster = cluster_residues[previous_cluster_index] + last_cluster
        
        # Remove the previous cluster as it's now part of the last one
        sorted_cluster_indices = sorted_cluster_indices[:-1]
        cluster_residues[sorted_cluster_indices[-1]] = last_cluster

        # Break if the concatenation meets the size
        if len(last_cluster) >= min_cluster_size:
            break

    if save_overlap_domains:
        save_overlap_domains_to_fasta(cluster_residues, sorted_cluster_indices, pdb_file, symmetry, min_cluster_size,min_cluster_size)

    if save_pdb_domains:
        save_pdb_clusters(structure, cluster_residues, pdb_file)

def save_overlap_domains_to_fasta(cluster_residues, sorted_cluster_indices, pdb_file, symmetry, min_cluster_size,min_ov_size):
    overlap_lengths, pairs, k, concatenated_sequence, domains_combined, overlap_length = [], [], 0, '', '', 0
    
    for i in range(len(sorted_cluster_indices) - 1):
        seq_records = []
        current_cluster_residue_list = cluster_residues[sorted_cluster_indices[i]]
        next_cluster_residue_list = cluster_residues[sorted_cluster_indices[i + 1]]

        current_cluster_sequence = ''.join(three_to_one(residue.get_resname()) for residue in current_cluster_residue_list)
        next_cluster_sequence = ''.join(three_to_one(residue.get_resname()) for residue in next_cluster_residue_list)
        concatenated_sequence_tmp = concatenated_sequence + current_cluster_sequence + next_cluster_sequence        
        overlap_length += len(next_cluster_sequence)

        if len(concatenated_sequence_tmp) >= min_cluster_size and overlap_length >= min_ov_size:
            domains_combined = domains_combined + f'{sorted_cluster_indices[i] + 1}-' + f'{sorted_cluster_indices[i + 1] + 1}' 
            for j in range(symmetry):
                seq_record = SeqRecord(Seq(concatenated_sequence_tmp), id=domains_combined, description='')
                seq_records.append(seq_record)
            overlap_lengths.append(overlap_length)
            pairs.append(domains_combined)
            domains_combined = ''
            fasta_filename = f'{pdb_file.split(".")[0]}_seq_{k}.fasta'
            SeqIO.write(seq_records, fasta_filename, "fasta")
            concatenated_sequence = ''
            overlap_length = 0
            k += 1
    
        else:
            concatenated_sequence = concatenated_sequence + current_cluster_sequence
            domains_combined = domains_combined + f'{sorted_cluster_indices[i] + 1}-'
        
    df_overhang = pd.DataFrame({'Concatenated Domains': pairs, 'Overhang Length': overlap_lengths})
    df_overhang.to_csv(f'{pdb_file.split(".")[0]}_overhangs.csv', index=False)

def save_pdb_clusters(structure, cluster_residues, pdb_file):
    # Get all chains in the structure
    chains = [chain for chain in structure.get_chains()]
    
    if len(chains) < 3:
        print("Warning: Expected three chains in the structure.")
        return
    
    first_chain = chains[0]
    second_chain = chains[1]
    third_chain = chains[2]
    
    # Get the number of residues in the first chain
    num_residues_in_chain = len(first_chain)

    for cluster_idx, cluster_residue_list in enumerate(cluster_residues):
        cluster_structure = structure.copy()

        # Extend the cluster_residue_list to include corresponding residue numbers from the second and third chains
        extended_cluster_residues = [residue.id[1] for residue in cluster_residue_list]  # Extract the residue numbers

        # Add residues from the second and third chains
        for residue_number in extended_cluster_residues.copy():
            extended_cluster_residues.append(residue_number + num_residues_in_chain)  # Second chain
            extended_cluster_residues.append(residue_number + 2 * num_residues_in_chain)  # Third chain

        class ResidueSelector(Select):
            def accept_residue(self, residue):
                # Check if the residue's number is part of the extended cluster residues
                return residue.id[1] in extended_cluster_residues

        # Save the structure for this cluster
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
