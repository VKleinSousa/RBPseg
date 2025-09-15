from Bio.PDB import PDBParser, PDBIO, Select
import pandas as pd
import numpy as np
import os
from .superimpose import calculate_per_residue_rmsd
from .rename_chains import rename_chains_spherical, rename_chains_hh
from .trim_chains import delete_overhang
from .merge_chains import merge_chains
def rename_chain_if_exists(structure, original_chain_id='A', new_chain_id='D'):
    """
    Rename the chain with ID `original_chain_id` to `new_chain_id` if it exists in the structure.
    """
    for model in structure:
        for chain in model:
            if chain.id == original_chain_id:
                chain.id = new_chain_id
                print(f"Renamed chain {original_chain_id} to {new_chain_id}")
    return structure


def remove_all_hydrogens(structure):
    # Iterate over all atoms in the structure and collect hydrogen atoms for removal
    atoms_to_remove = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.element == "H":  # Check if the atom is hydrogen
                        atoms_to_remove.append((residue, atom.get_id()))

    # Remove the collected hydrogen atoms
    for residue, atom_id in atoms_to_remove:
        residue.detach_child(atom_id)
    return(structure)
    
def remove_oxt_atom(structure):
    # Iterate over all residues in all chains and models to find and remove 'OXT'
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'OXT' in residue:
                    residue.detach_child('OXT')  # Remove the OXT atom
    return(structure)

def get_structure_name(pdb_file):
    # Extract the basename and remove the extension
    basename = os.path.basename(pdb_file)
    structure_name = os.path.splitext(basename)[0]
    return structure_name
    
def superimpose_and_merge(pdb_files,sequence_counts, overhang_size, superposition, name, overhang_list, change_bfactor, chain_mode=0):
    parser = PDBParser(QUIET=True)

    if len(pdb_files) > sequence_counts[0]:
        fixed_structure = None
        moving_structure = None
        rmsd = 100000
        
        ov_index = 0
        overhang_size = int(overhang_list.iloc[:, ov_index])
        
        for i in range(0, sequence_counts[0]):
            print('i is:',i)
            fixed_structure = parser.get_structure(get_structure_name(pdb_files[i]), pdb_files[i])
            print('Fixed_structure is', pdb_files[i])

            # Rename chain 'A' to 'D' if it exists
            fixed_structure = remove_all_hydrogens(fixed_structure)
            fixed_structure = remove_oxt_atom(fixed_structure)
            fixed_structure = rename_chain_if_exists(fixed_structure)
            
            print('Fixed structure pruned')

            for j in range(sequence_counts[0], sequence_counts[0]+sequence_counts[1]):
                moving_structure = parser.get_structure(get_structure_name(pdb_files[j]), pdb_files[j])
                print('moving_structure is', pdb_files[j])

                # Rename chain 'A' to 'D' if it exists
                moving_structure = remove_all_hydrogens(moving_structure)
                moving_structure = remove_oxt_atom(moving_structure)
                moving_structure = rename_chain_if_exists(moving_structure)
                print('Moving structure pruned')

                if superposition == 0:
                    print('Pair:',[i, j])
                    moving_structure_tmp, per_residue_rmsd, overall_rmsd = calculate_per_residue_rmsd(
                        fixed_structure, moving_structure, overhang_size
                    )
                    magnitudes = np.linalg.norm(per_residue_rmsd, axis=1)
                    best_point = np.argmin(magnitudes)
                    
                    if magnitudes[best_point] < rmsd:
                        pair = [i, j]
                        best_fixed_structure = parser.get_structure(get_structure_name(pdb_files[i]), pdb_files[i])
                        best_moving_structure_name = get_structure_name(pdb_files[i])+'+'+get_structure_name(pdb_files[j])
                        best_moving_structure = parser.get_structure(best_moving_structure_name, pdb_files[j])
                        print(f'New pair: {best_moving_structure_name} and interface rmsd: {magnitudes[best_point]}')
                        rmsd = magnitudes[best_point]

        print('Best Structural Pair:',pair)
        best_fixed_structure = parser.get_structure(get_structure_name(pdb_files[pair[0]]), pdb_files[pair[0]])
        best_moving_structure_name = get_structure_name(pdb_files[pair[1]])+'+'+get_structure_name(pdb_files[pair[1]])
        best_moving_structure = parser.get_structure(best_moving_structure_name, pdb_files[j])
        # Rename chains in the chosen structures
        fixed_structure = rename_chain_if_exists(best_fixed_structure)
        moving_structure = rename_chain_if_exists(best_moving_structure)
        
        # Remove hydrogens
        fixed_structure = remove_all_hydrogens(fixed_structure)
        moving_structure = remove_all_hydrogens(moving_structure)
        
        # Remove oxt
        fixed_structure = remove_oxt_atom(fixed_structure)
        moving_structure = remove_oxt_atom(moving_structure)

        
        print('Merging first two fractions...')
        moving_structure, per_residue_rmsd, overall_rmsd = calculate_per_residue_rmsd(
            fixed_structure, moving_structure, overhang_size
        )
        magnitudes = np.linalg.norm(per_residue_rmsd, axis=1)
        best_point = np.argmin(magnitudes)

        #per_chain_rmsd = pd.DataFrame(per_chain_rmsd)
        interface_rmsd = per_residue_rmsd[best_point]
        deleted_overhang = delete_overhang(fixed_structure, moving_structure, overhang_size, best_point)
        if chain_mode == 0:
            moving_structure, chain_pairs = rename_chains_spherical(fixed_structure, moving_structure, change_bfactor)
        elif chain_mode == 1:
            moving_structure, chain_pairs = rename_chains_hh(fixed_structure, moving_structure, interface_rmsd, change_bfactor)
            
        merge_chains(moving_structure, chain_pairs)
        
        fixed_structure = moving_structure
        fixed_structure = remove_oxt_atom(fixed_structure)
        

        for i in range(sequence_counts[0]+sequence_counts[1], len(pdb_files), sequence_counts[1]):
            #print(i)
            if i + sequence_counts[0] > len(pdb_files):
                print('Finished aligning and merging')
            else:
                print('Merging next segment')
                pair = []
                ov_index += 1
                
                for j in range(0, sequence_counts[0]):
                    current_structure = parser.get_structure(get_structure_name(pdb_files[i + j]), pdb_files[i + j])
                    print('Current structure:', get_structure_name(pdb_files[i + j]))

                    # Rename chain 'A' to 'D' if it exists
                    current_structure = rename_chain_if_exists(current_structure)
                    current_structure = remove_all_hydrogens(current_structure)
                    current_structure = remove_oxt_atom(current_structure)
                    
                    rmsd = 100000
                    moving_structure = current_structure
                    overhang_size = int(overhang_list.iloc[:, ov_index])
                    #print('Overhang size:', overhang_size)
                    
                    if superposition == 0:
                        print('Finding optimal merging point...')
                        moving_structure, per_residue_rmsd, overall_rmsd = calculate_per_residue_rmsd(
                            fixed_structure, moving_structure, overhang_size
                        )
                        magnitudes = np.linalg.norm(per_residue_rmsd, axis=1)
                        best_point = np.argmin(magnitudes)
                        #print('Best point for merging:',  best_point)
                        
                        if magnitudes[best_point] < rmsd:
                            pair = [j + i]
                            rmsd = magnitudes[best_point]
                    else:
                        result_superimposer = local_superimpose(fixed_structure, moving_structure, overhang_size)
                        print('Running local superimposition. RMSD not available in this mode')

                best_moving_structure_name = best_moving_structure_name +'_'+get_structure_name(pdb_files[pair[0]])
                moving_structure = parser.get_structure(best_moving_structure_name, pdb_files[pair[0]])

                # Rename chain in moving_structure if needed
                moving_structure = rename_chain_if_exists(moving_structure)
                moving_structure = remove_oxt_atom(moving_structure)
                moving_structure = remove_all_hydrogens(moving_structure)
                
                moving_structure, per_residue_rmsd, overall_rmsd= calculate_per_residue_rmsd(
                    fixed_structure, moving_structure, overhang_size
                )
                magnitudes = np.linalg.norm(per_residue_rmsd, axis=1)
                best_point = np.argmin(magnitudes)
                
                #per_chain_rmsd = pd.DataFrame(per_chain_rmsd)
                interface_rmsd = per_residue_rmsd[best_point]
                deleted_overhang = delete_overhang(fixed_structure, moving_structure, overhang_size, best_point)

                    
                if chain_mode == 0:
                    moving_structure, chain_pairs = rename_chains_spherical(fixed_structure, moving_structure, change_bfactor)
                elif chain_mode == 1:
                    moving_structure, chain_pairs = rename_chains_hh(fixed_structure, moving_structure, interface_rmsd, change_bfactor)
                else:
                    raise('Select a valid mode for the chain pairing. 0 : spherical 1: hierarchical')
                merge_chains(moving_structure, chain_pairs)
                fixed_structure = moving_structure
    else:
        moving_structure = parser.get_structure(get_structure_name(pdb_files[0]), pdb_files[0])

    pdb_io = PDBIO()
    pdb_io.set_structure(moving_structure)
    pdb_io.save(name)
    return moving_structure
    
if __name__ == "__main__":
    # Example usage
    pdb_files = ["path/to/pdb1.pdb", "path/to/pdb2.pdb", "path/to/pdb3.pdb"]  # Replace with actual PDB file paths
    overhang_size = 50
    superposition = 0
    name = "merged_output.pdb"
    overhang_list = pd.DataFrame([[50], [40], [30]])  # Replace with actual data
    change_bfactor = False

    final_structure = superimpose_and_merge(pdb_files, overhang_size, superposition, name, overhang_list, change_bfactor)
    print("Final structure saved as:", name)

