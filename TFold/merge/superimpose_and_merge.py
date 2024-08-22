from Bio.PDB import PDBParser, PDBIO
import pandas as pd
import numpy as np
from .superimpose import calculate_per_residue_rmsd
from .rename_chains import rename_chains_2
from .trim_chains import delete_overhang
from .merge_chains import merge_chains

def superimpose_and_merge(pdb_files, overhang_size, superposition, name, overhang_list, change_bfactor):
    parser = PDBParser(QUIET=True)
    
    if len(pdb_files) > 5:
        fixed_structure = None
        moving_structure = None
        rmsd = 100000
        
        ov_index = 0 
        overhang_size = int(overhang_list.iloc[:, ov_index])
        print(overhang_list)
        print(overhang_size)
        
        for i in range(0, 5):
            fixed_structure = parser.get_structure('fixed_structure', pdb_files[i])

            for j in range(5, 10):
                moving_structure = parser.get_structure('moving_structure', pdb_files[j])   
                
                if superposition == 0:
                    moving_structure, per_residue_rmsd, overall_rmsd, per_chain_rmsd = calculate_per_residue_rmsd(
                        fixed_structure, moving_structure, overhang_size
                    )
                    best_point = np.argmin(per_residue_rmsd)
                    per_chain_rmsd = pd.DataFrame(per_chain_rmsd)
                    interface_rmsd = per_chain_rmsd.iloc[:, [np.argmin(per_residue_rmsd)]]
                    
                    if np.mean(interface_rmsd) < rmsd:
                        pair = [i, j]
                        print('New pair and interface rmsd:', pair, np.mean(interface_rmsd))
                        rmsd = np.mean(interface_rmsd)

        fixed_structure = parser.get_structure('fixed_structure', pdb_files[pair[0]])
        moving_structure = parser.get_structure('moving_structure', pdb_files[pair[1]])
        moving_structure, per_residue_rmsd, overall_rmsd, per_chain_rmsd = calculate_per_residue_rmsd(
            fixed_structure, moving_structure, overhang_size
        )
        best_point = np.argmin(per_residue_rmsd)
        per_chain_rmsd = pd.DataFrame(per_chain_rmsd)
        interface_rmsd = per_chain_rmsd.iloc[:, [np.argmin(per_residue_rmsd)]]
        print(interface_rmsd , per_chain_rmsd)
        deleted_overhang = delete_overhang(fixed_structure, moving_structure, overhang_size, np.argmin(per_residue_rmsd))
        moving_structure, chain_pairs = rename_chains_2(fixed_structure, moving_structure, interface_rmsd, change_bfactor)
        merge_chains(moving_structure, chain_pairs)
        
        fixed_structure = moving_structure

        for i in range(10, len(pdb_files), 5):
            if i + 5 > len(pdb_files):
                print('Finished aligning and merging')
            else:
                print('Merging next segment')
                pair = []
                ov_index += 1 
                
                for j in range(0, 5):
                    current_structure = parser.get_structure('current_structure', pdb_files[i + j])
                    print('Current structure:', pdb_files[i + j])

                    rmsd = 100000
                    moving_structure = current_structure
                    overhang_size = int(overhang_list.iloc[:, ov_index])
                    print('Overhang size:', overhang_size)
                    
                    if superposition == 0:
                        moving_structure, per_residue_rmsd, overall_rmsd, per_chain_rmsd = calculate_per_residue_rmsd(
                            fixed_structure, moving_structure, overhang_size
                        )
                        best_point = np.argmin(per_residue_rmsd)
                        print('Best point:', len(per_residue_rmsd), best_point)
                        per_chain_rmsd = pd.DataFrame(per_chain_rmsd)
                        interface_rmsd = per_chain_rmsd.iloc[:, [np.argmin(per_residue_rmsd)]]
                        
                        if np.mean(interface_rmsd) < rmsd:
                            pair = [j + i]
                            rmsd = np.mean(interface_rmsd)
                    else:
                        result_superimposer = local_superimpose(fixed_structure, moving_structure, overhang_size)
                        print('Running local superimposition. RMSD not available in this mode')

                moving_structure = parser.get_structure('moving_structure', pdb_files[pair[0]])
                moving_structure, per_residue_rmsd, overall_rmsd, per_chain_rmsd = calculate_per_residue_rmsd(
                    fixed_structure, moving_structure, overhang_size
                )
                best_point = np.argmin(per_residue_rmsd)
                print('Best point:', len(per_residue_rmsd), best_point)
                per_chain_rmsd = pd.DataFrame(per_chain_rmsd)
                interface_rmsd = per_chain_rmsd.iloc[:, [np.argmin(per_residue_rmsd)]]

                deleted_overhang = delete_overhang(fixed_structure, moving_structure, overhang_size, np.argmin(per_residue_rmsd))
                moving_structure, chain_pairs = rename_chains_2(fixed_structure, moving_structure, interface_rmsd, change_bfactor)

                merge_chains(moving_structure, chain_pairs)
                fixed_structure = moving_structure
    else:
        moving_structure = parser.get_structure('current_structure', pdb_files[0])

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
