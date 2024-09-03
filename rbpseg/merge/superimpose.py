from Bio.PDB import PDBParser
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

def calculate_per_residue_rmsd(fixed_structure, moving_structure, overhang_size):
    # Perform the superimposition and obtain the transformed structure and RMSD
    print(overhang_size)
    print(fixed_structure)
    print(moving_structure)
    transformed_structure, overall_rmsd = superimpose_overhangs(fixed_structure, moving_structure, overhang_size)
    
    # Initialize a list to store per-residue RMSD for each chain
    per_chain_rmsd = []

    # Iterate through each chain in the moving structure
    for model_fixed, model_moving in zip(fixed_structure, transformed_structure):
        # Initialize a list to store per-residue RMSD for a single chain
        per_residue_rmsd_chain = []

        # Iterate through each chain_fixed and chain_moving
        for chain_fixed in model_fixed:
            for chain_moving in model_moving:
                residues_fixed = list(chain_fixed.get_residues())
                selected_residues_fixed = residues_fixed[-overhang_size:]
                residues_moving = list(chain_moving.get_residues())
                selected_residues_moving = residues_moving[:overhang_size]

                # Initialize a list to store per-residue RMSD for a single chain
                per_residue_rmsd_chain = []

                # Iterate through residues in the chain
                for residue_fixed, residue_moving in zip(selected_residues_fixed, selected_residues_moving):
                    # Initialize variables to calculate per-residue RMSD
                    residue_atoms_fixed = list(residue_fixed.get_atoms())
                    residue_atoms_moving = list(residue_moving.get_atoms())
                    n_atoms = len(residue_atoms_fixed)
                    total_squared_diff = 0.0

                    # Calculate the RMSD between corresponding atoms in the residue
                    for atom_fixed, atom_moving in zip(residue_atoms_fixed, residue_atoms_moving):
                        diff = atom_fixed.get_coord() - atom_moving.get_coord()
                        total_squared_diff += np.sum(diff ** 2)

                    # Calculate the RMSD for the residue
                    residue_rmsd = np.sqrt(total_squared_diff / n_atoms)
                    per_residue_rmsd_chain.append(residue_rmsd)

                # Append the per-residue RMSD for this chain to the list
                per_chain_rmsd.append(per_residue_rmsd_chain)

    # Calculate the average per-residue RMSD
    per_residue_rmsd = np.mean(per_chain_rmsd, axis=0)

    return transformed_structure, per_residue_rmsd, overall_rmsd, per_chain_rmsd
    
def superimpose_overhangs(fixed_structure, moving_structure, overhang_size):
    # Calculates the global superimposition of both models.
    last_fixed = []  # saves the C terminal residues of the fixed structure
    first_fixed = []  # saves the N terminal residues of the fixed structure
    
    for model in fixed_structure:
        for chain in model:
            residues = list(chain.get_residues())
            #print(residues.get_atoms())
            selected_residues_last = residues[-overhang_size:]
            selected_residues_first = residues[:overhang_size]
            for residue in selected_residues_last:
                last_fixed.extend(residue.get_atoms())
            for residue in selected_residues_first:
                first_fixed.extend(residue.get_atoms())

                    
    first_moving = []  # saves the N terminal residues of the moving structure
    last_moving = []  # saves the C terminal residues of the moving structure
    
    for model in moving_structure:
        for chain in model:
            residues = list(chain.get_residues())
            #print(residues.get_atoms())
            selected_residues_last = residues[-overhang_size:]
            selected_residues_first = residues[:overhang_size]
            for residue in selected_residues_first:
                first_moving.extend(residue.get_atoms())
            for residue in selected_residues_last:
                last_moving.extend(residue.get_atoms())

    # Calculate the centroids of the overhang atoms
    centroid_fixed = sum(atom.get_coord() for atom in last_fixed) / len(last_fixed)
    centroid_moving = sum(atom.get_coord() for atom in first_moving) / len(first_moving)
    centroid_fixed_first = sum(atom.get_coord() for atom in first_fixed) / len(first_fixed)
    centroid_moving_last = sum(atom.get_coord() for atom in last_moving) / len(last_moving)
    
    if centroid_fixed[0] * centroid_fixed[1] * centroid_fixed[2] * centroid_moving_last[0] * centroid_moving_last[1] * centroid_moving_last[2] < 0:
        # If models are pointing to different directions, flip the moving structure.
        for atom in moving_structure.get_atoms():
            atom.set_coord(-1 * atom.get_coord())
        # Recalculate the centroids of the moving structure
        centroid_moving = sum(atom.get_coord() for atom in first_moving) / len(first_moving)
        centroid_moving_last = sum(atom.get_coord() for atom in last_moving) / len(last_moving)
                    
    vector_fixed = np.array([atom.get_coord() for atom in last_fixed])
    vector_moving = np.array([atom.get_coord() for atom in first_moving])


    print(f"Fixed vector size: {vector_fixed.shape}")
    print(f"Moving vector size: {vector_moving.shape}")

    # Initialize the Superimposer with the selected atoms
    superimposer = SVDSuperimposer()
    superimposer.set(vector_fixed, vector_moving)
    superimposer.run()
    rot, tran = superimposer.get_rotran()
    
    for atom in moving_structure.get_atoms():
        if centroid_fixed[0] * centroid_fixed[1] * centroid_fixed[2] * centroid_moving_last[0] * centroid_moving_last[1] * centroid_moving_last[2] < 0:
            atom.set_coord(-1 * np.dot(atom.get_coord(), rot) + tran)
        else:
            atom.set_coord(np.dot(atom.get_coord(), rot) + tran)
    
    return moving_structure, superimposer.rms

if __name__ == "__main__":
    # Example usage
    parser = PDBParser(QUIET=True)
    fixed_structure = parser.get_structure("fixed_structure", "fixed.pdb")
    moving_structure = parser.get_structure("moving_structure", "moving.pdb")
    overhang_size = 10  # Example overhang size

    transformed_structure, per_residue_rmsd, overall_rmsd, per_chain_rmsd = calculate_per_residue_rmsd(fixed_structure, moving_structure, overhang_size)
    print("Overall RMSD:", overall_rmsd)
    print("Per-residue RMSD:", per_residue_rmsd)

