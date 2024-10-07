import os
import argparse
import pandas as pd
import re
from Bio.PDB import PDBParser, PDBIO
from .superimpose_and_merge import superimpose_and_merge
from .relax_model_openmm import relax_amber14_score  # Assuming this is an available module


def relax_model(merged_path,relaxed_merged_path,iter):
    relax_amber14_score(merged_path,relaxed_merged_path,iter=iter)

def load_overhang_file(overhang_file):
    """
    Load the overhang file and return it as a DataFrame.
    """
    if overhang_file:
        try:
            df = pd.read_csv(overhang_file)
            overhang_list = pd.DataFrame(df)
            overhang_list = overhang_list.transpose().iloc[1:]
            return overhang_list
        except pd.errors.EmptyDataError:
            print("Warning: Overhang file is empty.")
            return None
    return None

def apply_b_factors(merged_structure, relaxed_structure):
    """
    Apply B-factors from the merged structure to the relaxed structure.
    """
    bfactor_vector = [atom.get_bfactor() for model in merged_structure for atom in model.get_atoms()]
    bfactor_index = 0

    for model in relaxed_structure:
        for atom in model.get_atoms():
            if bfactor_index < len(bfactor_vector):
                atom.set_bfactor(bfactor_vector[bfactor_index])
                bfactor_index += 1
            else:
                print("Warning: B-factor vector is shorter than the number of atoms.")
                break

def main():
    # Create an ArgumentParser object to handle command-line arguments
    parser = argparse.ArgumentParser(description='Perform superimposition and merging.')

    # Add command-line flags for directory_path, overhang_size, and superimpose_function
    parser.add_argument('-d', '--directory', required=True, help='Path to the directory containing PDB files')
    parser.add_argument('-o', '--overhang', type=int, default=50, help='Overhang size (default: 50)')
    parser.add_argument('-of', '--overhang_file', type=str, help='Overhangs file')
    parser.add_argument('-f', '--function', type=int, default=0, help='Superimpose function: 0 to global, 1 to local (default: 1)')
    parser.add_argument('-n', '--save_name', type=str, default='merged.pdb', help='Name of the final file (default: merged.pdb)')
    parser.add_argument('-r', '--relax',  action='store_true', help='Add -r to run amber relaxation. ')
    parser.add_argument('-b', '--bfactor', action='store_true', help='Add -b to change the bfactor to zero in segment points')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Retrieve values from the parsed arguments
    directory_path = args.directory
    overhang_size = args.overhang
    superimpose_function = args.function
    change_bfactor = args.bfactor
    name = os.path.join(directory_path, args.save_name)
    overhang_file = args.overhang_file

    # Load the overhang sizes from the file, if provided
    overhang_list = load_overhang_file(overhang_file)

    # Get the list of PDB files from the specified directory
    file_list = os.listdir(directory_path)
    pdb_files = sorted([os.path.join(directory_path, filename) for filename in file_list if filename.endswith('.pdb')],
                       key=lambda x: tuple(map(int, re.findall(r'\d+', x))))
    print('Reordered PDB files:', pdb_files)

    # Check if there are multiple PDB files for alignment and merging
    if len(pdb_files) <= 1:
        print('No alignment and merging needed, only one model found.')
        return

    # Perform superimposition and merging
    merged_structure = superimpose_and_merge(pdb_files, overhang_size, superimpose_function, name, overhang_list, change_bfactor)

    # Optionally apply amber relaxation
    if args.relax == True:
        name_relax = os.path.join(directory_path, name + '_relax.pdb')
        relax_model(name, name_relax, 0)
        parser = PDBParser(QUIET=True)
        # Parse the current PDB file
        relax_structure = parser.get_structure('current_structure', name_relax)
        
        bfactor_vector = []

        # Extract B-factors from merged_structure and save to the vector
        for model1 in merged_structure:
            for atom1 in model1.get_atoms():
                bfactor_vector.append(atom1.get_bfactor())

        # Apply B-factors from the vector to relax_structure
        for model2 in relax_structure:
            for chain2 in model2:
                bfactor_index = 0
                for atom2 in chain2.get_atoms():
                    # Check if there are still B-factors in the vector
                    if bfactor_index < len(bfactor_vector):
                        # Set B-factor from the vector
                        atom2.set_bfactor(bfactor_vector[bfactor_index])
                        bfactor_index += 1
                    else:
                        # Handle the case where the vector is shorter than expected
                        print("Warning: B-factor vector is shorter than the number of atoms.")
                        break

        pdb_io = PDBIO()
        pdb_io.set_structure(relax_structure)
        pdb_io.save(name_relax)

if __name__ == "__main__":
    main()
