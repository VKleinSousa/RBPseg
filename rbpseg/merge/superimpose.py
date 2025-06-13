import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
import subprocess
import os
import tempfile


def calculate_per_residue_rmsd(fixed_structure, moving_structure, overhang_size, method="usalign", output_prefix="usalign"):
    """
    Calculate per-residue RMSD for two structures, including the transformed structure and overall RMSD.
    """
    print(f"Performing superimposition using method: {method}")
    
    # Perform superimposition
    transformed_structure, overall_rmsd, diff = superimpose_overhangs(
        fixed_structure, moving_structure, overhang_size, method=method, output_prefix=output_prefix
    )
    return transformed_structure, diff, overall_rmsd

def extract_atoms_from_pdb(pdb_file):
    """
    Extract all atoms from a PDB file.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        atoms (list): List of atoms in the PDB file.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    
    # Get the first chain
    first_chain = next(structure.get_chains())
    
    # Extract C-alpha atoms from the first chain
    calpha_atoms = [atom for residue in first_chain.get_residues()
                    if residue.get_resname() not in ['H_SO4', 'HOH']
                    for atom in residue if atom.get_name() == "CA"]
    
    return calpha_atoms

def superimpose_overhangs(fixed_structure, moving_structure, overhang_size, method="usalign", output_prefix="usalign"):
    """
    Superimpose two structures based on overhang regions.
    """
    if method == "usalign":
        # Create a temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            fixed_file = os.path.join(tmpdir, f"{output_prefix}_fixed.pdb")
            moving_file = os.path.join(tmpdir, f"{output_prefix}_moving.pdb")
            aligned_file_prefix = os.path.join(tmpdir, output_prefix)
            transform_matrix_file = os.path.join(tmpdir, f"{output_prefix}_matrix.txt")

            # Save overhangs to temporary PDB files
            extract_overhang_pdb(fixed_structure, overhang_size, fixed_file, position="C-terminal")
            extract_overhang_pdb(moving_structure, overhang_size, moving_file, position="N-terminal")

            # Run USalign and return the superimposed file
            aligned_pdb = run_usalign(fixed_file, moving_file, output_prefix=aligned_file_prefix)
            # Apply the transformation to the entire moving structure
            aligned_structure = apply_usalign_transformation(moving_structure, transform_matrix_file)
            
            fixed_atoms = extract_atoms_from_pdb(fixed_file)
            aligned_atoms = extract_overhang_pdb(aligned_structure, overhang_size, moving_file, position="N-terminal")
            aligned_atoms = extract_atoms_from_pdb(moving_file)
            
            rmsd, diff = calculate_rmsd(fixed_atoms, aligned_atoms)

            return aligned_structure, rmsd, diff


def extract_overhang_pdb(structure, overhang_size, output_file, position="C-terminal"):
    """
    Extract overhang residues and save them to a new PDB file.
    """
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, OverhangSelect(overhang_size, position))


def calculate_rmsd_first_chain(fixed_structure, aligned_structure):
    """
    Calculate the RMSD between the fixed atoms and the aligned moving atoms for the first chain.

    Parameters:
        fixed_structure (Structure): Biopython Structure object for the fixed structure.
        aligned_structure (Structure): Biopython Structure object for the aligned structure.

    Returns:
        rmsd (float): Root Mean Square Deviation between the two sets of atoms.
        diff (np.ndarray): Array of differences between the coordinates.
    """
    # Get the first chain from both structures
    fixed_chain = next(fixed_structure.get_chains())
    aligned_chain = next(aligned_structure.get_chains())

    # Get the atoms from the first chain
    fixed_atoms = [atom for atom in fixed_chain.get_atoms()]
    aligned_atoms = [atom for atom in aligned_chain.get_atoms()]

    # Ensure both chains have the same number of atoms
    if len(fixed_atoms) != len(aligned_atoms):
        raise ValueError("The fixed and aligned atoms must have the same number of atoms.")

    # Extract coordinates
    fixed_coords = np.array([atom.get_coord() for atom in fixed_atoms])
    aligned_coords = np.array([atom.get_coord() for atom in aligned_atoms])

    # Calculate RMSD
    diff = fixed_coords - aligned_coords
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd, diff
    
def calculate_rmsd(fixed_atoms, aligned_atoms):
    """
    Calculate the RMSD between the fixed atoms and the aligned moving atoms.

    Parameters:
        fixed_atoms (list): List of fixed atoms (C-alpha or other relevant atoms).
        aligned_atoms (list): List of aligned moving atoms (C-alpha or other relevant atoms).

    Returns:
        rmsd (float): Root Mean Square Deviation between the two sets of atoms.
        diff (np.ndarray): Array of differences between the coordinates.
    """
    # Ensure both overhangs have the same number of atoms
    if len(fixed_atoms) != len(aligned_atoms):
        raise ValueError("The fixed and aligned atoms must have the same number of atoms.")

    # Extract coordinates
    fixed_coords = np.array([atom.get_coord() for atom in fixed_atoms])
    aligned_coords = np.array([atom.get_coord() for atom in aligned_atoms])
    
    # Calculate RMSD
    diff = fixed_coords - aligned_coords
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

    return rmsd, diff
    
def run_usalign(fixed_file, moving_file, output_prefix="usalign_output"):
    """
    Run USalign in multimeric mode to superimpose two structures.
    """
    transform_matrix_file = f"{output_prefix}_matrix.txt"
    output_pdb = f"{output_prefix}_sup.pdb"
    cmd = f"USalign {moving_file} {fixed_file} -mm 1 -ter 1 -m {transform_matrix_file} -o {output_prefix}"
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(f"USalign failed:\n{result.stderr}")
    #rmsd = calculate_rmsd(fixed_file,output_pdb)
    return output_pdb #rmsd


def apply_usalign_transformation(structure, transform_matrix_file, output_pdb="transformed_structure.pdb"):
    """
    Apply the transformation matrix generated by USalign to a structure and save the transformed structure.

    Parameters:
        structure (Structure): The structure to transform.
        transform_matrix_file (str): Path to the file containing the transformation matrix.
        output_pdb (str): Path to save the transformed structure as a PDB file.
    """
    with open(transform_matrix_file, 'r') as f:
        lines = f.readlines()
    
    # Skip the header lines and extract the translation vector and rotation matrix
    try:
        data_lines = lines[2:5]  # Assuming the matrix and translation start from line 2
        translation = np.array([float(line.split()[1]) for line in data_lines])
        rotation = np.array([
            [float(x) for x in line.split()[2:]] for line in data_lines
        ])
    except ValueError as e:
        raise ValueError(f"Failed to parse transformation matrix: {e}")

    # Apply the transformation to all atoms in the structure
    for atom in structure.get_atoms():
        coord = np.array(atom.get_coord())
        transformed_coord = np.dot(rotation, coord) + translation
        atom.set_coord(transformed_coord)

    # Save the transformed structure to a PDB file
    io = PDBIO()
    io.set_structure(structure)

    return structure

class OverhangSelect(Select):
    """
    Selection class for extracting overhang residues.
    """
    def __init__(self, overhang_size, position="C-terminal"):
        self.overhang_size = overhang_size
        self.position = position

    def accept_residue(self, residue):
        chain_residues = list(residue.get_parent().get_residues())
        if self.position == "C-terminal":
            return residue in chain_residues[-self.overhang_size:]
        elif self.position == "N-terminal":
            return residue in chain_residues[:self.overhang_size]
        return False

# Main script
if __name__ == "__main__":
    parser = PDBParser(QUIET=True)
    fixed_structure = parser.get_structure("fixed", "fixed.pdb")
    moving_structure = parser.get_structure("moving", "moving.pdb")
    overhang_size = 10

    transformed_structure, per_residue_rmsd, overall_rmsd, per_chain_rmsd = calculate_per_residue_rmsd(
        fixed_structure, moving_structure, overhang_size, method="usalign"
    )
    print("Overall RMSD:", overall_rmsd)
    print("Per-residue RMSD:", per_residue_rmsd)
