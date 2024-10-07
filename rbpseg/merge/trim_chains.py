from Bio.PDB import PDBParser

def delete_overhang(fixed_structure, moving_structure, overhang_size, optimal_residue):
    """
    Delete the overhang residues from both the fixed and moving structures.
    """
    for chain_to_modify_fixed, chain_to_modify_moving in zip(fixed_structure[0], moving_structure[0]):
        residues_fixed = list(chain_to_modify_fixed.get_residues())
        segment_point_fixed = overhang_size - optimal_residue
        residues_moving = list(chain_to_modify_moving.get_residues())
        segment_point_moving = optimal_residue

        # Detach residues from the end of the fixed structure chain
        for residue in residues_fixed[-segment_point_fixed:]:
            chain_to_modify_fixed.detach_child(residue.id)
        
        # Detach residues from the beginning of the moving structure chain
        for residue in residues_moving[:segment_point_moving]:
            chain_to_modify_moving.detach_child(residue.id)

if __name__ == "__main__":
    # Example usage
    parser = PDBParser(QUIET=True)
    fixed_structure = parser.get_structure("fixed_structure", "fixed.pdb")
    moving_structure = parser.get_structure("moving_structure", "moving.pdb")

    overhang_size = 50  # Example overhang size
    optimal_residue = 10  # Example optimal residue position

    delete_overhang(fixed_structure, moving_structure, overhang_size, optimal_residue)
    print("Overhang residues deleted successfully.")
