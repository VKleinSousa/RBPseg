import os
import argparse
from Bio.PDB import MMCIFParser, PDBIO

def convert_cif_to_pdb(input_dir, output_dir):
    """
    Converts all .cif files in the input directory to .pdb format and saves them in the output directory.
    
    Args:
    - input_dir (str): Path to the directory containing .cif files.
    - output_dir (str): Path to the directory where .pdb files will be saved.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Initialize MMCIF parser and PDB writer
    parser = MMCIFParser(QUIET=True)
    pdb_writer = PDBIO()

    # Loop through all files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".cif"):
            cif_path = os.path.join(input_dir, file_name)
            pdb_path = os.path.join(output_dir, os.path.splitext(file_name)[0] + ".pdb")
            
            try:
                # Parse the CIF file
                structure = parser.get_structure(file_name, cif_path)
                # Write the structure to a PDB file
                pdb_writer.set_structure(structure)
                pdb_writer.save(pdb_path)
                print(f"Converted: {cif_path} -> {pdb_path}")
            except Exception as e:
                print(f"Error converting {cif_path}: {e}")

if __name__ == "__main__":
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Convert .cif files to .pdb format.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input directory containing .cif files.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory for .pdb files.")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the function with provided arguments
    convert_cif_to_pdb(args.input, args.output)

