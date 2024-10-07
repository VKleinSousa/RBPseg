import os
import shutil
import argparse

def prepare_files(input_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    # Iterate through each folder in the input path
    for folder_name in os.listdir(input_path):
        folder_path = os.path.join(input_path, folder_name)
        if os.path.isdir(folder_path):
            # Extract pdb_name and n from the folder name
            if "_seq_" in folder_name:
                pdb_name, seq_n = folder_name.split("_seq_")
                
                # Iterate through each .pdb file in the folder
                for pdb_file in os.listdir(folder_path):
                    if pdb_file.endswith(".pdb") and pdb_file.startswith("ranked_"):
                        # Extract the rank number (m) from the file name
                        rank_number = pdb_file.split("_")[1].split(".")[0]
                        
                        # Prepare the new file name
                        new_file_name = f"{seq_n}_{pdb_name}_{rank_number}.pdb"
                        
                        # Copy the file to the output directory with the new name
                        source_file = os.path.join(folder_path, pdb_file)
                        destination_file = os.path.join(output_path, new_file_name)
                        shutil.copyfile(source_file, destination_file)
                        
                        print(f"Copied {source_file} to {destination_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare PDB files for TFold merging.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input directory containing PDB folders.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output directory where the prepared files will be saved.')

    args = parser.parse_args()
    
    prepare_files(args.input, args.output)
