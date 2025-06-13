import os
import shutil
import argparse
from .convert_cif_to_pdb import convert_cif_to_pdb
from .copy_cif_files import copy_and_rename_cif_files

def prepare_files_af2(input_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    
    for folder_name in os.listdir(input_path):
        folder_path = os.path.join(input_path, folder_name)
        if os.path.isdir(folder_path):
            if "_seq_" in folder_name:
                pdb_name, seq_n = folder_name.split("_seq_")
                for pdb_file in os.listdir(folder_path):
                    if pdb_file.endswith(".pdb") and pdb_file.startswith("ranked_"):
                        rank_number = pdb_file.split("_")[1].split(".")[0]
                        new_file_name = f"seq_{seq_n}_{pdb_name}_{rank_number}.pdb"
                        source_file = os.path.join(folder_path, pdb_file)
                        destination_file = os.path.join(output_path, new_file_name)
                        shutil.copyfile(source_file, destination_file)
                        print(f"Copied {source_file} to {destination_file}")

def prepare_files_af3(input_path, output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    tmp_out = os.path.join(output_path, 'tmp')
    for directory in os.listdir(input_path):
        directory_path = os.path.join(input_path, directory)
        if os.path.isdir(directory_path):
            print(f"Processing directory: {directory_path}")
            copy_and_rename_cif_files(directory_path, tmp_out)
            convert_cif_to_pdb(tmp_out, output_path)
            for tmp_file in os.listdir(tmp_out):
                tmp_file_path = os.path.join(tmp_out, tmp_file)
                os.remove(tmp_file_path)
            os.rmdir(tmp_out)
        else:
            print(f"Skipping non-directory: {directory_path}")

def main():
    parser = argparse.ArgumentParser(description='Prepare PDB files for merging.')
    parser.add_argument('-i', '--input', required=True, help=' Path to the input directory containing alphafold prediction folders.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output directory where the prepared files will be saved.')
    parser.add_argument("-af3", "--af3", action='store_true', help='Add flag if files are coming from af3.')
    args = parser.parse_args()
    if args.af3:
        prepare_files_af3(args.input, args.output)
    else:
        prepare_files_af2(args.input, args.output)

if __name__ == "__main__":
    main()
