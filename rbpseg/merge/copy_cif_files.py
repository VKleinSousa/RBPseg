import os
import shutil
import argparse

def copy_and_rename_cif_files(source_dir, output_dir):
    """
    Copies .cif files from 'seed' subdirectories in source_dir to output_dir, 
    renaming them based on the parent directory structure.

    Args:
    - source_dir (str): The path to the source directory.
    - output_dir (str): The path to the output directory.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Get the base name of the source directory
    base_name = os.path.basename(os.path.normpath(source_dir))
    
    # Walk through the source directory
    for root, dirs, files in os.walk(source_dir):
        # Check if the current folder is a 'seed' folder
        if "seed" in os.path.basename(root):
            for file in files:
                if file.endswith(".cif"):
                    # Get the relative path from source_dir to the current directory
                    relative_path = os.path.relpath(root, source_dir)
                    # Build the new filename
                    new_name = f"{base_name}_seed_{relative_path.replace(os.sep, '_')}.cif"
                    # Source file path
                    src_file = os.path.join(root, file)
                    # Destination file path
                    dst_file = os.path.join(output_dir, new_name)
                    # Copy and rename the file
                    shutil.copy2(src_file, dst_file)
                    print(f"Copied and renamed: {src_file} -> {dst_file}")

if __name__ == "__main__":
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Copy and rename .cif files from seed subdirectories.")
    parser.add_argument("-s", "--source", required=True, help="Path to the source directory.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory.")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Call the function with provided arguments
    copy_and_rename_cif_files(args.source, args.output)

