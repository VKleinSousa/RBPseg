import os
import shutil
import re
import argparse

def copy_and_rename_cif_files(source_dir, output_dir):
    """
    Copies .cif files from subdirectories and renames them using:
    {normalized_seq_n}_{pdb_name}_{rank_number}.cif

    Ensures the lowest sequence index starts at 0, then increments by 1.
    """
    os.makedirs(output_dir, exist_ok=True)

    parent_dir = os.path.dirname(os.path.normpath(source_dir))

    # --- collect all seq_n values from sibling folders ---
    seq_nums = []
    for folder in os.listdir(parent_dir):
        match = re.match(r"(.+)_([0-9]+)$", folder)
        if match:
            seq_nums.append(int(match.group(2)))

    # build normalization map (lowest becomes 0, next 1, etc.)
    seq_nums_sorted = sorted(set(seq_nums))
    seq_map = {old: new for new, old in enumerate(seq_nums_sorted)}

    # --- parse current folder name ---
    folder_name = os.path.basename(os.path.normpath(source_dir))
    match = re.match(r"(.+)_([0-9]+)$", folder_name)
    if not match:
        raise ValueError(f"Folder name '{folder_name}' does not match expected pattern '<pdb_name>_<seq_n>'")

    pdb_name = match.group(1)
    seq_n = int(match.group(2))
    normalized_seq_n = seq_map[seq_n]

    # --- process CIF files ---
    cif_files = sorted([f for f in os.listdir(source_dir) if f.endswith(".cif")])
    for new_rank, file in enumerate(cif_files):
        src_path = os.path.join(source_dir, file)
        new_file_name = f"seq_{normalized_seq_n}_{pdb_name}_{new_rank}.cif"
        dst_path = os.path.join(output_dir, new_file_name)

        shutil.copy2(src_path, dst_path)
        print(f"Copied and renamed: {src_path} -> {dst_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Copy and rename .cif files, normalizing sequence indices.")
    parser.add_argument("-s", "--source", required=True, help="Path to the source directory.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output directory.")
    args = parser.parse_args()

    copy_and_rename_cif_files(args.source, args.output)

