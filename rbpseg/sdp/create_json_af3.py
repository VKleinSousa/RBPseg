#!/usr/bin/env python3
import os
import argparse
import json
from Bio import SeqIO
import string
import random

def create_json(fasta_file, output_file, rna_flag=None, dna_flag=None, ligand_flag=None, num_seeds=2):
    """
    Create a JSON file based on the provided protein FASTA file and optional flags.
    """
    sequences = []
    job_name = os.path.splitext(os.path.basename(fasta_file))[0]
    
    # Parse the protein FASTA file
    for idx, record in  enumerate(SeqIO.parse(fasta_file, "fasta")):
        protein_entry = {
            "protein": {
                "id": string.ascii_uppercase[idx],
                "sequence": str(record.seq)
            }
        }
        sequences.append(protein_entry)

    # Add RNA, DNA, and ligand entries if provided
    if rna_flag:
        sequences.append({"rna": {"id": "rna", "sequence": rna_flag}})
    if dna_flag:
        sequences.append({"dna": {"id": "dna", "sequence": dna_flag}})
    if ligand_flag:
        sequences.append({"ligand": {"id": "ligand", "name": ligand_flag}})
    
    # Generate random seeds
    model_seeds = [random.randint(0, 99999) for _ in range(num_seeds)]
    
    # Create the final JSON structure
    json_data = {
        "name": job_name,
        "modelSeeds": model_seeds,
        "sequences": sequences,
        "dialect": "alphafold3",  # Required
        "version": 1             # Required
    }
    
    # Write to output file
    with open(output_file, "w") as json_file:
        json.dump(json_data, json_file, indent=2)
    print(f"JSON file created at: {output_file}")


if __name__ == "__main__":
    # Define the argument parser
    parser = argparse.ArgumentParser(description="Create a JSON file for AlphaFold3 jobs.")
    parser.add_argument("-f", "--fasta", required=True, help="Protein FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Output JSON file.")
    parser.add_argument("--rna", help="Optional RNA sequence.")
    parser.add_argument("--dna", help="Optional DNA sequence.")
    parser.add_argument("--ligand", help="Optional ligand name.")
    parser.add_argument("--num-seeds", type=int, default=2, help="Number of random seeds to generate (default: 2).")
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Create the JSON file
    create_json(args.fasta, args.output, args.rna, args.dna, args.ligand, args.num_seeds)

