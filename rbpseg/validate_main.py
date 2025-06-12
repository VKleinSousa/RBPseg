import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, is_aa
from scipy.spatial.distance import cdist
import argparse


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def pclash(n_clashes, plddt, weight):
    return plddt / (1 + (weight * n_clashes))


def calculate_distances(res1, res2):
    return np.linalg.norm(res1 - res2)


def calculate_merge_score(n, distance, a=1, b=1):
    return a * (n / (n + 1)) + b * abs(distance - 3.8)


def calculate_pclash(structure, clash_threshold, weight, update_bfactor):
    model = structure[0]  # Assuming the first model
    results = {}
    chain_boundaries = []
    start_residue_index = 0

    chain_ca_atoms = {}  # Store C-alpha atoms for each chain

    # Collect C-alpha atoms for all chains
    for chain in model:
        ca_atoms = []
        for residue in chain:
            if is_aa(residue) and 'CA' in residue:
                ca_atoms.append(residue['CA'])
        if ca_atoms:
            chain_ca_atoms[chain.id] = ca_atoms

    chain_results = {}  # Ensure consistent use of dictionary

    # Iterate over chains to calculate clashes
    for chain_id, ca_atoms in chain_ca_atoms.items():
        chain_distances = cdist(
            np.array([atom.get_coord() for atom in ca_atoms]),
            np.array([atom.get_coord() for atom in ca_atoms])
        )

        chain_results[chain_id] = []  # Correctly use as dictionary
        for i, atom in enumerate(ca_atoms):
            residue = atom.get_parent()
            intra_clashes = np.sum((chain_distances[i] < clash_threshold) & (chain_distances[i] > 0))
            inter_clashes = 0

            # Calculate inter-chain clashes
            for other_chain_id, other_ca_atoms in chain_ca_atoms.items():
                if chain_id != other_chain_id:
                    inter_distances = cdist(
                        np.array([atom.get_coord()]),
                        np.array([other_atom.get_coord() for other_atom in other_ca_atoms])
                    )
                    inter_clashes += np.sum(inter_distances[0] < clash_threshold)

            # Calculate peptide bond distance
            if i < len(ca_atoms) - 1:
                next_atom = ca_atoms[i + 1]
                peptide_distance = calculate_distances(atom.get_coord(), next_atom.get_coord())
            else:
                peptide_distance = 0.0

            # Calculate merge score
            total_clashes = intra_clashes + inter_clashes
            merge_score = calculate_merge_score(total_clashes, peptide_distance)

            # Append results
            chain_results[chain_id].append(
                (
                    start_residue_index + i,
                    atom.get_bfactor(),  # pLDDT
                    pclash(total_clashes, atom.get_bfactor(), weight),
                    intra_clashes,
                    inter_clashes,
                    abs(peptide_distance - 3.8),
                    merge_score,
                )
            )

            if update_bfactor:
                atom.set_bfactor(merge_score)

        results[chain_id] = chain_results[chain_id]
        chain_boundaries.append(start_residue_index)
        start_residue_index += len(ca_atoms) + 10  # Add padding for visualization

    # Aggregated scores
    agg_merge_score = 0
    plddt_list = []
    for chain_id, chain_data in results.items():
        _, pLDDTs, _, _, _, _, merge_scores = zip(*chain_data)
        plddt_list.extend(pLDDTs)
        agg_merge_score += np.sum(np.array(pLDDTs) * np.array(merge_scores))

    # Normalize aggregated residuals
    total_residues = len(plddt_list)
    mean_plddt = np.mean(plddt_list)
    normalized_aggregated_residuals = agg_merge_score / (total_residues * mean_plddt)

    return structure, results, chain_boundaries, normalized_aggregated_residuals


def plot_pclash(results, chain_boundaries, save_path="pclash_plot.png"):
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), gridspec_kw={'height_ratios': [3, 1]})
    scatter_ax, bar_ax = axes

    x_offset = 0
    x_ticks = []

    for chain_id, chain_results in results.items():
        residue_indices, pLDDTs, pClashes, intra_clashes, inter_clashes, distances, merge_scores = zip(*[
            (res[0], res[1], res[2], res[3], res[4], res[5], res[6]) for res in chain_results
        ])
        x_values = np.arange(x_offset, x_offset + len(residue_indices))
        intra_diff = np.array(intra_clashes) > 0
        inter_diff = np.array(inter_clashes) > 0
        # Scatter plot for pLDDT and merge scores
        scatter_ax.scatter(x_values, pLDDTs, color='black', label='pLDDT' if chain_id == list(results.keys())[0] else "", s=10)
        scatter_ax.scatter(x_values, pClashes, color='skyblue', label='pClash' if chain_id == list(results.keys())[0] else "", s=10)
        scatter_ax.scatter(np.array(x_values)[intra_diff], np.array(pClashes)[intra_diff], color='orange', label='Intra-chain Clash' if chain_id == list(results.keys())[0] else "", s=15)
        scatter_ax.scatter(np.array(x_values)[inter_diff], np.array(pClashes)[inter_diff], color='red', label='Inter-chain Clash' if chain_id == list(results.keys())[0] else "", s=15)
        
        # Bar plot for merge scores
        bar_ax.bar(x_values, merge_scores, color='black', width=3.0, label="Merge Score" if chain_id == list(results.keys())[0] else "")

        # Annotate bars with residue numbers where merge score > 2
        for idx, (x, merge_score, residue_index) in enumerate(zip(x_values, merge_scores, residue_indices)):
            if merge_score > 2:  # Only annotate if merge score > 2
                bar_ax.text(x, merge_score + 0.1, str(residue_index), ha='center', va='bottom', fontsize=8, color='red')

        # Update x_offset and x_ticks for the next chain
        x_ticks.append(x_offset + len(residue_indices) // 2)
        x_offset += len(residue_indices) + 10  # Add padding between chains

    # Scatter plot formatting
    scatter_ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    scatter_ax.set_ylabel("Score", fontsize=12)
    scatter_ax.set_title("pLDDT and pClash Scores with Intra/Inter-chain Mismatches", fontsize=14, weight='bold')
    scatter_ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=10, frameon=False)

    # Bar chart formatting
    bar_ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    bar_ax.set_xlabel("Residue Index", fontsize=12)
    bar_ax.set_ylabel("Merge Score", fontsize=12)
    bar_ax.set_title("Merge Score by Residue", fontsize=14, weight='bold')

    # Shared x-axis ticks
    chain_labels = [f"Chain {chr(65+i)}" for i in range(len(chain_boundaries))]
    scatter_ax.set_xticks(x_ticks)
    scatter_ax.set_xticklabels(chain_labels, fontsize=10, rotation=0)
    bar_ax.set_xticks(x_ticks)
    bar_ax.set_xticklabels(chain_labels, fontsize=10, rotation=0)

    # Tight layout and save
    plt.tight_layout(rect=[0, 0.1, 1, 0.95])  # Adjust to make space for titles and legend
    plt.savefig(save_path, dpi=300)
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Calculate pClash and merge scores for a protein structure.")
    parser.add_argument("-p", "--pdb_file", required=True, help="Path to the PDB file.")
    parser.add_argument("-c", "--clash_threshold", type=float, default=2.0, help="Distance threshold for detecting clashes.")
    parser.add_argument("-w", "--weight", type=float, default=1.0, help="Weight factor for pClash calculation.")
    parser.add_argument("-u", "--update_bfactor", action="store_true", help="Update B-factor with merge scores.")
    parser.add_argument("-plot", "--plot", action='store_true', help="Generate plots.")
    parser.add_argument("-o", "--output_plot", default="pclash_plot.png", help="Filename for the output plot.")

    args = parser.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", args.pdb_file)

    structure, results, chain_boundaries, normalized_aggregated_residuals = calculate_pclash(
        structure, args.clash_threshold, args.weight, args.update_bfactor
    )

    print(f"Normalized Aggregated Residuals: {normalized_aggregated_residuals}")

    if args.plot:
        plot_pclash(results, chain_boundaries, args.output_plot)
        print(f"Plot saved as {args.output_plot}")


if __name__ == "__main__":
    main()

