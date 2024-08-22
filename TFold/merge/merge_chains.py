from Bio.PDB import PDBParser
import string

def merge_chains(moving_structure, chain_pairs):
    """
    Merge the chains in the moving structure based on the provided chain pairs.
    """
    number_chains = []  # Initialize an empty list

    for model in moving_structure:
        for chain in model:
            number_chains.append(chain)

    nmodels = len(number_chains) // 3
    alphabet = list(string.ascii_uppercase)[1:26]  # A-Z excluding the first letter

    for chain_pair in chain_pairs:
        chains_ = [None] * nmodels  # Initialize a list to store the chains
        for i, chain_id in enumerate(chain_pair):
            for model in moving_structure:
                for chain in model:
                    if chain.id == chain_id:
                        chains_[i] = chain
        
        chains_ = sorted(chains_, key=lambda x: x.id)
        
        if all(chain is not None for chain in chains_):
            len_chain_a = int(chains_[0].get_unpacked_list()[-1].id[1]) + 1
            k = 0  # Variable to manage the merging of residues

            for chain_obj in chains_[1:]:  # Loop over chain objects to merge
                k += 1
                for l, residue in enumerate(chain_obj.get_residues()):
                    old_id2 = list(residue.id)
                    old_id2[1] = len_chain_a + l + k * 1000
                    residue.id = tuple(old_id2)
                    chains_[0].add(residue)

            for j in range(1, nmodels):
                len_chain_a = len_chain_a + int(chains_[j].get_unpacked_list()[-1].id[1]) + 1 

            new_resnums1 = [i + 1 for i in range(len_chain_a)]

            for model in moving_structure:
                for chain in model:
                    for l, residue in enumerate(chains_[0].get_residues()):
                        res_id = list(residue.id)
                        res_id[1] = new_resnums1[l]
                        residue.id = tuple(res_id)

    # Delete unwanted chains (if needed)
    for i in alphabet[3:26]:      
        for model in moving_structure:
            for chain in model:
                if chain.id == i:
                    model.detach_child(chain.id)

if __name__ == "__main__":
    # Example usage
    parser = PDBParser(QUIET=True)
    moving_structure = parser.get_structure("moving_structure", "moving.pdb")

    chain_pairs = [('E', 'B'), ('F', 'B'), ('G', 'B')]  # Example chain pairs to merge

    merge_chains(moving_structure, chain_pairs)
    print("Chains merged successfully.")
