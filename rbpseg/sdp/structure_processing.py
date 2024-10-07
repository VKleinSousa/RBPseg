import numpy as np
from Bio.PDB import PDBParser, PDBIO
import logging

logger = logging.getLogger(__name__)

def load_structure(pdb_file, only_one_chain):
    logger.info(f"Loading structure from {pdb_file}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    if only_one_chain:
        # Extract chains and sort them alphabetically by their chain IDs
        chains = sorted(structure[0].get_chains(), key=lambda chain: chain.id)
        
        # Select the first chain in alphabetical order
        chain_id = chains[0].id if chains else None
        
        if chain_id:
            structure = extract_chain(structure, chain_id)
        else:
            logger.error("No chains found in the structure.")
            return None

    return structure
def preprocess_structure(structure, experimental_structure, normalize_plddt):
    if experimental_structure:
        structure = make_bfactor(structure)
    if normalize_plddt:
        structure = normalize_bfactor(structure)
    return structure

def normalize_bfactor(structure):
    try:
        max_bfactor = max(atom.get_bfactor() for atom in structure.get_atoms())
        min_bfactor = min(atom.get_bfactor() for atom in structure.get_atoms())

        for atom in structure.get_atoms():
            old_bfactor = atom.get_bfactor()
            normalized_bfactor = 100 * (old_bfactor - min_bfactor) / (max_bfactor - min_bfactor)
            atom.set_bfactor(normalized_bfactor)
        logger.info("B-factors normalized")
        return structure
    except Exception as e:
        logger.error(f"Error normalizing B-factors: {e}", exc_info=True)
        raise
        
def make_bfactor(structure):
    try:
        for atom in structure.get_atoms():
            atom.set_bfactor(100)
        logger.info("B-factors set")
        return structure
    except Exception as e:
        logger.error(f"Error setting B-factors to experimental structure: {e}", exc_info=True)
        raise
        
def extract_chain(structure, chain_id):
    from Bio import PDB
    try:
        extracted_structure = PDB.Structure.Structure('extracted')
        for model in structure:
            for chain in model:
                if chain.id == chain_id:
                    extracted_chain = chain.copy()
                    extracted_model = PDB.Model.Model(model.id)
                    extracted_model.add(extracted_chain)
                    extracted_structure.add(extracted_model)

        logger.info(f"Extracted chain {chain_id}")
        return extracted_structure
    except Exception as e:
        logger.error(f"Error extracting chain {chain_id}: {e}", exc_info=True)
        raise
