import numpy as np
from Bio.PDB import PDBParser, PDBIO
import logging

logger = logging.getLogger(__name__)

def load_structure(pdb_file, only_one_chain):
    logger.info(f"Loading structure from {pdb_file}")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    if only_one_chain:
        chain_id = 'A'  # Modify as needed
        structure = extract_chain(structure, chain_id)
    return structure

def preprocess_structure(structure, normalize_plddt):
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
