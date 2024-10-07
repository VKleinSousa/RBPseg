import numpy as np
from scipy.spatial.distance import cdist
import logging

logger = logging.getLogger(__name__)

def sdp(plddt1, plddt2, d, pair_distance_constant):
    if d == 0:
        return 1
    else:
        return 1 / (1 + np.exp(-1 * pair_distance_constant * (plddt1 + plddt2) / (d * d)))

def sdp_t(plddt1, plddt2, d, pair_distance_constant):
    if d == 0:
        return 1
    else:
        return 1 / (1 + np.exp(-1 * pair_distance_constant * (plddt1 * plddt2) / (d * d)))

def sdp_d(plddt1, plddt2, d, pair_distance_constant):
    if d == 0:
        return 1
    else:
        return 1 / (1 + np.exp(-1 * pair_distance_constant * (plddt1 + plddt2) / d))

def sdp_matrix_calculation(structure, pair_distance_constant, sdp_mode=0):
    try:
        # Filter out residues that don't have 'CA' atom
        residues = [residue for model in structure for chain in model for residue in chain if 'CA' in residue]
        
        num_residues = len(residues)
        sdp_matrix = np.zeros((num_residues, num_residues))

        # Extract coordinates and pLDDT values for all residues with 'CA' atoms
        coords = np.array([residue['CA'].get_coord() for residue in residues])
        plddts = np.array([residue['CA'].get_bfactor() for residue in residues])

        # Ensure coordinates array is 2D
        if coords.ndim != 2 or coords.shape[1] != 3:
            raise ValueError(f"Coords array must be 2D with shape (n_residues, 3), but got shape {coords.shape}")

        # Calculate the pairwise distances between residues
        distances = cdist(coords, coords)

        # Set a minimum threshold for distances to avoid division by zero
        epsilon = 1e-8
        distances = np.where(distances == 0, epsilon, distances)

        # Compute the SDP matrix based on the chosen mode
        if sdp_mode == 0:
            sdp_matrix = 1 / (1 + np.exp(-1 * pair_distance_constant * (plddts[:, None] + plddts[None, :]) / (distances ** 2)))
        elif sdp_mode == 1:
            sdp_matrix = 1 / (1 + np.exp(-1 * pair_distance_constant * (plddts[:, None] * plddts[None, :]) / (distances ** 2)))
        elif sdp_mode == 2:
            sdp_matrix = 1 / (1 + np.exp(-1 * pair_distance_constant * (plddts[:, None] + plddts[None, :]) / distances))

        # Set diagonal to 1
        np.fill_diagonal(sdp_matrix, 1)
        return sdp_matrix, np.mean(plddts), residues

    except KeyError as e:
        logger.error(f"Error accessing CA atom or B-factor for residues: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error during SDP matrix calculation: {e}", exc_info=True)
        raise


