import numpy as np
import pandas as pd
import random

def find_last_atom_position(fixed_structure, moving_structure):
    """
    Find the positions of the last atom in the fixed structure and the first atom in the moving structure.
    """
    last_atom_positions_fixed = []
    first_atom_positions_moving = []

    for model in fixed_structure:
        for chain in model:
            last_residue = list(chain.get_residues())[-1]
            last_atom = list(last_residue.get_atoms())[-1]
            last_atom_position = last_atom.get_coord()
            last_atom_positions_fixed.append(last_atom_position)

    for model in moving_structure:
        for chain in model:
            first_residue = list(chain.get_residues())[0]
            first_atom = list(first_residue.get_atoms())[0]
            first_atom_position = first_atom.get_coord()
            first_atom_positions_moving.append(first_atom_position)

    return np.array(last_atom_positions_fixed), np.array(first_atom_positions_moving)

def euclidean_distance(atom1_coord, atom2_coord):
    """
    Calculate the Euclidean distance between two atoms.
    """
    diff = atom1_coord - atom2_coord
    distance = np.sqrt(np.sum(diff ** 2))
    return distance

def find_triangle_center(vertices):
    """
    Calculate the centroid of the triangle.
    """
    center = np.mean(vertices, axis=0)
    return center

def calculate_triangle_height(vertices):
    """
    Calculate the height of a triangle given its vertices.
    """
    v1, v2, v3 = vertices
    base_length = np.linalg.norm(v2 - v1)
    hypotenuse = np.linalg.norm(v3 - v1)
    height = np.sqrt(hypotenuse ** 2 - (0.5 * base_length) ** 2)
    return height

def create_sphere(center, radius):
    """
    Create a sphere centered at 'center' with a radius 'radius'.
    """
    return {'center': center, 'radius': radius}

def check_line_sphere_intersection(vertex, fourth_point, sphere):
    """
    Check if a line (from vertex to fourth_point) intersects with a sphere.
    """
    line_direction = fourth_point - vertex
    vertex_to_center = sphere['center'] - vertex
    line_length = np.linalg.norm(line_direction)
    line_direction /= line_length

    t = np.dot(vertex_to_center, line_direction)
    closest_point = vertex + t * line_direction
    distance_to_center = np.linalg.norm(closest_point - sphere['center'])

    return distance_to_center < sphere['radius'], line_length

def rename_chains_2(fixed_structure, moving_structure, interface_rmsd, change_bfactor):
    """
    Rename chains in the moving structure based on the best matching pair with the fixed structure.
    """
    pairs = [('E', 'B'), ('F', 'B'), ('G', 'B'), ('E', 'C'), ('F', 'C'), ('G', 'C'), ('E', 'D'), ('F', 'D'), ('G', 'D')]
    index = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (1, 2), (0, 2), (1, 2), (2, 2)]
    fixed_index = [0, 0, 0, 1, 1, 1, 2, 2, 2]
    moving_index = [0, 1, 2, 0, 1, 2, 0, 1, 2]
    # Ensure interface_rmsd is a 1D array
    if isinstance(interface_rmsd, pd.DataFrame):
        # If it's a DataFrame with a single column
        interface_rmsd = interface_rmsd.squeeze()  # Convert to Series
    if isinstance(interface_rmsd, pd.Series):
        # Convert to NumPy array and flatten
        interface_rmsd = interface_rmsd.values.flatten()
    print(interface_rmsd)
    df = pd.DataFrame({
        'rms': interface_rmsd,
        'pair': pairs,
        'index': index,
        'moving_index': moving_index,
        'fixed_index': fixed_index
    })

    last_fixed, chains_fixed_positions = [], []
    first_moving, chains_moving_positions = [], []

    for model in fixed_structure:
        for chain in model:
            residues = list(chain.get_residues())
            last_residue = residues[-1]
            last_fixed.extend(last_residue.get_atoms())
            if change_bfactor:
                for atom in last_residue.get_atoms():
                    atom.set_bfactor(-1)
            pos_fixed = sum(atom.get_coord() for atom in last_fixed) / len(last_fixed)
            chains_fixed_positions.append(pos_fixed)
            last_fixed = []

    for model in moving_structure:
        for chain in model:
            residues = list(chain.get_residues())
            first_residue = residues[0]
            first_moving.extend(first_residue.get_atoms())
            if change_bfactor:
                for atom in first_residue.get_atoms():
                    atom.set_bfactor(-1)
            pos_moving = sum(atom.get_coord() for atom in first_moving) / len(first_moving)
            chains_moving_positions.append(pos_moving)
            first_moving = []

    vertices = np.array(chains_fixed_positions)
    vertices_up = np.array(chains_moving_positions)

    triangle_center = find_triangle_center(vertices_up)
    triangle_height = calculate_triangle_height(vertices_up)
    sphere_radius = triangle_height / 3.0
    sphere = create_sphere(triangle_center, sphere_radius)

    intersection_list = []
    connection_length_list = []
    for k in range(len(df['rms'])):
        i, j = df['index'].iloc[k]
        vertex = vertices[j]
        fourth_point = vertices_up[i]
        intersection, connection_length = check_line_sphere_intersection(fourth_point, vertex, sphere)
        intersection_list.append(intersection)
        connection_length_list.append(connection_length)

    df['connection_length'] = connection_length_list
    df['intersection'] = intersection_list
    df = df.sort_values('connection_length')

    chain_pairs = []
    select_pairs = [('E', 'B'), ('F', 'B'), ('G', 'B'), ('E', 'C'), ('F', 'C'), ('G', 'C'), ('E', 'D'), ('F', 'D'), ('G', 'D')]
    excluded_residues = set()
    k_2 = 0

    while k_2 < 2:
        filtered_df = df[df['pair'].apply(lambda x: len(set(x)) == len(x))]
        filtered_df = filtered_df[filtered_df['intersection'] == False]
        filtered_df['rms+connection_length'] = filtered_df['rms'] + filtered_df['connection_length']
        sorted_df = filtered_df.sort_values(by='rms+connection_length')

        if sorted_df.empty:
            lowest_rms_pair = random.choice(select_pairs)
            chain_pairs.append(lowest_rms_pair)
            excluded_residues.update(lowest_rms_pair)
            select_pairs = [pair for pair in select_pairs if not any(residue in pair for residue in excluded_residues)]
            if k_2 < 1:
                lowest_rms_pair = random.choice(select_pairs)
                chain_pairs.append(lowest_rms_pair)
                excluded_residues.update(lowest_rms_pair)
                select_pairs = [pair for pair in select_pairs if not any(residue in pair for residue in excluded_residues)]
                break
            else:
                break
        else:
            lowest_rms_pair = sorted_df.iloc[0]['pair']
            chain_pairs.append(lowest_rms_pair)
            k_2 += 1
            excluded_residues.update(lowest_rms_pair)

        select_pairs = [pair for pair in select_pairs if not any(residue in pair for residue in excluded_residues)]
        df = df[~df['pair'].apply(lambda x: any(residue in x for residue in excluded_residues))]

    if len(chain_pairs) < 3:
        chain_X = set(letter for pair in pairs for letter in pair)
        chain_Y = set(letter for pair in chain_pairs for letter in pair)
        missing_chains = sorted(list(chain_X - chain_Y), key=str.lower, reverse=True)
        chain_pairs = chain_pairs + [(missing_chains[0], missing_chains[1])]

    chain_pairs = sorted(chain_pairs, key=lambda x: x[0])
    
    for model in moving_structure:
        for i, chain in enumerate(model):
            if i < len(chain_pairs):
                chain.id = chain_pairs[i][0]

    for model in fixed_structure:
        for chain in model:
            moving_structure[0].add(chain.copy())

    return moving_structure, chain_pairs


if __name__ == "__main__":
    # Example usage
    parser = PDBParser(QUIET=True)
    fixed_structure = parser.get_structure("fixed_structure", "fixed.pdb")
    moving_structure = parser.get_structure("moving_structure", "moving.pdb")
    interface_rmsd = np.random.random(9)  # Example RMSD values
    change_bfactor = False  # Example b-factor modification flag

    moving_structure, chain_pairs = rename_chains_2(fixed_structure, moving_structure, interface_rmsd, change_bfactor)
    print("Chain pairs after renaming:", chain_pairs)
