import pandas as pd
import numpy as np
from sklearn.cluster import SpectralClustering
from sklearn.metrics import silhouette_score
from scipy.optimize import curve_fit
from sklearn.preprocessing import normalize
import argparse


def prepare_diagonal_matrix(mat, mode, max_classes):
    print("Preparing diagonal matrix...")
    pdb_codes = np.unique(np.concatenate([mat['PDBchain1'], mat['PDBchain2']]))
   
    if len(pdb_codes) < max_classes:
        print(f"Found {len(pdb_codes)} unique PDB codes. This is less than preset max_classes. max_classes readjusted to {len(pdb_codes)-1}.")
        max_classes = len(pdb_codes)-1
    else:
         print(f"Found {len(pdb_codes)} unique PDB codes.")
    mat_transformed = pd.DataFrame(0.0, index=pdb_codes, columns=pdb_codes)
    for i in range(len(mat)):
        if mode == 0:
            # take the mean TM score
            score = np.mean((mat['TM1'][i], mat['TM2'][i]))
        elif mode == 1:
            # take the max TM score
            score = max((mat['TM1'][i], mat['TM2'][i]))
        elif mode == 2:
            # take the min TM score
            score = min((mat['TM1'][i], mat['TM2'][i]))
        else:
            print('Error: mode invalid. Valid modes: 0, 1, 2.')
            return None  # Return None on error
        mat_transformed.loc[mat['PDBchain1'][i], mat['PDBchain2'][i]] = score
        mat_transformed.loc[mat['PDBchain2'][i], mat['PDBchain1'][i]] = score
    np.fill_diagonal(mat_transformed.values, 1)
    print("Diagonal matrix preparation complete.")
    return mat_transformed, max_classes


def calculate_cluster_means(identity_matrix, clusters, row_names):
    print("Calculating cluster means...")
    cluster_means = {}
    unique_clusters = set(clusters)
    
    mean_no_d = []
    for cluster_label in unique_clusters:
        # Get indices of members belonging to the current cluster
        cluster_indices = [i for i, label in enumerate(clusters) if label == cluster_label]
        
        # Filter rows and columns corresponding to the current cluster
        cluster_matrix = identity_matrix.iloc[cluster_indices, cluster_indices]
        
        # Calculate m_value excluding diagonal elements
        if len(cluster_indices) > 1:
            m_value = (cluster_matrix.sum(axis=1) - np.diag(cluster_matrix)) / (cluster_matrix.shape[1] - 1)
        else:
            m_value = np.array([1])  # Handle single-member clusters
            
        mean_no_d.append(list(m_value))
        mean_value = np.mean(m_value)  # Mean of all m_value entries
        
        if len(cluster_indices) > 1:
            std_value = np.std(m_value) / np.sqrt(len(m_value) - 1)  # Corrected standard error calculation
        else:
            std_value = 0  # No variance for single-member clusters
        
        # Store mean value along with cluster label
        cluster_means[cluster_label] = {
            'mean_value': np.round(mean_value, 2),
            'std_value': np.round(std_value, 2),
            'cluster_size': len(cluster_indices)
        }
    
    print("Cluster means calculation complete.")
    return cluster_means, mean_no_d


def calculate_sm(identity_matrix, max_clusters, min_clusters=3, n_rounds=3):
    print(f"Calculating spectral clustering metrics for {min_clusters} to {max_clusters} clusters...")
    results_df = []

    for i in range(min_clusters, max_clusters+1):
        print(f"Calculating for {i} clusters...")
        TM = []
        Sil = []
        sTM = []
        cluster_size = []
        cluster_max = []
        
        for j in range(n_rounds):
            #print(f"  Round {j+1}/{n_rounds}...")
            spectral_clustering = SpectralClustering(n_clusters=i, affinity='precomputed')
            clusters = spectral_clustering.fit_predict(identity_matrix)
            row_names = identity_matrix.index
            cluster_means, _ = calculate_cluster_means(identity_matrix, clusters, row_names)
            mean_values = [values['mean_value'] for values in cluster_means.values()]
            std_values = [values['std_value'] for values in cluster_means.values()]
            cluster_values = [values['cluster_size'] for values in cluster_means.values()]
            
            TM.append(np.mean(mean_values))
            sTM.append(np.mean(std_values))
            cluster_size.append(np.mean(cluster_values))
            cluster_max.append(max(cluster_values))
            Sil.append(silhouette_score(identity_matrix, clusters))
        
        results_df.append({
            'cluster_number': i,
            'TM': np.mean(TM),
            'TM_min': min(mean_values),
            'TM_max': max(mean_values),
            'sTM': np.mean(sTM),
            'cluster_size': np.mean(cluster_size),
            'cluster_max_size': np.mean(cluster_max),
            'Silhouette': np.mean(Sil)
        })
    
    print("Spectral clustering metrics calculation complete.")
    results_df = pd.DataFrame(results_df)
    
    metric = 1 / (1 - np.exp(-(results_df['TM_min'] + results_df['Silhouette']) / (results_df['cluster_max_size'] * results_df['sTM'])))
    results_df['SM'] = metric
    return results_df

def exponential_decay(x, a, b, c):
    return a * np.exp(-b * x)  + c 

def exponential_decay_derivative(x, a, b, c):
    return -a * b * np.exp(-b * x)

def exponential_decay_double_derivative(x, a, b, c):
    return a * b*b * np.exp(-b * x)

def fit_psm(results_df):
    print("Fitting exponential decay model...")
    #metric_array = metric.values
    #metric_2d = metric_array.reshape(-1, 1)
    #nmetric = normalize(metric_2d, axis=0).flatten()
    try:
        params, covariance = curve_fit(exponential_decay, results_df['cluster_number'], results_df['SM'])
        a, b, c = params
        
        fitted_values = exponential_decay(results_df['cluster_number'], a, b, c)
        residuals = results_df['SM'] - fitted_values
        SS_res = np.sum(residuals**2)
        mean_metric = np.mean(results_df['SM'])
        SS_tot = np.sum((results_df['SM'] - mean_metric)**2)
        R_squared = 1 - (SS_res / SS_tot)
        chi_squared = np.sum((residuals / np.std(results_df['SM']))**2)
    
        print(f'Exponential Decay fitted with chi_squared = {chi_squared} and R_squared = {R_squared}')
        return a, b, c
    except RuntimeError as e:
            print(f"Error during curve fitting: {e} /n error often happens when max_classes was too small.")
            return None, None, None


def calculate_optimal_epsilon(max_clusters, a, b, c, percentage=0.05, adaptive=True):
    print("Estimating optimal epsilon value...")
    second_derivatives = [abs(exponential_decay_double_derivative(i, a, b, c)) for i in range(2, max_clusters)]
    median_second_derivative = np.mean(second_derivatives)
    
    # Adaptive adjustment of percentage based on the scale of second derivatives
    if adaptive:
        if median_second_derivative < 1e-4:
            percentage *= 2  # Increase the percentage if derivatives are very small
        elif median_second_derivative > 1e-1:
            percentage /= 2  # Decrease the percentage if derivatives are large

    optimal_epsilon = percentage * median_second_derivative
    print(f"Estimated optimal epsilon: {optimal_epsilon}")
    return optimal_epsilon

def calculate_optimal_classes(max_clusters, a, b, c, epsilon=None, retry_limit=3):
    if epsilon is None:
        epsilon = calculate_optimal_epsilon(max_clusters, a, b, c)
    
    print("Calculating optimal number of clusters...")
    cluster_number = None
    retry_count = 0

    while retry_count < retry_limit:
        for i in range(2, max_clusters):
            ep = abs(exponential_decay_double_derivative(i, a, b, c))
            pSM = abs(exponential_decay(i, a, b, c))
            print(f'{i} classes: pSM={pSM}, deceleration={ep}')
            if ep < epsilon:
                print(f'Optimal number of classes is {i} with a pSM of {pSM}, and a deceleration of {ep}')
                cluster_number = i
                break
        
        if cluster_number is not None:
            break  # Exit the loop if an optimal cluster number is found
        
        # If no cluster number is found, retry with a larger epsilon
        retry_count += 1
        epsilon *= 10
        print(f"Retrying with a larger epsilon: {epsilon}")

    # If no optimal cluster number was found after retries, return the maximum number of clusters
    if cluster_number is None:
        cluster_number = max_clusters
        print(f"No optimal number of clusters found. Returning maximum number of clusters: {max_clusters}")
    
    return cluster_number


def main():
    # Create an ArgumentParser object to handle command-line arguments
    parser = argparse.ArgumentParser(description='Perform spectral clustering on TM score matrix.')

    # Add command-line flags for TM_file, max_classes, mode, epsilon, output, and dissimilarity
    parser.add_argument('-f', '--TM_file', required=True, help='Path to the TM score file')
    parser.add_argument('-c', '--max_classes', type=int, default=100, help='Max number of classes (Default: 100)')
    parser.add_argument('-m', '--mode', type=int, default=0, help='Modes: 0 (Default) takes the mean TM of both chains, 1 takes the maximum TM, 2 takes the minimum')
    parser.add_argument('-e', '--epsilon', type=float, default=None, help='Threshold for the second derivative in exponential decay to calculate the number of classes. If not provided, it will be estimated automatically.')
    parser.add_argument('-o', '--output', type=str, default='output', help='Output file path (default: output)')
    parser.add_argument('-d', '--dissimilarity', type=bool, default=False, help='If clustering based on dissimilarity')
    
    # Parse the command-line arguments
    args = parser.parse_args()

    # Retrieve values from the parsed arguments
    TM_path = args.TM_file
    max_classes = args.max_classes
    mode = args.mode
    epsilon = args.epsilon
    output = args.output
    dissimilarity = args.dissimilarity

    print(f"Reading TM score file: {TM_path}")
    mat = pd.read_csv(TM_path, sep="\t")

    # Prepare the identity matrix and adjust max_classes if needed
    identity_matrix, max_classes = prepare_diagonal_matrix(mat, mode, max_classes)
    
    if dissimilarity:
        print("Converting identity matrix to dissimilarity matrix...")
        identity_matrix = 1 - identity_matrix
    
    sm_df = calculate_sm(identity_matrix, max_classes)
   
    a, b, c = fit_psm(sm_df)
    # Estimate epsilon if not provided
    if epsilon is None:
        epsilon = calculate_optimal_epsilon(max_classes, a, b, c)
    cluster_number = calculate_optimal_classes(max_classes, a, b, c, epsilon)
    
    if cluster_number is not None:
        print(f"Performing spectral clustering with {cluster_number} clusters...")
        spectral_clustering = SpectralClustering(n_clusters=cluster_number, affinity='precomputed')
        clusters = spectral_clustering.fit_predict(identity_matrix)
        row_names = identity_matrix.index
        cluster_means, classes_mean_no_diagonal = calculate_cluster_means(identity_matrix, clusters, row_names)
        
        df_classes = pd.DataFrame({'Protein_name': row_names, 'Classes': clusters})
        df_classes.to_csv(f'{output}.classes.csv', index=False)
        print(f"Saved cluster classifications to {output}.classes.csv")
        
        means_df = pd.DataFrame(classes_mean_no_diagonal)
        means_df.to_csv(f'{output}.means.csv', index=False)
        print(f"Saved mean values to {output}.means.csv")
    else:
        print("Failed to determine an optimal number of clusters.")
    

if __name__ == "__main__":
    main()
