import scanpy as sc

adata = sc.read('/path/to/your/data.h5ad')

import numpy as np

def find_border_spots(adata, cluster_label):
    cluster_indices = np.where(adata.obs['clusters'] == cluster_label)[0]
    neighbors_matrix = adata.obsp['connectivities']  # Adjusted to access from `obsp`

    # Find border spots: those that have at least one neighbor not in the cluster
    border_spots = []
    for spot in cluster_indices:
        neighbors = neighbors_matrix[spot].nonzero()[1]
        if any(adata.obs['clusters'][neighbor] != cluster_label for neighbor in neighbors):
            border_spots.append(spot)
    
    return border_spots


border_spots = find_border_spots(adata, '4')

print("Number of border spots identified:", len(border_spots))
print("Sample border spots:", border_spots[:5])

import numpy as np
from scipy.spatial.distance import cdist

def find_distance_based_outside_borders(adata, border_spots, cluster_label, threshold):
    """ Find outside border spots based on distance rather than connectivities,
        ensuring these spots are not part of the original cluster. """
    coords = adata.obsm['spatial']
    cluster_indices = adata.obs['clusters'] != cluster_label  # Spots not in the original cluster
    distances = cdist(coords, coords)
    
    # Spots within 'threshold' distance are considered adjacent
    is_adjacent = distances <= threshold
    
    outside_borders = []
    for i in range(distances.shape[0]):
        # Spot must not be a border spot, must not be in the original cluster, and must be adjacent to any border spot
        if i not in border_spots and cluster_indices[i] and any(is_adjacent[i, j] for j in border_spots if j != i):
            outside_borders.append(i)
            
    return outside_borders

# Usage example:
#border_spots = find_border_spots(adata, 'specific_cluster_label')
outside_border_spots = find_distance_based_outside_borders(adata, border_spots, '4', threshold=455.05)

# You might need to adjust the 'threshold' depending on the scale of your spatial coordinates.

import numpy as np
from scipy.spatial.distance import cdist

def find_third_outside_borders(adata, outside_border_spots, border_spots, cluster_label, threshold):
    """Find third-level outside border spots based on distance, ensuring these spots are not part of the original cluster, border spots, or outside border spots."""
    coords = adata.obsm['spatial']
    # Combine conditions for not being a border spot, an outside border spot, or in the original cluster
    not_allowed_indices = np.concatenate((border_spots, outside_border_spots, np.where(adata.obs['clusters'] == cluster_label)[0]))

    distances = cdist(coords, coords)
    # Spots within 'threshold' distance are considered adjacent
    is_adjacent = distances <= threshold
    
    third_outside_borders = []
    for i in range(distances.shape[0]):
        # Spot must not be a border spot, an outside border spot, or in the original cluster, and must be adjacent to any outside border spot
        if i not in not_allowed_indices and any(is_adjacent[i, j] for j in outside_border_spots if j != i):
            third_outside_borders.append(i)

    return third_outside_borders

# Usage example:
third_outside_borders = find_third_outside_borders(adata, outside_border_spots, border_spots, '4', threshold=455.05)

# You might need to adjust the 'threshold' depending on the scale of your spatial coordinates.


#Use this script to find the threshold distance necessary
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

def plot_distance_distribution(adata):
    # Extract spatial coordinates
    coords = adata.obsm['spatial']
    
    # Calculate Euclidean distances between all pairs
    distances = cdist(coords, coords)
    
    # We replace zeros with np.inf to ignore self-comparisons in min calculation
    np.fill_diagonal(distances, np.inf)
    
    # Find the minimum distance to any other spot for each spot
    min_distances = np.min(distances, axis=1)

    # Plotting the distribution of minimum distances
    plt.figure(figsize=(10, 6))
    plt.hist(min_distances, bins=50, color='skyblue', alpha=0.7)
    plt.title('Distribution of Minimum Distances Between Spots')
    plt.xlabel('Distance')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

    # Print basic statistics to help decide on a threshold
    print("Mean minimum distance:", np.mean(min_distances))
    print("Median minimum distance:", np.median(min_distances))
    print("Minimum distance observed:", np.min(min_distances))

# Usage
plot_distance_distribution(adata)


# Reset all to 'other' first
adata.obs['group'] = 'other'

# Using iloc for direct integer location access
adata.obs.iloc[[int(spot) for spot in border_spots], adata.obs.columns.get_loc('group')] = 'border'
adata.obs.iloc[[int(spot) for spot in outside_border_spots], adata.obs.columns.get_loc('group')] = 'outside_border'
adata.obs.iloc[[int(spot) for spot in third_outside_borders], adata.obs.columns.get_loc('group')] = 'third_outside_border'

# Check the assignment
print(adata.obs['group'].value_counts())



import scanpy as sc

# Filter out 'other' from the analysis
adata_subset = adata[adata.obs['group'].isin(['4', 'border', 'outside_border'])].copy()

# Recompute the dendrogram for the current groups (if you are using it)
sc.tl.dendrogram(adata_subset, groupby='group')

# Run differential expression analysis
sc.tl.rank_genes_groups(adata_subset, groupby='group', groups=['4', 'border', 'outside_border'], method='wilcoxon')

# Visualizing the results with a heatmap
sc.pl.rank_genes_groups_heatmap(adata_subset, groups=['4', 'border', 'outside_border'], n_genes=20, show_gene_labels=True, cmap='viridis')
