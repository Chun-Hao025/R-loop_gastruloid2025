#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scanpy as sc
import pandas as pd
import scipy
import scipy.sparse as sp
import random
import networkx as nx
import scipy.stats as ss
from copy import deepcopy
from sklearn.preprocessing import normalize
from numpy.linalg import inv, pinv
from numpy.linalg.linalg import LinAlgError
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import find, csr_matrix, csc_matrix, lil_matrix, save_npz, load_npz
from scipy.stats import entropy
from scipy.sparse.csgraph import dijkstra
 



###github page:https://github.com/ShobiStassen/VIA/tree/master
from pyVIA.core import *
import pyVIA.datasets_via as datasets_via
#from core_working_ import *
import matplotlib.pyplot as plt
import pyVIA.core as via

import anndata
np.random.seed(4)
random.seed(4)

import gc
gc.collect()


cell_meta = pd.read_csv("/path/to/via/via_meta_new.csv", index_col=0)
harm=pd.read_csv("/path/to/via/via_harmony_new.csv", index_col=0)
umap=pd.read_csv("/path/to/via/via_umap_new.csv", index_col=0)

###get cluster results
###true_labels is not quit important in this step here since we just want to get cluster result through PARC method here.

true_labels_numeric= cell_meta['time_C'].tolist()
true_labels= cell_meta['ann_new_sub'].tolist()
harm=harm.to_numpy()
embedding=umap.to_numpy()
v1 = via.VIA(data=harm[:, 0:30], true_label=true_labels,
             cluster_graph_pruning=0.3,neighboring_terminal_states_threshold=3,small_pop=200,
             knn=40, too_big_factor=0.4,root_user=['Naive Pluripotent Cells'],dataset='group', random_seed=4, embedding = embedding,
            is_coarse=True, preserve_disconnected=True, pseudotime_threshold_TS=40, x_lazy=0.95,
         alpha_teleport=0.95, time_series=True, time_series_labels=true_labels_numeric, t_diff_step=3,
         edgebundle_pruning_twice=False, knn_sequential=30, knn_sequential_reverse=10,
             piegraph_arrow_head_width=0.15,piegraph_edgeweight_scalingfactor=2.5,edgebundle_pruning=0.3)

v1.run_VIA()
cell_meta['PARC_KNN40']=v1.labels
v1.to_csv("/path/to/via/via_meta_new1.csv", index=False) #Go back to seurat for annotation.

### Here, we aimed to test multiple parameter combinations for trajectory inference (specifically controlling the random walk step) to assess the reliability and robustness of predicted terminal clusters.
### We later observed that the random walk implementation from pecanpy (`from pecanpy import pecanpy as node2vec`) introduces slight variability in the results, even when setting a random seed through the VIA function. Specifically, the assignment `random_state = g.random_state` overrides the random seed set by VIA, causing subtle differences between runs of the same code.
### While we haven't yet identified a definitive solution within the source code to enforce complete reproducibility (pecanpy random walk), we found the results to be reasonably stable upon repeated runs. Moreover, by exploring multiple parameter combinations, we are still able to robustly evaluate the overall performance and consistency of our trajectory inference outcomes.

cell_meta = pd.read_csv("/path/to/via/via_meta_new1.csv", index_col=0) ## with annotation

# Define the criteria
criteria = cell_meta['treatment'] == 'control'

# Use np.where to get the locations (index positions) of all rows that match the criteria
matching_indices = np.where(criteria)[0]
cell_meta1=cell_meta.iloc[matching_indices]
true_labels1= cell_meta1['ann_KNN40'].tolist()
unique_elements1 = {val: idx for idx, val in enumerate(sorted(set(true_labels1)))}
labels_numeric1 = [unique_elements1[val] for val in true_labels1]
true_labels_numeric1= cell_meta1['time_C'].tolist()
harm1=harm[matching_indices,:]
umap1=embedding[matching_indices,:]

# Helper function to generate floating-point ranges
def frange(start, stop, step):
    while start <= stop:
        yield start
        start += step

# Define the parameter ranges
neighboring_thresholds = range(1, 7, 1)  # 1 to 6 with interval of 1
pseudotime_thresholds = range(40, 76, 5)  # 40 to 60 with interval of 5
x_lazy_vals = [round(x, 2) for x in frange(0.95, 0.99, 0.01)]  # 0.95 to 0.99 with interval of 0.01
alpha_teleport_vals = [round(x, 2) for x in frange(0.95, 0.99, 0.01)]  # 0.95 to 0.99 with interval of 0.01

# List to hold the results
results = []

# Function to run VIA for a single combination of parameters
def run_via_for_params(neighboring_threshold, pseudotime_threshold, x_lazy, alpha_teleport):
    # Initialize and run VIA with the given parameters
    v2 = via.VIA(data=harm1[:, 0:30], 
                 true_label=true_labels1, 
                 labels=labels_numeric1,
                 cluster_graph_pruning=0.15,
                 small_pop=200,
                 neighboring_terminal_states_threshold=neighboring_threshold,
                 knn=knn, 
                 too_big_factor=0.2,
                 root_user=['Naive Pluripotent Cells'],
                 dataset='group', 
                 random_seed=4, 
                 embedding=umap1, 
                 is_coarse=True, 
                 preserve_disconnected=True, 
                 pseudotime_threshold_TS=pseudotime_threshold, 
                 x_lazy=x_lazy,
                 alpha_teleport=alpha_teleport, 
                 time_series=True, 
                 time_series_labels=true_labels_numeric1, 
                 t_diff_step=2,
                 edgebundle_pruning_twice=False, 
                 knn_sequential=30,  
                 knn_sequential_reverse=10,
                 piegraph_arrow_head_width=0.15,
                 piegraph_edgeweight_scalingfactor=2.5,
                 edgebundle_pruning=0.15)
    
    # Run the VIA algorithm
    v2.run_VIA()

    # Return the parameter combination and the resulting terminal clusters
    return {
        'neighboring_threshold': neighboring_threshold,
        'pseudotime_threshold': pseudotime_threshold,
        'x_lazy': x_lazy,
        'alpha_teleport': alpha_teleport,
        'terminal_clusters': v2.terminal_clusters
    }
# Limit the number of workers to avoid exceeding memory limits
max_workers = 8  # Adjust this based on available memory (try 4, 8, etc.)

# Create a ProcessPoolExecutor for parallel execution
with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
    # Submit all parameter combinations for parallel execution
    future_to_params = {
        executor.submit(run_via_for_params, neighboring_threshold, pseudotime_threshold, x_lazy, alpha_teleport): 
        (neighboring_threshold, pseudotime_threshold, x_lazy, alpha_teleport)
        for neighboring_threshold in neighboring_thresholds
        for pseudotime_threshold in pseudotime_thresholds
        for x_lazy in x_lazy_vals
        for alpha_teleport in alpha_teleport_vals
    }
    
    # Collect the results as they complete
    for future in concurrent.futures.as_completed(future_to_params):
        try:
            result = future.result()
            results.append(result)
        except Exception as exc:
            params = future_to_params[future]
            print(f'Generated an exception with params {params}: {exc}')



# Convert the results to a DataFrame
results_df = pd.DataFrame(results)


# Save the results to a CSV or Excel file if needed
results_df.to_csv('/path/to/via_terminal_clusters_results.csv', index=False)

###Although pyVIA is nice to combine time information for trajectory inference, it only gave a cluster level pseudotime and would be unable to compare between control and Dox condition.
###Margaret was applied with the modification (time information included) to derive single cell level pseudotime



def convert_to_adjacency(n_augmented, d_augmented, n_cells):
    """
    Convert neighbor indices and distances from sequential_knn into adjacency matrices.
    
    Parameters:
    - n_augmented (ndarray): Neighbor indices of shape (n_cells, k).
    - d_augmented (ndarray): Distances corresponding to neighbor indices.
    - n_cells (int): Total number of cells.
    
    Returns:
    - adj_dist (csr_matrix): Distance-based adjacency matrix (n_cells x n_cells).
    - adj_conn (csr_matrix): Binary adjacency matrix (n_cells x n_cells).
    """

    row_idx = np.repeat(np.arange(n_cells), n_augmented.shape[1])  # Source nodes
    col_idx = n_augmented.flatten()  # Target nodes
    weights = d_augmented.flatten()  # Edge weights (distances)

    # Create sparse adjacency matrix with distances
    adj_dist = sp.csr_matrix((weights, (row_idx, col_idx)), shape=(n_cells, n_cells))

    # Convert to binary adjacency matrix (1 for connected, 0 for not connected)
    adj_conn = adj_dist.copy()
    adj_conn.data = np.ones_like(adj_conn.data)

    return adj_dist, adj_conn

def compute_undirected_cluster_connectivity(
    communities, adj, z_threshold=1.0, conn_threshold=None
):
    N = communities.shape[0]
    n_communities = np.unique(communities).shape[0]

    # Create cluster index
    clusters = {}
    for idx in np.unique(communities):
        cluster_idx = communities == idx
        clusters[idx] = cluster_idx

    undirected_cluster_connectivity = pd.DataFrame(
        np.zeros((n_communities, n_communities)),
        index=np.unique(communities),
        columns=np.unique(communities),
    )
    undirected_z_score = pd.DataFrame(
        np.zeros((n_communities, n_communities)),
        index=np.unique(communities),
        columns=np.unique(communities),
    )
    cluster_outgoing_edges = {}
    for i in np.unique(communities):
        cluster_i = clusters[i]

        # Compute the outgoing edges from the ith cluster
        adj_i = adj[cluster_i, :]
        adj_ii = adj_i[:, cluster_i]
        e_i = np.sum(adj_i) - np.sum(adj_ii)
        n_i = np.sum(cluster_i)
        cluster_outgoing_edges[i] = e_i

        for j in np.unique(communities):
            if i == j:
                continue
            # Compute the outgoing edges from the jth cluster
            cluster_j = clusters[j]
            adj_j = adj[cluster_j, :]
            adj_jj = adj_j[:, cluster_j]
            e_j = np.sum(adj_j) - np.sum(adj_jj)
            n_j = np.sum(cluster_j)

            # Compute the number of inter-edges from the ith to jth cluster
            adj_ij = adj_i[:, cluster_j]
            e_ij = np.sum(adj_ij)

            # Compute the number of inter-edges from the jth to ith cluster
            adj_ji = adj_j[:, cluster_i]
            e_ji = np.sum(adj_ji)
            e_sym = e_ij + e_ji

            # Compute the random assignment of edges from the ith to the jth
            # cluster under the PAGA binomial model
            e_sym_random = (e_i * n_j + e_j * n_i) / (N - 1)

            # Compute the cluster connectivity measure
            std_sym = (e_i * n_j * (N - n_j - 1) + e_j * n_i * (N - n_i - 1)) / (
                N - 1
            ) ** 2
            undirected_z_score.loc[i, j] = (e_sym - e_sym_random) / std_sym

            # Only add non-spurious edges based on a threshold
            undirected_cluster_connectivity.loc[i, j] = (e_sym - e_sym_random) / (
                e_i + e_j - e_sym_random
            )
            if conn_threshold is not None:
                if undirected_cluster_connectivity.loc[i, j] < conn_threshold:
                    undirected_cluster_connectivity.loc[i, j] = 0
            elif undirected_z_score.loc[i, j] < z_threshold:
                undirected_cluster_connectivity.loc[i, j] = 0
    return undirected_cluster_connectivity, undirected_z_score

def prune_network_edges(communities, adj_sc, adj_cluster):
    """Prune network edges using sparse slicing."""
    cluster_ids = np.unique(communities)

    # Create cluster index
    clusters = {idx: np.where(communities == idx)[0] for idx in cluster_ids}

    # Convert adjacency matrix to LIL format for efficient modification
    adj_sc_lil = lil_matrix(adj_sc)

    for c_idx in adj_cluster.index:
        cluster_i = clusters[c_idx]
        non_connected_clusters = adj_cluster.columns[adj_cluster.loc[c_idx, :] == 0]
        for nc_idx in non_connected_clusters:
            if nc_idx == c_idx:
                continue
            cluster_nc = clusters[nc_idx]

            # Remove edges between non-connected clusters
            adj_sc_lil[np.ix_(cluster_i, cluster_nc)] = 0

    return adj_sc_lil.tocsr()

def adjust_weights_with_time(adj_matrix, time_points, penalty, threshold):
    """
    Adjust adjacency matrix weights based on time points to enforce chronological order.
    """
    adj_matrix = adj_matrix.tocoo()  # Convert to COO for efficient indexing
    row, col, data = adj_matrix.row, adj_matrix.col, adj_matrix.data  # Extract COO format

    # Compute time differences only for existing edges
    time_diff = time_points[col] - time_points[row]

    # Ensure correct shape matching
    if time_diff.shape != data.shape:
        raise ValueError(f"Mismatch: time_diff has shape {time_diff.shape}, but adj_matrix.data has shape {data.shape}")

    # Identify edges violating time constraints
    penalty_mask_high = np.abs(time_diff) > threshold
    penalty_mask_small = (time_diff < 0) & (np.abs(time_diff) < threshold)

    # Apply penalties only to edges that exist
    data[penalty_mask_high] += penalty
    data[penalty_mask_small] += penalty / 2

    # Convert back to CSR format for efficient computation
    return csr_matrix((data, (row, col)), shape=adj_matrix.shape)

def compute_connectivity_graph(
    embeddings, communities, cluster_connectivities, mode="undirected"
):
    assert mode in ["directed", "undirected"]
    g = nx.Graph() if mode == "undirected" else nx.DiGraph()
    node_positions = {}
    cluster_ids = np.unique(communities)
    for i in cluster_ids:
        g.add_node(i)
        # determine the node pos for the cluster
        cluster_i = communities == i
        node_pos = np.mean(embeddings[cluster_i, :], axis=0)
        node_positions[i] = node_pos

    n_nodes = len(cluster_ids)
    for row_id, i in enumerate(cluster_ids):
        for col_id, j in enumerate(cluster_ids):
            if cluster_connectivities.loc[i, j] > 0:
                g.add_edge(
                    cluster_ids[row_id],
                    cluster_ids[col_id],
                    weight=cluster_connectivities.loc[i, j],
                )
    return g, node_positions

def compute_connectivity_graph(
    embeddings, communities, cluster_connectivities, mode="undirected"
):
    assert mode in ["directed", "undirected"]
    g = nx.Graph() if mode == "undirected" else nx.DiGraph()
    node_positions = {}
    cluster_ids = np.unique(communities)
    for i in cluster_ids:
        g.add_node(i)
        # determine the node pos for the cluster
        cluster_i = communities == i
        node_pos = np.mean(embeddings[cluster_i, :], axis=0)
        node_positions[i] = node_pos

    n_nodes = len(cluster_ids)
    for row_id, i in enumerate(cluster_ids):
        for col_id, j in enumerate(cluster_ids):
            if cluster_connectivities.loc[i, j] > 0:
                g.add_edge(
                    cluster_ids[row_id],
                    cluster_ids[col_id],
                    weight=cluster_connectivities.loc[i, j],
                )
    return g, node_positions

def compute_pseudotime_v3(
    ad,
    start_cell_ids,
    adj_dist,
    adj_cluster,
    time_points,  # New parameter for cell time points
    comm_key="metric_clusters",
    data_key="metric_embedding",
    prune_edges=True,  # Optional pruning
    penalty=None,  # High penalty for time-violating paths
    threshold=24,  # Fixed typo from "threhold"
):
    """Compute pseudotime with time constraint penalties to enforce chronological order."""
    
    communities = ad.obs[comm_key]
    cluster_ids = np.unique(communities)
    data = ad.obsm[data_key]
    
    # Create cluster index
    clusters = {idx: np.where(communities == idx)[0] for idx in cluster_ids}

    # Prune the initial adjacency matrix (optional)
    if prune_edges:
        adj_dist_pruned = prune_network_edges(communities, adj_dist, adj_cluster)
    else:
        adj_dist_pruned = adj_dist.copy()
    # Determine the max value in the adjacency matrix
    max_value = adj_dist_pruned.mean()
    
    # Use max_value if penalty is not provided
    penalty_value = penalty if penalty is not None else max_value

    # Add penalties for time-violating paths
    time_points_array = np.array(time_points)
    adj_dist_pruned = adjust_weights_with_time(
        adj_dist_pruned, time_points_array, penalty_value, threshold
    )

    # Compute initial pseudotime using Dijkstra
    obs_names = ad.obs_names.to_numpy()
    start_indices = [np.where(obs_names == s)[0][0] for s in start_cell_ids]
    pseudotime = dijkstra(
        adj_dist_pruned, directed=False, indices=start_indices, return_predecessors=False
    )

    # Convert to pandas Series
    pseudotime = pd.Series(np.min(pseudotime, axis=0), index=obs_names)

    # Update cluster-specific graphs
    for cluster_id, cluster_idx in clusters.items():
        # Extract the subgraph for this cluster
        adj_sc = adj_dist_pruned[cluster_idx, :][:, cluster_idx]
        cluster_pseudo = pseudotime.iloc[cluster_idx]

        # Identify the start cell within the cluster
        cluster_start_cell = cluster_pseudo.idxmin()
        start_idx = np.where(obs_names[cluster_idx] == cluster_start_cell)[0][0]

        # Reconnect the subgraph
        adj_sc = connect_graph(adj_sc, data[cluster_idx, :], start_idx)

        # Update the pruned adjacency matrix
        adj_dist_pruned[cluster_idx, :][:, cluster_idx] = adj_sc

    # Recompute pseudotime with updated graph
    pseudotime = dijkstra(
        adj_dist_pruned, directed=False, indices=start_indices, return_predecessors=False
    )
    pseudotime = pd.Series(np.min(pseudotime, axis=0), index=obs_names)

    # Ensure no infinite values
    pseudotime.replace([np.inf, -np.inf], np.nan, inplace=True)
    pseudotime.fillna(0, inplace=True)

    # Store pseudotime in AnnData object
    ad.obs["metric_pseudotime_v2"] = pseudotime
    return pseudotime

v3 = via.VIA(data=harm[:, 0:30], true_label=true_labels, labels=labels_numeric,
             cluster_graph_pruning=0.15,neighboring_terminal_states_threshold=3,small_pop=200,
             knn=40, too_big_factor=0.2,root_user=['Naive_Pluripotent_Cells'],dataset='group', random_seed=4, embedding = embedding,
            is_coarse=True, preserve_disconnected=True, pseudotime_threshold_TS=40, x_lazy=0.99,
         alpha_teleport=0.99, time_series=True, time_series_labels=true_labels_numeric, t_diff_step=3,
         edgebundle_pruning_twice=False, knn_sequential=30, knn_sequential_reverse=10,
             piegraph_arrow_head_width=0.15,piegraph_edgeweight_scalingfactor=2.5,edgebundle_pruning=0.15)
v3.knn_struct = via._construct_knn(v3.data, knn=v3.knn, distance=v3.distance, num_threads=v3.num_threads)
neighbors, distances = v3.knn_struct.knn_query(v3.data, k=v3.knn)
adjacency_augmented = None
adjacency = None
n_augmented, d_augmented = via.sequential_knn(v3.data, v3.time_series_labels, neighbors, distances,k_seq=v3.knn_sequential,k_reverse=v3.knn_sequential_reverse,
                                              num_threads=v3.num_threads, distance=v3.distance)
adjacency_augmented = v3._make_csrmatrix_noselfloop(n_augmented, d_augmented,time_series=v3.time_series,time_series_labels=v3.time_series_labels,t_diff_step=v3.t_diff_step)

# Example usage
n_cells = n_augmented.shape[0]
adj_dist1, adj_conn1 = convert_to_adjacency(n_augmented, d_augmented, n_cells)

com=np.array(labels_numeric)
un_connectivity1, un_z_score1 = compute_undirected_cluster_connectivity(com, adj_conn1, z_threshold=0.6)

ref=sc.read_h5ad("/path/to/via/ref_RH.h5ad")  ## convert/construct anndata object from seurat, see network analysis script for more detail. 
ref.obsm['X_umap']=embedding
ref.obs['label_numeric']=com

##### Define the criteria
criteria = (ref.obs['ann_new_noZX_KNN40'] == 'Naive_Pluripotent_Cells')

# Use np.where to get the locations (index positions) of all rows that match the criteria
matching_indices = np.where(criteria)[0]
start_ids=ref.obs_names[matching_indices]

G_undirected1, node_positions1 = compute_connectivity_graph(embedding, com, un_connectivity1)
adj_cluster1 = nx.to_pandas_adjacency(G_undirected1)

time=np.array(true_labels_numeric)

pseudotime3 = compute_pseudotime_v3(ref, start_ids, adj_dist1, adj_cluster1, time_points=time, comm_key="label_numeric", data_key="X_umap", prune_edges=True, penalty=None, threshold=24)

pseudotime_df = pd.DataFrame({
    "Condition": pseudotime3
})
pseudotime_df.to_csv("/path/to/via/pseu_condition.csv")

#the correlation between pseudotime and real timepoint
pseudotime_df['time']=cell_meta['time_C']
correlation = pseudotime_df['time'].corr(pseudotime_df['Condition'])
print(f'Correlation of condition3 with time {round(correlation * 100, 2)} %')

#cluster level pseudotime by average
average_expression = pseudotime_df.groupby('cluster')['Condition_3'].mean()
print(f'condtion_3: {np.array(average_expression)}')

###compare the pseudotime shift/change for each cell type between control and Dox
###There is subtle difference and it seems not improve too much comparing to cell type level pseudotime from pyVIA. However, we still used it for comprehensive map construction.
pseudotime_df['treatment']=cell_meta['treatment']
criteria = cell_meta['cell'] == 'RHKI'
pseudo1=pseudotime_df.loc[criteria]

# Create a FacetGrid (like facet_wrap(~cluster))
g = sns.FacetGrid(pseudo1, col="cluster", col_wrap=5, sharex=True, sharey=True)

# Map kdeplot to each facet (similar to ggplot2 facet_wrap)
g.map_dataframe(sns.kdeplot, x="Condition_3RT2", hue="treatment", fill=True, alpha=0.4)

# Adjust labels and titles
g.set_axis_labels("Feature Value", "Density")
g.set_titles(col_template="{col_name}")  # Titles for each facet
g.add_legend(title="Cell Type")
plt.show()

##scale pseudotime for the extreme value as it did in pyVIA
pseudotime=pseudotime_df['Condition']
very_high = np.mean(pseudotime) + 1.5 * np.std(pseudotime)
without_very_high_pt = [iii for iii in pseudotime if iii < 1000]
new_very_high = np.mean(without_very_high_pt) + np.std(without_very_high_pt)
# print('very high, and new very high', very_high, new_very_high)
new_hitting_times = [x if x < 1000 else 1000 for x in pseudotime]
hitting_times = np.asarray(new_hitting_times)
scaling_fac = 10 / max(hitting_times)
hitting_times = hitting_times * scaling_fac

ref.obs['pseudotime']=hitting_times

###For trajectory comparison, first, we constructed the trajectory map from all timepoints in control condition of our data plus publised data
###`t_diff=4` would equal to 34 hours difference for time adjanceyment matrix trimming for all timepoint data
### Parameters such as alpha_teleport and x_lazy were not quite important here. We would use full graph pyVIA with the potential edges to reconstruct the comprehensive maps.

##### Define the criteria
criteria = (cell_meta['treatment'] == 'control')

# Use np.where to get the locations (index positions) of all rows that match the criteria
matching_indices = np.where(criteria)[0]
cell_meta1=cell_meta.iloc[matching_indices]
true_labels1= cell_meta1['ann_KNN40'].tolist()
true_labels_numeric1= cell_meta1['time_C'].tolist()
unique_elements1 = {val: idx for idx, val in enumerate(sorted(set(true_labels1)))}
labels_numeric1 = [unique_elements1[val] for val in true_labels1]
harm1=harm[matching_indices,:]
umap1=embedding[matching_indices,:]

v4 = via.VIA(data=harm1[:, 0:30], true_label=true_labels1, labels=labels_numeric1,
             cluster_graph_pruning=0.3,neighboring_terminal_states_threshold=4,small_pop=200,
             knn=40, too_big_factor=0.2,root_user=['Naive_Pluripotent_Cells'],dataset='group', random_seed=4, embedding = umap1,
            is_coarse=True, preserve_disconnected=True, pseudotime_threshold_TS=40, x_lazy=0.90,
         alpha_teleport=0.99, time_series=True, time_series_labels=true_labels_numeric1, t_diff_step=3,
         edgebundle_pruning_twice=False, knn_sequential=30, knn_sequential_reverse=10,
             piegraph_arrow_head_width=0.15,piegraph_edgeweight_scalingfactor=2.5,edgebundle_pruning=0.3)

v4.run_VIA()

###`time_series=False` would set for RHKI control and Dox comparison since >1 would allow time interval larger then the reasonable time interval (time information would not be quite useful with just 3 timepoints)
##### RHKI control
criteria1 = (cell_meta['treatment'] == 'control') & (cell_meta['cell'] == 'RHKI')

# Use np.where to get the locations (index positions) of all rows that match the criteria
matching_indices1 = np.where(criteria1)[0]
cell_meta2=cell_meta.iloc[matching_indices1]
true_labels2= cell_meta2['ann_KNN40'].tolist()
true_labels_numeric2= cell_meta2['time_C'].tolist()
unique_elements2 = {val: idx for idx, val in enumerate(sorted(set(true_labels2)))}
labels_numeric2 = [unique_elements2[val] for val in true_labels2]
harm2=harm[matching_indices1,:]
umap2=embedding[matching_indices1,:]

v5 = via.VIA(data=harm2[:, 0:30], true_label=true_labels2, labels=labels_numeric2,
             cluster_graph_pruning=0.3,small_pop=200,neighboring_terminal_states_threshold=3,
             knn=40, too_big_factor=0.2,root_user=['Naive_Pluripotent_Cells'],dataset='group', random_seed=4, embedding = umap2,
            is_coarse=True, preserve_disconnected=False, preserve_disconnected_after_pruning=True, pseudotime_threshold_TS=40, x_lazy=0.94,
         alpha_teleport=0.94, time_series=False, time_series_labels=true_labels_numeric2, t_diff_step=2,
         edgebundle_pruning_twice=False, knn_sequential=30, knn_sequential_reverse=10,
             piegraph_arrow_head_width=0.15,piegraph_edgeweight_scalingfactor=2.5,edgebundle_pruning=0.3)
v5.run_VIA()

##### RHKI Dox
criteria2 = (cell_meta['treatment'] == 'Dox') & (cell_meta['cell'] == 'RHKI')

# Use np.where to get the locations (index positions) of all rows that match the criteria
matching_indices2 = np.where(criteria2)[0]
cell_meta3=cell_meta.iloc[matching_indices2]
true_labels3= cell_meta3['ann_KNN40'].tolist()
true_labels_numeric3= cell_meta3['time_C'].tolist()
unique_elements3 = {val: idx for idx, val in enumerate(sorted(set(true_labels3)))}
labels_numeric3 = [unique_elements3[val] for val in true_labels3]
harm3=harm[matching_indices2,:]
umap3=embedding[matching_indices2,:]

v6 = via.VIA(data=harm3[:, 0:30], true_label=true_labels3, labels=labels_numeric3,
             cluster_graph_pruning=0.3,small_pop=200,neighboring_terminal_states_threshold=3,
             knn=40, too_big_factor=0.2,root_user=['Naive_Pluripotent_Cells'],dataset='group', random_seed=4, embedding = umap3,
            is_coarse=True, preserve_disconnected=False, preserve_disconnected_after_pruning=True, pseudotime_threshold_TS=40, x_lazy=0.94,
         alpha_teleport=0.94, time_series=False, time_series_labels=true_labels_numeric3, t_diff_step=2,
         edgebundle_pruning_twice=False, knn_sequential=30, knn_sequential_reverse=10,
             piegraph_arrow_head_width=0.15,piegraph_edgeweight_scalingfactor=2.5,edgebundle_pruning=0.3)
v6.run_VIA()

##For simplify trajectory map in Fig. 5C (Sankey plot)
def plot_differentiation_flow2(
        via_object, idx: list = None, dpi=150, marker_lineages=[], label_node: list = [],
        do_log_flow: bool = True, fontsize: int = 8, alpha_factor: float = 0.9,
        majority_cluster_population_dict: dict = None, cmap_sankey='rainbow',
        title_str: str = 'Differentiation Flow', root_cluster_list: list = None,
        export_plot_path=None, export_matrix_path=None, export_format="pdf", rotate_plot=True
):
    import math
    import matplotlib.pyplot as plt
    from bokeh.io.export import export_svgs
    from bokeh.io.export import export_png
    from bokeh.io import output_file, save
    from holoviews import opts, dim
    import holoviews as hv
    from holoviews.plotting.util import process_cmap
    import plotly.graph_objects as go

    if len(marker_lineages) == 0:
        marker_lineages = via_object.terminal_clusters
    if root_cluster_list is None:
        root_cluster_list = via_object.root

    else:
        marker_lineages = [i for i in marker_lineages if i in via_object.labels]  # via_object.terminal_clusters]
    print(f'{datetime.now()}\tMarker_lineages: {marker_lineages}')

    # make the sankey node labels either using via_obect.true_label or the labels provided by the user
    print(f'{datetime.now()}\tStart dictionary modes')
    df_mode = pd.DataFrame()
    df_mode['cluster'] = via_object.labels
    # df_mode['celltype'] = pre_labels_celltype_df['fine'].tolist()#v0.true_label
    if len(label_node) > 0:
        df_mode['celltype'] = label_node  # v0.true_label
    else:
        df_mode['celltype'] = via_object.true_label
    majority_cluster_population_dict = df_mode.groupby(['cluster'])['celltype'].agg(
        lambda x: pd.Series.mode(x)[0])  # agg(pd.Series.mode would give all modes) #series
    majority_cluster_population_dict = majority_cluster_population_dict.to_dict()
    print(f'{datetime.now()}\tEnd dictionary modes')

    if idx is None: idx = np.arange(0, via_object.nsamples)
    # G = via_object.full_graph_shortpath
    n_original_comp, n_original_comp_labels = connected_components(via_object.csr_full_graph, directed=False)
    # G = via_object.full_graph_paths(via_object.data, n_original_comp)
    # knn_hnsw = _make_knn_embeddedspace(embedding)
    y_root = []
    x_root = []
    root1_list = []
    p1_sc_bp = np.nan_to_num(via_object.single_cell_bp[idx, :], nan=0.0, posinf=0.0, neginf=0.0)
    # row normalize
    row_sums = p1_sc_bp.sum(axis=1)
    p1_sc_bp = p1_sc_bp / row_sums[:,
                          np.newaxis]  # make rowsums a column vector where i'th entry is sum of i'th row in p1-sc-bp
    print(f'{datetime.now()}\tCheck sc pb {p1_sc_bp[0, :].sum()} ')

    p1_labels = np.asarray(via_object.labels)[idx]

    p1_cc = via_object.connected_comp_labels
    p1_sc_pt_markov = list(np.asarray(via_object.single_cell_pt_markov)[idx])
    X_data = via_object.data

    X_ds = X_data[idx, :]
    p_ds = hnswlib.Index(space='l2', dim=X_ds.shape[1])
    p_ds.init_index(max_elements=X_ds.shape[0], ef_construction=200, M=16)
    p_ds.add_items(X_ds)
    p_ds.set_ef(50)
    num_cluster = len(set(via_object.labels))
    G_orange = ig.Graph(n=num_cluster, edges=via_object.edgelist_maxout,
                        edge_attrs={'weight': via_object.edgeweights_maxout})
    for ii, r_i in enumerate(root_cluster_list):
        sankey_edges = []
        for fst_i in marker_lineages:

            path_orange = G_orange.get_shortest_paths(root_cluster_list[ii], to=fst_i)[0]
            # path_orange = G_orange.get_shortest_paths(3, to=fst_i)[0]
            # if the roots is in the same component as the terminal cluster, then print the path to output
            if len(path_orange) > 0:
                print(
                    f'{datetime.now()}\tCluster path on clustergraph starting from Root Cluster {root_cluster_list[ii]} to Terminal Cluster {fst_i}: {path_orange}')
                do_sankey = True
                if do_sankey:

                    import holoviews as hv
                    hv.extension('bokeh')
                    from bokeh.plotting import show
                    from holoviews import opts, dim

                    print(f"{datetime.now()}\tHoloviews for TC {fst_i}")

                    cluster_adjacency = via_object.cluster_adjacency
                    # row normalize
                    row_sums = cluster_adjacency.sum(axis=1)
                    cluster_adjacency_rownormed = cluster_adjacency / row_sums[:, np.newaxis]

                    for n_i in range(len(path_orange) - 1):
                        source = path_orange[n_i]
                        dest = path_orange[n_i + 1]
                        if n_i < len(path_orange) - 2:
                            if do_log_flow:

                                val_edge = round(math.log1p(cluster_adjacency_rownormed[source, dest]),
                                                 2)  # * cluster_population_dict[source] #  natural logarithm (base e) of 1 + x
                            else:
                                val_edge = round(cluster_adjacency_rownormed[source, dest], 2)
                                # print("clipping val edge")
                                # if val_edge > 0.5: val_edge = 0.5
                        else:
                            if dest in via_object.terminal_clusters:
                                ts_array_original = np.asarray(via_object.terminal_clusters)
                                loc_ts_current = np.where(ts_array_original == dest)[0][0]
                                print(f'dest {dest}, is at loc {loc_ts_current} on the bp_array')
                                if do_log_flow:
                                    val_edge = round(math.log1p(via_object.cluster_bp[source, loc_ts_current]),
                                                     2)  # * cluster_population_dict[source]
                                else:
                                    val_edge = round(via_object.cluster_bp[source, loc_ts_current], 2)
                                    # print("clipping val edge")
                                    # if val_edge > 0.5: val_edge = 0.5
                            else:
                                if do_log_flow:

                                    val_edge = round(math.log1p(cluster_adjacency_rownormed[source, dest]),
                                                     2)  # * cluster_population_dict[source] #  natural logarithm (base e) of 1 + x
                                else:
                                    val_edge = round(cluster_adjacency_rownormed[source, dest], 2)
                                    # print("clipping val edge")
                                    # val_edge = 0.5
                        # sankey_edges.append((majority_cluster_population_dict[source]+'_C'+str(source), majority_cluster_population_dict[dest]+'_C'+str(dest), val_edge))#, majority_cluster_population_dict[source],majority_cluster_population_dict[dest]))
                        sankey_edges.append((source, dest,
                                             val_edge))  # ,majority_cluster_population_dict[source]+'_C'+str(source),'magenta' ))

    
    # print(f'pre-final sankey set of edges and vals {len(sankey_edges)}, {sankey_edges}')
    source_dest = list(set(sankey_edges))
    # print(f'final sankey set of edges and vals {len(source_dest)}, {source_dest}')
    source_dest_df = pd.DataFrame(source_dest, columns=['Source', 'Dest', 'Count'])  # ,'Label','Color'])
    nodes_in_source_dest = list(set(set(source_dest_df.Source) | set(source_dest_df.Dest)))
    nodes_in_source_dest.sort()
    convert_old_to_new = {}
    convert_new_to_old = {}
    majority_newcluster_population_dict = {}
    for ei, ii in enumerate(nodes_in_source_dest):
        convert_old_to_new[ii] = ei
        convert_new_to_old[ei] = ii
        majority_newcluster_population_dict[ei] = majority_cluster_population_dict[ii]
    source_dest_new = []
    for tuple_ in source_dest:
        source_dest_new.append((convert_old_to_new[tuple_[0]], convert_old_to_new[tuple_[1]], tuple_[2]))
    # print('new source dest after reindexing', source_dest_new)
    # nodes = [majority_cluster_population_dict[i] for i in range(len(majority_cluster_population_dict))]
    # nodes = [majority_cluster_population_dict[i] for i in nodes_in_source_dest]
    nodes = [majority_newcluster_population_dict[key] + '_C' + str(convert_new_to_old[key]) for key in
                majority_newcluster_population_dict]
    # nodes = ['C' + str(convert_new_to_old[key]) for key in                 majority_newcluster_population_dict]
    # print('nodes', len(nodes), nodes,)
    # Export the result matrix if requested
    if export_matrix_path:
        print(f'{datetime.now()}\tExporting Sankey data matrix to {export_matrix_path}')
        source_dest_df.to_csv(export_matrix_path, index=False)
    # Build and export the Sankey plot
    # Prepare data for Plotly Sankey
    unique_nodes = sorted(set(source_dest_df['Source']).union(set(source_dest_df['Dest'])))
    node_labels = [majority_cluster_population_dict.get(node, f"Cluster {node}") for node in unique_nodes]
    node_map = {node: idx for idx, node in enumerate(unique_nodes)}

    # Map source and destination nodes to their indices
    source_dest_df['Mapped_Source'] = source_dest_df['Source'].map(node_map)
    source_dest_df['Mapped_Dest'] = source_dest_df['Dest'].map(node_map)

    # Create Plotly Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=40,  # Increase spacing between nodes
            thickness=30,
            line=dict(color="black", width=0.5),
            label=node_labels,
        ),
        link=dict(
            source=source_dest_df['Mapped_Source'],
            target=source_dest_df['Mapped_Dest'],
            value=source_dest_df['Count']
        )
    )])
    if rotate_plot:
        fig.update_layout(
            title_text=title_str,
            font=dict(size=20),
            xaxis=dict(autorange='reversed'),
            yaxis=dict(scaleanchor="x", scaleratio=1), 
        )
    else:
        fig.update_layout(title_text=title_str, font=dict(size=10))
    # Export plot
    if export_plot_path:
        print(f'{datetime.now()}\tExporting Sankey plot to {export_plot_path}')
        if export_plot_path.endswith('.svg') or export_plot_path.endswith('.pdf'):
            fig.write_image(export_plot_path, )
        else:
            print("Unsupported format. Only .svg or .pdf is supported for export.")

    return fig

plot_differentiation_flow2(via_object=v4,do_log_flow=True, root_cluster_list=[4], title_str="origin", rotate_plot=False,
                           export_plot_path="/path/to/via/differentiation_flow_control.svg",
                           export_matrix_path="/path/to/via/sankey_data_control.csv", export_format="svg")


pseudotime_df['condition_ST']=hitting_times
pseudotime_df1=pseudotime_df.loc[criteria] 
pseudotime_df2=pseudotime_df.loc[criteria1] 
pseudotime_df3=pseudotime_df.loc[criteria2] 
average_expression1 = pseudotime_df1.groupby('cluster')['condition_ST'].mean() ### cluster level pseudotime for all control data
average_expression2 = pseudotime_df2.groupby('cluster')['condition_ST'].mean() ### cluster level pseudotime for RHKI control
average_expression3 = pseudotime_df3.groupby('cluster')['condition_ST'].mean() ### cluster level pseudotime for RHKI Dox
print(average_expression1)
print(average_expression2)
print(average_expression3)

###For comprehensive map, only keep the potential edges based on all parameter combination trial above rather single parameter which only allow one most possible edge for each node.
###Based on these edges, recalculating weight as it did through pyVIA

#all control
edges_all=[(3, 0), (4, 3), (3,12), (6, 0), (7, 6), (11, 5), (10,15), (11, 10), (12, 6), 
           (13, 1), (13, 3), (13, 5), (13, 6), (14, 7), (14, 9), (15, 2), (15, 11), (6,1), (7,1), (7,8), (8,10)]
graph = ig.VertexClustering(v4.ig_full_graph, membership=v4.labels).cluster_graph(combine_edges='sum')
graph = recompute_weights(graph, Counter(v4.labels))  # returns csr matrix
from scipy.sparse.csgraph import minimum_spanning_tree
Tcsr = csr_mst(graph)
initial_links_n = len(graph.data)
n_comp, comp_labels = connected_components(csgraph=graph, directed=False, return_labels=True)
rows, cols = zip(*edges_all)
mask = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=graph.shape)
filtered_csr = graph.multiply(mask)
weights = filtered_csr.data / (np.std(filtered_csr.data))
edges = list(zip(*filtered_csr.nonzero()))
bias_weights_2_A = get_biased_weights(edges, weights,average_expression1,round=2)
n_clus = len(set(v4.labels))
temp_csr = csr_matrix((bias_weights_2_A, tuple(zip(*edges))),shape=(n_clus, n_clus))
temp_csr = temp_csr.transpose().todense() + temp_csr.todense()
temp_csr = np.tril(temp_csr, -1)  # elements along the main diagonal and above are set to zero
temp_csr = csr_matrix(temp_csr)
edgeweights_maxout_A = temp_csr.data
scale_factor = max(edgeweights_maxout_A) - min(edgeweights_maxout_A)
edgeweights_maxout_A = [((wi + .1) * 2.5 / scale_factor) + 0.1 for wi in edgeweights_maxout_A]

sources, targets = temp_csr.nonzero()
edgelist_maxout_A = list(zip(sources.tolist(), targets.tolist()))

#RHKI control
edges_all=[(3, 0), (4, 3), (3,12), (6, 0), (7, 6), (11, 5), (10,15), (11, 10), (12, 6), 
           (13, 1), (13, 3), (13, 5), (13, 6), (14, 7), (14, 9), (15, 2), (15, 11), (6,1), (7,1), (7,8), (8,10)]
graph = ig.VertexClustering(v5.ig_full_graph, membership=v5.labels).cluster_graph(combine_edges='sum')
graph = recompute_weights(graph, Counter(v5.labels))  # returns csr matrix
from scipy.sparse.csgraph import minimum_spanning_tree
Tcsr = csr_mst(graph)
initial_links_n = len(graph.data)
n_comp, comp_labels = connected_components(csgraph=graph, directed=False, return_labels=True)
rows, cols = zip(*edges_all)
mask = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=graph.shape)
filtered_csr = graph.multiply(mask)
weights = filtered_csr.data / (np.std(filtered_csr.data))
edges = list(zip(*filtered_csr.nonzero()))
bias_weights_2_C = get_biased_weights(edges, weights,average_expression2,round=2)
n_clus = len(set(v5.labels))
temp_csr = csr_matrix((bias_weights_2_C, tuple(zip(*edges))),shape=(n_clus, n_clus))
temp_csr = temp_csr.transpose().todense() + temp_csr.todense()
temp_csr = np.tril(temp_csr, -1)  # elements along the main diagonal and above are set to zero
temp_csr = csr_matrix(temp_csr)
edgeweights_maxout_C = temp_csr.data
scale_factor = max(edgeweights_maxout_C) - min(edgeweights_maxout_C)
edgeweights_maxout_C = [((wi + .1) * 2.5 / scale_factor) + 0.1 for wi in edgeweights_maxout_C]

sources, targets = temp_csr.nonzero()
edgelist_maxout_C = list(zip(sources.tolist(), targets.tolist()))

#RHKI Dox
edges_all=[(3, 0), (4, 3), (3,12) (6, 0), (7, 6), (11, 5), (10,15), (11, 10), (12, 6), 
           (13, 1), (13, 3), (13, 5), (13, 6), (14, 7), (14, 9), (15, 2), (15, 11), (6,1), (7,1), (7,8), (8,10)]
graph = ig.VertexClustering(v6.ig_full_graph, membership=v6.labels).cluster_graph(combine_edges='sum')
graph = recompute_weights(graph, Counter(v6.labels))  # returns csr matrix
from scipy.sparse.csgraph import minimum_spanning_tree
Tcsr = csr_mst(graph)
initial_links_n = len(graph.data)
n_comp, comp_labels = connected_components(csgraph=graph, directed=False, return_labels=True)
rows, cols = zip(*edges_all)
mask = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=graph.shape)
filtered_csr = graph.multiply(mask)
weights = filtered_csr.data / (np.std(filtered_csr.data))
edges = list(zip(*filtered_csr.nonzero()))
bias_weights_2_D = get_biased_weights(edges, weights,average_expression3,round=2)
n_clus = len(set(v6.labels))
temp_csr = csr_matrix((bias_weights_2_D, tuple(zip(*edges))),shape=(n_clus, n_clus))
temp_csr = temp_csr.transpose().todense() + temp_csr.todense()
temp_csr = np.tril(temp_csr, -1)  # elements along the main diagonal and above are set to zero
temp_csr = csr_matrix(temp_csr)
edgeweights_maxout_D = temp_csr.data
scale_factor = max(edgeweights_maxout_D) - min(edgeweights_maxout_D)
edgeweights_maxout_D = [((wi + .1) * 2.5 / scale_factor) + 0.1 for wi in edgeweights_maxout_D]

sources, targets = temp_csr.nonzero()
edgelist_maxout_D = list(zip(sources.tolist(), targets.tolist()))

###comprehensive map
import igraph as ig
import matplotlib.pyplot as plt
num_cluster = len(set(v4.labels))
graph1=ig.Graph(n=num_cluster, edges=edgelist_maxout_A, edge_attrs={'weight': edgeweights_maxout_A})
graph2=ig.Graph(n=num_cluster, edges=edgelist_maxout_C, edge_attrs={'weight': edgeweights_maxout_C})
graph3=ig.Graph(n=num_cluster, edges=edgelist_maxout_D, edge_attrs={'weight': edgeweights_maxout_D})
node_names = ['Caudal_Neuroectoderm1', 'Definitive_Endoderm_Guts', 'Endothelium', 'Epiblasts', 
              'Naive_Pluripotent_Cells', 'Nascent_Mesoderm', 'Neuromesodermal_Progenitors1', 
              'Neuromesodermal_Progenitors2', 'Neuromesodermal_Progenitors3', 'Neuron-like_Cells', 
              'Paraxial_Mesoderm_A', 'Paraxial_Mesoderm_B', 'Placodal_area', 'Primitive_Streak', 
              'Spinal_Cord', 'Splanchnic_Mesoderm']
graph1.vs['name'] = node_names[:len(graph1.vs)]


layout_fixed = graph1.layout('fruchterman_reingold')

weights1 = dict(zip(graph1.get_edgelist(), graph1.es['weight']))
edge_widths1 = [1+2*weights1[edge] for edge in graph1.get_edgelist()]
weights2 = dict(zip(graph2.get_edgelist(), graph2.es['weight']))
edge_widths2 = [1+2*weights2[edge] for edge in graph2.get_edgelist()]
weights3 = dict(zip(graph3.get_edgelist(), graph3.es['weight']))
edge_widths3 = [1+2*weights3[edge] for edge in graph3.get_edgelist()]

weights_diff = {}
for edge in graph1.get_edgelist():
    weight2 = weights2.get(edge, 0)
    weight3 = weights3.get(edge, 0)
    weights_diff[edge] = weight3 - weight2

# Step 8: Set edge colors and widths for the difference visualization
def compute_edge_colors_and_widths(all_edges, weights_diff):
    weight_diffs = np.array([abs(w) for w in weights_diff.values()])
    percentile = np.percentile(weight_diffs, 70)
    
    edge_colors = []
    edge_widths = []
    for edge in all_edges:
        if abs(weights_diff.get(edge, 0)) <= percentile:
            edge_colors.append('black')  # Within 25th percentile
        else:
            diff = weights_diff[edge]
            if diff == 0:
                edge_colors.append('wheat')  # Same weight
            elif diff > 0:
                edge_colors.append('skyblue')  # Higher in Graph1
            else:
                edge_colors.append('darkorange')  # Higher in Graph2
        edge_widths.append(2 + 5 * abs(weights_diff.get(edge, 0)))  # Scaled widths
    return edge_colors, edge_widths
edge_colors, edge_widths_C = compute_edge_colors_and_widths(graph1.get_edgelist(), weights_diff)

combine=ig.Graph(edges=graph1.get_edgelist())

fig, ax = plt.subplots(2,2, figsize=(16, 16))
ig.plot(
    graph1,
    layout=layout_fixed,
    target=ax[0,0],
    edge_width=edge_widths1,
    vertex_size=40,
    vertex_color="white",
    vertex_label=graph1.vs['name'],
    vertex_frame_color="black",
)

# Plot graph2
ig.plot(
    graph2,
    layout=layout_fixed,
    target=ax[0,1],
    edge_width=edge_widths2,
    vertex_size=40,
    vertex_color="white",
    vertex_label=graph1.vs['name'],
    vertex_frame_color="black",
)

# Plot graph3
ig.plot(
    graph3,
    layout=layout_fixed,
    target=ax[1,0],
    edge_width=edge_widths3,
    vertex_size=40,
    vertex_color="white",
    vertex_label=graph1.vs['name'],
    vertex_frame_color="black",
)

ig.plot(
    combine,
    layout=layout_fixed,
    target=ax[1,1],
    edge_width=edge_widths_C,
    edge_color=edge_colors,
    vertex_size=40,
    vertex_color="white",
    vertex_label=graph1.vs['name'],
    vertex_frame_color="black",
)



plt.tight_layout()
plt.savefig("/path/to/via/trajectory_difference.pdf", format="pdf", bbox_inches='tight')
plt.show()
plt.close()


