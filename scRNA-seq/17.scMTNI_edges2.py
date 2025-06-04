import pandas as pd
from pathlib import Path
import numpy as np
import os

### This script was used to recall the correlation coefficient from scMTNI analysis and we choose the max value as the represent value
def aggregate_edge_coeffs(network_list_file, 
                          consensus_edges_file, 
                          output_file):
    """
    1) Reads the list of per‐bootstrap network files (TF, target, coef).
    2) Computes, for each (TF, target), the coefficient with the largest absolute value
       *and* preserves its original sign.
    3) Merges that back onto the consensus‐edges table (which has at least TF, target, freq).
    4) Writes out consensus_edges + signed_coef as a 4th column.
    """

    # 1) load the list of network files
    paths = [Path(p) for p in open(network_list_file).read().splitlines() if p]
    
    # We'll collect each file's edges into one big DataFrame
    dfs = []
    for p in paths:
        df = pd.read_csv(p, sep='\t', header=None, names=['TF','target','coef'])
        dfs.append(df)
    all_edges = pd.concat(dfs, ignore_index=True)

    # 2) find the row per (TF,target) with max |coef|
    #    we do this by computing abs_coef, then sorting & dropping duplicates
    all_edges['abs_coef'] = all_edges['coef'].abs()
    best = (
        all_edges
        .sort_values(['TF','target','abs_coef'], ascending=[True,True,False])
        .drop_duplicates(subset=['TF','target'], keep='first')
        .loc[:, ['TF','target','coef']]
        .rename(columns={'coef':'best_coef'})
    )

    # 3) read the consensus‐edge table and merge
    cons = pd.read_csv(consensus_edges_file, sep='\t', header=None, names=['TF','target','coef'])
    merged = cons.merge(best, on=['TF','target'], how='left')
    # missing edges get NaN → fill with 0
    merged['best_coef'] = merged['best_coef'].fillna(0)

    # 4) write it out
    merged.to_csv(output_file, sep='\t', index=False, header=False)

# example usage:

clusters=pd.read_table('/path/to/GRN/scMTNI/cluster.txt', header=None)[0].tolist()
indir='/path/to/GRN/scMTNI/scMTNI_control' ### change to the all control, RHKI control or RHKI Dox directory
for cell in clusters:
    cluster_path=os.path.join(indir, f'control_N/{cell}')
    aggregate_edge_coeffs(
        network_list_file   = f'{cluster_path}/network_files.txt',
        consensus_edges_file= f'{cluster_path}/net_alledge.txt',
        output_file         = f'{cluster_path}/net_alledge_origcoef.txt'
    )