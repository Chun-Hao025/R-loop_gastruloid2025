#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scanpy as sc
import pandas as pd
import scipy
from math import sqrt
from scipy.spatial.distance import jensenshannon
import matplotlib.pyplot as plt
import seaborn as sns

import gc
from ctxcore.rnkdb import opendb, RankingDatabase  # for loading ranking databases
import os
import pandas as pd
import ast

random_seed = 0
np.random.seed(random_seed)



# Function to extract the first (target, weight) tuple from a cell value.
def extract_first_target_info(x):
    try:
        # Convert string representation to an actual list.
        parsed = ast.literal_eval(x)
        # Check if it's a list and non-empty.
        if isinstance(parsed, list) and len(parsed) > 0 and isinstance(parsed[0], tuple):
            target, weight = parsed[0]
            return pd.Series([target, weight])
        else:
            return pd.Series([None, None])
    except Exception:
        return pd.Series([None, None])


cell=['cluster0', 'cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5', 'cluster6', 'cluster7', 'cluster8',
      'cluster9', 'cluster10', 'cluster11', 'cluster12', 'cluster13', 'cluster14', 'cluster15', "total"]

indir='/path/to/GRN/'

for clu in cell:
    filename=os.path.join(indir, f'{clu}_reg1.csv')
    df = pd.read_csv(filename, header=[0, 1], index_col=[0, 1])
    df_long = df.reset_index()  # This makes the index columns ('TF', 'MotifID') regular columns.
    print(df_long.columns.tolist())
    df_long[['target', 'weight']] = df_long[('Enrichment', 'TargetGenes')].apply(extract_first_target_info)
    result = df_long[['TF', 'target', 'weight']]
    result.columns = [col[0] for col in result.columns.tolist()]
    result=result.drop_duplicates(subset=['TF', 'target', 'weight'])
    filename1=os.path.join(indir, f'{clu}_reg_TF_target.csv')
    result.to_csv(filename1)



def regulon_specificity_scores1(auc_mtx, cell_type_series):
    """
    Calculates the Regulon Specificity Scores (RSS). [doi: 10.1016/j.celrep.2018.10.045]

    :param auc_mtx: DataFrame with AUC values for all cells and regulons (n_cells x n_regulons).
    :param cell_type_series: Series with cell identifiers as index and cell type labels as values.
    :return: DataFrame with RSS values (cell type x regulon).
    """

    # Ensure inputs are valid pandas objects
    if not isinstance(auc_mtx, pd.DataFrame):
        raise TypeError("auc_mtx must be a pandas DataFrame")
    if not isinstance(cell_type_series, pd.Series):
        raise TypeError("cell_type_series must be a pandas Series")

    # Ensure alignment between cell_type_series and auc_mtx
    if not auc_mtx.index.equals(cell_type_series.index):
        raise ValueError("The index of auc_mtx must match the index of cell_type_series")

    # Get unique cell types and regulons
    cell_types = cell_type_series.unique()
    n_types = len(cell_types)
    regulons = auc_mtx.columns
    n_regulons = len(regulons)

    # Initialize the RSS values matrix
    rss_values = np.zeros((n_types, n_regulons))

    # Normalize AUC matrix once for efficiency
    auc_mtx_norm = auc_mtx.div(auc_mtx.sum(axis=0), axis=1)

    def rss(aucs, labels):
        """
        Calculate the RSS value for a single regulon and cell type.
        :param aucs: Normalized AUC values for a single regulon.
        :param labels: Binary array indicating cell type membership.
        :return: RSS value.
        """
        labels_norm = labels / labels.sum()
        return 1.0 - jensenshannon(aucs, labels_norm)

    # Compute RSS values
    for ridx, cell_type in enumerate(cell_types):
        # Create binary labels for the current cell type
        labels = (cell_type_series == cell_type).astype(int)

        for cidx, regulon_name in enumerate(regulons):
            # Compute RSS value for each regulon and cell type
            rss_values[ridx, cidx] = rss(auc_mtx_norm[regulon_name], labels)

    # Return the RSS values as a DataFrame
    return pd.DataFrame(data=rss_values, index=cell_types, columns=regulons)

cell_meta = pd.read_csv("/path/to/via/via_meta_new1.csv", index_col=0) ## with annotation

# Define the criteria
criteria = cell_meta['treatment'] == 'control'

# Use np.where to get the locations (index positions) of all rows that match the criteria
matching_indices = np.where(criteria)[0]
cell_meta1=cell_meta.iloc[matching_indices]
aux_mtx=pd.read_csv("/path/to/GRN/total_auc1.csv", index_col=0)
rss_cellType = regulon_specificity_scores1(aux_mtx, cell_meta1['ann_KNN40'] )
rss_cellType
rss_cellType.to_csv("/path/to/GRN/total_rss.csv", index_label=True)