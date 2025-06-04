#!/usr/bin/env python
# coding: utf-8

import os
import anndata
import sys
import warnings
from scipy.ndimage import convolve
import scipy.sparse as sp
from scipy import io
import gc
warnings.simplefilter("ignore", category=UserWarning)
import re
import numpy as np
import pandas as pd
import itertools
import scanpy as sc


os.chdir('/path/to/trajectory/')
###build up anndata object
X=io.mmread("ref_RH_count.mtx")
X=X.transpose().tocsr()
cell_meta=pd.read_csv("via_meta_new1.csv", low_memory=False, index_col=0)

with open("gene_name.csv", 'r') as f:
    gene_names = f.read().splitlines()

var=pd.DataFrame(data={}, index=gene_names)
velo = anndata.AnnData(X=X.astype(np.float32), obs=cell_meta, var=var)
velo.var.index.name='Gene'
velo.obs.index.name='CellID'

velo.write("/path/to/trajectory/ref_RH.h5ad", compression="gzip")

criteria = (velo.obs['treatment'] == "control")
velo_sub=velo[criteria,:]

# save a copy of the raw data
velo_sub.layers['counts']=velo_sub.X.copy()

# Total-count normalize (library-size correct) to 10,000 reads/cell
sc.pp.normalize_total(velo_sub)

# log transform the data.
sc.pp.log1p(velo_sub)

# identify highly variable genes.
sc.pp.highly_variable_genes(velo_sub, n_top_genes=5000, batch_key= 'replicate')


###prepare the input for SCENIC analysis to construct the initial network across all timepoint of control data.

celltype=['Caudal_Neuroectoderm1', 'Definitive_Endoderm_Guts', 'Endothelium', 'Epiblasts', 
              'Naive_Pluripotent_Cells', 'Nascent_Mesoderm', 'Neuromesodermal_Progenitors1', 
              'Neuromesodermal_Progenitors2', 'Neuromesodermal_Progenitors3', 'Neuron-like_Cells', 
              'Paraxial_Mesoderm_A', 'Paraxial_Mesoderm_B', 'Placodal_area', 'Primitive_Streak', 
              'Spinal_Cord', 'Splanchnic_Mesoderm']
cluster=["cluster0", "cluster1", "cluster2","cluster3","cluster4","cluster5","cluster6","cluster7","cluster8","cluster9","cluster10",
          "cluster11","cluster12","cluster13","cluster14","cluster15"]

for cell, clus in zip(celltype, cluster):
    criteria2 = (ref.obs['ann_KNN40'] == cell)
    tmp=ref[criteria2,:]
    tmp.write_loom(f'/path/to/GRN/{clus}.loom')
    del tmp
    gc.collect()

###total data to get the global variance for initial GRN construct through SCENIC
# keep only highly variable genes:
TF=pd.read_table('/path/to/GRN/TF_filter2.txt',header=None)
TF_genes = set(TF.iloc[:, 0])
keep_genes = velo_sub.var.index[velo_sub.var['highly_variable']]
combined_genes = TF_genes.union(set(keep_genes))
filtered_genes = [gene for gene in velo_sub.var_names if gene in combined_genes]
velo_sub = velo_sub[:, filtered_genes]

velo_sub.X=velo_sub.layers['counts']
velo_sub.write_loom(f'/path/to/GRN/total.loom')

