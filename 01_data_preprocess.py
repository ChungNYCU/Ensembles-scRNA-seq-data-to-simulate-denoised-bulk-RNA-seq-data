# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 19:17:28 2021

@author: chung
"""

import scanpy as sc
import os

directory = 'write/'

if not os.path.exists(directory):
    os.makedirs(directory)

adata = sc.read_csv('output/ensembled_data.csv', ',', first_column_names=True)
adata = adata.T
adata.var_names_make_unique()
results_file = 'write/ensembled_data.h5ad'

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# time = []
# for i in range(400): 
#     time.append('00')
# for i in range(400): 
#     time.append('01')
# for i in range(400): 
#     time.append('02')
# for i in range(400): 
#     time.append('04')
# for i in range(400): 
#     time.append('08')
# for i in range(400): 
#     time.append('18')
# typea = []
# for i in range(2400):
#     typea.append('Ensembled')

# adata.obs['time']=time
# adata.obs['type']=typea

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.umap(adata)
umap3 = sc.tl.umap(adata, n_components=3, copy=True)
# adata.write(results_file)