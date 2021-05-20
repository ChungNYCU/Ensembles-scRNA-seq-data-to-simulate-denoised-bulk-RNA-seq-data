# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 19:17:28 2021

@author: chung
"""

import scanpy as sc
import os
import ensembles as en

directory = 'write/'

if not os.path.exists(directory):
    os.makedirs(directory)

adata = en.get_adata('scRNA.csv', ',')

results_file = 'write/scRNA_origin.h5ad'

time = []
for i in range(400): 
    time.append('00')
for i in range(400): 
    time.append('01')
for i in range(400): 
    time.append('02')
for i in range(400): 
    time.append('04')
for i in range(400): 
    time.append('08')
for i in range(400): 
    time.append('18')

typea = []
for i in range(2400):
    typea.append('Ensembled')

en.add_labels(adata, 'time', time)
en.add_labels(adata, 'type', typea)


sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

sc.tl.pca(adata, svd_solver='arpack')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.umap(adata)
umap3 = sc.tl.umap(adata, n_components=3, copy=True)
adata.write(results_file)