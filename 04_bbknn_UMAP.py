# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 18:03:03 2020

@author: chung
"""
import scanpy as sc
import pandas as pd
import seaborn as sns
# from bbknn import bbknn

sc.settings.verbosity = 1             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(3, 3), facecolor='white')

adata_ref = sc.read_h5ad('write/ensembled_data.h5ad') #
adata = sc.read_h5ad('write/real_bulk.h5ad') #18 samples bulk dataset

sc.pl.umap(adata_ref, color='time',size=200)
sc.pl.umap(adata, color='time')

var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

sc.pp.pca(adata_ref)
sc.pp.neighbors(adata_ref)
sc.tl.umap(adata_ref)

# sc.tl.ingest(adata_ref, adata,  obs=None, embedding_method=('umap','pca'), labeling_method='knn', neighbors_key=None, inplace=True,)

adata.uns['time_colors'] = adata_ref.uns['time_colors']

adata_concat = adata_ref.concatenate(adata, batch_categories=['Ensembled', 'Real bulk'])

adata_concat.obs.type = adata_concat.obs.type.astype('category')
# adata_concat.obs.type.cat.reorder_categories(adata_ref.obs.type.cat.categories, inplace=True)  # fix category ordering
adata_concat.uns['time_colors'] = adata_ref.uns['time_colors']  # fix category colors

sc.pl.umap(adata_concat, color=['time', 'batch'])

# adata_concat = adata_ref.concatenate(adata, batch_categories=['Ensembled', 'Real bulk'])
# sc.tl.pca(adata_concat)
# sc.external.pp.bbknn(adata_concat, batch_key='batch')
# sc.tl.umap(adata_concat)
# sc.pl.umap(adata_concat, color=['batch', 'time'])
# sc.pl.umap(adata_concat, color=['time', 'batch'])