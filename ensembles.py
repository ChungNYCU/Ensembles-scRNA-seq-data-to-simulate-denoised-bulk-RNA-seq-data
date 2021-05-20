# -*- coding: utf-8 -*-
"""
Created on Tue May  4 12:26:37 2021

@author: chung
"""

import pandas as pd
import scanpy as sc
import copy
from scipy.stats import norm

def get_adata(filename):
    """Get processed data. Optional. 
    Can use scanpy for pre-processing and get adata.

    Keyword arguments:
    filename -- filename or filepath of csv or txt file
    """
    adata = sc.read_csv(filename, ',', first_column_names=True)
    adata = adata.T
    adata.var_names_make_unique()
    
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    sc.tl.pca(adata, svd_solver='arpack')
    
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    
    sc.tl.umap(adata)
    return adata

def co_gene(scRNA, BulkRNA) -> 'list':
    """Get common gene list between single-cell and Bulk data. Optional. 
    Esalier to data visualization and compare single-cell and bulk data.

    Keyword arguments:
    scRNA -- DataFrame of single-cell RNA-seq data, row = gene, column = sample.
    BulkRNA -- DataFrame of Bulk RNA-seq data, row = gene, column = sample.
    """
    scRNA = scRNA.index.tolist()
    BulkRNA = BulkRNA.index.tolist()
    co_gene_list = list(set(scRNA) & set(BulkRNA))
    return co_gene_list

def add_labels(adata, label_name:'str', labels:'list'):
    """Add labels to adata. Optional. 
    Esalier to data visualization and compare single-cell and bulk data.

    Keyword arguments:
    adata -- adata is processed by scanpy.
    label_name -- e.g. 'celltype'
    labels -- For every samples, need to have a corresponding label.
    """
    adata.obs[label_name]=labels

def get_nearest_neighbors(adata, sample:'int', k:'int') -> 'list':
    """Get nearest neighbors list.

    Keyword arguments:
    adata -- adata is processed by scanpy.
    sample -- The column index of dataset, correspond a sample.
    k -- Select k nearest neighbors.
    """
    nearest_neighbors = list(adata.obsp['distances'])[sample].toarray().flatten()
    knn_list=[]
    knn_list.append([sample,0])
    for j in range(len(nearest_neighbors)):
        if(nearest_neighbors[j]!=0):
            knn_list.append([j, nearest_neighbors[j]])
    knn_list.sort(key = lambda s: s[1])
    knn_list = knn_list[:k] 
    return knn_list

def get_sampling_p(knn_list, sampling_type:'int' = 0) -> 'list':
    """Get duplicate sampling probability for each neighbors.

    Keyword arguments:
    knn_list -- List of nearest neighbors.
    sampling_type -- 0 = Gaussian distribution sampling, 1 = Average sampling.
    """
    p=[]
    knn_copy = copy.deepcopy(knn_list)
    if(sampling_type==0):
        mean = 0
        variance = sum([(x[1]-mean)**2 for x in knn_copy])/len(knn_copy)
        std = variance ** 0.5
    
        for i in knn_copy:
            i.append(norm(mean, std).pdf(i[1]))
            
        normalize_value = 0
        for i in knn_copy:
            normalize_value += i[2]
        normalize_value = 1/normalize_value
        for i in knn_copy:
            i[2]*=normalize_value
        for i in knn_copy:
            p.append([i[0], i[2]])
    if(sampling_type==1):
        times=1/len(knn_copy)
        for i in knn_copy:
            p.append([i[0], times])
    return p

def resampling(df, times:'float') -> 'list':
    """Output = Sample * times/probability.

    Keyword arguments:
    df -- A sample column.
    times -- times or probability for duplicate sampling.
    """
    x = list(df.tolist())
    x[:] = [i * times for i in x]
    return x

def ensembles_data(df, knn_list, sample:'int', k:'int') -> 'list':
    """Ensembles a single-cell sample to ensembled sample (sample -> sample).

    Keyword arguments:
    df -- DataFrame of single-cell RNA-seq data, row = gene, column = sample.
    knn_list -- List of nearest neighbors.
    sample -- The column index of dataset, correspond a sample.
    k -- Select k nearest neighbors.
    """
    p = get_sampling_p(knn_list)
    ensembled_data = resampling(df.iloc[:, p[0][0]], 0) # 0=self
    for i in p:
        ensembled_data = list(map(lambda x,y: x + y, ensembled_data, resampling(df.iloc[:, i[0]], i[1])))
    return ensembled_data

def ensembles_df(df, adata, k):
    """Ensembles single-cell data to ensembled data (DataFrame -> DataFrame).

    Keyword arguments:
    df -- DataFrame of single-cell RNA-seq data, row = gene, column = sample.
    adata -- adata is processed by scanpy.
    k -- Select k nearest neighbors.
    """
    ensembled_df = pd.DataFrame()
    for i in range(len(df.columns)):
        knn_list = get_nearest_neighbors(adata, i, k)
        ensembled_df[i] = ensembles_data(df, knn_list, i, k)
    ensembled_df.index = df.index.values
    return ensembled_df

def ratio_(df, r):
    """Read count different between single-cell and bulk RNA-seq data, get ratio r through below function.
    r = (Library size (Bulk RNA) / sample size (Bulk RNA)) /  (Library size (scRNA) / sample size (scRNA))
    Ensembled data * ratio r can be viewed as Bulk RNA-seq data.

    Keyword arguments:
    df -- DataFrame of ensembled data, row = gene, column = sample.
    r -- Read count ratio between Bulk RNA-seq and scRNA-seq
    """
    df = df.multiply(r)
    return df