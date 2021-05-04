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

def co_gene(scRNA, Bulk) -> 'list':
    scRNA = scRNA.index.tolist()
    Bulk = Bulk.index.tolist()
    co_gene_list = list(set(scRNA) & set(Bulk))
    return co_gene_list

def add_labels(adata, label_name:'str', labels:'list'):
    adata.obs[label_name]=labels
    return adata

def get_nearest_neighbors(adata, sample:'int', k:'int') -> 'list':
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
    x = list(df.tolist())
    x[:] = [i * times for i in x]
    return x

def ensembles_data(df, knn_list, sample:'int', k:'int') -> 'list':
    p = get_sampling_p(knn_list)
    # ensembled_data = []*20
    ensembled_data = resampling(df.iloc[:, p[0][0]], 0) # 0=self
    for i in p:
        ensembled_data = list(map(lambda x,y: x + y, ensembled_data, resampling(df.iloc[:, i[0]], i[1])))
    return ensembled_data

def ensembles_df(df, adata, k):
    ensembled_df = pd.DataFrame()
    for i in range(len(df.columns)):
        print(i)
        knn_list = get_nearest_neighbors(adata, i, k)
        ensembled_df[i] = ensembles_data(df, knn_list, i, k)
    ensembled_df.index = df.index.values
    return ensembled_df

def ratio_(df, r):
    df = df.multiply(r)
    return df