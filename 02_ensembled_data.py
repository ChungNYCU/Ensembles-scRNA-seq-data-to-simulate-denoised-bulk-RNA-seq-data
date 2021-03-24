# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 19:54:22 2021

@author: chung
"""

import pandas as pd
import scanpy as sc

def resampling(df, times):
    x = list(df.tolist())
    x[:] = [i * times for i in x]
    return x

def resample_flow(i):
    Nlist = list(adata.obsp['connectivities'])[i].toarray().flatten() #neighbor list
    temp=[]#near neighbor list
    #select near neighbor
    for j in range(len(Nlist)):
        if(Nlist[j]!=0):
            temp.append([Nlist[j],j])
    temp.sort(reverse=True)
    temp = temp[:4]
            
    simu_bulk = resampling(scRNA_df[sample_list[i]],5) # 0=self
    for k in range(len(temp)):
        simu_bulk = list(map(lambda x,y: x + y, simu_bulk,resampling(scRNA_df[sample_list[temp[k][1]]], 5)))
        
    return simu_bulk

if __name__ == '__main__':
    file_path = 'write/scRNA_data.h5ad'
    origin_scRNA_path = 'dataset/GSE141834_scRNAseq_seuratV3_normalized.txt'

    adata = sc.read_h5ad(file_path)
    scRNA_df = pd.read_csv(origin_scRNA_path, sep='\t', index_col=0)
    
    sample_list = list(adata.obs.index)
    gene_list = adata.var.index.tolist()
    
    df_simu_bulk = pd.DataFrame()
    for i in range(len(sample_list)):
        #str(sample_list[i][4:6]) this string is depend on dataset
        df_simu_bulk['ensembled_Bulk_'+str(i)+'_'+str(sample_list[i][4:6])] = resample_flow(i) 
    
    df_simu_bulk.index = scRNA_df.index.values
    df_simu_bulk.to_csv('output/ensembled_data.csv')