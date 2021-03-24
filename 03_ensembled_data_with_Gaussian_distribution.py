# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 21:09:52 2020

@author: chung
"""

# import numpy as np
# from sklearn.datasets import load_iris, load_digits
# from sklearn.model_selection import train_test_split
# import matplotlib.pyplot as plt
# import seaborn as sns
import pandas as pd
import umap
import scanpy as sc
from scipy.stats import norm


def resampling(df, times):
    x = list(df.tolist())
    x[:] = [i * times * 100 for i in x]
    return x

def resample_flow(count):
    Nlist = list(adata.obsp['distances'])[count].toarray().flatten() #neighbor list
    temp=[]#near neighbor list
    #select near neighbor
    
    for j in range(len(Nlist)):
        if(Nlist[j]!=0):
            temp.append([j, Nlist[j]*0.2])
    temp.sort(key = lambda s: s[1])
    temp = temp[:5] 
    mean = 0
    temp.append([count,mean])
    
    for i in temp:
        i.append(norm(mean, sigmas[count]).pdf(i[1]))
    normalize_value = 0
    for i in temp:
        normalize_value += i[2]
    normalize_value = 1/normalize_value
    for i in temp:
        if(i[2]<0.000001):
            i[2]=0
        i[2]*=normalize_value

            
    simu_bulk = resampling(scRNA_df[sample_list[temp[-1][0]]], temp[-1][2]) # 0=self
    for k in range(len(temp)-1):
        simu_bulk = list(map(lambda x,y: x + y, simu_bulk,resampling(scRNA_df[sample_list[temp[k][0]]], temp[k][2])))

    return simu_bulk

# def normal_distribution(samples_distance, sigma):
#     for i in range(len(samples_distance)):
#         rv = norm()
#         samples_distance[i][0]=
        
#     return samples_distance


if __name__ == '__main__':
    file_path = 'write/scRNA_origin.h5ad'
    origin_scRNA_path = 'dataset/GSE141834_scRNAseq_seuratV3_normalized.txt'

    adata = sc.read_h5ad(file_path)
    scRNA_df = pd.read_csv(origin_scRNA_path, sep='\t', index_col=0)
    
    sample_list = list(adata.obs.index)
    gene_list = adata.var.index.tolist()
    
    inputDF = scRNA_df.T
    
    reducer = umap.UMAP(random_state=42)
    reducer.fit(inputDF)
    sigmas = reducer._sigmas
    
    df_simu_bulk = pd.DataFrame()
    
    for count in range(len(sample_list)):
        #str(sample_list[i][4:6]) this string is depend on dataset
        print(sample_list[count])
        df_simu_bulk['ensembled_Bulk_'+str(count)+'_'+str(sample_list[count][4:6])] = resample_flow(count) 
    
    df_simu_bulk.index = scRNA_df.index.values
    df_simu_bulk.to_csv('output/ensembled_data_origin_v2.csv')




# embedding = reducer.transform(inputDF)
# # Verify that the result of calling transform is
# # idenitical to accessing the embedding_ attribute
# assert(np.all(embedding == reducer.embedding_))
# embedding.shape

# color = []
# color.append(0)
# for i in range(1,3022):
#     color.append(1)

# plt.scatter(embedding[:, 0], embedding[:, 1], c=color, cmap='Spectral', s=1)
# plt.gca().set_aspect('equal', 'datalim')
# plt.colorbar(boundaries=np.arange(11)-0.5).set_ticks(np.arange(10))
# plt.title('scRNA dataset', fontsize=12)




