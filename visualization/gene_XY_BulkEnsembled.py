#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 16:22:40 2021

@author: chung
"""

import seaborn as sns
import pandas as pd

def plot(x,y):
    a=bulk_df.iloc[:,x].sort_index()
    b=ensembled_df.iloc[:,y].sort_index()
    c=scRNA_df.iloc[:,y].sort_index()
    datatype=pd.Series(['ensembled']*len(b)+['scRNA']*len(c))
    
    a=pd.concat([a,a])
    b=pd.concat([b,c])
    datatype.index=a.index
    
    frames = [a, b, datatype]
    result = pd.concat(frames, axis=1)
    result = result[result[a.name] >= 5] 
    result = result[result[a.name] <= 18] 
    # result = result[result[0] <= 150]
    
    result = result.rename(columns={a.name: "Bulk RNA", 0: "ensembled and scRNA", 1:"Data Type"})

    sns.scatterplot(data=result, x="Bulk RNA", y="ensembled and scRNA", hue='Data Type', s=1)
    sns_plot=sns.scatterplot(data=result, x="Bulk RNA", y="ensembled and scRNA", hue='Data Type', s=1)
    sns_plot.figure.savefig('gene_gene_expression3.png', dpi=600)
    

scRNA_filename = 'GSE141834_scRNAseq_seuratV3_normalized.txt'
bulk_filename = 'GSE141834_bulkRNAseq_normalized_counts.txt'
ensembled_filename = 'ensembled_data_origin_v2_k.csv'

scRNA_df = pd.read_csv('dataset/' + scRNA_filename, sep='\t')
bulk_df = pd.read_csv('dataset/' + bulk_filename, sep='\t')
ensembled_df = pd.read_csv('dataset/' + ensembled_filename, sep=',',index_col=0)

scRNA_gene_name = scRNA_df.index.tolist()
RNA_gene_name = bulk_df.index.tolist()
co_gene_name = list(set(scRNA_gene_name) & set(RNA_gene_name))

del_scRNA_gene_name=list(set(scRNA_gene_name) - set(co_gene_name))
del_RNA_gene_name=list(set(RNA_gene_name) - set(co_gene_name))

scRNA_df.drop(del_scRNA_gene_name, inplace=True)
bulk_df.drop(del_RNA_gene_name, inplace=True)
ensembled_df.drop(del_scRNA_gene_name, inplace=True)
#x<18 y<2399
plot(17,2399)




# , hue="time", style="time"