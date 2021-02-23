# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:26:20 2021

@author: chung
"""
import pandas as pd
import os

scRNA_df = pd.read_csv('GSE141834_scRNAseq_seuratV3_normalized.txt', sep='\t')
RNA_df = pd.read_csv('GSE141834_bulkRNAseq_normalized_counts.txt', sep='\t')

scRNA_gene_name = scRNA_df.index.tolist()
RNA_gene_name = RNA_df.index.tolist()
co_gene_name = list(set(scRNA_gene_name) & set(RNA_gene_name))

del_scRNA_gene_name=list(set(scRNA_gene_name) - set(co_gene_name))
del_RNA_gene_name=list(set(RNA_gene_name) - set(co_gene_name))

scRNA_df.drop(del_scRNA_gene_name, inplace=True)
RNA_df.drop(del_RNA_gene_name, inplace=True)

directory = 'output/'

if not os.path.exists(directory):
    os.makedirs(directory)

scRNA_df.to_csv('output/scRNA_co_gene_filtered.csv')
RNA_df.to_csv('output/RNA_co_gene_filtered.csv')