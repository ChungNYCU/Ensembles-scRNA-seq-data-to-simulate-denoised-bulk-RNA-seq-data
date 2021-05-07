# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:26:20 2021

@author: chung
"""
import pandas as pd
import os
import ensembles as en

scRNA_filename = 'GSE141834_scRNAseq_seuratV3_normalized.txt'
bulkRNA_filename = 'GSE141834_bulkRNAseq_normalized_counts.txt'

scRNA_df = pd.read_csv('dataset/' + scRNA_filename, sep='\t')
RNA_df = pd.read_csv('dataset/' + bulkRNA_filename, sep='\t')

co_gene_name = en.get_common_gene(scRNA_df, RNA_df)

RNA_df=RNA_df.loc[co_gene_name]
scRNA_df=scRNA_df.loc[co_gene_name]

directory = 'output/'

if not os.path.exists(directory):
    os.makedirs(directory)

scRNA_df.to_csv('output/scRNA_co_gene_filtered.csv')
RNA_df.to_csv('output/RNA_co_gene_filtered.csv')