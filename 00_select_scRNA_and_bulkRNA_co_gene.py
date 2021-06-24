# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 17:26:20 2021

@author: chung
"""
import pandas as pd
import os
import ensembles as en

scRNA_filename = 'scRNA.csv'
bulkRNA_filename = 'RNA.csv'

scRNA_df = pd.read_csv(scRNA_filename, sep=',', index_col=0)
RNA_df = pd.read_csv(bulkRNA_filename, sep=',', index_col=0)

co_gene_name = en.get_common_gene(scRNA_df, RNA_df)

RNA_df=RNA_df.loc[co_gene_name]
scRNA_df=scRNA_df.loc[co_gene_name]

directory = 'output/'

if not os.path.exists(directory):
    os.makedirs(directory)

scRNA_df.to_csv('output/scRNA_co_gene_filtered.csv')
RNA_df.to_csv('output/RNA_co_gene_filtered.csv')