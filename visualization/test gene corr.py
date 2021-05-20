#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 09:37:32 2021

@author: chung
"""
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

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

a=scRNA_df
a=a.sort_index()
a=a.iloc[:50]
a=a.T
# a=a.iloc[:500]

b=ensembled_df
b=b.sort_index()
b=b.iloc[:50]
b=b.T
# b=b.iloc[:500]

c=bulk_df
c=c.sort_index()
c=c.iloc[:50]
c=c.T

# b=b.iloc[:500]

sc_heatmap=a.corr()
en_heatmap=b.corr()
bu_heatmap=c.corr()
# sns.heatmap(sc_heatmap)

sns_plot=sns.heatmap(sc_heatmap)
sns_plot.figure.savefig('sc_heatmap_gene_gene.png', dpi=600)
sns_plot=sns.heatmap(en_heatmap)
sns_plot.figure.savefig('en_heatmap_gene_gene.png', dpi=600)
sns_plot=sns.heatmap(bu_heatmap)
sns_plot.figure.savefig('bu_heatmap_gene_gene.png', dpi=600)
