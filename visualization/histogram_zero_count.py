#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 15:45:11 2021

@author: chung
"""
import seaborn as sns
import pandas as pd
scRNA_filename = 'GSE141834_scRNAseq_seuratV3_normalized.txt'
ensembled_filename = 'ensembled_data_origin_v2_k.csv'


scRNA_df = pd.read_csv('dataset/' + scRNA_filename, sep='\t')
ensembled_df = pd.read_csv('dataset/' + ensembled_filename, sep=',',index_col=0)

scRNA_zero_count=scRNA_df.isin([0]).sum().tolist()
ensembled_zero_count=ensembled_df.isin([0]).sum().tolist()
df = pd.DataFrame()
a=scRNA_zero_count+ensembled_zero_count
b=['scRNA']*2400+['ensembled']*2400
df['Missing value']=a
df['Data type']=b

sns.histplot(data=df, x='Missing value', hue='Data type')

sns_plot=sns.histplot(data=df, x='Missing value', hue='Data type')
sns_plot.figure.savefig('missing value histogram.png', dpi=600)