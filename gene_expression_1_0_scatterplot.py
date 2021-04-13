#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 17:14:40 2021

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

scRNA_not_zero = [len(scRNA_df)-x for x in scRNA_zero_count]
ensembled_not_zero = [len(ensembled_df)-x for x in ensembled_zero_count]

res = scRNA_not_zero+ensembled_not_zero
data_type=['scRNA']*2400+['ensembled']*2400
Y=[]
for i in range(2400):
    Y.append(i)
Y=Y+Y


df=pd.DataFrame()
df['Gene Expression']=res
df['Data type']=data_type
df['INDEX']=Y

sns.scatterplot(data=df, x='INDEX', y='Gene Expression', hue='Data type', s=5)

sns_plot=sns.scatterplot(data=df, x='INDEX', y='Gene Expression', hue='Data type', s=5)
sns_plot.figure.savefig('gene 1 0 scatter.png', dpi=600)

df = pd.DataFrame()
a=scRNA_not_zero+ensembled_not_zero
b=['scRNA']*2400+['ensembled']*2400
df['Expression']=a
df['Data type']=b

sns.histplot(data=df, x='Expression', hue='Data type')

sns_plot=sns.scatterplot(data=df, x='INDEX', y='Gene Expression', hue='Data type', s=5)
sns_plot.figure.savefig('gene 1 0 scatter.png', dpi=600)