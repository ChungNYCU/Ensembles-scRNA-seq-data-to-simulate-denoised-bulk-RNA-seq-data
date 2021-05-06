#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 15:45:11 2021

@author: chung
"""
import seaborn as sns
import pandas as pd
import random
scRNA_filename = 'GSE141834_scRNAseq_seuratV3_normalized.txt'
bulk_filename = 'GSE141834_bulkRNAseq_normalized_counts.txt'
ensembled_filename = 'ensembled_data_origin_v2_k.csv'


scRNA_df = pd.read_csv('dataset/' + scRNA_filename, sep='\t')
bulk_df = pd.read_csv('dataset/' + bulk_filename, sep='\t')
ensembled_df = pd.read_csv('dataset/' + ensembled_filename, sep=',',index_col=0)

scRNA_read_count=[]
bulk_read_count=[]
ensembled_read_count=[]

for i in scRNA_df.columns:
    scRNA_read_count.append(scRNA_df[i].sum())
    
for i in bulk_df.columns:
    bulk_read_count.append(bulk_df[i].sum())
    
for i in ensembled_df.columns:
    ensembled_read_count.append(ensembled_df[i].sum())
    
sample_of_random_scRNA = random.sample(scRNA_read_count, 18)
sample_of_random_ensemble = random.sample(ensembled_read_count, 72)
    

    

df = pd.DataFrame()
a=sample_of_random_scRNA+sample_of_random_ensemble+bulk_read_count
b=['scRNA']*18+['ensembled']*72+['Bulk']*18
df['Read Count']=a
df['Data type']=b

sns.histplot(data=df, x='Read Count', hue='Data type')

sns_plot=sns.histplot(data=df, x='Read Count', hue='Data type')
sns_plot.figure.savefig('read count histogram.png', dpi=600)
