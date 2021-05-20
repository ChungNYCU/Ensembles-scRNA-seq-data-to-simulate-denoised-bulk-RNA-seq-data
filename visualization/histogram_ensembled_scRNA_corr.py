#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 15:45:11 2021

@author: chung
"""
import seaborn as sns
import pandas as pd
import numpy as np

filename = 'BCO_18_all.csv'
data = pd.read_csv('dataset/' + filename, sep=',', index_col=0)

data = data.values
data = data.T
corr = np.corrcoef(data)
temp=corr[0:18]
temp=temp.T
temp=temp[18:]
res=[]
for i in temp:
    res.append(sum(i)/len(i))
# res=np.reshape(res, (2, 18))
df = pd.DataFrame()
df['INDEX']=res
df['Data type']=['ensembled']*18+['scRNA']*18
sns.histplot(data=df, x='INDEX', hue='Data type')


# scRNA_zero_count=scRNA_df.isin([0]).sum().tolist()
# ensembled_zero_count=ensembled_df.isin([0]).sum().tolist()
# df = pd.DataFrame()
# a=scRNA_zero_count+ensembled_zero_count
# b=['scRNA']*2400+['ensembled']*2400
# df['Missing value']=a
# df['Data type']=b

# sns.histplot(data=df, x='Missing value', hue='Data type')

# sns_plot=sns.histplot(data=df, x='Missing value', hue='Data type')
# sns_plot.figure.savefig('missing value histogram.png', dpi=600)