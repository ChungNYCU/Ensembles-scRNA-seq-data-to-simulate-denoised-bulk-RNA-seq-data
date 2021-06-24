# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 19:54:22 2021

@author: chung
"""

import pandas as pd
import scanpy as sc
import ensembles as en

if __name__ == '__main__':
    file_path = 'write/scRNA_origin.h5ad'

    adata = sc.read_h5ad(file_path)
    scRNA_df = pd.read_csv('scRNA.csv', sep=',', index_col=0)
    RNA_df = pd.read_csv('RNA.csv', sep=',', index_col=0)
    
    ensembled = en.ensembles_df(scRNA_df, adata, 5)
    r = en.get_ratio_r(scRNA_df, RNA_df)
    ensembled = en.ratio_(ensembled, r)