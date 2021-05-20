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
    origin_scRNA_path = 'dataset/GSE141834_scRNAseq_seuratV3_normalized.txt'

    adata = sc.read_h5ad(file_path)
    scRNA_df = pd.read_csv(origin_scRNA_path, sep='\t', index_col=0)
    
    ensembled = en.ensembles_df(scRNA_df, adata, 5)
    ensembled = en.ratio_(ensembled, 25)