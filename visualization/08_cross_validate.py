# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 14:23:51 2021

@author: chung
"""
from sklearn import datasets, linear_model
from sklearn.model_selection import cross_validate
from sklearn.metrics import make_scorer
from sklearn.metrics import confusion_matrix
from sklearn.svm import LinearSVC
import pandas as pd

# diabetes = datasets.load_diabetes()
ensembled = pd.read_csv('dataset/ensembled_data_origin_v2.csv', sep=',', index_col=0)
bulk = pd.read_csv('dataset/GSE141834_bulkRNAseq_normalized_counts.txt', sep='\t')
X = bulk
X = X.T
y = []
for i in range(3):
    y.append('00')
    y.append('01')
    y.append('02')
    y.append('04')
    y.append('08')
    y.append('18')
y.sort()
# bulk = pd.read_csv('dataset/GSE141834_bulkRNAseq_normalized_counts.txt', sep='\t')

lasso = linear_model.Lasso()

cv_results = cross_validate(lasso, X, y)
sorted(cv_results.keys())

cv_results['test_score']

scores = cross_validate(lasso, X, y, cv=3,
                        scoring=('r2', 'neg_mean_squared_error'),
                        return_train_score=True)
print(scores['test_neg_mean_squared_error'])

print(scores['train_r2'])