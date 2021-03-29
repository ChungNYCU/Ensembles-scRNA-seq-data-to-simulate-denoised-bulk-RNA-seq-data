# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_validate
from sklearn.metrics import make_scorer
from sklearn.metrics import confusion_matrix
import pandas as pd

ensembled = pd.read_csv('dataset/ensembled_data_origin_v2.csv', sep=',', index_col=0)
bulk = pd.read_csv('dataset/GSE141834_scRNAseq_seuratV3_normalized.txt', sep='\t')
X = bulk

# scRNA_gene_name = ensembled.index.tolist()
# RNA_gene_name = bulk.index.tolist()
# co_gene_name = list(set(scRNA_gene_name) & set(RNA_gene_name))

# del_scRNA_gene_name=list(set(scRNA_gene_name) - set(co_gene_name))
# del_RNA_gene_name=list(set(RNA_gene_name) - set(co_gene_name))

# ensembled.drop(del_scRNA_gene_name, inplace=True)
# bulk.drop(del_RNA_gene_name, inplace=True)

X = X.T
y = []
for i in range(400):
    y.append('00')
    y.append('01')
    y.append('02')
    y.append('04')
    y.append('08')
    y.append('18')
y.sort()

clf = LogisticRegression(random_state=0).fit(X, y)

# clf.score(X, y)

cv_results = cross_validate(clf, X, y, cv=3)
sorted(cv_results.keys())

cv_results['test_score']

scores = cross_validate(clf, X, y, cv=3,
                        scoring=('r2', 'neg_mean_squared_error'),
                        return_train_score=True)
print(scores['test_neg_mean_squared_error'])

print(scores['train_r2'])
