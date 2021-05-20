# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 08:59:36 2021

@author: chung
"""
from sklearn.svm import LinearSVC
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_validate
import pandas as pd

ensembled = pd.read_csv('dataset/ensembled_data.csv', sep=',', index_col=0)
bulk = pd.read_csv('dataset/GSE141834_bulkRNAseq_normalized_counts.txt', sep='\t')
X = ensembled

X = X.T
# bulk = bulk.T
y = []
for i in range(400):
    y.append('00')
    y.append('01')
    y.append('02')
    y.append('04')
    y.append('08')
    y.append('18')
y.sort()
clf = make_pipeline(StandardScaler(), LinearSVC(random_state=0, tol=1e-5))

cv_results = cross_validate(clf, X, y, cv=3)
sorted(cv_results.keys())

cv_results['test_score']

scores = cross_validate(clf, X, y, cv=3,
                        scoring=('r2', 'neg_mean_squared_error'),
                        return_train_score=True)
print(scores['test_neg_mean_squared_error'])

print(scores['train_r2'])

# ---------------------------------------------------
# clf.fit(X, y)

# print(clf.score(X, y, sample_weight=None))

# print(clf.named_steps['linearsvc'].coef_)
# print(clf.named_steps['linearsvc'].intercept_)
# c=clf.decision_function(X)

# for i in range(7):
#     b=bulk.iloc[i]
#     print(clf.predict([b]))
    
# import seaborn as sns 
# import pandas as pd
# import matplotlib.pyplot as plt

# ax = sns.heatmap(c)
# ax.invert_yaxis()
# print(c)
# plt.show()