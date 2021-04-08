import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import ElasticNet
import os

os.chdir('C:\\GitHub\\')
files={'Leukemia cells': 'Leukemia CancerRX exp.xlsx',
       'NB cells': 'NB CancerRX exp.xlsx',
       'NB patients': 'Vergsteeg exp.xlsx',
       'AML patients': 'Verhaal exp.xlsx'}
genesets={'RTK': 'RTK genes.xlsx'}
dataset='NB cells'
geneset='RTK'
HDAC=['HDAC1', 'HDAC2', 'HDAC3', 'HDAC6', 'HDAC8']
df= pd.read_excel(files[dataset])
df_genes= pd.read_excel(genesets[geneset])
data_genes= set(dict.fromkeys(df.Gene.dropna()))
input_genes= set(dict.fromkeys(df_genes.Gene.dropna()))
df.set_index('Gene', inplace=True)
genes=list(input_genes & data_genes)

"""
#Code used for penatly parameters determination
X=df.loc[HDAC].transpose()
y=df.transpose()[gene]
elastic=ElasticNet(normalize=True, tol=0.005)
search=GridSearchCV(estimator=elastic,param_grid={'alpha':[0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1],'l1_ratio':[.1,.2,.4,.6,.8]},scoring='neg_mean_squared_error',n_jobs=1,refit=True,cv=5)
search=GridSearchCV(estimator=elastic,param_grid={'alpha':[0.0001, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 1],'l1_ratio':[.1,.2,.4,.6,.8]},scoring='r2',n_jobs=1,refit=True,cv=3)
search.fit(X,y)
print(search.best_params_)
print(abs(search.best_score_))
"""

param={'NB cells': [0.01, 0.2], 'Leukemia cells': [0.001, 0.8], 'NB patients': [0.01, 0.4], 'AML patients': [0.001, 0.8]}

elastic=ElasticNet(normalize=True,alpha=param[dataset][0],l1_ratio=param[dataset][1], selection='random')
df_elastic=pd.DataFrame()
for gene in genes:
    X=df.loc[HDAC].transpose()
    y=df.transpose()[gene]
    elastic.fit(X,y)
    coef_dict_baseline = {}
    for coef, feat in zip(elastic.coef_,X.columns):
        coef_dict_baseline[feat] = coef
    df_temp = pd.DataFrame.from_dict(coef_dict_baseline, orient='index')
    df_elastic = pd.concat([df_elastic, df_temp], axis=1)
df_elastic.columns= genes
df_elastic.to_excel('Elastic_'+dataset+'_'+geneset+'_HDAC.xlsx')
