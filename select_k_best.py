# -*- coding:utf-8 -*-

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest,SelectFromModel
from sklearn.feature_selection import chi2,f_classif
from sklearn.ensemble import GradientBoostingClassifier,RandomForestClassifier
from xgboost import XGBClassifier
from utils import write_feature



posX = pd.read_csv('./human_cell_line/data/ZCPseKNC/human_pos_ZCPseKNC_all_csv.txt',header=None)
negX = pd.read_csv('./human_cell_line/data/ZCPseKNC/human_neg_ZCPseKNC_all_csv.txt',header=None)
# print(posX.columns[np.where(np.isnan(posX))[1]])
posX = np.array(posX)
negX = np.array(negX)
X = np.concatenate((posX,negX),axis=0)
# print(X.shape)

posy = np.ones((len(posX)))
negy = np.zeros((len(negX)))
negy[0:len(negX)] = -1
y = np.concatenate((posy,negy))

# feature_importance = XGBClassifier().fit(X,y).feature_importances_
# feature_importance_array = np.array(feature_importance)
# sort_index = list(np.argsort(feature_importance_array))
# for i in range(-1,-11,-1):
#     index = sort_index[i]
#     print('feature_name:'+feature_name[index]+'         '+'importance:'+str(feature_importance[index]))

X_new = SelectKBest(f_classif,k=200).fit_transform(X,y)
# print(X_new.shape)
posX_new = X_new[0:len(posX)]
negX_new = X_new[len(posX):len(posX)+len(negX)]
# print(posX_new.shape)
pos_outPath = './human_cell_line/data/ZCPseKNC/human_pos_ZCPseKNC_top200_csv.txt'
neg_outPath = './human_cell_line/data/ZCPseKNC/human_neg_ZCPseKNC_top200_csv.txt'
write_feature(posX_new, pos_outPath, '+1', 'csv')
write_feature(negX_new, neg_outPath, '-1', 'csv')
