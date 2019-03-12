# -*- coding:utf-8 -*-

########################################################################
# 对初始合并的特征做特征选择， 同时记录重要性最高的10个特征
########################################################################

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest,SelectFromModel
from sklearn.ensemble import GradientBoostingClassifier,RandomForestClassifier
from xgboost import XGBClassifier
from utils import write_feature

from Zcurve_PseKNC import Zcurve_PseKNC

feature_name = []
posX = []
negX = []
parameters = ((3,0.1,3),(3,0.1,5),(3,0.2,3),(3,0.5,3),(3,0.6,3),
       (3,0.8,3),(3,0.9,1),(3,0.9,3),(3,0.9,6),(3,0.9,8))
for k,w,lamb in parameters:
    feature_name.extend(Zcurve_PseKNC('./Archaea/archaea_pos80.txt',k,w,lamb,flag=1))
    posX.extend(Zcurve_PseKNC('./Archaea/archaea_pos80.txt',k,w,lamb))
    negX.extend(Zcurve_PseKNC('./Archaea/archaea_neg80.txt',k,w,lamb))
# posX = pd.read_csv('./Archaea/archaea_ZPseKNC/archaea_ZPseKNC_pos_top10_csv.txt',header=None)
# negX = pd.read_csv('./Archaea/archaea_ZPseKNC/archaea_ZPseKNC_neg_top10_csv.txt',header=None)

posX = np.array(posX)
negX = np.array(negX)
X = np.concatenate((posX,negX),axis=0)

posy = np.ones((518))
negy = np.zeros((1072))
negy[0:1072] = -1
y = np.concatenate((posy,negy))

feature_importance = XGBClassifier().fit(X,y).feature_importances_
feature_importance_array = np.array(feature_importance)
sort_index = list(np.argsort(feature_importance_array))
for i in range(-1,-11,-1):
    index = sort_index[i]
    print('feature_name:'+feature_name[index]+'         '+'importance:'+str(feature_importance[index]))

X_new = SelectFromModel(XGBClassifier()).fit_transform(X,y)
# print(X_new.shape)
posX_new = X_new[0:518]
negX_new = X_new[518:1590]
# print(posX_new.shape)
pos_outPath = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_pos_xgb_top10_svm.txt'
neg_outPath = './Archaea/archaea_ZPseKNC/archaea_ZPseKNC_neg_xgb_top10_svm.txt'
write_feature(posX_new, pos_outPath, '+1', 'svm')
write_feature(negX_new, neg_outPath, '-1', 'svm')
