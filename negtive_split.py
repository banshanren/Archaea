# -*- coding:utf-8 -*-

import pandas as pd
import numpy as np 
from utils import write_feature

negX = pd.read_csv('./Archaea/archaea_ZPseKNC/archaea_ZPseKNC_neg_XGBC_top10_csv.txt',header=None)
step = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]
split_size = []
for i in step:
    split_size.append(int(518*i))
# print(split_size)

for s in split_size:
    index  = np.random.choice(1072,size=s,replace=False)	
    negX_new = negX.loc[index]
    negX_new = negX_new.values
    # print(negX_new.shape)
    neg_outPath = './Archaea/archaea_ZPseKNC/archaea_negtive_split/archaea_ZPseKNC_neg_XGBC_top10_'+str(s)+'_csv.txt'
    write_feature(negX_new, neg_outPath, '-1', 'csv')

# negX = pd.read_csv('./Archaea/archaea_ZPseKNC/archaea_ZPseKNC_neg_XGBC_top10_csv.txt',header=None)

# inc = int(518*0.1)      #增量
# # print(inc)

# cnt = int((1072-518)/inc)          #做的次数
# for i in range(cnt+1):
#     negX_new = negX.loc[:517+inc*i]               #闭区间
#     negX_new = negX_new.values
#     # print(negX_new.shape)
#     s = str(518+inc*i)
#     neg_outPath = './Archaea/archaea_ZPseKNC/archaea_negtive_split/archaea_ZPseKNC_neg_XGBC_top10_'+str(s)+'_svm.txt'
#     write_feature(negX_new, neg_outPath, '-1', 'svm')
