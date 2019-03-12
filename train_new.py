# -*- coding: utf-8 -*-

# from libsvm.checkdata import check_data
# from libsvm.python.svmutil import *
# from libsvm.python.plotroc import *
import math
import numpy as np
import pandas as pd
from sklearn import svm
from sklearn.model_selection import KFold,StratifiedKFold,LeaveOneOut
from sklearn.utils import check_random_state
from imblearn.over_sampling import SMOTE
from sklearn.metrics import roc_auc_score, f1_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.externals import joblib
from xgboost import XGBClassifier


def evaluation(y_test, y_pred, y_score):
    '''
    计算所有的评价指标
    参数：
        y_test: 真实的样本标签
        y_pred: 预测的样本标签
        y_score: 预测样本的概率值
    中间变量：
        TP: true positive rate
        TN: true negative rate
        FP: false positive rate
        FN: false negative rate
    '''
    TP = 0.0
    TN = 0.0
    FP = 0.0
    FN = 0.0
    for i in range(len(y_test)):
        if y_test[i] == 1.0 and y_pred[i] == 1.0:
            TP += 1.0
        elif y_test[i] == 1.0 and y_pred[i] == -1.0:
            FN += 1.0
        elif y_test[i] == -1.0 and y_pred[i] == 1.0:
            FP += 1.0
        elif y_test[i] == -1.0 and y_pred[i] == -1.0:
            TN += 1.0
    try:
        SN = TP / (TP + FN)
    except ZeroDivisionError:
        SN = 0.0
    try:
        SP = TN / (FP + TN)
    except ZeroDivisionError:
        SP = 0.0
    try:
        ACC = (TP + TN) / (TP + TN + FP + FN)
    except ZeroDivisionError:
        ACC = 0.0
    try:
        MCC = (TP * TN - FP * FN) / math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    except ZeroDivisionError:
        MCC = 0.0
    try:
        # AUC = plot_roc(deci_value, origin_labels, output, title, roc, roc_data_file)
        AUC = roc_auc_score(y_test, y_score)
    except ZeroDivisionError:
        AUC = 0.0
    fscore = f1_score(y_test, y_pred)
    return ACC, MCC, AUC, SN, SP, fscore
    
def cross_validation(X, y, clf, option = '5', smote=False, flag=0, ml='svc'):
    '''
    对训练数据进行交叉验证
    参数：
        X: 样本数据的特征向量
        y: 样本数据的标签
        option: 交叉验证类型，包括5折，10折交叉以及jackknife验证
        smote: 设置是否在交叉验证过程进行smote
        ml: 分类器类型，包括linear SVC和kernel SVC等
    '''
    if option == '5' or option == '10':
        option = int(option)
        skf = StratifiedKFold(n_splits=option, shuffle=True)  
        cv_split = list(skf.split(X, y))
    elif option == 'j':
        loo = LeaveOneOut()
        cv_split = list(loo.split(X, y))
    else:
        print('error cv option!')
        return -1
    y_score_all = [0.0 for i in range(len(y))]
    y_pred_all = [0.0 for i in range(len(y))]
    for i, (train_index, test_index) in enumerate(cv_split):
        X_train = X[train_index]
        y_train = y[train_index]
        if smote:
            # estimator = svm.SVC(class_weight='balanced', random_state=check_random_state(None), kernel='linear')
            # estimator = svm.SVC(class_weight='balanced', random_state=check_random_state(None))
            X_train, y_train = SMOTE(kind='svm').fit_sample(X_train, y_train)
            # X_train, y_train = Smote(sampling_rate=2).fit_sample(X_train, y_train)
        X_test = X[test_index]
        y_test = y[test_index]

        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        # y_score = clf.decision_function(X_test)
        if ml == 'lsvc' or ml == 'svc':
            y_score = clf.decision_function(X_test)
        else:
            y_score = clf.predict_proba(X_test)[:, 1]
        
        for j in range(len(test_index)):
            y_pred_all[test_index[j]] = y_pred[j]
            y_score_all[test_index[j]] = y_score[j]

    ACC, MCC, AUC, Sn, Sp, fscore = evaluation(y, y_pred_all, y_score_all)
    if flag == 0:
        return ACC, MCC, AUC, Sn, Sp
    else:
        return ACC, MCC, AUC, Sn, Sp, fscore



def para_optimization(X, y, ml = 'linear', eval='AUC'):
    '''
    参数优化，选取最优的参数用于后续分类
    参数：
        X: 样本数据的特征向量
        y: 样本数据的标签
        ml: 分类器类型，包括包括linear SVC和kernel SVC
        eva: 设置用于筛选最优参数的评价指标
    '''
    best_c = 0
    best_g = 0
    best_AUC = 0
    best_MCC = 0
    best_ACC = 0
    if ml == 'linear':
        for c in range(-5,10, 1):
            lsvc = svm.LinearSVC(C=2**c)
            # lsvc = svm.SVC(C=2**c, kernel='linear')
            print('current c:', c)
            ACC, MCC, AUC, Sn, Sp = cross_validation(X, y, lsvc, option = '5')
            if eval == 'AUC':
                if AUC > best_AUC:
                    best_AUC = AUC
                    best_c = c
            elif eval == 'MCC':
                if MCC > best_MCC:
                    best_MCC = MCC
                    best_c = c
            elif eval == 'ACC':
                if ACC > best_ACC:
                    best_ACC = ACC
                    best_c = c
            else:
                print('wrong evaluation matrix!')
                return -1

    elif ml == 'kernel':
        for c in range(-5, 10, 2):
            for g in range(-10, 5, 2):
                svc = svm.SVC(C=2**c, gamma=2**g)
                ACC, MCC, AUC, Sn, Sp = cross_validation(X, y, svc, option = '5')
                print('current c and g:', c, g)
                print('current score:', ACC, MCC, AUC, Sn, Sp)
                if eval == 'AUC':
                    if AUC > best_AUC:
                        best_AUC = AUC
                        best_c = c
                        best_g = g
                elif eval == 'MCC':
                    if MCC > best_MCC:
                        best_MCC = MCC
                        best_c = c
                        best_g = g
                elif eval == 'ACC':
                    if ACC > best_ACC:
                        best_ACC = ACC
                        best_c = c
                        best_g = g
                else:
                    print('wrong evaluation matrix!')
                    return -1
    else:
        print('wrong option!')
        return -1
    return best_c, best_g

if __name__ == "__main__":
    # pos_data = pd.read_csv('./data/Zcurve/archaea_pos_Zcurve_k_5_csv.txt', header=None)
    # neg_data = pd.read_csv('./data/Zcurve/archaea_neg_Zcurve_k_5_csv.txt', header=None)
    pos_data = pd.read_csv('./data/archaea_pos_ZCPseKNC_k3_w0.6_l3_csv.txt', header=None)
    neg_data = pd.read_csv('./data/archaea_neg_ZCPseKNC_k3_w0.6_l3_csv.txt', header=None)
    pos_len = pos_data.shape[0]
    neg_len = neg_data.shape[0]

    data = pd.concat([pos_data, neg_data], axis=0, ignore_index=True)

    X = data.values
    y = np.array([1.0 for i in range(pos_len)] + [-1.0 for i in range(neg_len)])

    # print('parameter optimization start')
    best_c, best_g = para_optimization(X, y, ml='linear', eval='AUC')
    print('best c', 2**best_c)
    # best_c, best_g = para_optimization(X, y, ml='kernel', eval='AUC')
    # print('best c and best g:', 2**best_c, 2**best_g)

    lsvc = svm.LinearSVC(C=2**best_c)
    # lsvc = svm.SVC(C=2**best_c, random_state=check_random_state(None), kernel='linear')
    # svc = svm.SVC(C=2**best_c, gamma=2**best_g, random_state=check_random_state(None), probability=True)
    # xgb = XGBClassifier(max_depth=5, n_estimators=200, min_child_weight=1)
    # rfc = RandomForestClassifier(n_estimators=50)

    # clf = svc.fit(X,y)
    # joblib.dump(clf, './model/archaea_train_model.m')

    print('start cross validation')
    # ACC, MCC, AUC, Sn, Sp, y_score_all = cross_validation(X, y, lsvc, option = '5', smote=True, flag = 1)     # 5 fold
    # print(ACC, MCC, AUC, Sn, Sp)
    # with open('./human_cell_line/prob_value/human_FV45_top800_lsvc_LinearSVMSmote_v5_prob_value.txt','w+') as f:
    #     for i in range(len(y)):
    #         f.write(str(y[i]) + '\t' + str(y_score_all[i]) + '\n')

    skf = StratifiedKFold(n_splits=5, shuffle=True)  # 5 fold
    skf_split = list(skf.split(X, y))
    loo = LeaveOneOut()                              
    loo_split = list(loo.split(X, y))
    y_score_all = [0.0 for i in range(pos_len + neg_len)]
    y_pred_all = [0.0 for i in range(pos_len + neg_len)]
    for i, (train_index, test_index) in enumerate(skf_split):
        # print("Fold", i)
        X_train = X[train_index]
        y_train = y[train_index]
        # estimator = svm.SVC(class_weight='balanced', random_state=check_random_state(None), kernel='linear')
        # estimator = svm.SVC(class_weight='balanced', random_state=check_random_state(None))
        # X_train, y_train = SMOTE(kind = 'svm').fit_sample(X_train, y_train)
        # # X_train, y_train = Smote(sampling_rate=2).fit_sample(X_train, y_train)
        X_test = X[test_index]
        y_test = y[test_index]
        lsvc.fit(X_train, y_train)
        y_pred = lsvc.predict(X_test)
        # y_score = lsvc.predict_proba(X_test)[:, 1]
        y_score = lsvc.decision_function(X_test)
        for j in range(len(test_index)):
            y_pred_all[test_index[j]] = y_pred[j]
            y_score_all[test_index[j]] = y_score[j]
    ACC, MCC, AUC, Sn, Sp, Fscore = evaluation(y, y_pred_all, y_score_all)
    print(ACC, MCC, AUC, Sn, Sp)
    # with open('./human_cell_line/prob_value/human_ZCPseKNC_top800_lsvc_v5_prob_value.txt','w+') as f:
    #     for i in range(len(y)):
    #         f.write(str(y[i]) + '\t' + str(y_score_all[i]) + '\n')