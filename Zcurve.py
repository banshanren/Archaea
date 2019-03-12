# -*- coding:utf-8 -*-

########################################################################
# 计算Z curve特征，同时保存中间结果即kmer频率值以及特征向量的属性名，
# 用于ZCPseKNC特征的计算
########################################################################

import operator

from utils import write_feature,generate_kmer


# 矩阵合并
def merge_dict(*dict_args):
    dict_merged = {}
    for item in dict_args:
        dict_merged.update(item)
    return dict_merged

def Zcurve_feature_one(seq):
    '''
    计算kemr为单核苷酸的Z curve特征
    参数：
        seq: 由ATCG组成的DNA序列
    '''
    featureVector = []
    A=[0,0,0]
    G=[0,0,0]
    C=[0,0,0]
    T=[0,0,0]
    x=[0,0,0]
    y=[0,0,0]
    z=[0,0,0]
    n=[0,0,0]
    for i in range(len(seq)):
        if seq[i] == 'A':
            j = i%3
            A[j] += 1
            n[j] += 1
        elif seq[i] == 'G':
            j = i%3
            G[j] += 1
            n[j] += 1
        elif seq[i] == 'C':
            j = i%3
            C[j] += 1
            n[j] += 1
        elif seq[i] == 'T':
            j = i%3
            T[j] += 1
            n[j] += 1

    for i in range(3):
        A[i] = A[i]/n[i]
        G[i] = G[i]/n[i]
        C[i] = C[i]/n[i]
        T[i] = T[i]/n[i]
    for i in range(3):
        x[i] = (A[i]+G[i])-(C[i]+T[i])
        y[i] = (A[i]+C[i])-(G[i]+T[i])
        z[i] = (A[i]+T[i])-(G[i]+C[i])
    featureVector.extend(x)
    featureVector.extend(y)
    featureVector.extend(z)
    return featureVector


def Zcurve_feature_general(seq,k):
    '''
    计算可变kemr的Z curve特征
    参数：
        seq: 由ATCG组成的DNA序列
        k: kmer片段的长度
    '''
    featureVector = []
    XA = []
    XG = []
    XC = []
    XT = []
    x=[]
    y=[]
    z=[]
    feature_name = []
    x_name = []
    y_name = []
    z_name = []
    X = generate_kmer(k-1)
    for i in range(3):                         #initialization with 0
        XA.append({(pre_kmer+'A'):0 for pre_kmer in X})
        XG.append({(pre_kmer+'G'):0 for pre_kmer in X})
        XC.append({(pre_kmer+'C'):0 for pre_kmer in X})
        XT.append({(pre_kmer+'T'):0 for pre_kmer in X})
    # XT = {(i+'T'):0 for i in X}
    # XT = [XT]*3                           # this is an incrrect way because XT is just reference of list
    L = len(seq)
    for i in range(k-1,L):
        j = i%3
        kmer = seq[(i-k+1):(i+1)]                  # i:j include i but doesn't include j
        if seq[i] == 'A':                          # count the amount of codon specific kmer in a sequence
            XA[j][kmer] += 1                       
        elif seq[i] == 'G':
            XG[j][kmer] += 1
        elif seq[i] == 'C':
            XC[j][kmer] += 1
        elif seq[i] == 'T':
            XT[j][kmer] += 1
    XA_new = [sorted(i.items(),key=operator.itemgetter(0)) for i in XA]          # convert dic to list, sort in alpha order which is convenient in later calculation
    XG_new = [sorted(i.items(),key=operator.itemgetter(0)) for i in XG]
    XC_new = [sorted(i.items(),key=operator.itemgetter(0)) for i in XC]
    XT_new = [sorted(i.items(),key=operator.itemgetter(0)) for i in XT]
    dict_len = len(XA_new[0])
    for i in range(3):
        for j in range(dict_len):
            N_ij = XA_new[i][j][1]+XG_new[i][j][1]+XC_new[i][j][1]+XT_new[i][j][1]      # total count
            if N_ij==0:
                A_ij = 0                   
                G_ij = 0
                C_ij = 0
                T_ij = 0
            else:
                A_ij = XA_new[i][j][1]/float(N_ij)         # calculate the frequency of codon specific kmer    
                G_ij = XG_new[i][j][1]/float(N_ij)
                C_ij = XC_new[i][j][1]/float(N_ij)
                T_ij = XT_new[i][j][1]/float(N_ij)
                XA[i][XA_new[i][j][0]] = A_ij                
                XG[i][XG_new[i][j][0]] = G_ij
                XC[i][XC_new[i][j][0]] = C_ij
                XT[i][XT_new[i][j][0]] = T_ij
            x.append((A_ij+G_ij)-(C_ij+T_ij))            # calculate the Z curve variable
            y.append((A_ij+C_ij)-(G_ij+T_ij))
            z.append((A_ij+T_ij)-(G_ij+C_ij))
            codon = str(i+1)                             #codon position
            kmer_ = XA_new[i][j][0][:-1]                  # k-1 mer
            x_name.append('x_'+codon+'_'+kmer_)          # save freture name
            y_name.append('y_'+codon+'_'+kmer_)
            z_name.append('z_'+codon+'_'+kmer_)
    featureVector.extend(x)
    featureVector.extend(y)
    featureVector.extend(z)
    feature_name.extend(x_name)
    feature_name.extend(y_name)
    feature_name.extend(z_name)
    p_kmer = []
    for i in range(3):
        # dict_merged = dict(XA[i].items()+XG[i].items()+XC[i].items()+XT[i].items())
        dict_merged = merge_dict(XA[i],XG[i],XC[i],XT[i])
        p_kmer.append(dict_merged)
    return featureVector,feature_name,p_kmer


def Zcurve(fileName,k):
    '''
    对数据集文件的每条序列样本计算Z curve特征
    参数：
        fileName：样本数据的文件名，需要包含文件路径，要求数据格式为fasta格式
        k: kmer片段的长度
    '''
    if k>5 or k<1:
        print('invalide k value!')
        return -1
    feature_list = []
    with open(fileName) as f:
        for line in f.readlines():
            if line[0] == '>':
                continue
            else:
                seq = line.strip().upper()
                featureVector,feature_name,p_kmer = Zcurve_feature_general(seq,k)
                feature_list.append(featureVector)
    return feature_list

def Zcurve_one(fileName):
    '''
    对数据集文件的每条序列样本计算Z curve特征
    参数：
        fileName：样本数据的文件名，需要包含文件路径，要求数据格式为fasta格式
    '''
    feature_list = []
    with open(fileName) as f:
        for line in f.readlines():
            if line[0] == '>':
                continue
            else:
                seq = line.strip().upper()
                featureVector = Zcurve_feature_one(seq)
                feature_list.append(featureVector)
    return feature_list

if __name__ == '__main__':

    # gc_a = Zcurve_one('./Archaea/archaea_pos80.txt')
    # write_feature(gc_a, './Archaea/Zcurve/archaea_pos_Zcurve_k_1_csv.txt', '+1', 'csv')
    # gc_na = Zcurve_one('./Archaea/archaea_neg80.txt')
    # write_feature(gc_na, './Archaea/Zcurve/archaea_neg_Zcurve_k_1_csv.txt', '-1', 'csv')

    # for k in (1,2,3,4,5):
    #     gc_a = Zcurve('./human_cell_line/data/human_pos.txt',k)
    #     write_feature(gc_a, './human_cell_line/data/Zcurve/human_pos_Zcurve'+'_k_'+str(k)+'_csv.txt', '+1', 'csv')
    #     gc_na = Zcurve('./human_cell_line/data/human_neg.txt',k)
    #     write_feature(gc_na, './human_cell_line/data/Zcurve/human_neg_Zcurve'+'_k_'+str(k)+'_csv.txt', '-1', 'csv')

    gc_a_all = Zcurve('./human_cell_line/data/human_pos.txt',1)
    print('gc_a_all', len(gc_a_all),len(gc_a_all[0]))
    gc_na_all = Zcurve('./human_cell_line/data/human_neg.txt',1)
    for k in (2,3):
        gc_a = Zcurve('./human_cell_line/data/human_pos.txt',k)
        print(k)
        print('gc_a', len(gc_a),len(gc_a[0]))
        for i in range(len(gc_a)):
            gc_a_all[i].extend(gc_a[i])
        print('gc_a_all', len(gc_a_all),len(gc_a_all[0]))
        gc_na = Zcurve('./human_cell_line/data/human_neg.txt',k)
        for i in range(len(gc_na)):
            gc_na_all[i].extend(gc_na[i])
    write_feature(gc_a_all, './human_cell_line/data/Zcurve/human_pos_Zcurve_all_3_csv.txt', '+1', 'csv')
    write_feature(gc_na_all, './human_cell_line/data/Zcurve/human_neg_Zcurve_all_3_csv.txt', '-1', 'csv')