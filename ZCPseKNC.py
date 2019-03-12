# -*- coding:utf-8 -*-

########################################################################
# 计算ZCPseKNC特征，首先计算Z curve特征作为初始局部特征，之后利用Z curve中间
# 结果即kmer频率值来计算序列初始全局特征。最后进行归一化来得到最终的特征。
########################################################################

import operator

from utils import write_feature, generate_kmer

# 矩阵合并
def merge_dict(*dict_args):
    dict_merged = {}
    for item in dict_args:
        dict_merged.update(item)
    return dict_merged


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



def interval_kmer_composition(seq,k,gap,p_kmer):
    '''
    利用Z curve中间结果计算全局特征
    参数：
        seq: 由ATCG组成的DNA序列
        k: kmer片段的长度
        gap: 全局变量的间隔
        p_kmer：kmer片段对应的p值
    '''
    p_sum = 0
    L = len(seq)
    for i in range(k-1,L-gap):
        kmer1 = seq[(i-k+1):(i+1)]
        kmer2 = seq[(i-k+1+gap):(i+1+gap)]
        j1 = i%3
        j2 = (i+gap)%3
        p_sum += (p_kmer[j1][kmer1]-p_kmer[j2][kmer2])**2
    p_average = p_sum/(L-gap-k+1)
    return p_average

def ZCPseKNC_feature(seq,k,w,lamb):
    '''
    计算ZCPseKNC特征
    参数：
        seq: 由ATCG组成的DNA序列
        k: kmer片段的长度
        w: 全局变量的权重
        lamb: 计算全局变量的最大间隔
    '''
    feature_combine = []
    feature_lambda = []
    feature_name_new = []
    feature_zcurve,feature_name,p_kmer= Zcurve_feature_general(seq,k)

    for name in feature_name:
        feature_name_new.append('k'+str(k)+'_'+'w'+str(w)+'_'+'lamb'+str(lamb)+'_'+name)

    for gap in range(1,lamb+1):
        gap_feature = interval_kmer_composition(seq,k,gap,p_kmer)
        feature_lambda.append(gap_feature)
        feature_name_new.append('k'+str(k)+'_'+'w'+str(w)+'_'+'lamb'+str(lamb)+'_'+'l'+str(gap))
    
    #calculate the normalized feature
    feature_sum = 0
    for feature_value in feature_zcurve:
        feature_sum += feature_value
    for feature_value in feature_lambda:
        feature_sum += feature_value*w
    for feature_value in feature_zcurve:
        feature_value_normal = feature_value/feature_sum
        feature_combine.append(feature_value_normal)
    for feature_value in feature_lambda:
        feature_value_normal = feature_value*w/feature_sum
        feature_combine.append(feature_value_normal)

    return feature_combine,feature_name_new


def ZCPseKNC(fileName,k,w,lamb,flag=False):
    '''
    对数据集文件的每条序列样本计算ZCPseKNC特征
    参数：
        fileName：样本数据的文件名，需要包含文件路径，要求数据格式为fasta格式
        k: kmer片段的长度
        w: 全局变量的权重
        lamb: 计算全局变量的最大间隔
        flag: flag=False输出特征向量, flag=True输出特征名
    '''
    feature_list = []
    if lamb > 10 or lamb<0:
        print('invalide lambda value!')
        return -1
    with open(fileName,'r') as f:
        for line in f.readlines():
            if line[0] == '>':
                continue
            else:
                line = line.strip().upper()
                featureVector,feature_name = ZCPseKNC_feature(line,k,w,lamb)
            feature_list.append(featureVector)
    if flag==False:
        return feature_list
    else:
        return feature_name


if __name__ == '__main__':
    # for k in (1,2,3,4,5):
    #     for w in (0.1,0.3,0.5,0.7,0.9):
    #         for lamb in (1,2,3,4,5):
    #             gc_a = ZCPseKNC('./Archaea/archaea_pos80.txt',k,w,lamb)
    #             write_feature(gc_a, './Archaea/ZCPseKNC/pos/archaea_ZCPseKNC_pos'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '+1', 'csv')
    #             gc_na = ZCPseKNC('./Archaea/archaea_neg80.txt',k,w,lamb)
    #             write_feature(gc_na, './Archaea/ZCPseKNC/neg/archaea_ZCPseKNC_neg'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '-1', 'csv')

    # for k in (1,2,3,4,5):
    #     for w in (0.1,0.5,0.9):
    #         lamb = 5
    #         gc_a = ZCPseKNC('./human_cell_line/data/human_pos.txt',k,w,lamb)
    #         write_feature(gc_a, './human_cell_line/data/ZCPseKNC/ZCPseKNC_pos'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '+1', 'csv')
    #         gc_na = ZCPseKNC('./human_cell_line/data/human_neg.txt',k,w,lamb)
    #         write_feature(gc_na, './human_cell_line/data/ZCPseKNC/ZCPseKNC_neg'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '-1', 'csv')

    # for k in (1,2,3):
    #         w = 0.3
    #         lamb = 5
    #         gc_a = ZCPseKNC('./human_cell_line/data/human_pos.txt',k,w,lamb)
    #         write_feature(gc_a, './human_cell_line/data/ZCPseKNC_k123/pos/ZCPseKNC_pos'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '+1', 'csv')
    #         gc_na = ZCPseKNC('./human_cell_line/data/human_neg.txt',k,w,lamb)
    #         write_feature(gc_na, './human_cell_line/data/ZCPseKNC_k123/neg/ZCPseKNC_neg'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '-1', 'csv')

    k = 3
    w = 0.3
    lamb = 5
    gc_a = ZCPseKNC('./human_cell_line/data/human_pos.txt',k,w,lamb)
    write_feature(gc_a, './human_cell_line/data/human_pos_ZCPseKNC'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '+1', 'csv')
    gc_na = ZCPseKNC('./human_cell_line/data/human_neg.txt',k,w,lamb)
    write_feature(gc_na, './human_cell_line/data/human_neg_ZCPseKNC'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'_csv.txt', '-1', 'csv')
 