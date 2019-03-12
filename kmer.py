# -*- coding:utf-8 -*-

########################################################################
# 计算最简单的kmer特征，保存特征结果
########################################################################

import operator
from utils import write_feature, generate_kmer


def kmer_feature(seq, k):
    '''
    计算一条序列的kmer特征
    参数：
        seq: 由ATCG组成的DNA序列
        k: kmer片段的长度
    '''
    kmer_list = generate_kmer(k)
    kmer_freq = {kmer:0 for kmer in kmer_list}
    L = len(seq)
    for i in range(L-k+1):
        kmer = seq[i:i+k]
        kmer_freq[kmer] += 1
    freq_list = list(kmer_freq.values())
    for i in range(len(freq_list)):
        freq_list[i] = freq_list[i]/(L-k+1)
    return freq_list

def Kmer(fileName, k):
    '''
    对文件的每条序列样本计算kmer特征
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
                featureVector = kmer_feature(seq,k)
                feature_list.append(featureVector)
    return feature_list

if __name__ == '__main__':

    for k in [1,2,3,4,5]:
        gc_a = Kmer('./human_cell_line/data/human_pos.txt', k)
        write_feature(gc_a, './human_cell_line/data/kmer_me/human_pos_kmer'+'_k_'+str(k)+'_csv.txt', '+1', 'csv')
        gc_na = Kmer('./human_cell_line/data/human_neg.txt', k)
        write_feature(gc_na, './human_cell_line/data/kmer_me/human_neg_kmer'+'_k_'+str(k)+'_csv.txt', '-1', 'csv')