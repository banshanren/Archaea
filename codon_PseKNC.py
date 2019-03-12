# -*- coding:utf-8 -*-

from Zcurve import Zcurve_feature_general
from utils import write_feature

def interval_kmer_composition(seq,k,gap,p_kmer):
    '''
        compute lambda feature
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

def codon_PseKNC_feature(seq,k,w,lamb):
    feature_combine = []
    feature_lambda = []
    feature_name_new = []
    feature_codon = []
    feature_zcurve,feature_name,p_kmer= Zcurve_feature_general(seq,k)
    for i in range(3):
        feature_codon.extend(list(p_kmer[i].values()))

    for name in feature_name:
        feature_name_new.append('k'+str(k)+'_'+'w'+str(w)+'_'+'lamb'+str(lamb)+'_'+name)

    for gap in range(1,lamb+1):
        gap_feature = interval_kmer_composition(seq,k,gap,p_kmer)
        feature_lambda.append(gap_feature)
        feature_name_new.append('k'+str(k)+'_'+'w'+str(w)+'_'+'lamb'+str(lamb)+'_'+'l'+str(gap))
    
    #calculate the normalized di
    feature_sum = 0
    for feature_value in feature_codon:
        feature_sum += feature_value
    for feature_value in feature_lambda:
        feature_sum += feature_value*w
    for feature_value in feature_codon:
        feature_value_normal = feature_value/feature_sum
        feature_combine.append(feature_value_normal)
    for feature_value in feature_lambda:
        feature_value_normal = feature_value*w/feature_sum
        feature_combine.append(feature_value_normal)

    return feature_combine,feature_name_new

def codon_PseKNC(fileName,k,w,lamb,flag=0):
    '''
    k: kmer
    w: weight
    lamb: lambda
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
                line = line.strip()
                featureVector,feature_name = codon_PseKNC_feature(line,k,w,lamb)
            feature_list.append(featureVector)
    if flag==0:
        return feature_list
    else:
        return feature_name

if __name__ == '__main__':
    for k in (1,2,3):
        for w in (0.1,0.5,0.9):
            for lamb in (1,3,5):
                gc_a = codon_PseKNC('./Archaea/archaea_pos80.txt',k,w,lamb)
                write_feature(gc_a, './Archaea/codon_PseKNC/codon_PseKNC_pos_svm'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'.txt', '+1', 'svm')
                gc_na = codon_PseKNC('./Archaea/archaea_neg80.txt',k,w,lamb)
                write_feature(gc_na, './Archaea/codon_PseKNC/codon_PseKNC_neg_svm'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'.txt', '-1', 'svm')

    # k=3
    # w=0.6
    # lamb=3
    # gc_a = Zcurve_PseKNC('./Archaea/archaea_pos80.txt',k,w,lamb)
    # # write_feature(gc_a, './tmp_Zcurve+PseKNC_pos_svm'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'.txt', '+1', 'svm')
    # # gc_na = Zcurve_PseKNC('./Archaea/archaea_neg80.txt',k,w,lamb)
    # # write_feature(gc_na, './tmp_Zcurve+PseKNC_neg_svm'+'_k_'+str(k)+'_w_'+str(w)+'_lambda_'+str(lamb)+'.txt', '-1', 'svm')
