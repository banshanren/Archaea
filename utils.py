# -*- coding:utf-8 -*-
import numpy as np
import os

def splice(inputFilePath,outputFilePath):
    '''
    The origin file splite a sequence into several lines, this function just splices them into one line.
    '''
    with open(inputFilePath) as inFile:
        with open(outputFilePath,'w') as outFile:
            seq = ''
            for line in inFile.readlines():
                if line[0] == '>':
                    if seq != '':
                        outFile.write(seq+'\n')
                    outFile.write(line)
                    seq = ''
                else:
                    seq += line.strip()
            outFile.write(seq+'\n')

def get_ID(filePath):
    '''
    get gene ID form dna seq
    '''
    IDList = []
    with open(filePath) as f:
        for line in f.readlines():
            if line[0] == '>':
                IDList.append(line)
    return IDList

def select_seq(inputFilePath,outputFilePath,IDList):
    '''
    select seq from origin pro seq file based on ID 
    '''
    with open(outputFilePath,'w') as outFile:
        with open(inputFilePath) as inFile:
            seqID = ''
            seq = ''
            flag = 0
            for line in inFile.readlines():
                if line[0] == '>':
                    if seq!='' and flag==1:
                        outFile.write(seqID+'\n')
                        outFile.write(seq+'\n')
                        flag = 0
                    if line in IDList:
                        flag = 1
                    seqID = line.strip()
                    seq = ''
                else:
                    seq += line.strip()
            if seq!='' and flag==1:
                outFile.write(seqID+'\n')
                outFile.write(seq+'\n')

def protein_pre(inputFilePath,outputFilePath):
    '''
    remove protein sequences which contain unstandard amino acids
        @date:2018/06/26
    ''' 
    aaList = 'ACDEFGHIKLMNPQRSTVWY'                       
    with open(outputFilePath,'w') as outFile:
        with open(inputFilePath) as inFile:
            seqID = ''
            seq = ''
            for line in inFile.readlines():
                if line[0] == '>':
                    if seq!='':
                        flag = 0
                        for i in seq:
                            if i not in aaList:
                                flag = 1
                                break
                        if flag == 0:
                            outFile.write(seqID+'\n')
                            outFile.write(seq+'\n')
                    seqID = line.strip()
                    seq = ''
                else:
                    seq += line.strip()
            if seq!='':
                flag = 0
                for i in seq:
                    if i not in aaList:
                        flag = 1
                        break
                if flag == 0:
                    outFile.write(seqID+'\n')
                    outFile.write(seq+'\n')

def rare_AA(inputFilePath):
    '''
    calculate the frequence of rare amino acid in protein seq 
    '''
    aaList = 'ACDEFGHIKLMNPQRSTVWY'
    with open(inputFilePath) as inFile:
        for line in inFile.readlines():
            if line[0] != '>':
                t = 0
                n = len(line.strip())
                for i in line.strip():
                    if i not in aaList:
                        t += 1
                print(t)

def TMHMM_pre(inputFilePath,outputFilePath,label):
    '''
    Extract TMHMM feature to be svm formate
    '''
    with open(outputFilePath,'w') as outFile:
        with open(inputFilePath) as inFile:
            for line in inFile.readlines():
                line = line.strip().split('\t')
                tmp = line[4].split('=')
                predTHM = tmp[1]
                outFile.write(label+' '+'1:'+predTHM+'\n')

def normalization(inputPosFilePath,inputNegFilePath,outputPosFilePath,outputNegFilePath):
    '''
    Normalize feature
    '''
    posList = []
    negList = []
    allList = []
    with open(inputPosFilePath) as inFilePos:
        with open(inputNegFilePath) as inFileNeg:
            PosLines = inFilePos.readlines()
            PosLen = len(PosLines)
            NegLines = inFileNeg.readlines()
            NegLen = len(NegLines)
            for line in PosLines:
                val = line.strip().split(':')
                allList.append(int(val[1]))
            for line in NegLines:
                val = line.strip().split(':')
                allList.append(int(val[1]))
    maxVal = max(allList)
    minVal = min(allList)
    for i in range(len(allList)):
        allList[i] = (allList[i]-minVal)/(maxVal-minVal)
    print(allList)


def csv2svm(inPath,outPath,label,title=False):
    '''
    transfer csv formate to libsvm formate
    '''
    with open(outPath,'w') as fo:
        with open(inPath) as fi:
            if title==True:
                line=fi.readline()
            for line in fi.readlines():
                line = line.strip().split(',')
                fo.write(label)
                for i in range(len(line)):
                    fo.write(' '+str(i+1)+':'+line[i])
                fo.write('\n')

def svm2csv(inPath,outPath,title=False):
    '''
    transfer libsvm formate to csv formate
    '''
    with open(outPath,'w') as fo:
        with open(inPath) as fi:
            if title==True:
                line=fi.readline()
            for line in fi.readlines():
                line = line.strip().split(' ')
                fo.write(line[1].split(':')[1])
                for i in range(2,len(line)):
                    fo.write(','+line[i].split(':')[1])
                fo.write('\n')

def combine_feature_v1(fList, outPath, formation='svm'):
    '''
    combine feature from two or more feature files
    '''
    with open(outPath,'w+') as fo:
        if formation == 'svm':
            for line1 in fList[0].readlines():
                line1 = line1.strip()
                line1_list = line1.split(' ')
                l1 = int(line1_list[-1].split(':')[0])
                # print('------------------',l1)
                for i in range(1,len(fList)):
                    line2 = fList[i].readline().strip().split(' ')
                    l2 = len(line2)
                    for i in range(1, l2):
                        list = line2[i].split(':')
                        line1 = line1 + ' ' + str(l1 + i) + ":" + list[1]
                    l1 = l1+l2-1
                fo.write(line1 + '\n')
        elif formation == 'csv':
            for line1 in fList[0].readlines():
                line1 = line1.strip()
                for i in range(1,len(fList)):
                    line2 = fList[i].readline().strip()
                    line1 = line1 + ',' + line2
                fo.write(line1+ '\n')
        else:
            print('wrong formation')
            exit()

def combine_feature(inPath, fileNameList, outPath, formation='svm'):
    '''
    combine feature from two or more feature files
    '''
    fList=[]
    with open(outPath,'w+') as fo:
        for fileName in fileNameList:
            fList.append(open(inPath+fileName))
        if formation == 'svm':
            for line1 in fList[0].readlines():
                line1 = line1.strip()
                line1_list = line1.split(' ')
                l1 = int(line1_list[-1].split(':')[0])
                # print('------------------',l1)
                for i in range(1,len(fList)):
                    line2 = fList[i].readline().strip().split(' ')
                    l2 = len(line2)
                    for i in range(1, l2):
                        list = line2[i].split(':')
                        line1 = line1 + ' ' + str(l1 + i) + ":" + list[1]
                    l1 = l1+l2-1
                fo.write(line1 + '\n')
        elif formation == 'csv':
            for line1 in fList[0].readlines():
                line1 = line1.strip()
                for i in range(1,len(fList)):
                    line2 = fList[i].readline().strip()
                    line1 = line1 + ',' + line2
                fo.write(line1+ '\n')
        else:
            print('wrong formation')
            exit()

def write_feature(featureList, outPath, label, formation='svm'):
    with open(outPath, 'w') as f:
        if formation=='svm':
            for featureVector in featureList:
                f.write(label+' ')
                f.write(str(1)+':'+str(featureVector[0]))
                for i in range(1,len(featureVector)):
                    f.write(' '+str(i+1)+':'+str(featureVector[i]))
                f.write('\n')
        elif formation=='csv':
            for featureVector in featureList:
                f.write(str(featureVector[0]))
                for i in range(1,len(featureVector)):
                    f.write(','+str(featureVector[i]))
                f.write('\n')
        else:
            print('formation error')

def generate_kmer(k):
    '''
    generate kmers for DNA seq
    '''
    base = ['A','G','C','T']
    kmer_list = ['']
    for i in range(k):
        new_kmer_list = []
        for kmer in kmer_list:
            for c in base:
                new_kmer_list.append(kmer+c)
        kmer_list = new_kmer_list
    return kmer_list

if __name__ == "__main__":
    csv2svm('./human_cell_line/data/ZCPseKNC_k123/human_pos_ZCPseKNC_k123_csv.txt', './human_cell_line/data/ZCPseKNC_k123/human_pos_ZCPseKNC_k123_svm.txt', label='+1')
    csv2svm('./human_cell_line/data/ZCPseKNC_k123/human_neg_ZCPseKNC_k123_csv.txt', './human_cell_line/data/ZCPseKNC_k123/human_neg_ZCPseKNC_k123_svm.txt', label='-1')
    # svm2csv('./human_cell_line/data/Kmer/human_pos_svm_Kmer_k_4.txt', './human_cell_line/data/Kmer/human_pos_Kmer_k_4_csv.txt')
    # svm2csv('./human_cell_line/data/Kmer/human_neg_svm_Kmer_k_4.txt', './human_cell_line/data/Kmer/human_neg_Kmer_k_4_csv.txt')

    # in_pos = './human_cell_line/data/ZCPseKNC_k123/pos/'
    # file_name_pos = os.listdir(in_pos)
    # out_pos = './human_cell_line/data/ZCPseKNC_k123/human_pos_ZCPseKNC_k123_csv.txt'
    # combine_feature(in_pos, file_name_pos, out_pos, formation='csv')
    # in_neg = './human_cell_line/data/ZCPseKNC_k123/neg/'
    # file_name_neg = os.listdir(in_neg)
    # out_neg = './human_cell_line/data/ZCPseKNC_k123/human_neg_ZCPseKNC_k123_csv.txt'
    # combine_feature(in_neg, file_name_neg, out_neg, formation='csv')