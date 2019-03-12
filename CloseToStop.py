#coding:utf-8

# import numpy as np
# import pandas as pd
# import matplotlib as mpl
# import matplotlib.pyplot as plt


CTSlist = ['AAG','GAG','CAG','TGG','TCG','TTG','TAA','TAC','TAT',
                'AAA','GAA','CAA','TGA','TCA','TTA','TAG','TAC','TAT',
                'AGA','GGA','CGA','TAA','TCA','TTA','TGG','TGC','TGT']

def Frequency(seq):
    cts = 0
    total = 0
    for i in range(0,len(seq),3):
        if seq[i:i+3] in CTSlist:                                  
            cts += 1
        total += 1
    return cts/total
        
def CloseToStop(fileName):
    freqList = []
    with open(fileName) as f:
        for line in f.readlines():
            if line[0] == '>':
                continue
            else:
                freqList.append(Frequency(line))
    return freqList

def writeFeature(featureList, label):
    with open('./ctsfeature'+label+'.txt', 'w') as f:
        for i in featureList:
            f.write(label+' ')
            f.write('1:'+str(i))
            f.write('\n')

cts_a = CloseToStop('./archaea_pos80.txt')
cts_na = CloseToStop('./archaea_neg80.txt')

writeFeature(cts_a, '+1')
writeFeature(cts_na, '-1')
