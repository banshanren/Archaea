#coding:utf-8

# import numpy as np
# import pandas as pd
# import matplotlib as mpl
# import matplotlib.pyplot as plt

def Frequency(seq):
    gc = 0
    total = 0
    for i in seq:
        if i == 'G' or i == 'C':                                   #??
            gc += 1
        total += 1
    return gc/total
        
def GCcontent(fileName):
    freqList = []
    with open(fileName) as f:
        for line in f.readlines():
            if line[0] == '>':
                continue
            else:
                freqList.append(Frequency(line))
    return freqList

def writeFeature(featureList, label):
    with open('./GCfeature'+label+'.txt', 'w') as f:
        for i in featureList:
            f.write(label+' ')
            f.write('1:'+str(i))
            f.write('\n')

gc_a = GCcontent('./archaea_pos.txt')
gc_na = GCcontent('./archaea_neg.txt')

writeFeature(gc_a, '+1')
writeFeature(gc_na, '-1')

# with open('./GCfeature_pos.txt','w') as f:
#     for i in gc_a:
#         f.write('+1'+' ')
#         f.write('1:'+str(i))
#         f.write('\n')

# with open('./GCfeature_neg.txt','w') as f:
#     for i in gc_na:
#         f.write('-1'+' ')
#         f.write('1:'+str(i))
#         f.write('\n')



# gcArray_a = np.array(GCcontent('./archaea_pos.txt'))
# gcArray_na = np.array(GCcontent('./archaea_neg.txt'))

# plt.hist(gcArray_a,range=(0.2,0.6), bins=40, label='positive', edgecolor='k', alpha = 0.7)
# plt.hist(gcArray_na,range=(0.2,0.6), bins=40, label='negtive', edgecolor='k', alpha = 0.6)
# plt.title('GC frequency histgram')
# plt.xlabel('frequency')
# plt.ylabel('sample number')

# plt.legend()
# plt.show()                                 
