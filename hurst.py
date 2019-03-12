#-*-conding:utf-8-*-
def cha_to_num(c):
    if c=='A':
        return '1'
    elif c=='G':
        return '2'
    elif c=='C':
        return '3'
    elif c=='T':
        return '4'

def trans_seq_to_num(inPath,outPath,codon):
    if codon==0:
        step=1
        start=0
    elif codon>=1 and codon<=3:
        step=3
        start=codon-1
    else:
        print("paramerer error")
        return -1
    with open(outPath,'w') as fo:
        with open(inPath) as fi:
            for line in fi.readlines():
                asNum = ''
                if line[0]=='>':
                    continue
                else:
                    line = line.strip().upper()
                    for i in range(start,len(line),step):
                        asNum += cha_to_num(line[i])
                    fo.write(asNum+'\n')

trans_seq_to_num('./data/archaea_pos80.txt','./data/hurst/archaea_pos80_num.txt',0)
trans_seq_to_num('./data/archaea_pos80.txt','./data/hurst/archaea_pos80_num1.txt',1)
trans_seq_to_num('./data/archaea_pos80.txt','./data/hurst/archaea_pos80_num2.txt',2)
trans_seq_to_num('./data/archaea_pos80.txt','./data/hurst/archaea_pos80_num3.txt',3)
trans_seq_to_num('./data/archaea_neg80.txt','./data/hurst/archaea_neg80_num.txt',0)
trans_seq_to_num('./data/archaea_neg80.txt','./data/hurst/archaea_neg80_num1.txt',1)
trans_seq_to_num('./data/archaea_neg80.txt','./data/hurst/archaea_neg80_num2.txt',2)
trans_seq_to_num('./data/archaea_neg80.txt','./data/hurst/archaea_neg80_num3.txt',3)
