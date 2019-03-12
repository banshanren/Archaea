with open('./archaea_neg80_fixed_len.txt','w') as fneg_new:
    with open("./archaea_neg80.txt") as fneg:
        lines = fneg.readlines()
        flag = 10000
        for i in range(len(lines)):
            if lines[i][0] == '>':
                fneg_new.write(lines[i])
            else:
                fneg_new.write(lines[i][0:95]+'\n')
        #         print(len(lines[i]))
        #         flag = min(flag,len(lines[i]))
        # print(flag)