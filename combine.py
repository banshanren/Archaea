from sys import argv

formation = 'svm'

if argv[len(argv)-2] == '-f':
    formation = argv[len(argv)-1]

point = 0
for i in range(1, len(argv)):
    if argv[i] == '-o':
        point = i
if point == 0:
    print("not set output file")
else:
    f = []
    fout = open(argv[point+1], 'w')
    for i in range(1, point):
        f.append(open(argv[i], 'r')) 
    if formation == 'svm':
        for line1 in f[0].readlines():
            line1 = line1.strip(' \n')
            line1_list = line1.split(' ')
            l1 = int(line1_list[-1].split(':')[0])
            # print('------------------',l1)
            for i in range(1,point-1):
                line2 = f[i].readline().strip(' \n').split(' ')
                l2 = len(line2)
                for i in range(1, l2):
                    list = line2[i].split(':')
                    line1 = line1 + ' ' + str(l1 + i) + ":" + list[1]
                l1 = l1+l2-1
            fout.write(line1 + '\n')
    elif formation == 'tab':
        for line1 in f[0].readlines():
            line1 = line1.strip(' \n')
            for i in range(1,point-1):
                line2 = f[i].readline().strip(' \n')
                line1 = line1 + ' ' + line2
            fout.write(line1 + '\n')
    else:
        print('wrong formation')
        exit()
    for i in range(0,point-3):
        f[i].close()
    fout.close()


