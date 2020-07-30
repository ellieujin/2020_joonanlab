import os,sys

f = sys.argv[1]
with open(f) as fh:
        for l in fh:
                if l[0] == '#':
                    sys.stdout.write(l)
                elif l[0] != '#':
                        info = l.split('\t')
                        QUAL = info[5]
                        FILTER = info[6]
                        num = [10,11,12,13,17]
                        for i in num:
                            if info[i].split(':')[0] == '1/1':
                                for j in range(9,19):
                                    if j == 18 : break
                                    if j != i:
                                        if info[j].split(':')[0] == '1/1': break
                                if j==18 and FILTER == 'PASS' and QUAL == '999' and info[2].split(':')[0] != 'MantaBND' : sys.stdout.write(l)

