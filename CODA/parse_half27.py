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
                        num = [13,15,16,19,24,27,28,29,31,32]
                        for i in num:
                            if info[i].split(':')[0] == '0/1':
                                for j in range(9,35): 
                                    if j != i:
                                        if info[j].split(':')[0] == '0/1' or info[j].split(':')[0] == '1/1': break
                                if j==34 and FILTER == 'PASS' and QUAL == '999': sys.stdout.write(l)
