import os,sys

f = sys.argv[1]
with open(f) as fh:
        for l in fh:
                if l[0] == '#':
                    sys.stdout.write(l)
                if l[0] != '#':
                        info = l.split('\t')
                        QUAL = info[5]
                        FILTER = info[6]
                        num = [10, 14, 16,17,19,22,23,25,29,32,33]
                        for i in num:
                            if info[i].split(':')[0] == '0/1':
                                for j in range(9,34):
                                    if j != i:
                                        if info[j].split(':')[0] == '0/1' or info[j].split(':')[0] == '1/1': break
                                if j==33 and FILTER == 'PASS' and QUAL == '999': sys.stdout.write(l)

