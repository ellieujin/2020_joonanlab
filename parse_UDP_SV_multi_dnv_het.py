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
                        if info[13].split(':')[0] == '0/1' and info[9].split(':')[0] == '0/0' and info[15].split(':')[0] == '0/0' and FILTER == 'PASS' and QUAL == '999' and info[2].split(':')[0] != 'MantaBND' : sys.stdout.write(l)
                        elif info[10].split(':')[0] == '0/1' and info[16].split(':')[0] == '0/0' and info[14].split(':')[0] == '0/0' and FILTER == 'PASS' and QUAL == '999' and info[2].split(':')[0] != 'MantaBND' : sys.stdout.write(l)

