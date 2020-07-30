import os,sys

f = sys.argv[1]
with open(f) as fh:
        for l in fh:
                if l[0] == '#':
                    sys.stdout.write(l)
                elif l[0] != '#':
                    info = l.split('\t')
                    CHROM = info[0]
                    REF = info[3]
                    if REF != 'T' and CHROM != 'chr1': sys.stdout.write(l)

