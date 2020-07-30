import sys, os, glob
list_sample = glob.glob('/home/titan/mosdepth/bam/*.bam')
sample = []
for s in list_sample:
    sample.append(s.split('/')[5])
print(sample)
for p in sample:
    print("samtools index " + p)
    os.system("samtools index " + p)



