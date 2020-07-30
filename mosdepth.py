import sys, os, glob, re
list_sample = glob.glob('/home/titan/mosdepth/bam/*.bam')
sample = []
for s in list_sample:
    sample.append(s.split('/')[5])
print(sample)
sampleName = []
for n in sample:
    name = re.sub("_sorted.bam", "", n)
    sampleName.append(name)
print(sampleName)
for p in sampleName:
    print("time mosdepth --by /home/titan/resources/interval_list/exome.bed.gz " + p + " /home/titan/mosdepth/bam/" + p + "_sorted.bam")
    os.system("time mosdepth --by /home/titan/resources/interval_list/exome.bed.gz " + p + " /home/titan/mosdepth/bam/" + p + "_sorted.bam")


