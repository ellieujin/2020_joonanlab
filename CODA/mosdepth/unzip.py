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
for a in sampleName:
    print("time gzip -d /home/titan/mosdepth/" + a + "/*.gz")
    os.system("time gzip -d /home/titan/mosdepth/" + a + "/*.gz")


