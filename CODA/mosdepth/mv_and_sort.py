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
    print("time mkdir " + p) #run this cmd at mosdepth, not bam
    os.system("time mkdir " + p)
for a in sampleName:
    print("time mv /home/titan/mosdepth/bam/" + a + ".* /home/titan/mosdepth/" + a)
    os.system("time mv /home/titan/mosdepth/bam/" + a + ".* /home/titan/mosdepth/" + a)


