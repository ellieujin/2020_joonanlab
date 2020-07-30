import sys, os, glob
list_sample = glob.glob('/home/titan/mosdepth/bam/bambam/*.bam')
cmd = []
for s in list_sample:
    cmd.append(" --bam " + s)
all_cmd = ""
all_cmd = all_cmd.join(cmd)
print(all_cmd)
print("time ~/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py " + all_cmd + "  --referenceFasta /home/titan/resources/gatk_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --runDir /home/titan/manta/half_samples_200528/ --callRegions /home/titan/resources/interval_list/exome.bed.gz --exome")
os.system("time ~/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py " + all_cmd + "  --referenceFasta /home/titan/resources/gatk_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --runDir /home/titan/manta/half_samples_200528/ --callRegions /home/titan/resources/interval_list/exome.bed.gz --exome")
