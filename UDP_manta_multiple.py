import sys, os, glob
list_sample = glob.glob('/titan2/UDP/*.bam')
cmd = []
for s in list_sample:
    cmd.append(" --bam " + s)
all_cmd = ""
all_cmd = all_cmd.join(cmd)
print(all_cmd)
print("time ~/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py " + all_cmd + "  --referenceFasta /home/titan/resources/gatk_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --runDir /titan2/UDP/UDP_SV_200710/")
os.system("time ~/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py " + all_cmd + " --referenceFasta /home/titan/resources/gatk_bundle/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta --runDir /titan2/UDP/UDP_SV_200710/")
