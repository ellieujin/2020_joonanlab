from difflib import context_diff

with open('vep_CODA_DNV_SV_200630_GRCH37.vcf', 'r') as f1:
    with open('vep_CODA_DNV_SV_200630_hg19.vcf', 'r') as f2:
        diff = context_diff(f1.readlines(), f2.readlines(), fromfile='f1', tofile='f2')
        for line in diff:
            print(line)
