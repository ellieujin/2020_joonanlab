library(tidyverse)
library(bedr)
library(dplyr)
library(ggplot2)
library(cowplot)

csq_name = c('Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'MANE', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'SOURCE', 'GENE_PHENO', 'NEAREST', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'SV_overlap_AF', 'SV_overlap_PC', 'SV_overlap_name', 'gnomADg', 'gnomADg_AF')

## import data
hom = read_tsv("df_csq_hom.tsv")



## make table
make_table <- function(file__){
  n_variant = nrow(file__)/length(csq_name)
  tb = tibble()
  for (i in seq(n_variant)){
    temp = slice(file__, (73*(i-1)+1): (73*i)) %>% 
      select(variant, s, csq) %>%
      mutate(CSQ = csq_name) %>% 
      spread(CSQ, csq)
    tb = bind_rows(tb, temp)
  }
  return(tb)
}

hom = make_table(hom)


###삭제가 필요한 컬럼들
CSQ_t <- function(file__){
  file__ = file__ %>% select(-AA_AF, -AF, -AFR_AF, -AMR_AF, -EA_AF, -EAS_AF, -EUR_AF, 
                             -gnomAD_AF, -gnomAD_AFR_AF,-gnomAD_AMR_AF,-gnomAD_ASJ_AF,
                             -gnomAD_EAS_AF, -gnomAD_FIN_AF, -gnomAD_NFE_AF, -gnomAD_OTH_AF,
                             -gnomAD_SAS_AF, -gnomADg, MAX_AF,	-MAX_AF_POPS, -SAS_AF,
                             -Allele, -HGVS_OFFSET, -SOMATIC)
  return(file__)
}

hom = CSQ_t(hom)

hom = hom %>% filter(BIOTYPE == 'protein_coding')

###########################################################################################
####Feature_type 컬럼에서 RegulatoryFeature로 표시된 entry의 경우, 해당 유전자가 표기되는 대신 ENSR이라는 Ensembl의 regulatory(Feature)가 표시되어 있음. 따라서 NEAREST에 등장한 transcript ID를 가지고 유전자 이름을 붙이는 작업을 해야할 필요가 있음. 
##간단하게 하는 방법: VEP annotation에 사용된 GENCODE 버젼을 확인한 후, GENCODE 웹사이트에서 file을 다운로드.
##다음 gist를 이용하여 테이블을 만들어줌. https://gist.github.com/joonan30/3e3815c5937388c0365b3701db89a0f1
##이후 merge를 이용하여 gene name을 합쳐줌.

library(rtracklayer)
d <- rtracklayer::import('../../Resources/gencode/gencode.v32.annotation.gtf.gz') %>% 
  as.data.frame() %>% 
  filter(type == 'transcript') %>%
  select(transcript_id, gene_name)

split_ <- function(x){
  x = strsplit(x, "\\.")[[1]][1]
  return(x)
}
d['transcript_id'] = apply(d['transcript_id'], 1, split_)


transc_gene <- function(file__){
  file__rest = file__ %>% filter(Feature_type!='RegulatoryFeature'| is.na(Feature_type))
  file__rf = file__ %>% filter(Feature_type=='RegulatoryFeature') %>% left_join(d, by =c('NEAREST'='transcript_id')) %>% mutate(SYMBOL=gene_name) %>% select(-gene_name)
  file__ = rbind(file__rest, file__rf) %>% arrange(s, variant)
  return(file__)
}


snv = transc_gene(snv)
indel = transc_gene(indel)

write_csv(snv, "HQvariants_consequence_snv_het_p.csv")
write_csv(indel, "HQvariants_consequence_indel_het_p.csv")

