library(tidyverse)
library(bedr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
options(stringsAsFactors = F)

# read vcf
dnv_het <- read.vcf("/titan2/UDP_SV/output/vep_UDP_SV_dnv_het_200715.vcf")
rare_hom <- read.vcf("/titan2/UDP_SV/output/vep_UDP_SV_rare_hom_200715.vcf")
single_het <- read.vcf("/titan2/UDP_SV/output/vep_UDP_SV_single_het_200715.vcf")
single_hom <- read.vcf("/titan2/UDP_SV/output/vep_UDP_SV_single_hom_200715.vcf")


# remove multi allelic
# gnomad AF == 0 or NA  --> these are filtered by hail, but nothing filtered...

# filter only chr 1~22, X, Y
# filter only DEL, DUP?  / remove BND(translocation), INS
# coding region
dnv_het$vcf = dnv_het$vcf %>%                               # 188 -> 49(rm BND) -> 32(rm INS) -> 22(protein coding only)
  filter(!str_detect(ID, regex('MantaBND'))) %>%
  filter(!str_detect(ID, regex('MantaINS'))) %>%
  filter(str_detect(INFO, regex('protein_coding'))) %>%
  filter(!str_detect(CHROM, fixed('v1')))
rare_hom$vcf = rare_hom$vcf %>%                             # 675 -> 628 -> 422 -> 186
  filter(!str_detect(ID, regex('MantaBND'))) %>%
  filter(!str_detect(ID, regex('MantaINS'))) %>%
  filter(str_detect(INFO, regex('protein_coding'))) %>%
  filter(!str_detect(CHROM, fixed('v1')))
single_het$vcf = single_het$vcf %>%                             # 297 -> 173 -> 113 -> 58
  filter(!str_detect(ID, regex('MantaBND'))) %>%
  filter(!str_detect(ID, regex('MantaINS'))) %>%
  filter(str_detect(INFO, regex('protein_coding'))) %>%
  filter(!str_detect(CHROM, fixed('v1')))
single_hom$vcf = single_hom$vcf %>%                             # 1080 -> 1014 -> 619 -> 275
  filter(!str_detect(ID, regex('MantaBND'))) %>%
  filter(!str_detect(ID, regex('MantaINS'))) %>%
  filter(str_detect(INFO, regex('protein_coding'))) %>%
  filter(!str_detect(CHROM, fixed('v1')))

# make table by samples
# single het 5, 9, 12
#select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_5-1`, -`wgs_9-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

single_het_5 = single_het
single_het_5$vcf = single_het$vcf %>%  # 58 -> 27
  filter(str_detect(`wgs_5-1`, regex('0/1'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_9-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

single_het_9 = single_het
single_het_9$vcf = single_het$vcf %>%  # 58 -> 13
  filter(str_detect(`wgs_9-1`, regex('0/1'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_5-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

single_het_12 = single_het
single_het_12$vcf = single_het$vcf %>%  # 58 -> 18
  filter(str_detect(`wgs_12-1`, regex('0/1'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_5-1`, -`wgs_9-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`)

# single hom

single_hom_5 = single_hom
single_hom_5$vcf = single_hom$vcf %>%  # 275 -> 69
  filter(str_detect(`wgs_5-1`, regex('1/1'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_9-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

single_hom_9 = single_hom
single_hom_9$vcf = single_hom$vcf %>%  # 275 -> 62
  filter(str_detect(`wgs_9-1`, regex('1/1'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_5-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

single_hom_12 = single_hom
single_hom_12$vcf = single_hom$vcf %>%  # 275 -> 65
  filter(str_detect(`wgs_12-1`, regex('1/1'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_5-1`, -`wgs_9-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`)

# dnv het 3, 8

dnv_het_3 = dnv_het
dnv_het_3$vcf = dnv_het$vcf %>%   # 22 -> 13
  filter(str_detect(`wgs_3-1`, regex('0/1'))) %>%
  filter(str_detect(`wgs_3-2`, regex('0/0'))) %>%
  filter(str_detect(`wgs_3-3`, regex('0/0'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_5-1`, -`wgs_9-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

dnv_het_8 = dnv_het
dnv_het_8$vcf = dnv_het$vcf %>%  # 22 -> 9
  filter(str_detect(`wgs_8-1`, regex('0/1'))) %>%
  filter(str_detect(`wgs_8-2`, regex('0/0'))) %>%
  filter(str_detect(`wgs_8-3`, regex('0/0'))) %>%
  select(-`wgs_3-2`, -`wgs_5-1`, -`wgs_9-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

# rare hom
rare_hom_3 = rare_hom
rare_hom_3$vcf = rare_hom$vcf %>%   # 186 -> 103
  filter(str_detect(`wgs_3-1`, regex('1/1'))) %>%
  filter(str_detect(`wgs_3-2`, regex('0/1'))) %>%
  filter(str_detect(`wgs_3-3`, regex('0/1'))) %>%
  select(-`wgs_3-2`, -`wgs_8-1`, -`wgs_5-1`, -`wgs_9-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)

rare_hom_8 = rare_hom
rare_hom_8$vcf = rare_hom$vcf %>%  # 186 -> 85
  filter(str_detect(`wgs_8-1`, regex('1/1'))) %>%
  filter(str_detect(`wgs_8-2`, regex('0/1'))) %>%
  filter(str_detect(`wgs_8-3`, regex('0/1'))) %>%
  select(-`wgs_3-2`, -`wgs_5-1`, -`wgs_9-1`, -`wgs_3-1`, -`wgs_8-3`,  -`wgs_3-3`,  -`wgs_8-2`, -`wgs_12-1`)


# write table to do hail
write.vcf(dnv_het_3, "/home/titan/Hail/UDP/temp_dnv_het_3.vcf")
write.vcf(dnv_het_8, "/home/titan/Hail/UDP/temp_dnv_het_8.vcf")
write.vcf(rare_hom_3, "/home/titan/Hail/UDP/temp_rare_hom_3.vcf")
write.vcf(rare_hom_8, "/home/titan/Hail/UDP/temp_rare_hom_8.vcf")

write.vcf(single_het_5, "/home/titan/Hail/UDP/temp_single_het_5.vcf")
write.vcf(single_het_9, "/home/titan/Hail/UDP/temp_single_het_9.vcf")
write.vcf(single_het_12, "/home/titan/Hail/UDP/temp_single_het_12.vcf")
write.vcf(single_hom_5, "/home/titan/Hail/UDP/temp_single_hom_5.vcf")
write.vcf(single_hom_9, "/home/titan/Hail/UDP/temp_single_hom_9.vcf")
write.vcf(single_hom_12, "/home/titan/Hail/UDP/temp_single_hom_12.vcf")

############################################## read result table from hail ########################################################

csq_single_het_5 = read_tsv("/titan2/UDP_SV/df_csq_single_het_5.tsv") # 68766
csq_single_het_9 = read_tsv("/titan2/UDP_SV/df_csq_single_het_9.tsv") # 28178
csq_single_het_12 = read_tsv("/titan2/UDP_SV/df_csq_single_het_12.tsv") # 17228

csq_single_hom_5 = read_tsv("/titan2/UDP_SV/df_csq_single_hom_5.tsv") # 30222
csq_single_hom_9 = read_tsv("/titan2/UDP_SV/df_csq_single_hom_9.tsv") # 31974
csq_single_hom_12 = read_tsv("/titan2/UDP_SV/df_csq_single_hom_12.tsv") # 34164

csq_dnv_het_3_1 = read_tsv("/titan2/UDP_SV/df_csq_dnv_het_3_1.tsv") #### 264260
csq_dnv_het_3_2 = read_tsv("/titan2/UDP_SV/df_csq_dnv_het_3_2.tsv") ####
csq_dnv_het_3_3 = read_tsv("/titan2/UDP_SV/df_csq_dnv_het_3_3.tsv") #### 9417

csq_dnv_het_8 = read_tsv("/titan2/UDP_SV/df_csq_dnv_het_8.tsv") # 2920

csq_rare_hom_3 = read_tsv("/titan2/UDP_SV/df_csq_rare_hom_3.tsv")  # 45041
csq_rare_hom_8 = read_tsv("/titan2/UDP_SV/df_csq_rare_hom_8.tsv")  # 257033


# make table (spread)
csq_name = c('Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'MANE', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'SOURCE', 'GENE_PHENO', 'NEAREST', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'SV_overlap_AF', 'SV_overlap_PC', 'SV_overlap_name', 'gnomADg', 'gnomADg_AF')

make_table <- function(file__){
  n_variant = nrow(file__)/length(csq_name)
  tb = tibble()
  for (i in seq(n_variant)){
    temp = slice(file__, (73*(i-1)+1): (73*i)) %>% 
      mutate(END = unlist(strsplit(info, split = "\"END\":"))[2]) %>%
      mutate(END = unlist(strsplit(END, split = ",\"CIPOS\":"))[1]) %>%
      select(variant, s, csq, GQ, PL, PR, SR, END) %>%
      mutate(CSQ = csq_name) %>% 
      spread(CSQ, csq)
    tb = bind_rows(tb, temp)
  }
  return(tb)
}

ccsq_single_het_5 = make_table(csq_single_het_5) # 942, 80
ccsq_single_het_9 = make_table(csq_single_het_9) # 386
ccsq_single_het_12 = make_table(csq_single_het_12) # 236

ccsq_single_hom_5 = make_table(csq_single_hom_5) #414
ccsq_single_hom_9 = make_table(csq_single_hom_9) #438
ccsq_single_hom_12 = make_table(csq_single_hom_12) #468,

ccsq_dnv_het_3_1 = make_table(csq_dnv_het_3_1)  #### 3620
# ccsq_dnv_het_3_2 = make_table(csq_dnv_het_3_2)  ####
ccsq_dnv_het_3_3 = make_table(csq_dnv_het_3_3)  #### 129
ccsq_dnv_het_8 = make_table(csq_dnv_het_8) #40

ccsq_rare_hom_3 = make_table(csq_rare_hom_3)  # 617
ccsq_rare_hom_8 = make_table(csq_rare_hom_8) #3521


###삭제가 필요한 컬럼들
CSQ_t <- function(file__){
  file__ = file__ %>% select(-AA_AF, -AF, -AFR_AF, -AMR_AF, -EA_AF, -EAS_AF, -EUR_AF, 
                             -gnomAD_AF, -gnomAD_AFR_AF,-gnomAD_AMR_AF,-gnomAD_ASJ_AF,
                             -gnomAD_EAS_AF, -gnomAD_FIN_AF, -gnomAD_NFE_AF, -gnomAD_OTH_AF,
                             -gnomAD_SAS_AF, -gnomADg, MAX_AF,	-MAX_AF_POPS, -SAS_AF,
                             -Allele, -HGVS_OFFSET, -SOMATIC, -DISTANCE, -miRNA, -MOTIF_NAME, 
                             -MOTIF_POS, -MOTIF_SCORE_CHANGE, -PolyPhen, -SIFT)
  file__ = file__[,c("variant","END","s","Consequence","SYMBOL","SYMBOL_SOURCE" ,"GQ","PL","PR","SR",        
                     "Amino_acids","APPRIS","BIOTYPE","CANONICAL","CCDS","cDNA_position","CDS_position",
                     "CLIN_SIG","Codons","DOMAINS","ENSP","Existing_variation","EXON","Feature","Feature_type",
                     "FLAGS","Gene","GENE_PHENO", "gnomADg_AF","HGNC_ID","HGVSc", "HGVSp","HIGH_INF_POS",
                     "IMPACT", "INTRON","MANE", "MAX_AF", "NEAREST","PHENO","Protein_position","PUBMED",
                     "SOURCE","STRAND", "SV_overlap_AF","SV_overlap_name","SV_overlap_PC", "SWISSPROT","TREMBL",          
                     "TSL","UNIPARC","VARIANT_CLASS")]
  return(file__)
}


ccsq_single_het_5 = CSQ_t(ccsq_single_het_5) # only column decreases to 51
ccsq_single_het_9 = CSQ_t(ccsq_single_het_9)
ccsq_single_het_12 = CSQ_t(ccsq_single_het_12)

ccsq_single_hom_5 = CSQ_t(ccsq_single_hom_5)
ccsq_single_hom_9 = CSQ_t(ccsq_single_hom_9)
ccsq_single_hom_12 = CSQ_t(ccsq_single_hom_12)

ccsq_dnv_het_3_1 = CSQ_t(ccsq_dnv_het_3_1) ####
# ccsq_dnv_het_3_2 = CSQ_t(ccsq_dnv_het_3_2) ####
ccsq_dnv_het_3_3 = CSQ_t(ccsq_dnv_het_3_3) ####
ccsq_dnv_het_8 = CSQ_t(ccsq_dnv_het_8)

ccsq_rare_hom_3 = CSQ_t(ccsq_rare_hom_3)  # 617,53
ccsq_rare_hom_8 = CSQ_t(ccsq_rare_hom_8)


####    1.    ######################## protein conding only + SV AF == NA ####################
# gnomad SV 1, 2 .... reported -> remove?

##SV_overlap_AF=frequency of overlapping SV
##SV_overlap_PC=percent of input SV covered by reference SV
##SV_overlap_name=name of overlapping SV

filter_pc_na <- function(file__){
  file__ = file__ %>% 
    filter(BIOTYPE == 'protein_coding') %>%
    filter((is.na(SV_overlap_AF)) == 1)
  return(file__)
}

ccsq_single_het_5 = filter_pc_na(ccsq_single_het_5) # 5 ,53
ccsq_single_het_9 = filter_pc_na(ccsq_single_het_9) # 5, 53
ccsq_single_het_12 = filter_pc_na(ccsq_single_het_12) # 26,53

ccsq_single_hom_5 = filter_pc_na(ccsq_single_hom_5) # 44,53
ccsq_single_hom_9 = filter_pc_na(ccsq_single_hom_9) # 44,53
ccsq_single_hom_12 = filter_pc_na(ccsq_single_hom_12) # 60,53

ccsq_dnv_het_3_1 = filter_pc_na(ccsq_dnv_het_3_1) #### 63
ccsq_dnv_het_3_2 = filter_pc_na(ccsq_dnv_het_3_2) ####
ccsq_dnv_het_3_3 = filter_pc_na(ccsq_dnv_het_3_3) #### 13

ccsq_dnv_het_8 = filter_pc_na(ccsq_dnv_het_8) # 9, 53

ccsq_rare_hom_3 = filter_pc_na(ccsq_rare_hom_3)  # 78, 53
ccsq_rare_hom_8 = filter_pc_na(ccsq_rare_hom_8) #93,53



write_csv(ccsq_single_het_5, "/home/titan/Hail/UDP/output/single_het_5.csv")
write_csv(ccsq_single_het_9, "/home/titan/Hail/UDP/output/single_het_9.csv")
write_csv(ccsq_single_het_12, "/home/titan/Hail/UDP/output/single_het_12.csv")

write_csv(ccsq_single_hom_5, "/home/titan/Hail/UDP/output/single_hom_5.csv")
write_csv(ccsq_single_hom_9, "/home/titan/Hail/UDP/output/single_hom_9.csv")
write_csv(ccsq_single_hom_12, "/home/titan/Hail/UDP/output/single_hom_12.csv")

write_csv(ccsq_dnv_het_3_1, "/home/titan/Hail/UDP/output/dnv_het_3_1.csv")
#write_csv(ccsq_dnv_het_3_2, "/home/titan/Hail/UDP/output/dnv_het_3_2.csv")
write_csv(ccsq_dnv_het_3_3, "/home/titan/Hail/UDP/output/dnv_het_3_3.csv")
write_csv(ccsq_dnv_het_8, "/home/titan/Hail/UDP/output/dnv_het_8.csv")

write_csv(ccsq_rare_hom_3, "/home/titan/Hail/UDP/output/rare_hom_3.csv")
write_csv(ccsq_rare_hom_8, "/home/titan/Hail/UDP/output/rare_hom_8.csv")


####    2.    ######################## GQ == 999  +  coding consequence + biotype == protein coding ####################

filter_pc_gq <- function(file__){
  file__ = file__ %>% 
    filter(!str_detect(Consequence, fixed('TF'))) %>%
    #filter(str_detect(Consequence, fixed('coding_sequence_variant')) | str_detect(Consequence, fixed('transcript_ablation'))) %>%
    filter(GQ > 900)
  return(file__)
}

table(ccsq_single_het_5$Consequence)
table(ccsq_single_het_5$Consequence, ccsq_single_het_5$IMPACT)

hist_5<-ccsq_single_het_5 %>% 
  #filter(str_detect(Consequence, fixed('coding_sequence_variant')) | str_detect(Consequence, fixed('transcript_ablation'))) %>%
  filter(!str_detect(Consequence, fixed('regulatory'))) %>%
  filter(!str_detect(Consequence, fixed('non_coding'))) %>%
  filter(!str_detect(Consequence, fixed('TF'))) %>%
  filter(!str_detect(Consequence, fixed('truncation'))) %>%
           ggplot() + 
           geom_violin(aes(Consequence,GQ)) + 
           ggtitle("single het 5") + 
           theme(axis.text.x=element_text(angle=20, hjust=1))

hist_12<-ccsq_single_het_12 %>% 
  #filter(str_detect(Consequence, fixed('coding_sequence_variant')) | str_detect(Consequence, fixed('transcript_ablation'))) %>%
  filter(!str_detect(Consequence, fixed('regulatory'))) %>%
  filter(!str_detect(Consequence, fixed('non_coding'))) %>%
  filter(!str_detect(Consequence, fixed('TF'))) %>%
  filter(!str_detect(Consequence, fixed('truncation'))) %>%
  ggplot() + 
  geom_violin(aes(Consequence,GQ)) + 
  ggtitle("single het 5") + 
  theme(axis.text.x=element_text(angle=20, hjust=1))

hist_8<-ccsq_rare_hom_8 %>% 
  #filter(str_detect(Consequence, fixed('coding_sequence_variant')) | str_detect(Consequence, fixed('transcript_ablation'))) %>%
  filter(!str_detect(Consequence, fixed('regulatory'))) %>%
  filter(!str_detect(Consequence, fixed('non_coding'))) %>%
  filter(!str_detect(Consequence, fixed('TF'))) %>%
  filter(!str_detect(Consequence, fixed('truncation'))) %>%
  ggplot() + 
  geom_violin(aes(Consequence,GQ)) + 
  ggtitle("single het 5") + 
  theme(axis.text.x=element_text(angle=20, hjust=1))

hist_5
hist_12
hist_8

ggsave("/home/titan/hist_5.png", hist_5)
ggsave("/home/titan/hist_8.png", hist_8)
ggsave("/home/titan/hist_12.png", hist_12)

cccsq_single_het_5 = filter_pc_gq(ccsq_single_het_5) # 5 ,53
cccsq_single_het_9 = filter_pc_gq(ccsq_single_het_9) # 5, 53
cccsq_single_het_12 = filter_pc_gq(ccsq_single_het_12) # 26,53

cccsq_single_hom_5 = filter_pc_gq(ccsq_single_hom_5) # 44,53
cccsq_single_hom_9 = filter_pc_gq(ccsq_single_hom_9) # 44,53
cccsq_single_hom_12 = filter_pc_gq(ccsq_single_hom_12) # 60,53

cccsq_dnv_het_3_1 = filter_pc_gq(ccsq_dnv_het_3_1) #### 63
#cccsq_dnv_het_3_2 = filter_pc_gq(ccsq_dnv_het_3_2) ####
cccsq_dnv_het_3_3 = filter_pc_gq(ccsq_dnv_het_3_3) #### 13

cccsq_dnv_het_8 = filter_pc_gq(ccsq_dnv_het_8) # 9, 53

cccsq_rare_hom_3 = filter_pc_gq(ccsq_rare_hom_3)  # 78, 53
cccsq_rare_hom_8 = filter_pc_gq(ccsq_rare_hom_8) #93,53


write_csv(cccsq_single_het_5, "/home/titan/Hail/UDP/output/gq_HQvariants_consequence_single_het_5.csv")
write_csv(cccsq_single_het_9, "/home/titan/Hail/UDP/output/gq_HQvariants_consequence_single_het_9.csv")
write_csv(cccsq_single_het_12, "/home/titan/Hail/UDP/output/gq_HQvariants_consequence_single_het_12.csv")


## gnomad AF annotate

#### wrtie csv
write_csv(cccsq_single_het_5, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_single_het_5.csv")
write_csv(cccsq_single_het_9, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_single_het_9.csv")
write_csv(cccsq_single_het_12, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_single_het_12.csv")

write_csv(cccsq_single_hom_5, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_single_hom_5.csv")
write_csv(cccsq_single_hom_9, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_single_hom_9.csv")
write_csv(cccsq_single_hom_12, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_single_hom_12.csv")

write_csv(cccsq_dnv_het_3_1, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_dnv_het_3_1.csv")
#write_csv(cccsq_dnv_het_3_2, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_dnv_het_3_2.csv")
write_csv(cccsq_dnv_het_3_3, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_dnv_het_3_3.csv")
write_csv(cccsq_dnv_het_8, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_dnv_het_8.csv")

write_csv(cccsq_rare_hom_3, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_rare_hom_3.csv")
write_csv(cccsq_rare_hom_8, "/home/titan/Hail/UDP/output/200728_HQvariants_consequence_rare_hom_8.csv")


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





#'CHROM',  'POS',     'ID',      'REF',    'ALT',     'QUAL',    'FILTER',  'INFO',    'FORMAT',  'wgs_3-2', 'wgs_8-1', 'wgs_5-1', 'wgs_9-1', 'wgs_3-1', wgs_8-3',  'wgs_3-3',  'wgs_8-2', 'wgs_12-1'
csq_name = c('Allele', 'Consequence', 'IMPACT', 'SYMBOL', 'Gene', 'Feature_type', 'Feature', 'BIOTYPE', 'EXON', 'INTRON', 'HGVSc', 'HGVSp', 'cDNA_position', 'CDS_position', 'Protein_position', 'Amino_acids', 'Codons', 'Existing_variation', 'DISTANCE', 'STRAND', 'FLAGS', 'VARIANT_CLASS', 'SYMBOL_SOURCE', 'HGNC_ID', 'CANONICAL', 'MANE', 'TSL', 'APPRIS', 'CCDS', 'ENSP', 'SWISSPROT', 'TREMBL', 'UNIPARC', 'SOURCE', 'GENE_PHENO', 'NEAREST', 'SIFT', 'PolyPhen', 'DOMAINS', 'miRNA', 'HGVS_OFFSET', 'AF', 'AFR_AF', 'AMR_AF', 'EAS_AF', 'EUR_AF', 'SAS_AF', 'AA_AF', 'EA_AF', 'gnomAD_AF', 'gnomAD_AFR_AF', 'gnomAD_AMR_AF', 'gnomAD_ASJ_AF', 'gnomAD_EAS_AF', 'gnomAD_FIN_AF', 'gnomAD_NFE_AF', 'gnomAD_OTH_AF', 'gnomAD_SAS_AF', 'MAX_AF', 'MAX_AF_POPS', 'CLIN_SIG', 'SOMATIC', 'PHENO', 'PUBMED', 'MOTIF_NAME', 'MOTIF_POS', 'HIGH_INF_POS', 'MOTIF_SCORE_CHANGE', 'SV_overlap_AF', 'SV_overlap_PC', 'SV_overlap_name', 'gnomADg', 'gnomADg_AF')
CIEND: array<str>, 
CIPOS: array<str>, 
CHR2: str, 
END: int32, 
MAPQ: int32, 
RE: int32, 
IMPRECISE: bool, 
PRECISE: bool, 
SVLEN: int32, 
SVMETHOD: str, 
SVTYPE: str, 
SUPP_VEC: str, 
SUPP: str, 
STRANDS: str, 
CSQ: array<str>, 
gnomADg: array<str>, 
gnomADg_AF: array<str>
