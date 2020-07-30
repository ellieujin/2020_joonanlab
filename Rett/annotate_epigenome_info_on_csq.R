library(bedr)
library(dplyr)
library(ggplot2)
library(cowplot)

###############################################################################################
################################### preliminary work ##########################################
###############################################################################################
SNV_broad = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_SNV_broad_narrow_compared.tsv', sep = '\t', header = T)
INDEL_broad = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_INDEL_broad_narrow_compared.tsv', sep = '\t', header = T)

# remove redundant columns
SNV_broad <- SNV_broad[,c(1,2,52:90)]  # this number selects only broad peaks
INDEL_broad <- INDEL_broad[,c(1,2,52:90)]

# make variants unique
SNV_broad = SNV_broad[!duplicated(SNV_broad$variant),]
INDEL_broad = INDEL_broad[!duplicated(INDEL_broad$variant),]

# mutate a column for each variant's total peak
SNV_broad <- SNV_broad %>% mutate(total = rowSums(SNV_broad[,c(3:41)]))
INDEL_broad <- INDEL_broad %>% mutate(total = rowSums(INDEL_broad[,c(3:41)]))

SNV_narrow = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_SNV_broad_narrow_compared.tsv', sep = '\t', header = T)
INDEL_narrow = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_INDEL_broad_narrow_compared.tsv', sep = '\t', header = T)

# remove redundant columns
SNV_narrow <- SNV_narrow[,c(1,2,91:129)]  # this number selects only broad peaks
INDEL_narrow <- INDEL_narrow[,c(1,2,91:129)]

# make variants unique
SNV_narrow = SNV_narrow[!duplicated(SNV_narrow$variant),]
INDEL_narrow = INDEL_narrow[!duplicated(INDEL_narrow$variant),]

# mutate a column for each variant's total peak
SNV_narrow <- SNV_narrow %>% mutate(total = rowSums(SNV_narrow[,c(3:41)]))
INDEL_narrow <- INDEL_narrow %>% mutate(total = rowSums(INDEL_narrow[,c(3:41)]))

###############################################################################################
################################### annotating epigenome on csq ###############################
###############################################################################################

SNV = read.table('/home/titan/Hail/rett/tables/20200213_rett_snv_csq_list.csv', sep = ',', header = T)
INDEL = read.table('/home/titan/Hail/rett/tables/20200213_rett_indel_csq_list.csv', sep = ',', header = T)

# select only more than 10 peaks or female fetal brain peak (making data frame to merge at csq)

SNV_temp = SNV_broad %>% filter(total > 10 | isRE_DHS_Fetal_Brain_Female == 1) %>% select(c(variant, total, isRE_DHS_Fetal_Brain_Female))
SNV_narrow = SNV_narrow %>% select(c(variant, total, narrow_isRE_DHS_Fetal_Brain_Female))
SNV_temp = merge(SNV_temp, SNV_narrow, by = 'variant', all.x = T)
names(SNV_temp)[names(SNV_temp) == 'total.x'] <- 'total_broad_peak'
names(SNV_temp)[names(SNV_temp) == 'total.y'] <- 'total_narrow_peak'

INDEL_temp = INDEL_broad %>% filter(total > 0 | isRE_DHS_Fetal_Brain_Female == 1) %>% select(c(variant, total, isRE_DHS_Fetal_Brain_Female))
INDEL_narrow = INDEL_narrow %>% select(c(variant, total, narrow_isRE_DHS_Fetal_Brain_Female))
INDEL_temp = merge(INDEL_temp, INDEL_narrow, by = 'variant', all.x = T)
names(INDEL_temp)[names(INDEL_temp) == 'total.x'] <- 'total_broad_peak'
names(INDEL_temp)[names(INDEL_temp) == 'total.y'] <- 'total_narrow_peak'

# merge epigenome information on csq file

SNV_final = merge(SNV_temp, SNV, by = 'variant', all.x = T)
SNV_final = SNV_final[!duplicated(SNV_final$Gene),]

INDEL_final = merge(INDEL_temp, INDEL, by = 'variant', all.x = T)
INDEL_final = INDEL_final[!duplicated(INDEL_final$Gene),]

# write table with final data

write.table(SNV_final, "/home/titan/Hail/rett/tables/20200224_rett_snv_epigenome_csq_list.tsv", sep = '\t', row.names = F)
write.table(INDEL_final, "/home/titan/Hail/rett/tables/20200224_rett_indel_epigenome_csq_list.tsv", sep = '\t', row.names = F)
