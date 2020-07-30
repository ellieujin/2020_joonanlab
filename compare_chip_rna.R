##  This script is to evaluate the correlation of H3K27ac and H3K27me3 with RNA-seq data
## Set the run
setwd("~/Dropbox/NeuroDevEpi/data_20190507/")
options(stringsAsFactors = F)
library(bedr)
library(dplyr)
library(ggplot2)
source('function_PermutationNeuroDevEpi.R')
library(rtracklayer) ###
library(cowplot)
## Set the run parameter
peak.cutoff <- 12
tss_type = 'longTX_upstreamOnly'
tss_dist = 2000

date = '20191217'
run_tag = paste('ExpPeak_',
                'cutoff_', peak.cutoff, '_',
                'tss_', tss_type, '_',
                'dist_', tss_dist, 
                sep='')

ages = c('E15', 'E18', 'P2')

###################################################
### PART 1 - data munging                       ###
### rerun from here with different peak cutoffs ###
###################################################
print(paste("Running permutation analysis for ", TF, "  /peak.cutoff=", peak.cutoff, sep=""))
print("Part 1 - data wrangling")

## Update the columns
cols_chipseq <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 
                  'signalValue', 'peak_pval', 'qValue', 'peak')



## RNA expression
f_rna = 'table.Siavash_20191004_logCPM_Matrix.rowMeans.txt'
e = read.delim(f_rna) %>% select(-gene_id) %>% dplyr::rename(E15='wt_e15', E18='wt_e18', P2='wt_p2') #rename at two different packages

## Read in H3K27ac for matching age  (H3K27ac=predict where the transcription occurs)
list_h3k27ac = list()
list_h3k27ac[['E15']] = "Sfd-Cx-E15-N2_bc4_GCACTA_L008_R1_001_peaks.narrowPeak"
list_h3k27ac[['E18']] = "Sfd-Cx-E18-N2_BC10_GGAGAA_L002_R1_001_peaks.narrowPeak"
list_h3k27ac[['P2']] = "P2-H3K27Ac_S51_L008_R1_001_peaks.narrowPeak"


## Read in H3K27me3 for matching age  (H3K27me3=where the transcription is suppressed)
list_h3k27me3 = list()
list_h3k27me3[['E15']] = "Sfd-Cx-E15-N2_bc5_ACCTCA_L008_R1_001_peaks.narrowPeak"
list_h3k27me3[['E18']] = "Sfd-Cx-E18-N1_H3K27me3_AGCATG_L001_R1_001_peaks.narrowPeak"
list_h3k27me3[['P2']] = "P2-H3K27me3_S52_L008_R1_001_peaks.narrowPeak"


##########change age to look for a different age#################
age<-"E15"

h3k27me3.tf <- read.table(list_h3k27me3[[age]], header=F)
colnames(h3k27me3.tf) <- cols_chipseq
h3k27me3.tf <- subset(h3k27me3.tf, h3k27me3.tf[,"qValue"]>= peak.cutoff) ##select qValue > 12
## table(h3k27me3.tf$chrom) ->lots of chrUn~~
h3k27me3.tf = h3k27me3.tf[h3k27me3.tf$chrom %in% c(paste('chr', seq(1:22), sep=''), 'chrX', 'chrY'),]

## TSS defined by Gencode Transcripts (TSS=transcription start site)
f_tss = paste('tssMatrix_',tss_type,'_',tss_dist,'.M23.20191217.txt',sep='')
print (paste('TSS matrix', f_tss, sep=': '))
tss.merged <- read.delim(f_tss)
tss.merged$isPLI995 = ifelse( tss.merged$isPLI995==1, ifelse( tss.merged$isASC_FDR01==1 | tss.merged$isDDD299==1, 0, 1), 0)
      ##??????what does PLI995 and all the other columns mean???


## assign expression status (present or absent)
tss.merged <- merge(tss.merged, e, by = 'gene_name', all.x=T)

## subset columns to create a bed file
tss.bed = tss.merged[,c('chrom', 'start', 'end', 'gene_name')]

#################################################
### PART 2 - assign features for each mm9 TSS ###
#################################################
print("Part 2 - promoter and distal peak mapping")
print(paste("tss type = ", tss_type, "  &  tss distance = ", tss_dist, sep=""))
# for each annotated promoter, report
# promoter overlap
# nearest non-promoter


## Assign H3K27ac enrichment
d_k27ac = list()
for (age in ages){
  print (age)
  h3k27ac.tf <- read.table(list_h3k27ac[[age]], header=F)
  colnames(h3k27ac.tf) <- cols_chipseq
  h3k27ac.tf <- subset(h3k27ac.tf, h3k27ac.tf[,"qValue"]>= peak.cutoff)
  h3k27ac.tf = h3k27ac.tf[h3k27ac.tf$chrom %in% c(paste('chr', seq(1:22), sep=''), 'chrX', 'chrY'),]

  peak.bed <- h3k27ac.tf %>% 
    select(1,2,3) %>%
    dplyr::rename(chr = 1, start = 2, end = 3) %>% 
    mutate(ID = paste(chr, start, end, sep='_'))

  ## assign H3K27ac binding
  d_k27ac[[age]] <- assign.peak(tss.data = tss.merged,     #TSS file with RNA expression information
                                tss.bed = tss.bed,         #TSS file in a bed format with gene name (contains promoter regions)
                                peak.bed = peak.bed,       #Chip-Seq peaks by TF
                                peak.name = paste('K27ac', age, sep='_'))
}


## Assign H3K27me3 enrichment
d_k27me3 = list()
for (age in ages){
  print (age)
  h3k27me3.tf <- read.table(list_h3k27me3[[age]], header=F)
  colnames(h3k27me3.tf) <- cols_chipseq
  h3k27me3.tf <- subset(h3k27me3.tf, h3k27me3.tf[,"qValue"]>= peak.cutoff)
  h3k27me3.tf = h3k27me3.tf[h3k27me3.tf$chrom %in% c(paste('chr', seq(1:22), sep=''), 'chrX', 'chrY'),]
  
  peak.bed <- h3k27me3.tf %>% select(1,2,3) %>% 
    dplyr::rename(chr = 1, start = 2, end = 3) %>%
    mutate(ID = paste(chr, start, end, sep='_')) 
  
  ## assign H3K27ac binding
  d_k27me3[[age]] <- assign.peak(tss.data = tss.merged, 
                                 tss.bed = tss.bed, 
                                 peak.bed = peak.bed, 
                                 peak.name = paste('K27me3', age, sep='_'))
}

############################################################
### PART 3 - making plot to compare ChIP-seq and RNA-seq ###
############################################################
print("Part 3 - making plot to compare ChIP-seq and RNA-seq ")


###################################################### H3K27ac - E15 ########################################################

# compare gene expression with promoter peaks and no peaks  ->   promoters with peaks have more expression than non-promoters
p1<-d_k27ac$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27ac_E15.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  ggplot(aes(isPromoter, E15)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +  # order of discrete data
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27ac Promoter") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# only promoter gene expression by gene type
d_k27ac$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27ac_E15.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  filter(isPromoter == 'promoter') %>%
  ggplot(aes(isPromoter, E15)) + geom_point(aes(color=gene_type)) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27ac Promoter by Gene Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# gene expression by mutated(PTV) or not
d_k27ac$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27ac_E15.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  mutate(isMutation = ifelse(!is.na(mut.ptv), 'PTV', 'Not Mutated')) %>%
  ggplot(aes(isPromoter, E15)) + geom_point(aes(color=isMutation)) +
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27ac Promoter by PTV") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# ---------------------------------------------promoter

# compare gene expression with enhancer peaks and no peaks  ->  enhancers with peaks have slightly more expression than non-enhancers
p2<-d_k27ac$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isEnhancer = ifelse(K27ac_E15.distal != '.', 'K27ac_enhancer', 'No_K27ac_enhancer')) %>%
  ggplot(aes(isEnhancer, E15)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_enhancer', 'No_K27ac_enhancer')) + 
  labs(y="Gene Expression") +
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27ac Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# compare gene expression with either promoter and enhancer peaks and non
p3<-d_k27ac$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27ac_E15.promoter != '.' & K27ac_E15.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter_and_Enhancer, E15)) + geom_boxplot() + 
  scale_x_discrete(limits = c('either', 'one or neither')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in Either K27ac Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))


# -----------------------------------------enhancer & either promoter and enhancer

H27ac_E15_cowp<-plot_grid(p1, p2, p3)
ggsave("/home/titan/Compare_Chip_RNA/plot/H27ac_E15_cowp_20200107.png",H27ac_E15_cowp)

H27ac_E15_facet<-d_k27ac$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27ac_E15.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  mutate(isEnhancer = ifelse(K27ac_E15.distal != '.', 'K27ac_enhancer', 'No_K27ac_enhancer')) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27ac_E15.promoter != '.' & K27ac_E15.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter, E15)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +
  facet_wrap(. ~ isEnhancer) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27ac Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))
ggsave("/home/titan/Compare_Chip_RNA/plot/H27ac_E15_facet_20200107.png",H27ac_E15_facet)


######################################################## H3K27me3 - E15 #############################################################

##me3 means the gene expression is suppressed ->  peaks should have less gene expression / 1 means it is me3 peak at promoter
d_k27me3$E15 %>% mutate(isPromoter = ifelse(K27me3_E15.promoter != '.', '1', '0')) %>% 
  ggplot(aes(isPromoter, E15)) + geom_boxplot()

# compare gene expression with promoter peaks and no peaks  ->   promoters with peaks have more expression than non-promoters
d_k27me3$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27me3_E15.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  ggplot(aes(isPromoter, E15)) + geom_jitter() + 
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +  # order of discrete data
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27me3 Promoter") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# only promoter gene expression by gene type
d_k27me3$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27me3_E15.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  filter(isPromoter == 'K27me3_promoter') %>%
  ggplot(aes(isPromoter, E15)) + geom_jitter(aes(color=gene_type)) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27me3 Promoter by Gene Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# gene expression by mutated(PTV) or not
d_k27me3$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27me3_E15.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  mutate(isMutation = ifelse(!is.na(mut.ptv), 'PTV', 'Not Mutated')) %>%
  ggplot(aes(isPromoter, E15)) + geom_jitter(aes(color=isMutation)) +
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27me3 Promoter by PTV") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# ---------------------------------------------promoter

# compare gene expression with enhancer peaks and no peaks  ->  enhancers with peaks have slightly more expression than non-enhancers
p5<-d_k27me3$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isEnhancer = ifelse(K27me3_E15.distal != '.', 'K27me3_enhancer', 'No_K27me3_enhancer')) %>%
  ggplot(aes(isEnhancer, E15)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_enhancer', 'No_K27me3_enhancer')) + 
  labs(y="Gene Expression") +
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27me3 Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# compare gene expression with either promoter and enhancer peaks and non
p6<-d_k27me3$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27me3_E15.promoter != '.' & K27me3_E15.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter_and_Enhancer, E15)) + geom_boxplot() + 
  scale_x_discrete(limits = c('either', 'one or neither')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in Either K27me3 Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))


# -----------------------------------------enhancer & either promoter and enhancer

H27me3_E15_cowp<-plot_grid(p4, p5, p6)
ggsave("/home/titan/Compare_Chip_RNA/plot/H27me3_E15_cowp_20200107.png",H27me3_E15_cowp)


H27me3_E15_facet<-d_k27me3$E15 %>% 
  subset(!is.na(E15)) %>%
  mutate(isPromoter = ifelse(K27me3_E15.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  mutate(isEnhancer = ifelse(K27me3_E15.distal != '.', 'K27me3_enhancer', 'No_K27me3_enhancer')) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27me3_E15.promoter != '.' & K27me3_E15.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter, E15)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +
  facet_wrap(. ~ isEnhancer) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E15 \n in K27me3 Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))
ggsave("/home/titan/Compare_Chip_RNA/plot/H27me3_E15_facet_20200107.png",H27me3_E15_facet)

###################################################### H3K27ac - E18 ########################################################

# compare gene expression with promoter peaks and no peaks  ->   promoters with peaks have more expression than non-promoters
p7<-d_k27ac$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27ac_E18.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  ggplot(aes(isPromoter, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +  # order of discrete data
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27ac Promoter") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# only promoter gene expression by gene type
d_k27ac$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27ac_E18.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  filter(isPromoter == 'promoter') %>%
  ggplot(aes(isPromoter, E18)) + geom_point(aes(color=gene_type)) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27ac Promoter by Gene Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# gene expression by mutated(PTV) or not
d_k27ac$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27ac_E18.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  mutate(isMutation = ifelse(!is.na(mut.ptv), 'PTV', 'Not Mutated')) %>%
  ggplot(aes(isPromoter, E18)) + geom_point(aes(color=isMutation)) +
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27ac Promoter by PTV") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# ---------------------------------------------promoter

# compare gene expression with enhancer peaks and no peaks  ->  enhancers with peaks have slightly more expression than non-enhancers
p8<-d_k27ac$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isEnhancer = ifelse(K27ac_E18.distal != '.', 'K27ac_enhancer', 'No_K27ac_enhancer')) %>%
  ggplot(aes(isEnhancer, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_enhancer', 'No_K27ac_enhancer')) + 
  labs(y="Gene Expression") +
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27ac Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# compare gene expression with either promoter and enhancer peaks and non
p9<-d_k27ac$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27ac_E18.promoter != '.' & K27ac_E18.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter_and_Enhancer, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('either', 'one or neither')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in Either K27ac Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))


# -----------------------------------------enhancer & either promoter and enhancer

H27ac_E18_cowp<-plot_grid(p7, p8, p9)
ggsave("/home/titan/Compare_Chip_RNA/plot/H27ac_E18_cowp_20200107.png",H27ac_E18_cowp)

H27ac_E18_facet<-d_k27ac$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27ac_E18.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  mutate(isEnhancer = ifelse(K27ac_E18.distal != '.', 'K27ac_enhancer', 'No_K27ac_enhancer')) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27ac_E18.promoter != '.' & K27ac_E18.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +
  facet_wrap(. ~ isEnhancer) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27ac Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))
ggsave("/home/titan/Compare_Chip_RNA/plot/H27ac_E18_facet_20200107.png",H27ac_E18_facet)


######################################################## H3K27me3 - E18 #############################################################

##me3 means the gene expression is suppressed ->  peaks should have less gene expression / 1 means it is me3 peak at promoter
d_k27me3$E18 %>% mutate(isPromoter = ifelse(K27me3_E18.promoter != '.', '1', '0')) %>% 
  ggplot(aes(isPromoter, E18)) + geom_boxplot()

# compare gene expression with promoter peaks and no peaks  ->   promoters with peaks have more expression than non-promoters
p10<-d_k27me3$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27me3_E18.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  ggplot(aes(isPromoter, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +  # order of discrete data
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27me3 Promoter") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# only promoter gene expression by gene type
d_k27me3$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27me3_E18.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  filter(isPromoter == 'K27me3_promoter') %>%
  ggplot(aes(isPromoter, E18)) + geom_jitter(aes(color=gene_type)) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27me3 Promoter by Gene Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# gene expression by mutated(PTV) or not
d_k27me3$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27me3_E18.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  mutate(isMutation = ifelse(!is.na(mut.ptv), 'PTV', 'Not Mutated')) %>%
  ggplot(aes(isPromoter, E18)) + geom_jitter(aes(color=isMutation)) +
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27me3 Promoter by PTV") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# ---------------------------------------------promoter

# compare gene expression with enhancer peaks and no peaks  ->  enhancers with peaks have slightly more expression than non-enhancers
p11<-d_k27me3$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isEnhancer = ifelse(K27me3_E18.distal != '.', 'K27me3_enhancer', 'No_K27me3_enhancer')) %>%
  ggplot(aes(isEnhancer, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_enhancer', 'No_K27me3_enhancer')) + 
  labs(y="Gene Expression") +
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27me3 Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# compare gene expression with either promoter and enhancer peaks and non
p12<-d_k27me3$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27me3_E18.promoter != '.' & K27me3_E18.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter_and_Enhancer, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('either', 'one or neither')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in Either K27me3 Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))


# -----------------------------------------enhancer & either promoter and enhancer

H27me3_E18_cowp<-plot_grid(p10, p11, p12)
ggsave("/home/titan/Compare_Chip_RNA/plot/H27me3_E18_cowp_20200107.png",H27me3_E18_cowp)


H27me3_E18_facet<-d_k27me3$E18 %>% 
  subset(!is.na(E18)) %>%
  mutate(isPromoter = ifelse(K27me3_E18.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  mutate(isEnhancer = ifelse(K27me3_E18.distal != '.', 'K27me3_enhancer', 'No_K27me3_enhancer')) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27me3_E18.promoter != '.' & K27me3_E18.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter, E18)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +
  facet_wrap(. ~ isEnhancer) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage E18 \n in K27me3 Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))
ggsave("/home/titan/Compare_Chip_RNA/plot/H27me3_E18_facet_20200107.png",H27me3_E18_facet)

###################################################### H3K27ac - P2 ########################################################

# compare gene expression with promoter peaks and no peaks  ->   promoters with peaks have more expression than non-promoters
p13<-d_k27ac$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27ac_P2.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  ggplot(aes(isPromoter, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +  # order of discrete data
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27ac Promoter") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# only promoter gene expression by gene type
d_k27ac$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27ac_P2.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  filter(isPromoter == 'promoter') %>%
  ggplot(aes(isPromoter, P2)) + geom_point(aes(color=gene_type)) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27ac Promoter by Gene Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# gene expression by mutated(PTV) or not
d_k27ac$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27ac_P2.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  mutate(isMutation = ifelse(!is.na(mut.ptv), 'PTV', 'Not Mutated')) %>%
  ggplot(aes(isPromoter, P2)) + geom_point(aes(color=isMutation)) +
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27ac Promoter by PTV") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# ---------------------------------------------promoter

# compare gene expression with enhancer peaks and no peaks  ->  enhancers with peaks have slightly more expression than non-enhancers
p14<-d_k27ac$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isEnhancer = ifelse(K27ac_P2.distal != '.', 'K27ac_enhancer', 'No_K27ac_enhancer')) %>%
  ggplot(aes(isEnhancer, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_enhancer', 'No_K27ac_enhancer')) + 
  labs(y="Gene Expression") +
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27ac Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# compare gene expression with either promoter and enhancer peaks and non
p15<-d_k27ac$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27ac_P2.promoter != '.' & K27ac_P2.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter_and_Enhancer, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('either', 'one or neither')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in Either K27ac Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))


# -----------------------------------------enhancer & either promoter and enhancer

H27ac_P2_cowp<-plot_grid(p13, p14, p15)
ggsave("/home/titan/Compare_Chip_RNA/plot/H27ac_P2_cowp_20200107.png",H27ac_P2_cowp)

H27ac_P2_facet<-d_k27ac$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27ac_P2.promoter != '.', 'K27ac_promoter', 'No_K27ac_promoter')) %>% 
  mutate(isEnhancer = ifelse(K27ac_P2.distal != '.', 'K27ac_enhancer', 'No_K27ac_enhancer')) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27ac_P2.promoter != '.' & K27ac_P2.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27ac_promoter', 'No_K27ac_promoter')) +
  facet_wrap(. ~ isEnhancer) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27ac Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))
ggsave("/home/titan/Compare_Chip_RNA/plot/H27ac_P2_facet_20200107.png",H27ac_P2_facet)


######################################################## H3K27me3 - P2 #############################################################

##me3 means the gene expression is suppressed ->  peaks should have less gene expression / 1 means it is me3 peak at promoter
d_k27me3$P2 %>% mutate(isPromoter = ifelse(K27me3_P2.promoter != '.', '1', '0')) %>% 
  ggplot(aes(isPromoter, P2)) + geom_boxplot()

# compare gene expression with promoter peaks and no peaks  ->   promoters with peaks have more expression than non-promoters
p16<-d_k27me3$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27me3_P2.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  ggplot(aes(isPromoter, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +  # order of discrete data
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27me3 Promoter") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# only promoter gene expression by gene type
d_k27me3$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27me3_P2.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  filter(isPromoter == 'K27me3_promoter') %>%
  ggplot(aes(isPromoter, P2)) + geom_jitter(aes(color=gene_type)) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27me3 Promoter by Gene Type") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# gene expression by mutated(PTV) or not
d_k27me3$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27me3_P2.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  mutate(isMutation = ifelse(!is.na(mut.ptv), 'PTV', 'Not Mutated')) %>%
  ggplot(aes(isPromoter, P2)) + geom_jitter(aes(color=isMutation)) +
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27me3 Promoter by PTV") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))

# ---------------------------------------------promoter

# compare gene expression with enhancer peaks and no peaks  ->  enhancers with peaks have slightly more expression than non-enhancers
p17<-d_k27me3$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isEnhancer = ifelse(K27me3_P2.distal != '.', 'K27me3_enhancer', 'No_K27me3_enhancer')) %>%
  ggplot(aes(isEnhancer, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_enhancer', 'No_K27me3_enhancer')) + 
  labs(y="Gene Expression") +
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27me3 Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))

# compare gene expression with either promoter and enhancer peaks and non
p18<-d_k27me3$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27me3_P2.promoter != '.' & K27me3_P2.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter_and_Enhancer, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('either', 'one or neither')) +
  labs(y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in Either K27me3 Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.title = element_text(size = 12), axis.text.x = element_text(size = 10, face = "bold"))


# -----------------------------------------enhancer & either promoter and enhancer

H27me3_P2_cowp<-plot_grid(p16, p17, p18)
ggsave("/home/titan/Compare_Chip_RNA/plot/H27me3_P2_cowp_20200107.png",H27me3_P2_cowp)


H27me3_P2_facet<-d_k27me3$P2 %>% 
  subset(!is.na(P2)) %>%
  mutate(isPromoter = ifelse(K27me3_P2.promoter != '.', 'K27me3_promoter', 'No_K27me3_promoter')) %>% 
  mutate(isEnhancer = ifelse(K27me3_P2.distal != '.', 'K27me3_enhancer', 'No_K27me3_enhancer')) %>%
  mutate(isPromoter_and_Enhancer = ifelse(K27me3_P2.promoter != '.' & K27me3_P2.distal != '.', 'either', 'one or neither')) %>%
  ggplot(aes(isPromoter, P2)) + geom_boxplot() + 
  scale_x_discrete(limits = c('K27me3_promoter', 'No_K27me3_promoter')) +
  facet_wrap(. ~ isEnhancer) +
  labs(x="",y="Gene Expression") + 
  ggtitle("mm9 Gene Expression at Stage P2 \n in K27me3 Promoter and Enhancer") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 15), axis.text.x = element_text(size = 10, face = "bold"))
ggsave("/home/titan/Compare_Chip_RNA/plot/H27me3_P2_facet_20200107.png",H27me3_P2_facet)

