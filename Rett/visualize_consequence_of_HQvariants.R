library(bedr)
library(dplyr)
library(ggplot2)
library(cowplot)

####################################################################################################
################################### only with broad peaks ##########################################
####################################################################################################

SNV_broad = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_SNV_broad_narrow_compared.tsv', sep = '\t', header = T)
INDEL_broad = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_INDEL_broad_narrow_compared.tsv', sep = '\t', header = T)

##### 1. making unique variants

# remove redundant columns
SNV_broad <- SNV_broad[,c(1,2,52:90)]  # this number selects only broad peaks
INDEL_broad <- INDEL_broad[,c(1,2,52:90)]

# make variants unique
SNV_broad = SNV_broad[!duplicated(SNV_broad$variant),]
INDEL_broad = INDEL_broad[!duplicated(INDEL_broad$variant),]

##### 2. mutate a column for each variant's total peak
SNV_broad <- SNV_broad %>% mutate(total = rowSums(SNV_broad[,c(3:41)]))
INDEL_broad <- INDEL_broad %>% mutate(total = rowSums(INDEL_broad[,c(3:41)]))

# making plots for this question(to look for each variant's total peak)
plot_variant_SNV <- SNV_broad %>% 
  filter(total > 10) %>%
  ggplot(aes(x = reorder(variant,-total), total)) + 
  geom_point() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  #geom_label(aes(label = total)) +
  ggtitle("SNVs with more than \n 10 DHS broad peaks", subtitle = "(39 Epigenome RoadMap tissues in total)") + xlab("Single Nucleotide Variant (SNV)") + ylab("Number of tissues \n with DHS peaks")

plot_variant_SNV

plot_variant_SNV_l <- SNV_broad %>% 
  filter(total <= 10 & total >= 5) %>%
  ggplot(aes(x = reorder(variant,-total), total)) + 
  geom_point() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.subtitle = element_text(size = 15, hjust = 0.5)) +
  #geom_label(aes(label = total)) +
  ggtitle("SNVs with 5~10 DHS peaks", subtitle = "(39 Epigenome RoadMap tissues in total)") + xlab("Single Nucleotide Variant (SNV)") + ylab("Number of tissues \n with DHS peaks")

plot_variant_SNV_l

plot_variant_SNV_ll <- SNV_broad %>% 
  filter(total < 5 & total > 0) %>%
  ggplot(aes(x = reorder(variant,-total), total)) + 
  geom_point() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        plot.subtitle = element_text(size = 15, hjust = 0.5)) +
  #geom_label(aes(label = total)) +
  ggtitle("SNVs with 1~4 DHS peaks", subtitle = "(39 Epigenome RoadMap tissues in total)") + xlab("Single Nucleotide Variant (SNV)") + ylab("Number of tissues \n with DHS peaks")

plot_variant_SNV_ll

plot_variant_INDEL <- INDEL_broad %>% 
  filter(total > 0) %>%
  ggplot(aes(x = reorder(variant,-total), total)) + 
  geom_point() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  #geom_label(aes(label = total)) +
  ggtitle("Indels with more than \n 0 DHS broad peak", subtitle = "(39 Epigenome RoadMap tissues in total)") + xlab("Indel") + ylab("Number of tissues \n with DHS peaks")

plot_variant_INDEL

ggsave("~/Hail/rett/plots/plot_variant_SNV.png", plot_variant_SNV)
ggsave("~/Hail/rett/plots/plot_variant_SNV_l.png", plot_variant_SNV_l)
ggsave("~/Hail/rett/plots/plot_variant_SNV_ll.png", plot_variant_SNV_ll)
ggsave("~/Hail/rett/plots/plot_variant_INDEL.png", plot_variant_INDEL)

#### 3. make a new data frame to analyze which tissue has the most DNVs with peaks
peakSums_byTissue_SNV_broad = data.frame(tissue = names(SNV_broad[,c(3:41)]), peakSums=colSums(SNV_broad[,c(3:41)]))
rownames(peakSums_byTissue_SNV_broad) <- NULL

peakSums_byTissue_INDEL_broad = data.frame(tissue = names(INDEL_broad[,c(3:41)]), peakSums=colSums(INDEL_broad[,c(3:41)]))
rownames(peakSums_byTissue_INDEL_broad) <- NULL

# rename the columns for legibility

peakSums_byTissue_SNV_broad$tissue <- gsub('isRE_DHS_','', peakSums_byTissue_SNV_broad$tissue)
peakSums_byTissue_SNV_broad$tissue <- gsub('_',' ', peakSums_byTissue_SNV_broad$tissue)

peakSums_byTissue_INDEL_broad$tissue <- gsub('isRE_DHS_','', peakSums_byTissue_INDEL_broad$tissue)
peakSums_byTissue_INDEL_broad$tissue <- gsub('_',' ', peakSums_byTissue_INDEL_broad$tissue)

peakSums_byTissue_SNV_broad$tissue[peakSums_byTissue_SNV_broad$tissue == 'iPS DF 6 9 Cells'] <- 'iPS DF 6.9 Cells'
peakSums_byTissue_SNV_broad$tissue[peakSums_byTissue_SNV_broad$tissue == 'iPS DF 19 11 Cells'] <- 'iPS DF 19.11 Cells'
peakSums_byTissue_SNV_broad$tissue[peakSums_byTissue_SNV_broad$tissue == 'Primary hematopoietic stem cells G CSF mobilized Male'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Male'
peakSums_byTissue_SNV_broad$tissue[peakSums_byTissue_SNV_broad$tissue == 'Primary hematopoietic stem cells G CSF mobilized Female'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Female'
peakSums_byTissue_SNV_broad$tissue[peakSums_byTissue_SNV_broad$tissue == 'Breast variant Human Mammary Epithelial Cells'] <- 'Breast variant Human Mammary Epithelial Cells (vHMEC)'

peakSums_byTissue_INDEL_broad$tissue[peakSums_byTissue_INDEL_broad$tissue == 'iPS DF 6 9 Cells'] <- 'iPS DF 6.9 Cells'
peakSums_byTissue_INDEL_broad$tissue[peakSums_byTissue_INDEL_broad$tissue == 'iPS DF 19 11 Cells'] <- 'iPS DF 19.11 Cells'
peakSums_byTissue_INDEL_broad$tissue[peakSums_byTissue_INDEL_broad$tissue == 'Primary hematopoietic stem cells G CSF mobilized Male'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Male'
peakSums_byTissue_INDEL_broad$tissue[peakSums_byTissue_INDEL_broad$tissue == 'Primary hematopoietic stem cells G CSF mobilized Female'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Female'
peakSums_byTissue_INDEL_broad$tissue[peakSums_byTissue_INDEL_broad$tissue == 'Breast variant Human Mammary Epithelial Cells'] <- 'Breast variant Human Mammary Epithelial Cells (vHMEC)'

# making tissue group 
# use comment.char = "" to read '#' as a str (R usually ignores infos after #)
# stringsAsFactors = F to prevent changing characters to factors
Roadmap = read.table('~/Downloads/Roadmap.csv', sep = ',', header = T, comment.char = "", stringsAsFactors = F)
names(Roadmap)[names(Roadmap) == "Standardized.Epigenome.name"] <- c("tissue")

peakSums_byTissue_SNV_broad = merge(peakSums_byTissue_SNV_broad, Roadmap, by = 'tissue', all.x = T)
peakSums_byTissue_INDEL_broad = merge(peakSums_byTissue_INDEL_broad, Roadmap, by = 'tissue', all.x = T)

# tissue names do not match where 'peripheral' is included
peakSums_byTissue_SNV_broad$GROUP[peakSums_byTissue_SNV_broad$tissue == 'Primary monocytes from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_SNV_broad$COLOR[peakSums_byTissue_SNV_broad$tissue == 'Primary monocytes from peripheral blood'] = '#678C69'
peakSums_byTissue_SNV_broad$GROUP[peakSums_byTissue_SNV_broad$tissue == 'Primary Natural Killer cells from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_SNV_broad$COLOR[peakSums_byTissue_SNV_broad$tissue == 'Primary Natural Killer cells from peripheral blood'] = '#678C69'
peakSums_byTissue_SNV_broad$GROUP[peakSums_byTissue_SNV_broad$tissue == 'Primary T cells from peripheral blood'] = 'Blood & T-cell'
peakSums_byTissue_SNV_broad$COLOR[peakSums_byTissue_SNV_broad$tissue == 'Primary T cells from peripheral blood'] = '#55A354'

peakSums_byTissue_INDEL_broad$GROUP[peakSums_byTissue_INDEL_broad$tissue == 'Primary monocytes from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_INDEL_broad$COLOR[peakSums_byTissue_INDEL_broad$tissue == 'Primary monocytes from peripheral blood'] = '#678C69'
peakSums_byTissue_INDEL_broad$GROUP[peakSums_byTissue_INDEL_broad$tissue == 'Primary Natural Killer cells from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_INDEL_broad$COLOR[peakSums_byTissue_INDEL_broad$tissue == 'Primary Natural Killer cells from peripheral blood'] = '#678C69'
peakSums_byTissue_INDEL_broad$GROUP[peakSums_byTissue_INDEL_broad$tissue == 'Primary T cells from peripheral blood'] = 'Blood & T-cell'
peakSums_byTissue_INDEL_broad$COLOR[peakSums_byTissue_INDEL_broad$tissue == 'Primary T cells from peripheral blood'] = '#55A354'

# too many groups -> less readable -> unify groups
peakSums_byTissue_SNV_broad$GROUP[peakSums_byTissue_SNV_broad$GROUP %in% c('IMR90', 'Thymus', 'Heart')] = 'Other'
peakSums_byTissue_SNV_broad$GROUP[peakSums_byTissue_SNV_broad$GROUP %in% c('ES-deriv', 'ESC')] = 'ESC & ESC-derived'
peakSums_byTissue_SNV_broad$GROUP[peakSums_byTissue_SNV_broad$GROUP %in% c('HSC & B-cell', 'Blood & T-cell')] = 'Blood'

peakSums_byTissue_INDEL_broad$GROUP[peakSums_byTissue_INDEL_broad$GROUP %in% c('IMR90', 'Thymus', 'Heart')] = 'Other'
peakSums_byTissue_INDEL_broad$GROUP[peakSums_byTissue_INDEL_broad$GROUP %in% c('ES-deriv', 'ESC')] = 'ESC & ESC-derived'
peakSums_byTissue_INDEL_broad$GROUP[peakSums_byTissue_INDEL_broad$GROUP %in% c('HSC & B-cell', 'Blood & T-cell')] = 'Blood'

	
# making plots for this question
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_tissue_SNV <- peakSums_byTissue_SNV_broad %>%
  ggplot(aes(x = reorder(tissue,peakSums), peakSums, color = GROUP)) + 
  geom_point() +
  theme(plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position = 'none') +
  scale_colour_manual(values=cbbPalette) +
  coord_flip() +
  #geom_label(aes(label = peakSums)) +
  ggtitle("Number of SNVs with DHS broad peak \n in Epigenome RoadMap Tissues") + xlab("Epigenome RoadMap Tissues") + ylab("Number of SNVs with DHS peak")

plot_tissue_SNV

plot_tissue_INDEL <- peakSums_byTissue_INDEL_broad %>%
  ggplot(aes(x = reorder(tissue, peakSums), peakSums, color = GROUP)) + 
  geom_point() +
  theme(plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position = 'none') +
  scale_colour_manual(values=cbbPalette) +
  coord_flip() +
  #geom_label(aes(label = peakSums)) +
  ggtitle("Number of Indels with DHS broad peak \n in Epigenome RoadMap Tissues") + xlab("Epigenome RoadMap Tissues") + ylab("Number of Indels with DHS peak")

plot_tissue_INDEL

ggsave("~/Hail/rett/plots/plot_tissue_SNV.png", plot_tissue_SNV)
ggsave("~/Hail/rett/plots/plot_tissue_INDEL.png", plot_tissue_INDEL)


#### 4. How many DNVs and peaks are in each sample?

sample = c('TWGS1_1', 'TWGS2_1', 'TWGS3_1', 'TWGS4_1', 'TWGS5_1', 'TWGS6_1', 'TWGS7_1', 'TWGS8_1', 'TWGS9_1')
DNV_SNV_all = c() # the number of SNVs per sample
DNV_SNV_peak = c() # the number of SNVs with peak per sample
DNV_indel_all = c()
DNV_indel_peak = c()

bySample_SNV <- function(samp, i){
  DNV_all = SNV_broad %>% filter(s == samp) %>% count() %>% as.integer()
  DNV_peak = SNV_broad %>% filter(s == samp & total != 0) %>% count() %>% as.integer()
  DNV = list(DNV_all, DNV_peak) # return list to return multiple objects
  return(DNV)
}

for(i in 1:9){
  DNV_SNV_all[i] = bySample_SNV(sample[i])[1]
  DNV_SNV_peak[i] = bySample_SNV(sample[i])[2]
}

bySample_indel <- function(samp, i){
  DNV_all = INDEL_broad %>% filter(s == samp) %>% count() %>% as.integer()
  DNV_peak = INDEL_broad %>% filter(s == samp & total != 0) %>% count() %>% as.integer()
  DNV = list(DNV_all, DNV_peak)
  return(DNV)
}

for(i in 1:9){
  DNV_indel_all[i] = bySample_indel(sample[i])[1]
  DNV_indel_peak[i] = bySample_indel(sample[i])[2]
}

# making list to vectors to make data frame
DNV_SNV_all = unlist(DNV_SNV_all) 
DNV_SNV_peak = unlist(DNV_SNV_peak) 
DNV_indel_all = unlist(DNV_indel_all)
DNV_indel_peak = unlist(DNV_indel_peak)

# making data frame
sums_bySample_SNV_broad <- data.frame(sample, DNV_SNV_all, DNV_SNV_peak, stringsAsFactors = F)
sums_bySample_indel_broad <- data.frame(sample, DNV_indel_all, DNV_indel_peak, stringsAsFactors = F)

# making plots for this question
plot_sample_SNV <- sums_bySample_SNV_broad %>%
  ggplot() +
  geom_point(aes(sample, DNV_SNV_all), col = "red") +
  geom_point(aes(sample, DNV_SNV_peak), col = "blue") +
  geom_text(aes('TWGS8_1', 40, label = "SNV"), col = "red") +
  geom_text(aes('TWGS8_1', 20, label ="SNV with DHS broad peak"), col = "blue") +
  ggtitle("Number of SNVs and SNVs with DHS broad peak \n in each sample") + xlab("Sample") + ylab("Count") +
  theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

plot_sample_SNV

plot_sample_INDEL <- sums_bySample_indel_broad %>%
  ggplot() +
  geom_point(aes(sample, DNV_indel_all), col = "red") +
  geom_point(aes(sample, DNV_indel_peak), col = "blue") +
  geom_text(aes('TWGS8_1', 6.5, label = "Indel"), col = "red") +
  geom_text(aes('TWGS8_1', 6, label ="Indel with DHS broad peak"), col = "blue") +
  ggtitle("Number of Indels and Indels with DHS broad peak \n in each sample") + xlab("Sample") + ylab("Count") +
  theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

plot_sample_INDEL

ggsave("~/Hail/rett/plots/plot_sample_SNV.png", plot_sample_SNV)
ggsave("~/Hail/rett/plots/plot_sample_INDEL.png", plot_sample_INDEL)



####################################################################################################
################################### only with narrow peaks #########################################
####################################################################################################

SNV_narrow = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_SNV_broad_narrow_compared.tsv', sep = '\t', header = T)
INDEL_narrow = read.table('~/Hail/rett/tables/rett_HQvariants_consequence_INDEL_broad_narrow_compared.tsv', sep = '\t', header = T)

##### 1. making unique variants

# remove redundant columns
SNV_narrow <- SNV_narrow[,c(1,2,91:129)]  # this number selects only broad peaks
INDEL_narrow <- INDEL_narrow[,c(1,2,91:129)]

# make variants unique
SNV_narrow = SNV_narrow[!duplicated(SNV_narrow$variant),]
INDEL_narrow = INDEL_narrow[!duplicated(INDEL_narrow$variant),]

##### 2. mutate a column for each variant's total peak
SNV_narrow <- SNV_narrow %>% mutate(total = rowSums(SNV_narrow[,c(3:41)]))
INDEL_narrow <- INDEL_narrow %>% mutate(total = rowSums(INDEL_narrow[,c(3:41)]))

# making plots for this question(to look for each variant's total peak)
plot_variant_SNV_narrow <- SNV_narrow %>% 
  filter(total > 10) %>%
  ggplot(aes(x = reorder(variant,-total), total)) + 
  geom_point() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  scale_y_continuous(breaks = seq(10,20,2)) +
  #geom_label(aes(label = total)) +
  ggtitle("SNVs with more than \n 10 DHS narrow peaks", subtitle = "(39 Epigenome RoadMap tissues in total)") + xlab("Single Nucleotide Variant (SNV)") + ylab("")

plot_variant_SNV_narrow

plot_variant_INDEL_narrow <- INDEL_narrow %>% 
  filter(total > 0) %>%
  ggplot(aes(x = reorder(variant,-total), total)) + 
  geom_point() +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  scale_y_continuous(breaks = seq(0,15,2)) +
  #geom_label(aes(label = total)) +
  ggtitle("Indels with more than \n 0 DHS narrow peak", subtitle = "(39 Epigenome RoadMap tissues in total)") + xlab("Indel") + ylab("")

plot_variant_INDEL_narrow

ggsave("~/Hail/rett/plots/plot_variant_SNV_narrow.png", plot_variant_SNV_narrow)
ggsave("~/Hail/rett/plots/plot_variant_INDEL_narrow.png", plot_variant_INDEL_narrow)

#### 3. make a new data frame to analyze which tissue has the most DNVs with peaks
peakSums_byTissue_SNV_narrow = data.frame(tissue = names(SNV_narrow[,c(3:41)]), peakSums=colSums(SNV_narrow[,c(3:41)]))
rownames(peakSums_byTissue_SNV_narrow) <- NULL

peakSums_byTissue_INDEL_narrow = data.frame(tissue = names(INDEL_narrow[,c(3:41)]), peakSums=colSums(INDEL_narrow[,c(3:41)]))
rownames(peakSums_byTissue_INDEL_narrow) <- NULL

# rename the columns for legibility

peakSums_byTissue_SNV_narrow$tissue <- gsub('isRE_DHS_','', peakSums_byTissue_SNV_narrow$tissue)
peakSums_byTissue_SNV_narrow$tissue <- gsub('_',' ', peakSums_byTissue_SNV_narrow$tissue)
peakSums_byTissue_SNV_narrow$tissue <- gsub('narrow ','', peakSums_byTissue_SNV_narrow$tissue)

peakSums_byTissue_INDEL_narrow$tissue <- gsub('isRE_DHS_','', peakSums_byTissue_INDEL_narrow$tissue)
peakSums_byTissue_INDEL_narrow$tissue <- gsub('_',' ', peakSums_byTissue_INDEL_narrow$tissue)
peakSums_byTissue_INDEL_narrow$tissue <- gsub('narrow ','', peakSums_byTissue_INDEL_narrow$tissue)

peakSums_byTissue_SNV_narrow$tissue[peakSums_byTissue_SNV_narrow$tissue == 'iPS DF 6 9 Cells'] <- 'iPS DF 6.9 Cells'
peakSums_byTissue_SNV_narrow$tissue[peakSums_byTissue_SNV_narrow$tissue == 'iPS DF 19 11 Cells'] <- 'iPS DF 19.11 Cells'
peakSums_byTissue_SNV_narrow$tissue[peakSums_byTissue_SNV_narrow$tissue == 'Primary hematopoietic stem cells G CSF mobilized Male'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Male'
peakSums_byTissue_SNV_narrow$tissue[peakSums_byTissue_SNV_narrow$tissue == 'Primary hematopoietic stem cells G CSF mobilized Female'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Female'
peakSums_byTissue_SNV_narrow$tissue[peakSums_byTissue_SNV_narrow$tissue == 'Breast variant Human Mammary Epithelial Cells'] <- 'Breast variant Human Mammary Epithelial Cells (vHMEC)'

peakSums_byTissue_INDEL_narrow$tissue[peakSums_byTissue_INDEL_narrow$tissue == 'iPS DF 6 9 Cells'] <- 'iPS DF 6.9 Cells'
peakSums_byTissue_INDEL_narrow$tissue[peakSums_byTissue_INDEL_narrow$tissue == 'iPS DF 19 11 Cells'] <- 'iPS DF 19.11 Cells'
peakSums_byTissue_INDEL_narrow$tissue[peakSums_byTissue_INDEL_narrow$tissue == 'Primary hematopoietic stem cells G CSF mobilized Male'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Male'
peakSums_byTissue_INDEL_narrow$tissue[peakSums_byTissue_INDEL_narrow$tissue == 'Primary hematopoietic stem cells G CSF mobilized Female'] <- 'Primary hematopoietic stem cells G-CSF-mobilized Female'
peakSums_byTissue_INDEL_narrow$tissue[peakSums_byTissue_INDEL_narrow$tissue == 'Breast variant Human Mammary Epithelial Cells'] <- 'Breast variant Human Mammary Epithelial Cells (vHMEC)'

# making tissue group 
# use comment.char = "" to read '#' as a str (R usually ignores infos after #)
# stringsAsFactors = F to prevent changing characters to factors
Roadmap = read.table('~/Downloads/Roadmap.csv', sep = ',', header = T, comment.char = "", stringsAsFactors = F)
names(Roadmap)[names(Roadmap) == "Standardized.Epigenome.name"] <- c("tissue")

peakSums_byTissue_SNV_narrow = merge(peakSums_byTissue_SNV_narrow, Roadmap, by = 'tissue', all.x = T)
peakSums_byTissue_INDEL_narrow = merge(peakSums_byTissue_INDEL_narrow, Roadmap, by = 'tissue', all.x = T)

# tissue names do not match where 'peripheral' is included
peakSums_byTissue_SNV_narrow$GROUP[peakSums_byTissue_SNV_narrow$tissue == 'Primary monocytes from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_SNV_narrow$COLOR[peakSums_byTissue_SNV_narrow$tissue == 'Primary monocytes from peripheral blood'] = '#678C69'
peakSums_byTissue_SNV_narrow$GROUP[peakSums_byTissue_SNV_narrow$tissue == 'Primary Natural Killer cells from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_SNV_narrow$COLOR[peakSums_byTissue_SNV_narrow$tissue == 'Primary Natural Killer cells from peripheral blood'] = '#678C69'
peakSums_byTissue_SNV_narrow$GROUP[peakSums_byTissue_SNV_narrow$tissue == 'Primary T cells from peripheral blood'] = 'Blood & T-cell'
peakSums_byTissue_SNV_narrow$COLOR[peakSums_byTissue_SNV_narrow$tissue == 'Primary T cells from peripheral blood'] = '#55A354'

peakSums_byTissue_INDEL_narrow$GROUP[peakSums_byTissue_INDEL_narrow$tissue == 'Primary monocytes from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_INDEL_narrow$COLOR[peakSums_byTissue_INDEL_narrow$tissue == 'Primary monocytes from peripheral blood'] = '#678C69'
peakSums_byTissue_INDEL_narrow$GROUP[peakSums_byTissue_INDEL_narrow$tissue == 'Primary Natural Killer cells from peripheral blood'] = 'HSC & B-cell'
peakSums_byTissue_INDEL_narrow$COLOR[peakSums_byTissue_INDEL_narrow$tissue == 'Primary Natural Killer cells from peripheral blood'] = '#678C69'
peakSums_byTissue_INDEL_narrow$GROUP[peakSums_byTissue_INDEL_narrow$tissue == 'Primary T cells from peripheral blood'] = 'Blood & T-cell'
peakSums_byTissue_INDEL_narrow$COLOR[peakSums_byTissue_INDEL_narrow$tissue == 'Primary T cells from peripheral blood'] = '#55A354'

# too many groups -> less readable -> unify groups
peakSums_byTissue_SNV_narrow$GROUP[peakSums_byTissue_SNV_narrow$GROUP %in% c('IMR90', 'Thymus', 'Heart')] = 'Other'
peakSums_byTissue_SNV_narrow$GROUP[peakSums_byTissue_SNV_narrow$GROUP %in% c('ES-deriv', 'ESC')] = 'ESC & ESC-derived'
peakSums_byTissue_SNV_narrow$GROUP[peakSums_byTissue_SNV_narrow$GROUP %in% c('HSC & B-cell', 'Blood & T-cell')] = 'Blood'

peakSums_byTissue_INDEL_narrow$GROUP[peakSums_byTissue_INDEL_narrow$GROUP %in% c('IMR90', 'Thymus', 'Heart')] = 'Other'
peakSums_byTissue_INDEL_narrow$GROUP[peakSums_byTissue_INDEL_narrow$GROUP %in% c('ES-deriv', 'ESC')] = 'ESC & ESC-derived'
peakSums_byTissue_INDEL_narrow$GROUP[peakSums_byTissue_INDEL_narrow$GROUP %in% c('HSC & B-cell', 'Blood & T-cell')] = 'Blood'


# making plots for this question
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_tissue_SNV_narrow <- peakSums_byTissue_SNV_narrow %>%
  ggplot(aes(x = reorder(tissue, peakSums), peakSums, color = GROUP)) + 
  geom_point() +
  theme(plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  scale_colour_manual(values=cbbPalette) +
  coord_flip() +
  #geom_label(aes(label = peakSums)) +
  ggtitle("Number of SNVs with DHS narrow peak \n in Epigenome RoadMap Tissues") + xlab("") + ylab("Number of SNVs with DHS peak")

plot_tissue_SNV_narrow

plot_tissue_INDEL_narrow <- peakSums_byTissue_INDEL_narrow %>%
  ggplot(aes(x = reorder(tissue, peakSums), peakSums, color = GROUP)) + 
  geom_point() +
  theme(plot.title = element_text(size=20, hjust = 0.5),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) +
  scale_colour_manual(values=cbbPalette) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0,2,1)) +
  #geom_label(aes(label = peakSums)) +
  ggtitle("Number of Indels with DHS narrow peak \n in Epigenome RoadMap Tissues") + xlab("") + ylab("Number of Indels with DHS peak")

plot_tissue_INDEL_narrow

ggsave("~/Hail/rett/plots/plot_tissue_SNV_narrow.png", plot_tissue_SNV_narrow)
ggsave("~/Hail/rett/plots/plot_tissue_INDEL_narrow.png", plot_tissue_INDEL_narrow)


#### 4. How many DNVs and peaks are in each sample?

sample = c('TWGS1_1', 'TWGS2_1', 'TWGS3_1', 'TWGS4_1', 'TWGS5_1', 'TWGS6_1', 'TWGS7_1', 'TWGS8_1', 'TWGS9_1')
DNV_SNV_all_narrow = c() # the number of SNVs per sample
DNV_SNV_peak_narrow = c() # the number of SNVs with peak per sample
DNV_indel_all_narrow = c()
DNV_indel_peak_narrow = c()

bySample_SNV_narrow <- function(samp, i){
  DNV_all_narrow = SNV_narrow %>% filter(s == samp) %>% count() %>% as.integer()
  DNV_peak_narrow = SNV_narrow %>% filter(s == samp & total != 0) %>% count() %>% as.integer()
  DNV = list(DNV_all_narrow, DNV_peak_narrow) # return list to return multiple objects
  return(DNV)
}

for(i in 1:9){
  DNV_SNV_all_narrow[i] = bySample_SNV_narrow(sample[i])[1]
  DNV_SNV_peak_narrow[i] = bySample_SNV_narrow(sample[i])[2]
}

bySample_indel_narrow <- function(samp, i){
  DNV_all_narrow = INDEL_narrow %>% filter(s == samp) %>% count() %>% as.integer()
  DNV_peak_narrow = INDEL_narrow %>% filter(s == samp & total != 0) %>% count() %>% as.integer()
  DNV = list(DNV_all_narrow, DNV_peak_narrow)
  return(DNV)
}

for(i in 1:9){
  DNV_indel_all_narrow[i] = bySample_indel_narrow(sample[i])[1]
  DNV_indel_peak_narrow[i] = bySample_indel_narrow(sample[i])[2]
}

# making list to vectors to make data frame
DNV_SNV_all_narrow = unlist(DNV_SNV_all_narrow) 
DNV_SNV_peak_narrow = unlist(DNV_SNV_peak_narrow) 
DNV_indel_all_narrow = unlist(DNV_indel_all_narrow)
DNV_indel_peak_narrow = unlist(DNV_indel_peak_narrow)

# making data frame
sums_bySample_SNV_narrow <- data.frame(sample, DNV_SNV_all_narrow, DNV_SNV_peak_narrow, stringsAsFactors = F)
sums_bySample_indel_narrow <- data.frame(sample, DNV_indel_all_narrow, DNV_indel_peak_narrow, stringsAsFactors = F)

# making plots for this question
plot_sample_SNV_narrow <- sums_bySample_SNV_narrow %>%
  ggplot() +
  geom_point(aes(sample, DNV_SNV_all_narrow), col = "red") +
  geom_point(aes(sample, DNV_SNV_peak_narrow), col = "blue") +
  geom_text(aes('TWGS8_1', 40, label = "SNV"), col = "red") +
  geom_text(aes('TWGS8_1', 20, label ="SNV with DHS narrow peak"), col = "blue") +
  ggtitle("Number of SNVs and SNVs with DHS narrow peak \n in each sample") + xlab("Sample") + ylab("Count") +
  theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

plot_sample_SNV_narrow

plot_sample_INDEL_narrow <- sums_bySample_indel_narrow %>%
  ggplot() +
  geom_point(aes(sample, DNV_indel_all_narrow), col = "red") +
  geom_point(aes(sample, DNV_indel_peak_narrow), col = "blue") +
  geom_text(aes('TWGS8_1', 6.5, label = "Indel"), col = "red") +
  geom_text(aes('TWGS8_1', 6, label ="Indel with DHS narrow peak"), col = "blue") +
  ggtitle("Number of Indels and Indels with DHS narrow peak \n in each sample") + xlab("Sample") + ylab("Count") +
  theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

plot_sample_INDEL_narrow

ggsave("~/Hail/rett/plots/plot_sample_SNV_narrow.png", plot_sample_SNV_narrow)
ggsave("~/Hail/rett/plots/plot_sample_INDEL_narrow.png", plot_sample_INDEL_narrow)




####################################################################################################
################################### merge two plots to compare #####################################
####################################################################################################

plot_variant_SNV_grid = plot_grid(plot_variant_SNV, plot_variant_SNV_narrow)
plot_variant_INDEL_grid = plot_grid(plot_variant_INDEL, plot_variant_INDEL_narrow)
plot_grid(plot_tissue_SNV,plot_tissue_SNV_narrow) # failed
plot_grid(plot_tissue_INDEL,plot_tissue_INDEL_narrow) # failed
plot_grid(plot_sample_SNV, plot_sample_SNV_narrow, ncol = 1)

plot_sample_SNV_grid <- sums_bySample_SNV_broad %>% 
  mutate(DNV_SNV_peak_narrow = sums_bySample_SNV_narrow$DNV_SNV_peak_narrow) %>%
  ggplot() +
  geom_point(aes(sample, DNV_SNV_all), col = "red") +
  geom_point(aes(sample, DNV_SNV_peak), col = "blue") +
  geom_point(aes(sample, DNV_SNV_peak_narrow), col = "green") +
  geom_text(aes('TWGS8_1', 40, label = "SNV"), col = "red") +
  geom_text(aes('TWGS8_1', 30, label ="SNV with DHS broad peak"), col = "blue") +
  geom_text(aes('TWGS8_1', 20, label ="SNV with DHS narrow peak"), col = "green") +
  ggtitle("Number of SNVs and SNVs with DHS peaks \n in each sample") + xlab("Sample") + ylab("Count") +
  theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))

plot_sample_INDEL_grid <- sums_bySample_indel_broad %>% 
  mutate(DNV_indel_peak_narrow = sums_bySample_indel_narrow$DNV_indel_peak_narrow) %>%
  ggplot() +
  geom_point(aes(sample, DNV_indel_all), col = "red") +
  geom_point(aes(sample, DNV_indel_peak), col = "blue") +
  geom_point(aes(sample, DNV_indel_peak_narrow), col = "green") +
  geom_text(aes('TWGS8_1', 8, label = "Indel"), col = "red") +
  geom_text(aes('TWGS8_1', 7.5, label ="Indel with DHS broad peak"), col = "blue") +
  geom_text(aes('TWGS8_1', 7, label ="Indel with DHS narrow peak"), col = "green") +
  ggtitle("Number of Indels and Indels with DHS peaks \n in each sample") + xlab("Sample") + ylab("Count") +
  theme(plot.title = element_text(size=25, face="bold", hjust = 0.5),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
  
ggsave("~/Hail/rett/plots/plot_variant_SNV_grid.png", plot_variant_SNV_grid)
ggsave("~/Hail/rett/plots/plot_variant_INDEL_grid.png", plot_variant_INDEL_grid)

ggsave("~/Hail/rett/plots/plot_tissue_SNV_for_grid.png", plot_tissue_SNV)
ggsave("~/Hail/rett/plots/plot_tissue_SNV_narrow_for_grid.png", plot_tissue_SNV_narrow)
ggsave("~/Hail/rett/plots/plot_tissue_INDEL_narrow_for_grid.png", plot_tissue_INDEL_narrow)
ggsave("~/Hail/rett/plots/plot_tissue_INDEL_for_grid.png", plot_tissue_INDEL)

ggsave("~/Hail/rett/plots/plot_sample_SNV_grid.png", plot_sample_SNV_grid)
ggsave("~/Hail/rett/plots/plot_sample_INDEL_grid.png", plot_sample_INDEL_grid)
