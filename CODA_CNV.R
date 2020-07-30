library(bedr)
library(dplyr)
library(ggplot2)
library(cowplot)

fam_id_0_PM <- read.table("/home/titan/mosdepth/JEW-PM/JEW-PM-with-bed.regions.bed", 
                          header = F, sep = "\t", stringsAsFactors = T, quote = "")
colnames(fam_id_0_PM) <- c("chr", "start", "end", "average_depth")
fam_id_0_M <- read.table("/home/titan/mosdepth/JEW-M/JEW-M-with-bed.regions.bed", 
                         header = F, sep = "\t", stringsAsFactors = F, quote = "")
colnames(fam_id_0_M) <- c("chr", "start", "end", "average_depth")
fam_id_0_P <- read.table("/home/titan/mosdepth/JEW-P/JEW-P-with-bed.regions.bed", 
                         header = F, sep = "\t", stringsAsFactors = F, quote = "")
colnames(fam_id_0_P) <- c("chr", "start", "end", "average_depth")


p1<-fam_id_0_PM %>%
  ggplot(aes(x=start + (end-start)/2, y=average_depth)) +
    geom_point() +
    geom_segment(aes(x=start, y=average_depth , xend=end, yend=average_depth)) +
    facet_wrap(~chr, ncol=6) +
    xlab("Position") + ylab("Average Depth") + ggtitle("CODA Family0 - PM") +
    theme_bw() +
    theme(axis.text.x = element_blank())
    
p1
ggsave("/home/titan/mosdepth/CODA_mosdepth_plot_JEW_PM_all_chromosomes.png", p1)

p2<-fam_id_0_P %>%
  ggplot(aes(x=start + (end-start)/2, y=average_depth)) +
  geom_point() +
  geom_segment(aes(x=start, y=average_depth , xend=end, yend=average_depth)) +
  facet_wrap(~chr, ncol=6) +
  xlab("Position") + ylab("Average Depth") + ggtitle("CODA Family0 - P") +
  theme_bw() +
  theme(axis.text.x = element_blank())
p2
ggsave("/home/titan/mosdepth/CODA_mosdepth_plot_JEW_P_all_chromosomes.png", p2)

p3<-fam_id_0_M %>%
  ggplot(aes(x=start + (end-start)/2, y=average_depth)) +
  geom_point() +
  geom_segment(aes(x=start, y=average_depth , xend=end, yend=average_depth)) +
  facet_wrap(~chr, ncol=6) +
  xlab("Position") + ylab("Average Depth") + ggtitle("CODA Family0 - M") +
  theme_bw() +
  theme(axis.text.x = element_blank())
p3
ggsave("/home/titan/mosdepth/CODA_mosdepth_plot_JEW_M_all_chromosomes.png", p3)


#### from this kind of plot, we can get information from all chromosomes, 
#### but it is not easy to compare each individuals
# Q. Why does coverage exist in chrY in mother?

## Let's look at chromosomes one by one.
pp1<-fam_id_0_PM %>%
  filter(fam_id_0_PM$chr == "chr1") %>%
  head(n=10000) %>%
  ggplot(aes(x=start + (end-start)/2, y=average_depth)) +
  geom_point() +
  geom_segment(aes(x=start, y=average_depth , xend=end, yend=average_depth)) +
  xlab("Position") + ylab("Average Depth") +
  ggtitle("CODA Family0 PM Chr1") +
  theme_bw()

pp2<-fam_id_0_P %>%
  filter(fam_id_0_P$chr == "chr1") %>%
  head(n=10000) %>%
  ggplot(aes(x=start + (end-start)/2, y=average_depth)) +
  geom_point() +
  geom_segment(aes(x=start, y=average_depth , xend=end, yend=average_depth)) +
  xlab("Position") + ylab("Average Depth") +
  ggtitle("CODA Family0 P Chr1") +
  theme_bw()

pp3<-fam_id_0_M %>%
  filter(fam_id_0_M$chr == "chr1") %>%
  head(n=10000) %>%
  ggplot(aes(x=start + (end-start)/2, y=average_depth)) +
  geom_point() +
  geom_segment(aes(x=start, y=average_depth , xend=end, yend=average_depth)) +
  xlab("Position") + ylab("Average Depth") +
  ggtitle("CODA Family0 M Chr1") +
  theme_bw()

pp4<-plot_grid(pp1, pp2, pp3)
ggsave("/home/titan/mosdepth/CODA_mosdepth_plot_Family0_Chr1.png", pp4)

#### Can't make decision through just watch and compare... there must be exact criterion..

########################################################################################
################ New Phase: Make Z score column to each sample #########################
########################################################################################

############################## 1. Make function to read all the sample data

# Match sample name and family id
sample_list = data.frame(list.dirs("/home/titan/mosdepth", full.names = F), stringsAsFactors = F)
sample_list = data.frame(sample_list[-c(1,2),])
colnames(sample_list) = "name"
family <- c("fam000_M", "fam000_P", "fam000_PM", 
            "fam001_M", "fam001_PM", 
            "fam002_M", "fam002_P", "fam002_PM", 
            "fam004_M", "fam004_P", "fam004_PM", 
            "fam005_P", "fam005_PM", 
            "fam006_M", "fam006_P", "fam006_PM",
            "fam007_M", "fam007_P", "fam007_PM", 
            "fam008_M", "fam008_PM",
            "fam009_M", "fam009_PM",
            "fam010_M", "fam010_PM",
            "fam011_M", "fam011_PM",
            "fam012_M", "fam012_PM",
            "fam013_M", "fam013_P", "fam013_PM",
            "fam014_M", "fam014_P", "fam014_PM",
            "fam015_PM", 
            "fam016_M", "fam016_PM",
            "fam003_M", "fam003_PM", "fam003_P")
sample_list["family"] <- family

# read.table function
read_sample <- function(name, id){
  id <- read.table(paste0("/home/titan/mosdepth/", name, "/", name, ".regions.bed"), 
                   header = F, sep = "\t", stringsAsFactors = F, quote = "", 
                   col.names = c("chr", "start", "end", "coverage"))
  return(data.frame(id))
}
for (i in 1:nrow(sample_list)){
  assign(family[i], read_sample(sample_list$name[i], sample_list$family[i]))
}

################### 2. Make function to make Z score column to each data
calc_Zscore <- function(sample){         # edit this function to change 'mean' to 'median'
  sample["Zscore"] <- (sample$coverage-mean(sample$coverage))/mean(sample$coverage)
  return(sample)
}
for(i in 1:nrow(sample_list)){
  data.frame(family[i]) = calc_Zscore(data.frame(family[i]))
}
fam013_M = calc_Zscore(fam013_M)
fam013_P = calc_Zscore(fam013_P)
fam013_PM = calc_Zscore(fam013_PM)

################################ read SNU013 family manta vcf ####
fam013_manta <- read.table("/home/titan/manta/SNU013-CB_family/results/variants/diploidSV.vcf", 
                          header = F, sep = "\t", stringsAsFactors = T, quote = "", 
                          col.names = c('CHROM',	'POS', 'ID',	'REF',	'ALT', 'QUAL',	'FILTER',	'INFO',	'FORMAT',	'SNU013-CB-M',	'SNU013-CB-P',	'SNU013-CB-PM'))

chr12	9311717	MantaBND:132:0:1:0:0:0:0	G	G]chr12:9420150]	25	PASS	SVTYPE=BND;MATEID=MantaBND:132:0:1:0:0:0:1;IMPRECISE;CIPOS=-89,89;BND_DEPTH=1409;MATE_BND_DEPTH=733	GT:FT:GQ:PL:PR	0/0:HomRef;MinGQ:5:47,0,999:197,30	0/0:HomRef:347:0,297,999:224,18	0/1:PASS:23:73,0,999:214,34

chr12	9420150	MantaBND:132:0:1:0:0:0:1	G	G]chr12:9311717]	25	PASS	SVTYPE=BND;MATEID=MantaBND:132:0:1:0:0:0:0;IMPRECISE;CIPOS=-178,179;BND_DEPTH=733;MATE_BND_DEPTH=1409	GT:FT:GQ:PL:PR	0/0:HomRef;MinGQ:5:47,0,999:197,30	0/0:HomRef:347:0,297,999:224,18	0/1:PASS:23:73,0,999:214,34

chr12	40505017	MantaDUP:TANDEM:161:0:2:0:0:0	A	<DUP:TANDEM>	51	PASS	END=40510965;SVTYPE=DUP;SVLEN=5948;IMPRECISE;CIPOS=-165,166;CIEND=-116,117	GT:FT:GQ:PL:PR	0/0:HomRef:119:0,69,999:71,8	0/0:HomRef:223:0,173,999:80,4	0/1:PASS:51:101,0,999:67,17

chrX	7843816	MantaBND:12:0:1:0:0:0:0	G	G]chrX:8170138]	52	PASS	SVTYPE=BND;MATEID=MantaBND:12:0:1:0:0:0:1;IMPRECISE;CIPOS=-76,76;BND_DEPTH=531;MATE_BND_DEPTH=220	GT:FT:GQ:PL:PR	0/0:HomRef;MinGQ:9:41,0,734:49,10	0/0:HomRef:58:0,8,329:21,3	0/1:PASS:51:101,0,397:27,10

chrX	8170138	MantaBND:12:0:1:0:0:0:1	C	C]chrX:7843816]	52	PASS	SVTYPE=BND;MATEID=MantaBND:12:0:1:0:0:0:0;IMPRECISE;CIPOS=-125,125;BND_DEPTH=220;MATE_BND_DEPTH=531	GT:FT:GQ:PL:PR	0/0:HomRef;MinGQ:9:41,0,734:49,10	0/0:HomRef:58:0,8,329:21,3	0/1:PASS:51:101,0,397:27,10

temp_fam013_M <- fam013_M[,-4]
temp <- merge(fam013_M, fam013_P, by = c('chr', 'start', 'end'), all = T )
fam013 <- merge(temp, fam013_PM, by = c('chr', 'start', 'end'), all = T )
fam013 <- fam013[,-c(4,6,8)]
names(fam013)[4] = "Zscore_M"
names(fam013)[5] = "Zscore_P"
names(fam013)[6] = "Zscore_PM"

p5 <- fam013 %>%
  filter(fam013$chr == "chr2") %>%
  filter(start > 178420000 & start < 178460000) %>%
  ggplot() +
  geom_point(aes(start,Zscore_M), col = "blue") +
  geom_point(aes(end,Zscore_M), col = "blue") +
  geom_point(aes(start,Zscore_P), col = "red") +
  geom_point(aes(end,Zscore_P), col = "red") +
  geom_point(aes(start,Zscore_PM)) +
  geom_point(aes(end,Zscore_PM)) +
  ggtitle("CODA fam013 Maternal Deletion") + xlab("chr2 exons") + ylab("Z score") +
  geom_text(aes(178460000, 1.5), label = "M", col="red") +
  geom_text(aes(178460000, 1.25), label = "P", col="blue") +
  geom_text(aes(178460000, 1), label = "PM") +
  geom_vline(xintercept = 178432254, col = "yellow") +
  geom_vline(xintercept = 178450965, col = "yellow")
ggsave("~/manta/SNU013-CB_family/CODA_manta_fam013_maternal_deletion.png", p5)

      fam013 %>%
  filter(fam013$chr == "chr13") %>%
  filter(start > 21145000 & start < 21165000) %>%
  ggplot() +
  geom_point(aes(start,Zscore_M), col = "blue") +
  geom_point(aes(end,Zscore_M), col = "blue") +
  geom_point(aes(start,Zscore_P), col = "red") +
  geom_point(aes(end,Zscore_P), col = "red") +
  geom_point(aes(start,Zscore_PM)) +
  geom_point(aes(end,Zscore_PM)) +
  ggtitle("CODA fam013 paternal deletion") + xlab("chr13 start position") + ylab("Z score") +
  geom_text(aes(21165000, 1), label = "M", col="red") +
  geom_text(aes(21165000, 0.75), label = "P", col="blue") +
  geom_text(aes(21165000, 0.5), label = "PM") +
  geom_vline(xintercept = 21155811, col = "yellow") +
  geom_vline(xintercept = 21157921, col = "yellow")
ggsave("~/manta/SNU013-CB_family/CODA_manta_fam013_CNV3.png", p3)

p3 <- fam013 %>%
  filter(fam013$chr == "chrX") %>%
  filter(start > 7843000 & start < 8171000) %>%
  ggplot() +
  geom_point(aes(start,Zscore_M), col = "blue") +
  geom_point(aes(start,Zscore_P), col = "red") +
  geom_point(aes(start,Zscore_PM)) +
  ggtitle("CODA fam013 de novo CNV3") + xlab("chrX start position") + ylab("Z score") +
  geom_text(aes(8200000, 1), label = "M", col="red") +
  geom_text(aes(8200000, 0.75), label = "P", col="blue") +
  geom_text(aes(8200000, 0.5), label = "PM") +
  geom_vline(xintercept = 7843816, col = "yellow") +
  geom_vline(xintercept = 8170138, col = "yellow")
ggsave("~/manta/SNU013-CB_family/CODA_manta_fam013_CNV3.png", p3)

p2 <- fam013 %>%
  filter(fam013$chr == "chr12") %>%
  filter(start > 40505000 & start < 40512000) %>%
  ggplot() +
  geom_point(aes(start,Zscore_M), col = "blue") +
  geom_point(aes(start,Zscore_P), col = "red") +
  geom_point(aes(start,Zscore_PM)) +
  ggtitle("CODA fam013 de novo CNV2") + xlab("chr12 start position") + ylab("Z score") +
  geom_text(aes(40512000, 1.5), label = "M", col="red") +
  geom_text(aes(40512000, 1.25), label = "P", col="blue") +
  geom_text(aes(40512000, 1), label = "PM") +
  geom_vline(xintercept = 40505017, col = "yellow") +
  geom_vline(xintercept = 40510965, col = "yellow")
ggsave("~/manta/SNU013-CB_family/CODA_manta_fam013_CNV2.png", p2)

p1 <- fam013 %>%
  filter(fam013$chr == "chr12") %>%
  filter(start > 9300000 & start < 9500000) %>%
  ggplot() +
  geom_point(aes(start,Zscore_M), col = "blue") +
  geom_point(aes(start,Zscore_P), col = "red") +
  geom_point(aes(start,Zscore_PM)) +
  ggtitle("CODA fam013 de novo CNV1") + xlab("chr12 start position") + ylab("Z score") +
  geom_text(aes(9449000, 8), label = "M", col="red") +
  geom_text(aes(9449000, 7.5), label = "P", col="blue") +
  geom_text(aes(9449000, 7), label = "PM") +
  geom_vline(xintercept = 9311717, col = "yellow") +
  geom_vline(xintercept = 9420150, col = "yellow")
ggsave("~/manta/SNU013-CB_family/CODA_manta_fam013_CNV1.png", p1)
  

bind_cols()
####################### 3. Merge Z score columns
d_samples = list()
for (sample in family){
  print(sample)
  #d_samples[[sample]] = read.table(paste0("/home/titan/mosdepth/", sample_list$name[sample_list$family == sample], 
   #                                       "/", sample_list$name[sample_list$family == sample], ".regions.bed"), 
    #                               header = F, sep = "\t", stringsAsFactors = F, quote = "", 
     #                              col.names = c("chr", "start", "end", "coverage"))
  d_samples[[sample]]['Zscore'] = (d_samples[[sample]]['coverage'] - mean(d_samples[[sample]]['coverage'])) / mean(d_samples[[sample]]['coverage'])
   
}
    
list_h3k27ac = list()
list_h3k27ac[['E15']] = "Sfd-Cx-E15-N2_bc4_GCACTA_L008_R1_001_peaks.narrowPeak"
list_h3k27ac[['E18']] = "Sfd-Cx-E18-N2_BC10_GGAGAA_L002_R1_001_peaks.narrowPeak"
list_h3k27ac[['P2']] = "P2-H3K27Ac_S51_L008_R1_001_peaks.narrowPeak"
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











