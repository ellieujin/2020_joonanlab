library(bedr)
library(dplyr)
library(ggplot2)
library(cowplot)

CODA_1 <- read.table("/home/titan/manta/half_samples_200527/results/variants/diploidSV.vcf", 
                           header = F, sep = "\t", stringsAsFactors = T, quote = "", 
                           col.names = c('CHROM',	'POS', 'ID',	'REF',	'ALT', 'QUAL',	'FILTER',	'INFO',	'FORMAT',	
                                         'SNU008-JJE-M',	'JEW-P',	'SNU004-LSM-M',	'SNU009-AJW-M',	'SNU008-JJE-PM',	
                                         'SNU001-YSW-M',	'SNU001-YSW-PM',	'SNUB-F0061-1LJW',	'SNUB-F0061-1JYJ',	
                                         'SNU007-KYJ-P',	'SNU002-YSM-PM',	'SNU004-LSM-P',	'JEW-M',	'SNU006-KHW-M',	
                                         'SNU005-KHS-P',	'SNU007-KYJ-PM',	'SNU006-KHW-P',	'SNU002-YSM-M',	'SNU006-KHW-PM',	
                                         'SNU009-AJW-PM',	'JEW-PM',	'SNUB-F0061-2LSB',	'SNU004-LSM-PM',	'SNU005-KHS-PM',	
                                         'SNU002-YSM-P',	'SNU007-KYJ-M'))

CODA_2 <- read.table("/home/titan/manta/half_samples_200528/results/variants/diploidSV.vcf", 
                     header = F, sep = "\t", stringsAsFactors = T, quote = "", 
                     col.names = c('CHROM',	'POS', 'ID',	'REF',	'ALT', 'QUAL',	'FILTER',	'INFO',	'FORMAT',	
                                   'SNU017-KTY-M', 'SNU013-CB-PM', 'SNU020-JEJ-P', 'SNU012-PJY-M', 'SNU013-CB-M', 
                                   'SNU019-PHE-PM', 'SNU017-KTY-P', 'SNU014-PTJ-PM', 'SNU016-JCW-PM', 'SNU019-PHE-M', 
                                   'SNU015-MSO-PM', 'SNU011-KHR-M', 'SNU014-PTJ-P', 'SNU012-PJY-PM', 'SNU017-KTY-PM', 
                                   'SNU016-JCW-M', 'SNU011-KHR-PM', 'SNU013-CB-P', 'SNU014-PTJ-M', 'SNU018-CYJ-M', 
                                   'SNU018-CYJ-PM', 'SNU019-PHE-P', 'SNU010-KSH-M', 'SNU020-JEJ-PM', 'SNU010-KSH-PM'))

CODA <- merge(x=CODA_1, y=CODA_2, by = c('CHROM',	'POS', 'ID',	'REF',	'ALT', 'QUAL',	'FILTER',	'INFO',	'FORMAT'), all =T)


DNV_1 <- read.table("/home/titan/manta/CODA_DNV_CNV_200527.txt",
                    header = F, sep = "\t", stringsAsFactors = F, quote = "",
                    comment.char = "")

DNV_2 <- read.table("/home/titan/manta/CODA_DNV_CNV_200528.txt",
                    header = F, sep = "\t", stringsAsFactors = F, quote = "",
                    comment.char = "")
  #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SNU017-KTY-M SNU013-CB-PM SNU020-JEJ-P SNU012-PJY-M SNU013-CB-M SNU019-PHE-PM SNU017-KTY-P SNU014-PTJ-PM SNU016-JCW-PM SNU019-PHE-M SNU015-MSO-PM SNU011-KHR-M SNU014-PTJ-P SNU012-PJY-PM SNU017-KTY-PM SNU016-JCW-M SNU011-KHR-PM SNU013-CB-P SNU014-PTJ-M SNU018-CYJ-M SNU018-CYJ-PM SNU019-PHE-P SNU010-KSH-M SNU020-JEJ-PM SNU010-KSH-PM
