## function to take gene set and assign promoter and neaerst distal
assign.peak <- function(tss.data, tss.bed, peak.bed, peak.name) {
  # tss.data = tss.merged
  # tss.bed = tss.bed
  # peak.bed = peak.bed
  # peak.name = TF
  # tss.data is the TSS file with RNA expression information.
  # tss.bed is the TSS file in a bed format with gene name (contains promoter regions).
  # peak.bed is the Chip-Seq peaks by TF.
  # peak.name
  
  # generate sorted bed for tss and peak
  tss.bed.sort <- bedr.sort.region(tss.bed, verbose=F)  ##from chr1 to 22 and from front region to rear region
  

  # generate tss overlap
  peak.bed.sort = bedr.sort.region(peak.bed[is.valid.region(peak.bed, verbose=F), ], verbose=F)  ##from chr1 to 22 and from front region to rear region
  
  
  # find a overlapped promoter regions with Chip-Seq peaks  #########standard : tss -> find which tss have peaks
  promoter.overlap <- bedr(
    input = list(a = tss.bed.sort, b = peak.bed.sort), 
    method = "intersect", 
    params = "-wao",
    verbose = F
  )
  
  promoter.overlap <- promoter.overlap[,c(1:4,8)] # 8 is overlapped region; 9 is overlap size
  promoter.overlap <- unique(promoter.overlap)   # there are lots of same rows -> 'unique' makes these unique
  colnames(promoter.overlap) <- c("chr", "start", "end", "name", paste(peak.name, ".promoter", sep=""))
  
  # generate nearest non-overlapping peak   #########standard : peak(ChIP-seq) -> find which peak is not promoter
  # first perform reverse overlap to identify peaks with no promoter
  peak.promoter.overlap <- bedr(
    input = list(a = peak.bed.sort, b = tss.bed.sort),  # change a and b -> means change the standard
    method = "intersect", 
    params = "-wao", 
    verbose=F
  )
  
  peak.promoter.overlap <- subset(peak.promoter.overlap, peak.promoter.overlap[,9]=='0')  # select V9 is 0 -> non overlapping peaks
  peak.promoter.overlap <- bedr.sort.region(peak.promoter.overlap[,1:4], verbose=FALSE)
  
  # non-overlapping peaks will be used for find distal peaks (enhancers)  ## peaks for non-promoter is distal peaks = enhancers
  #########standard : tss -> find among tss which is not promoter
  distal.overlap <- bedr(
    input = list(a = tss.bed.sort, b = peak.promoter.overlap), 
    method = "closest",    # intersect will certainly results in null -> so find the closest promoter ! that will indicate that this peak means the enhancer of that promoter
    params = "-d -io", 
    verbose=F
  )
  
  distal.overlap <- distal.overlap[,c(1:4,8,9)]
  colnames(distal.overlap) <- c("chr", "start", "end", "name", paste(peak.name, ".distal", sep=""), 
                                paste(peak.name, ".distal.distance", sep=""))
  overlap.merge <- merge(promoter.overlap, 
                         subset(distal.overlap, select=-c(chr,start,end)), by='name')
  final.data <- merge(tss.data, 
                      subset(overlap.merge, select=-c(chr,start,end)), by.x='gene_name', by.y = 'name')
  final.data <- unique(final.data)
  
  return(final.data)
}

## function to perform general permutation test and generate histogram output
geneset.perm.test <- function(binary.score.reference, binary.score.test, iterations, plot.name) {
  count.criteria <- vector(length=iterations)
  test.size <- length(binary.score.test)
  for (index in 1:iterations) {
    count.criteria[index] <- sum(sample(binary.score.reference, test.size, replace = F))
  }
  z <- (sum(binary.score.test) - mean(count.criteria))/sd(count.criteria)
  p <- 1 - pnorm(z) # one-sided
  return(list(count.criteria, z, p))
}

## function to perform permutation test and generate histogram output with added step for matched expression distribution
expression.matched.permutation <- function(binary.score.reference, binary.score.test, iterations, plot.name, exp.all.bin, exp.test.bin) {
  bin.value <- unique(exp.test.bin)
  sampled.elements <- vector("list", length=iterations)
  for (bin.index in 1:length(bin.value)) {
    # first generate test bin proportions for matching
    bin.count <- length(which(exp.test.bin==bin.value[bin.index])) #asd
    bin.set <- which(as.character(exp.all.bin)==as.character(bin.value[bin.index]))
    
    # next select elements for bin across iterations
    for (iteration.index in 1:iterations) {
      add.elements <- sample(bin.set, bin.count)
      sampled.elements[[iteration.index]] <-c(sampled.elements[[iteration.index]], add.elements)
    }
  }
  count.criteria <- vector(length=iterations)
  for (iteration.index in 1:iterations) {
    count.criteria[iteration.index] <- sum(binary.score.reference[sampled.elements[[iteration.index]]])
  }
  z <- (sum(binary.score.test) - mean(count.criteria))/sd(count.criteria)
  p <- 1 - pnorm(z) # one-sided
  return(list(count.criteria,z,p))
}


# Function: plot heatmap
plotHeatmap <- function(gsA, gsB, orderA, orderB, title){
  # Run fisher test
  df = run_geneEnrich(gsA, gsB)
  
  # Update the plot setting
  df = df[df$setA != 'g_bg' & df$setB != 'g_bg',] # Remove tests over grey and background
  df$fisher_padj <- p.adjust(df$fisher_p, method = 'bonferroni', n = nrow(df))
  
  df1 = df
  df1$setA_names <- factor(df1$setA, levels=orderA)
  df1$setB_names <- factor(df1$setB, levels=rev(orderB))
  
  df1$OR <- ifelse(df1$fisher_OR<1/8, 1/8, ifelse(df1$fisher_OR >=8, 8, df1$fisher_OR ))
  
  # Plot
  ylabels = c('C25 Dividing RG', 
              'C10 RG', 
              'C7 Dividing IPCs', 
              'C17 IPCs', 
              'C13 NN', 
              'C18 NN', 
              'C3 Early ExN', 
              'C2 Early/Late ExN', 
              'C16 Early/Late ExN', 
              'C8 MGE NPCs', 
              'C11 MGE RG', 
              'C23 MGE-der NN', 
              'C1 Striatal InN', 
              'C6 MGE-der InN', 
              'C15 CGE-der InN', 
              'C20 Choroid', 
              'C9 Endothelial', 
              'C4 OPC/Astro', 
              'C19 Microglia')
  
  p <- ggplot(df1, aes(setA_names, setB_names)) + 
    geom_tile(fill='white', size=0.1, color='grey', show.legend = F) +
    geom_point(aes(fill = log2(OR)), 
               shape=22, size=1.5, color='black', show.legend = T, stroke=0.1) +
    geom_point(data= df1[df1$fisher_padj > 0.05 & df1$fisher_p <= 0.05,],
               aes(fill = log2(OR)), 
               shape=22, size=4, color='black', show.legend = F, stroke=0.1) +
    geom_tile(data= df1[df1$fisher_padj<= 0.05,],
              aes(fill = log2(OR)),
              size=0.1, color='black', show.legend = F) +
    geom_point(data= df1[df1$fisher_padj<= 0.05,], 
               shape=8, size=1, color='black', show.legend = F) +
    scale_fill_gradient2(low="#6A90CA", mid = 'white', high="#CD2836", 
                         limits=c(log2(1/8),log2(8))) +
    labs(title = title, x='', y='') +
    coord_equal() +
    theme_minimal(base_size=7) + 
    scale_y_discrete(labels=rev(ylabels)) + 
    theme(axis.ticks=element_blank(), 
          panel.border=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x=element_text(size=7, angle=320, hjust = 0, colour="grey50") ) 
  
  out = list(df,p)
  return(out)
}

# Function: Run Fisher
run_geneEnrich <- function(genesetA, genesetB){
  # genesetA: Gene list of interest (e.g. single cell cluster genes, eGenes)
  # genesetB: Gene list for testing (e.g. gene set for autism neurobiology)
  # genesetA_bg: Background gene list of interest
  # genesetB_bg: Background gene list of a cell type
  
  # Find the background gene list
  genesetA_bg = sort(unique(unlist(genesetA)))
  genesetB_bg = sort(unique(unlist(genesetB)))
  
  print (length(genesetA_bg))
  print (length(genesetB_bg))
  
  # Run fisher exact test
  cols = c('setA', 'setB', 
           'setA_size', 'setB_size', 
           'setA_size_background', 'setB_size_background', 
           'setA_size_hits', 
           'fisher_OR', 'fisher_p', 'overlaps')
  res = as.data.frame(matrix(nrow=0, ncol=length(cols)))
  for (genesetB_1 in names(genesetB)){
    print (genesetB_1)
    genesetB1 = genesetB[[genesetB_1]]
    genesetB_bg1 = genesetB_bg[!(genesetB_bg %in% genesetB1)]
    for (genesetA_1 in names(genesetA)){
      genesetA1 = genesetA[[genesetA_1]]
      genesetA_bg1 = genesetA_bg[!(genesetA_bg %in% genesetA1)]
      
      # create a matrix
      mat = matrix( c( sum(genesetB1 %in% genesetA1), sum(genesetA_bg1 %in% genesetB1),
                       sum(genesetB_bg1 %in% genesetA1), sum(genesetA_bg1 %in% genesetB_bg1)), 
                    ncol=2   )
      
      fis = fisher.test(mat)
      
      out = as.data.frame(matrix(c(genesetA_1, genesetB_1,
                                   length(genesetA1), length(genesetB1), 
                                   length(genesetA_bg), length(genesetB_bg), 
                                   sum(genesetB1 %in% genesetA1), 
                                   fis$estimate, fis$p.value, 
                                   paste( intersect(genesetB1, genesetA1), collapse=',' )) , ncol=length(cols)))
      colnames(out) = cols
      res = rbind.data.frame(res, out)
    }
  }
  
  # update the type
  res$fisher_OR <- as.numeric(res$fisher_OR)
  res$fisher_p <- as.numeric(res$fisher_p)
  
  return(res)
}

