#make bedfiles of nucleotide resolution
  remove(list = ls())
  hg19_sizes <- read.table("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/hg19.sizes", sep = "\t", header = F, stringsAsFactors = FALSE)
  chromosomes <- c(seq(1:22),"X","Y")
  for (i in 1:length(chromosomes)){
    print(i)
    chr <- paste0("chr",chromosomes[i])
    start <- 1
    end <- hg19_sizes[i,2] #must end on this nucleotide (cannot go past)
    mod <- end-start
    resolution <- 3 #sample every 4 nucleotides.
    
    remainder <- end%%resolution
    if (remainder == 0){end <- end - 1} #end value is in limbo position, new end is one bp before. 
    if (remainder == 1){end <- end - 2} #end value is in the right position, new end is two bp before.
    if (remainder == 2){} #end value is in the left position, keep it there.
    
    #rowdim <- length(seq(start,end,resolution)) #slow
    rowdim <- floor(end/3) + end%%3
    genemat <- matrix(0,rowdim,3) #large matrix
    genemat[,1] <- chr
    genemat[,2] <- seq(start,end,resolution) #cannot have overlap! cannot be next to each other! need at least gap of 1 between end[N] and start[N+1], difference of 2
    genemat[,3] <- as.numeric(genemat[,2])+1 #each is 1 bp wide
    save(list = ls(all.names = TRUE), file = paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",chr,"_bedfile.RData"), envir = .GlobalEnv)
    remove(list = "genemat")
  }
  
