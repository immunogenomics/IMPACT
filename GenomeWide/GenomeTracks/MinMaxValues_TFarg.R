args = commandArgs(trailingOnly=TRUE)
TFs <- args
for (j in 1:length(TFs)){
 minmat <- matrix(0,24,1) #chromosomes by TF
 maxmat <- matrix(0,24,1)
 TF_to_plot <- TFs[j]
 chromosomes <- c(seq(1:22),"X","Y")
 for (i in 1:length(chromosomes)){
    	chr <- chromosomes[i]
	print(i)
	dat <- read.table(paste0(TF_to_plot,"_IMPACT_predictions_chr",chr,".txt.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)
  	minmat[i,j] <- min(dat$V1, na.rm = T)
	maxmat[i,j] <- max(dat$V1, na.rm = T)
	print(minmat[i,j])
	print(maxmat[i,j])
	}

write.table(minmat, paste0("MinChrWidePredictionsperTF_",TF_to_plot,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) 
write.table(maxmat, paste0("MaxChrWidePredictionsperTF_",TF_to_plot,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) 
}

