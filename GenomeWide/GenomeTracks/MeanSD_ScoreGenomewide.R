TFs <- c("Tbet","Gata3","Stat3","Foxp3")
for (j in 1:length(TFs)){
 meanmat <- matrix(0,24,1) #chromosomes by TF
 sdmat <- matrix(0,24,1)
 TF <- TFs[j]
 chromosomes <- c(seq(1:22),"X","Y")
 for (i in 1:length(chromosomes)){
    	chr <- chromosomes[i]
	print(i)
	dat <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",TF,"_IMPACT_scaled_predictions_chr",chr,".txt.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)
  	meanmat[i,1] <- mean(as.numeric(dat$x), na.rm = T)
	sdmat[i,1] <- sd(as.numeric(dat$x), na.rm = T)
	print(meanmat[i,1])
	print(sdmat[i,1])
	}

print(meanmat)
write.table(meanmat, paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/MeanChrWidePredictionsperTF_",TF,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) 
write.table(sdmat, paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/StDevChrWidePredictionsperTF_",TF,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
}
