TFs <- c("Tbet","Gata3","Stat3","Foxp3")
for (j in 1:length(TFs)){
 minmat <- matrix(0,24,1)
 maxmat <- matrix(0,24,1) 
 TF <- TFs[j]
 chromosomes <- c(seq(1:22),"X","Y")
 for (i in 1:length(chromosomes)){
    	chr <- chromosomes[i]
	print(i)
	dat <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",TF,"_IMPACT_predictions_chr",chr,".txt.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)
  	minmat[i,1] <- min(dat$V1, na.rm = T)
	maxmat[i,1] <- max(dat$V1, na.rm = T)
	print(minmat[i,1])
	print(maxmat[i,1])
	}

write.table(minmat, paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/MinChrWidePredictionsperTF_",TF,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) 
write.table(maxmat, paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/MaxChrWidePredictionsperTF_",TF,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) 
}


for (j in 1:length(TFs)){
 TF <- TFs[j]
 minmat <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/MinChrWidePredictionsperTF_",TF,".txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
 maxmat <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/MaxChrWidePredictionsperTF_",TF,".txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
 chromosomes <- c(seq(1:22),"X","Y")
 for (i in 1:length(chromosomes)){
        chr <- chromosomes[i]
        print(i)
        mymin <- as.numeric(minmat[i,1])
        mymax <- as.numeric(maxmat[i,1])
        dat <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",TF,"_IMPACT_predictions_chr",chr,".txt.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)
        dat_scaled <- (dat$V1 - mymin)/(mymax-mymin)

        gz1 <- gzfile(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",TF,"_IMPACT_scaled_predictions_chr",chr,".txt.gz"), "w")
        write.table(dat_scaled, gz1, sep = "\t", quote = F, row.names = FALSE)
        close(gz1)
        }
}




