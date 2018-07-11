args = commandArgs(trailingOnly=TRUE)

TFs <- args

for (j in 1:length(TFs)){
 TF <- TFs[j]
 chromosomes <- c(seq(1:22),"X","Y")
 minmat <- read.table(paste0("MinChrWidePredictionsperTF_",TF,".txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
 maxmat <- read.table(paste0("MaxChrWidePredictionsperTF_",TF,".txt"),sep = "\t", header = F, stringsAsFactors = FALSE)


 for (i in 1:length(chromosomes)){
    	chr <- chromosomes[i]
	print(i)
	mymin <- as.numeric(minmat[i,j])
	mymax <- as.numeric(maxmat[i,j])
	dat <- read.table(paste0(TF,"_IMPACT_predictions_chr",chr,".txt.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)
	dat_scaled <- (dat$V1 - mymin)/(mymax-mymin)

	gz1 <- gzfile(paste0(TF,"_IMPACT_scaled_predictions_chr",chr,".txt.gz"), "w")
        write.table(dat_scaled, gz1, sep = "\t", quote = F, row.names = FALSE)
        close(gz1)
	}
}

