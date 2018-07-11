args = commandArgs(trailingOnly=TRUE)

TFs <- args
for (tf in 1:length(TFs)){

  resolution <- 3
  for (i in 1:22){
    print(i)
    chr <- paste0("chr",i) 
    pred <- read.table(paste0("../GenomeTracks/",TFs[tf],"_IMPACT_scaled_predictions_",chr,".txt.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)


	#s-LDSC
	annot <- read.table(paste0("EUR.",i,".annot.gz"), sep = "\t", header = T, stringsAsFactors = FALSE) 
	myrow <- floor(annot$BP/resolution)
  	annot$Chromatin <- pred[myrow,1]
  	gz1 <- gzfile(paste0(TFs[tf],"_GenomewideTrack_IMPACT.",i,".annot.gz"), "w") #no per chromosome scaling. 
  	write.table(annot, gz1, sep = "\t", quote = F, row.names = FALSE)
  	close(gz1)


	#annot <- read.table(paste0("EAS.",i,".annot.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)
        #myrow <- floor(annot$BP/resolution)
        #annot$Chromatin <- pred[myrow,1]
        #gz1 <- gzfile(paste0(TFs[tf],"_GenomewideTrack_IMPACT_EAS.",i,".annot.gz"), "w") #no per chromosome scaling. 
        #write.table(annot, gz1, sep = "\t", quote = F, row.names = FALSE)
        #close(gz1)
}

}


