#do this for each TF pair 
TFpairmat <- matrix(0,6,2)
TFpairmat[1,] <- c(1,2)
TFpairmat[2,] <- c(1,3)
TFpairmat[3,] <- c(1,4)
TFpairmat[4,] <- c(2,3)
TFpairmat[5,] <- c(2,4)
TFpairmat[6,] <- c(3,4)

TFs <- c("Tbet","Gata3","Stat3","Foxp3") 
thresholds <- c(0,seq(0.01,0.1,0.01),0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)

for (j in 1:nrow(TFpairmat)){

print(paste0("Working on TF pair: ",TFs[TFpairmat[j,1]], " and ", TFs[TFpairmat[j,2]]))

datmat <- matrix(0,22,length(thresholds))
 for (i in 1:22){
    print(paste0("chr",i))
    #for (tf in 1:length(TFs)){
    for (tf in 1:2){
    mytf <- TFs[TFpairmat[j,tf]] #was TFs[tf] below
    annot <- read.table(paste0("tiffany/IMPACT_manuscript/GenomeWide/ldsc/",mytf,"_GenomewideTrack_IMPACT.",i,".annot.gz"), sep = "\t", header = T, stringsAsFactors = FALSE)
    assign(paste0("annot",tf), annot)
    }

    for (thr in 1:length(thresholds)){

	for (tf in 1:2){ #was for (tf in 1:length(TFs)){
	annot <- get(paste0("annot",tf))
	w <- which(annot$Chromatin >= thresholds[thr])
	assign(paste0("w",tf), w)
	}
   	#overlap <- intersect(w4,intersect(w3,intersect(w1, w2)))
	overlap <- intersect(w1, w2)
	print(length(overlap))
	datmat[i,thr] <- length(overlap)
   }
}

write.table(datmat,paste0("tiffany/IMPACT_manuscript/GenomeWide/ldsc/FindHowManySNPsareCommonlyAnnotated_",TFs[TFpairmat[j,1]],"vs",TFs[TFpairmat[j,2]],".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

sum_over_chr <- rowSums(t(datmat))
sum_over_chr <- sum_over_chr/sum_over_chr[1]
sum_over_chr <- cbind(thresholds,sum_over_chr)
write.table(sum_over_chr,paste0("tiffany/IMPACT_manuscript/GenomeWide/ldsc/FindHowManySNPsareCommonlyAnnotated_ProprtionoverallChrs_",TFs[TFpairmat[j,1]],"vs",TFs[TFpairmat[j,2]],".txt"),sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

}
