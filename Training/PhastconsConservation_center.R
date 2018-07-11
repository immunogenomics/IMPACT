args = commandArgs(trailingOnly=TRUE)

#phastcons46way <- read.table("../../../phastCons46wayPrimates.txt.gz", sep = "\t", header = F, stringsAsFactors = FALSE)
#phastcons46way_short <- phastcons46way[,c(2,3,4,13)] #chr, start, end, score
load("phastCons.RData") 

TFs <- args[1]
numTFs <- length(TFs)
for (mastertf in 1:numTFs){
        print(mastertf)
	pos_bed <- read.table(paste0("train_test_positive_bed_",TFs[mastertf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
	neg_bed <- read.table(paste0("train_test_negative_bed_",TFs[mastertf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
	traintest_toannotate <- rbind(pos_bed,neg_bed)

#Phastcons goes here 
AvgConservationScores <- numeric()
for (i in 1:nrow(traintest_toannotate)){
  print(i)
  phastcons46way_chr <- phastcons46way_short[which(phastcons46way_short[,1] == traintest_toannotate[i,1]),]
  w1 <- which(phastcons46way_chr[,2] <= traintest_toannotate[i,2])
  w2 <- which(phastcons46way_chr[,3] >= traintest_toannotate[i,2])
  iw1 <- intersect(w1,w2)
  w1 <- which(phastcons46way_chr[,2] <= traintest_toannotate[i,3])
  w2 <- which(phastcons46way_chr[,3] >= traintest_toannotate[i,3])
  iw2 <- intersect(w1,w2)
  #weighted average of conservation score
  #each region is 1024 bp, so we'd only ever overlap 2!
  iw <- union(iw1,iw2)
  #find correct proportion of overlap
  if (length(iw) == 2){
    proportion_2nd <- (traintest_toannotate[i,3]-phastcons46way_chr[iw[2],2])/(traintest_toannotate[i,3]-traintest_toannotate[i,2])
    proportion_1st <- 1-proportion_2nd
    AvgConservationScore <- phastcons46way_chr[iw[2],4]*proportion_2nd + phastcons46way_chr[iw[1],4]*proportion_1st
  }
  if (length(iw) == 1){
    AvgConservationScore <- phastcons46way_chr[iw,4]
  }
  AvgConservationScores <- c(AvgConservationScores, AvgConservationScore)
}

write.table(AvgConservationScores, paste0("ExtraFeatures_update_",TFs[mastertf],"only_center.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
}




