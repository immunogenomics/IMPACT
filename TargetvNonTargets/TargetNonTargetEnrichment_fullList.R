#look for target enrichment over non targets. 

#get files with ChIP targets 
Tbet_targets <- read.table("/tiffany/IMPACT_manuscript/TargetvNonTargets/Tbet_ChIP_target_list_pcgenes_fix.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
Gata3_targets <- read.table("/tiffany/IMPACT_manuscript/TargetvNonTargets/Gata3_ChIP_target_list_pcgenes_fix.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
Stat3_targets <- read.table("/tiffany/IMPACT_manuscript/TargetvNonTargets/Stat3_ChIP_target_list_pcgenes_fix.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
Foxp3_targets <- read.table("/tiffany/IMPACT_manuscript/TargetvNonTargets/Foxp3_ChIP_target_list_pcgenes_fix.txt", sep = "\t", header = F, stringsAsFactors = FALSE)

allgenes <- read.table("/tiffany/IMPACT_manuscript/TargetvNonTargets/allgenes.txt", sep = "\t", header = F, stringsAsFactors = FALSE) #can't reduce, because genes will combine with one another. 
w <- which(nchar(allgenes$V1)>5) #remove hap1/hap2 genes
allgenes <- allgenes[-w,]
allgenenames <- allgenes[,4] #not unique names 
uniquegenenames <- unique(allgenenames) #redefine nontarget genes as genes that aren't in the target list. 

TFs <- c("Tbet","Gata3","Stat3","Foxp3")

#extract only for gene. select maybe all random target and all non target genes
regregion <- 20000
#target_mat <- matrix(0,1000,1)
#nontarget_mat <- matrix(0,1000,1)
for (i in 1:length(TFs)){
  #target_predictions:
  tf_targets <- get(paste0(TFs[i],"_targets"))
  tf_targets <- tf_targets$V1
  m <- match(tf_targets, uniquegenenames)
  m <- m[is.na(m) == F]
  possible_nontargets <- uniquegenenames[-m] #remove if there is a match 
   
  target_mat <- matrix(0,length(tf_targets),1)
  nontarget_mat <- matrix(0,length(possible_nontargets),1)
  for (tar in 1:nrow(target_mat)){
    print(tar)
    gene <- tf_targets[tar]
    w <- which(allgenenames == gene) #don't want unique list here
    mystart <- min(allgenes$V2[w])
    myend <- max(allgenes$V3[w])
    mychr <- allgenes$V1[w[1]]
    
    myend <- myend + regregion
    mystart <- mystart - regregion
    mystart <- max(0,mystart)
    maxrange <- myend - mystart
    pred <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",TFs[i],"_IMPACT_scaled_predictions_",mychr,".txt.gz"), skip = (floor(mystart/3)), nrows = ceiling(maxrange/3), header = T) #added 1 to skip because we want to skip the header row. but skip until the row before the actual start, so take the floor instead. 
    target_mat[tar,1]<- mean(pred[,1])
    print(mean(pred[,1]))
    }

   for (tar in 1:nrow(nontarget_mat)){      
    print(tar)
    gene <- possible_nontargets[tar]
    w <- which(allgenenames == gene)
    mystart <- min(allgenes$V2[w])
    myend <- max(allgenes$V3[w])
    mychr <- allgenes$V1[w[1]]
    
    myend <- myend + regregion
    mystart <- mystart - regregion
    mystart <- max(0,mystart)
    maxrange <- myend - mystart
    pred <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",TFs[i],"_IMPACT_scaled_predictions_",mychr,".txt.gz"), skip = (floor(mystart/3)), nrows = ceiling(maxrange/3), header = T) #added 1 to skip because we want to skip the header row. but skip until the row before the actual start, so take the floor instead. 
    nontarget_mat[tar,1]<- mean(pred[,1])
    print(mean(pred[,1]))
   }
  write.table(target_mat, paste0("/tiffany/IMPACT_manuscript/TargetvNonTargets/TargetGenePredictionsPerTF_",TFs[i],".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
  write.table(nontarget_mat, paste0("/tiffany/IMPACT_manuscript/TargetvNonTargets/NonTargetGenePredictionsPerTF_",TFs[i],".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) 
}

