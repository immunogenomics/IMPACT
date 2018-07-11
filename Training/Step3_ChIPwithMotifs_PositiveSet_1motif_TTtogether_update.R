library(GenomicRanges)

train_df <- read.table("/tiffany/IMPACT_manuscript/Training/train_df.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
instances <- read.table("/tiffany/IMPACT_manuscript/Training/TFs.findinstances.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
MotifLibrary <- read.table("/tiffany/IMPACT_manuscript/Training/MotifThresholdList.txt", sep = "\t", header = F, stringsAsFactors = FALSE)

colnames(train_df) <- c("seqnames","starts","ends","names","scores","strands")
#select ChIP peaks with a matching motif under them.
train_df_full <- train_df
m <- match(train_df[,4], instances[,1]) 
train_df <- train_df_full[which(is.na(m) == F),] #get matrix of only peaks with motif instance, now recenter
#add window around peaks 
#windowsize <- 380
train_df_newpeaks <- matrix(0,nrow(train_df),6) #dummy var number of rows 
for (i in 1:nrow(train_df)){
  print(paste0("row ", i, " of ", nrow(train_df)))
  w <- which(instances$PositionID == train_df[i,4]) #match peakID, gauranteed to have a match 
  #only choose 1 motif, take the one with the top score
  w1 <- w[match(max(instances$MotifScore[w]), instances$MotifScore[w])]
  if (instances$Strand[w1] == "+"){
          centerpeak <- instances$Offset[w1] + floor(0.5*(nchar(instances$Sequence[w1]))) #this is relative to start 
        }else{
          centerpeak <- instances$Offset[w1] - floor(0.5*(nchar(instances$Sequence[w1])))
        }
	newstart <- train_df$starts[i]+centerpeak #previous version of impact -windowsize/2 #establish left bound of peaks
        newend <- newstart + 1 #previous version of impact +windowsize #establish right bound of peaks 
        train_df_newpeaks[i,1] <- as.character(train_df$seqnames[i])
        train_df_newpeaks[i,2] <- newstart
        train_df_newpeaks[i,3] <- newend
        train_df_newpeaks[i,4] <- as.character(train_df$names[i]) #preserve old information
        train_df_newpeaks[i,5] <- as.character(train_df$scores[i])
        train_df_newpeaks[i,6] <- as.character(train_df$strands[i])
}

write.table(train_df_newpeaks, "train_df_newpeaks.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

train_df_newpeaks_realdf <- data.frame(seqnames=train_df_newpeaks[,1],
                                       starts=as.numeric(train_df_newpeaks[,2]),
                                       ends=as.numeric(train_df_newpeaks[,3]),
                                       names=train_df_newpeaks[,5],
                                       scores=0,
                                       strands="+")
#shuffle regions

s <- sample(seq(1,nrow(train_df_newpeaks_realdf),1),nrow(train_df_newpeaks_realdf),replace = F)

train_df_newpeaks_realdf <- train_df_newpeaks_realdf[s,]
#put training and test together! 
#use partition function later to make it easier to sample many times! 
train_test_peaks <- unique(train_df_newpeaks_realdf)
colnames(train_test_peaks) <- c("seqnames","starts","ends","names","scores","strands")

#clean up 
if (length(which(train_test_peaks$start < 0)) > 0){
  train_test_peaks$start[which(train_test_peaks$start < 0)] <- 0
}

t <- train_test_peaks[,3]-train_test_peaks[,2]
w <- which(t %% 2 != 0) #odd
if (length(w) > 0){train_test_peaks[w,3] <- train_test_peaks[w,3] + 1}

w1 <- which(train_test_peaks$chr == "chrM")
w2 <- which(nchar(train_test_peaks$chr) > 5)
if (length(c(w1,w2)) > 0){train_test_peaks <- train_test_peaks[-c(w1,w2),]}

write.table(train_test_peaks, "/tiffany/IMPACT_manuscript/Training/train_test_positive_bed.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)








