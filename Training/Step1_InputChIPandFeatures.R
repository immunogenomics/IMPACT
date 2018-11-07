args = commandArgs(trailingOnly=TRUE) #creates string vector with arguments

#1: tf 2: cell-type 3: num input files

motif <- paste0(args[1],"_",args[2])

#Input ChIP
chip <- matrix(0,2,3)
for (i in 1:as.numeric(args[3])){
  inputchip <- read.table(paste0("Input_",args[1],"_",i,".bed.gz"), header = F, stringsAsFactors = FALSE)
  chip <- rbind(chip, inputchip[,1:3])
}
chip <- chip[-c(1,2),]
bed <- unique(cbind(chip[,1:3], 0,motif,"+"))
bed[,4] <- paste0("peak",seq(1,nrow(bed),1))

#reduce regions to non-overlapping. 
library(GenomicRanges)
colnames(bed) <- c('chr','start','end','id','score','strand') #score may be number of reads
bed_Granges <- with(bed, GRanges(chr, IRanges(start,end), strand, score, id = id))
bed_Granges <- reduce(bed_Granges)
bed <- cbind(as.character(seqnames(bed_Granges)), start(bed_Granges), end(bed_Granges), paste0("peak",seq(1,length(bed_Granges),1)), "+")
write.table(bed, "train_df.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

pwm <- read.table(paste0("Motif_",motif,".txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
write.table(pwm, "TF.1.PWM.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

#Input Motifs
MotifThreshvals <- pwm[which(pwm[,4]=="PH"),3]
MotifThreshnames <- pwm[which(pwm[,4]=="PH"),2]
MotifThreshinfo <- cbind(MotifThreshnames, MotifThreshvals)

write.table(pwm, "AllTFs_PWM.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
write.table(MotifThreshinfo, "MotifThresholdList.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)













