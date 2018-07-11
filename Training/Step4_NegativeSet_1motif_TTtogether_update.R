library(GenomicRanges)

train_df_newpeaks <- read.table("/tiffany/IMPACT_manuscript/Training/train_df_newpeaks.txt", sep = "\t", header = F, stringsAsFactors = FALSE)

MotifThreshold <- read.table("/tiffany/IMPACT_manuscript/Training/MotifThresholdList.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
motifs <- MotifThreshold[1,1]
instscan <- read.table("/tiffany/IMPACT_manuscript/Training/scanMotifsgenomewide.1.15000.sort.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
instscan <- instscan[order(instscan[,6], decreasing = T),]   
bed <- read.table("/tiffany/IMPACT_manuscript/Training/train_df.txt", sep = "\t", header = F, stringsAsFactors = FALSE) #old version: ChIPbed1.txt
colnames(bed) <- c('chr','start','end','id','score','strand') #score may be number of reads
Granges <- with(bed, GRanges(chr, IRanges(start,end), strand, score, id = id))
assign(paste0(motifs,"_Granges"),Granges)

allinstscan <- instscan
sapply_instscan <- sapply(1:nrow(allinstscan), function(x) strsplit(allinstscan[x,1], "-")[[1]][1])
allinstscan[,1] <- sapply_instscan
allinstscan[,6] <- allinstscan[,6]/MotifThreshold[1,2]

w <- which(allinstscan[,6] < 1) 
if(length(w)>0){allinstscan <- allinstscan[-w,]}

motifmatrix <- matrix(0,nrow(MotifThreshold),3)
#col1: motif name with cell type (1 per row)
#col2: degenerate motif name (used in scanMotifsGenome.pl), get with table(instances_scan[,1])
#col3: Granges files
motifmatrix[,1] <- motifs 
motifmatrix[,2] <- motifs
motifmatrix[,3] <- paste0(motifs,"_Granges")

#windowsize <- 380
print(paste0("working on ",motifmatrix[1,1]))
  w_mcp <- which(train_df_newpeaks[,5] == motifmatrix[1,2]) #motifcentered peaks
  motifcentered_specificTF <- train_df_newpeaks[w_mcp,]
  colnames(motifcentered_specificTF) <- c('chr','start','end','id','score','strand') 
  Granges_obj_mc <- with(motifcentered_specificTF, GRanges(chr, IRanges(start,end), strand, score, id = id))
  
  w_inst <- which(allinstscan[,1] == motifmatrix[1,2]) #instances pertaining to specific motif
  instances_specificTF <- allinstscan[w_inst,]
  inst_bed <- instances_specificTF[,c(2,3,4)]
  inst_bed <- cbind(inst_bed, motifmatrix[1,1], instances_specificTF[,6], "+") 
  #region as wide as a motif 
  motif_width <- as.numeric(inst_bed[1,3]) - as.numeric(inst_bed[,2])
  inst_bed[,2] <- as.numeric(inst_bed[,2]) + floor(motif_width/2)
  inst_bed[,3] <- as.numeric(inst_bed[,2]) + 1

  #inst_bed[,3] <- as.numeric(inst_bed[,3])#previous version of impact: +(windowsize-10)/2
  #inst_bed[,2] <- as.numeric(inst_bed[,2])#previous version of impact: -(windowsize-10)/2
  inst_bed <- unique(inst_bed)
  colnames(inst_bed) <- c('chr','start','end','id','score','strand') 
  inst_bed_Granges <- with(inst_bed, GRanges(chr, IRanges(start,end), strand, score, id = id)) #instances Granges obj
  
  int_1 <- intersect(inst_bed_Granges, get(motifmatrix[1,3]))
  int_2 <- intersect(inst_bed_Granges, Granges_obj_mc)
  
  ms_1 <- match(start(int_1)-1, inst_bed$start)
  me_1 <- match(end(int_1), inst_bed$end)
  a <- c(ms_1,me_1)
  a <- unique(a[is.na(a) == F])
  
  ms_2 <- match(start(int_2)-1, inst_bed$start)
  me_2 <- match(end(int_2), inst_bed$end)
  b <- c(ms_2,me_2)
  b <- unique(b[is.na(b) == F])
  
  ab <- unique(c(a,b))
  if (length(ab) > 0){noChIP <- inst_bed[-ab,] #good negative peaks
  }else{
    noChIP <- inst_bed
  }
write.table(noChIP, "/tiffany/IMPACT_manuscript/Training/noChIP.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
colnames(noChIP) <- c("seqnames","starts","ends","names","scores","strands")

#shuffle 
s <- sample(seq(1,nrow(noChIP),1),nrow(noChIP),replace = F)
noChIP <- noChIP[s,]
noChIP <- unique(noChIP)

##clean up 
if (length(which(noChIP$start < 0))){
  noChIP$start[which(noChIP$start < 0)] <- 0
}

t <- noChIP[,3]-noChIP[,2]
w <- which(t %% 2 != 0) #odd
if (length(w) > 0){noChIP[w,3] <- noChIP[w,3] + 1}

w1 <- which(noChIP$chr == "chrM")
w2 <- which(nchar(noChIP$chr) > 5)
if (length(c(w1,w2)) > 0){noChIP <- noChIP[-c(w1,w2),]}

write.table(noChIP, "/tiffany/IMPACT_manuscript/Training/train_test_negative_bed.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)












