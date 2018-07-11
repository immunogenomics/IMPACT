##define targets by genes with Tbet in their regulatory regions.
#look at tbet input (combination of 3 chip-seqs)
#get gene info for 20,000 PC snps and check each for intersection 
#if intersection, target, if no intersection, not a target
library(GenomicRanges)

TFs <- c("Tbet","Gata3","Stat3","Foxp3")
for (j in 1:length(TFs)){
  
TF_to_plot <- TFs[j]
TF_chip <- read.table(paste0("/tiffany/IMPACT_manuscript/TargetvNonTargets/",TF_to_plot,"_ChIPseq_Inputs.bed.gz"), sep = "\t", header = F, stringsAsFactors = FALSE)
TF_chip <- cbind(TF_chip[,1:3],0,0,"+")
colnames(TF_chip) <- c('chr','start','end','id','score','strand')

# allgenes <- read.table("allgenesfromUCSC.txt.gz", sep = "\t", header = F, stringsAsFactors = FALSE)
# w <- which((as.numeric(allgenes$V4)-as.numeric(allgenes$V3)) == 0) #remove all 0 width
# allgenes <- allgenes[-w,]
# allgenes2 <- cbind(allgenes[,1], allgenes[,3],allgenes[,4],allgenes[,5],0,allgenes[,2])
# write.table(allgenes2, "allgenes.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

allgenes <- read.table(paste0("/tiffany/IMPACT_manuscript/TargetvNonTargets/allgenes.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #can't reduce, because genes will combine with one another. 

pc_genes <- unique(allgenes$V4)
target_list <- numeric()
nontarget_list <- numeric()
regulatoryregion <- 20000 #20 kb- martha
for (i in 1:length(pc_genes)){ #20431
  print(i)
  mypcgene <- allgenes[which(allgenes$V4 == pc_genes[i]),]
  #add regulatory region, 20 kb up and downstream.
  start <- min(mypcgene[,2]-regulatoryregion)
  end <- max(mypcgene[,3]+regulatoryregion)
  #colnames(mypcgene) <- c('chr','start','end','id','score','strand')
  #mypcgene_granges <- with(mypcgene, GRanges(chr, IRanges(start,end), strand, score, id = id))
  #r <- intersect(mypcgene_granges, TF_chip_granges) #isn't picking up, like MNase- if first region is larger than the others, it returns 0.  
  TF_chip_chr <- TF_chip[(TF_chip[,1]==mypcgene[,1]),]
  w1 <- which(TF_chip_chr[,2]>= start)
  w2 <- which(TF_chip_chr[,3]<= end)
  r <- intersect(w1,w2)
  if (length(r) == 0){
    #print("nontarget")
    #not a target 
    nontarget_list <- c(nontarget_list, pc_genes[i])
  }else{
    #print("target")
    #a target 
    target_list <- c(target_list, pc_genes[i])
  }
} 
write.table(target_list, paste0("/tiffany/IMPACT_manuscript/TargetvNonTargets/",TFs[j],"_ChIP_target_list_pcgenes_fix.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
#write.table(nontarget_list, paste0("/tiffany/IMPACT_manuscript/TargetvNonTargets/",TFs[j],"_ChIP_nontarget_list_pcgenes_fix.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
}





























# #from target lists, intersect with 200 top STRING interactions. how many of 200 also have chip?  
# #from nontarget list, exclude any genes also in the STRING list. 
# 
# #remove # from first line
# string_TBX21 <- read.table("/Users/amariuta/Documents/SRLAB/CD4RegMap/TBX21_string_200interactions.tsv", quote = "", header = T, row.names = NULL, stringsAsFactors = FALSE)  
# string_GATA3 <- read.table("/Users/amariuta/Documents/SRLAB/CD4RegMap/GATA3_200interactions.tsv", quote = "", header = T, row.names = NULL, stringsAsFactors = FALSE)  
# string_STAT3 <- read.table("/Users/amariuta/Documents/SRLAB/CD4RegMap/STAT3_200interactions.tsv", quote = "", header = T, row.names = NULL, stringsAsFactors = FALSE)  
# string_FOXP3 <- read.table("/Users/amariuta/Documents/SRLAB/CD4RegMap/FOXP3_200interactions.tsv", quote = "", header = T, row.names = NULL, stringsAsFactors = FALSE)  
# 
# interactors_TBX21 <- unique(c(string_TBX21$node1[which(string_TBX21$node2 == "TBX21")], string_TBX21$node2[which(string_TBX21$node1 == "TBX21")])) #152, up to 200, wont show if there aren't more, these are the other TFs in the Th1 program we are trying to define. 
# interactors_GATA3 <- unique(c(string_GATA3$node1[which(string_GATA3$node2 == "GATA3")], string_GATA3$node2[which(string_GATA3$node1 == "GATA3")])) #180
# interactors_STAT3 <- unique(c(string_STAT3$node1[which(string_STAT3$node2 == "STAT3")], string_STAT3$node2[which(string_STAT3$node1 == "STAT3")])) #200
# interactors_FOXP3 <- unique(c(string_FOXP3$node1[which(string_FOXP3$node2 == "FOXP3")], string_FOXP3$node2[which(string_FOXP3$node1 == "FOXP3")])) #200
# 
# 
# 
# #remove all chromosomes that are chr6_apd_hap1, chr6_cox_hap2, etc 
# 
# allgenes <- read.table(paste0("/Users/amariuta/allgenes2.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #can't reduce, because genes will combine with one another. 
# w <- which(nchar(allgenes$V1)>5) #remove hap1/hap2 genes
# allgenes <- allgenes[-w,]
# 
# #TFs <- c("Tbet","Gata3","Stat3","Foxp3")
# TFs <- "Foxp3"
# #Proteins <- c("TBX21","GATA3","STAT3","FOXP3")
# Proteins <- "FOXP3"
# for (j in 1:length(TFs)){
#   targets <- read.table(paste0("/Users/amariuta/Documents/SRLAB/CD4RegMap/",TFs[j],"_ChIP_target_list_pcgenes_fix.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #can't reduce, because genes will combine with one another. 
#   interactors <- get(paste0("interactors_",Proteins[j]))
#   confident_targetgenes <- intersect(targets$V1, interactors)
#   
#   #get gene coordinate info. match to col 4 of allgenes
#   nucleotideres_bedfile_allgenes <- matrix(0,2,6) #initiate
#   for (ct in 1:length(confident_targetgenes)){
#     print(ct)
#     w <- which(allgenes[,4] == confident_targetgenes[ct])
#     mypcgene <- allgenes[w,]
#     chr <- unique(mypcgene[,1]) #enforce that there is just 1 unique value
#     start <- min(mypcgene[,2]-regulatoryregion)
#     end <- max(mypcgene[,3]+regulatoryregion)
#     mod <- end-start
#     resolution <- 10 #sample every 10 nucleotides.
#       while(mod%%resolution != 0){
#         start <- start-1
#         mod <- end-start
#       }
#       rowdim <- length(seq(start,end,resolution))
#       genemat <- matrix(0,rowdim,6)
#       genemat[,1] <- chr
#       genemat[,2] <- seq(start,end,resolution) #cannot have overlap! cannot be next to each other! need at least gap of 1 between end[N] and start[N+1], difference of 2
#       genemat[,3] <- as.numeric(genemat[,2])+2 #each is 2 bp wide
#       genemat[,4] <- paste0("peak",1:nrow(genemat),"_",confident_targetgenes[ct],"_","target") #need to be unique for homer, include gene and target/nontarget status
#       genemat[,5] <- 0
#       genemat[,6] <- "+"
#       nucleotideres_bedfile_allgenes <- rbind(nucleotideres_bedfile_allgenes, genemat)
#     }
#     nucleotideres_bedfile_allgenes <- nucleotideres_bedfile_allgenes[-c(1,2),]
#     print(nrow(nucleotideres_bedfile_allgenes))
#     write.table(nucleotideres_bedfile_allgenes, paste0("/Users/amariuta/Documents/SRLAB/CD4RegMap/target.",TFs[j],".lowres",resolution,".20kbeitherside.200STRINGintChIP.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) #200 string interactors intersected with ChIP
#     nucleotideres_bedfile_allgenes[,4] <- paste0("peak",1:nrow(nucleotideres_bedfile_allgenes))
#     write.table(nucleotideres_bedfile_allgenes, paste0("/Users/amariuta/Documents/SRLAB/CD4RegMap/target.",TFs[j],".lowres",resolution,".20kbeitherside.200STRINGintChIP.annotate.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
#     
# 
#   nontargets <- read.table(paste0("/Users/amariuta/Documents/SRLAB/CD4RegMap/",TFs[j],"_ChIP_nontarget_list_pcgenes_fix.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #can't reduce, because genes will combine with one another. 
#   nontargets <- nontargets$V1
#   #do any match to STRING? remove these 
#   removethese <- intersect(interactors, nontargets) #remove these genes
#   m <- match(removethese,nontargets)
#   if (length(m)>0){nontargets <- nontargets[-m]}
#   #get gene coordinate info. match to col 4 of allgenes
#   #then sample 200 of them. 
#   nontargets <- nontargets[sample(1:length(nontargets),size = 200, replace = F)]
#   nucleotideres_bedfile_allgenes <- matrix(0,2,6) #initiate
#   for (ct in 1:length(nontargets)){
#     print(ct)
#     w <- which(allgenes[,4] == nontargets[ct])
#     if (length(w)>0){ #otherwise just skip and go onto the next gene. 
#       mypcgene <- allgenes[w,]
#       chr <- unique(mypcgene[,1]) #enforce that there is just 1 unique value, skip if gene is on multiple chromosomes. 
#       if (length(chr)==1){
#         start <- min(mypcgene[,2]-regulatoryregion)
#         end <- max(mypcgene[,3]+regulatoryregion)
#         mod <- end-start
#         resolution <- 10 #sample every 10 nucleotides.
#         while(mod%%resolution != 0){
#           start <- start-1
#           mod <- end-start
#         }
#         rowdim <- length(seq(start,end,resolution))
#         genemat <- matrix(0,rowdim,6)
#         genemat[,1] <- chr
#         genemat[,2] <- seq(start,end,resolution) #cannot have overlap! cannot be next to each other! need at least gap of 1 between end[N] and start[N+1], difference of 2
#         genemat[,3] <- as.numeric(genemat[,2])+2 #each is 2 bp wide
#         genemat[,4] <- paste0("peak",1:nrow(genemat),"_",nontargets[ct],"_","nontarget") #need to be unique for homer, include gene and target/nontarget status
#         genemat[,5] <- 0
#         genemat[,6] <- "+"
#         nucleotideres_bedfile_allgenes <- rbind(nucleotideres_bedfile_allgenes, genemat)
#       }
#     }
#   }
#   nucleotideres_bedfile_allgenes <- nucleotideres_bedfile_allgenes[-c(1,2),]
#   print(nrow(nucleotideres_bedfile_allgenes))
#   write.table(nucleotideres_bedfile_allgenes, paste0("/Users/amariuta/Documents/SRLAB/CD4RegMap/nontarget.",TFs[j],".lowres",resolution,".20kbeitherside.200STRINGintChIP.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE) #200 string interactors intersected with ChIP
#   nucleotideres_bedfile_allgenes[,4] <- paste0("peak",1:nrow(nucleotideres_bedfile_allgenes))
#   write.table(nucleotideres_bedfile_allgenes, paste0("/Users/amariuta/Documents/SRLAB/CD4RegMap/nontarget.",TFs[j],".lowres",resolution,".20kbeitherside.200STRINGintChIP.annotate.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
#   
# }
# 



