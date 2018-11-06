remove(list = ls()) 

library(GenomicRanges)
#datasets
celltype_directories <- c("Tcell","FetalBrain","OtherCellTypes","CellTypeNonSpecific")
count <- 0
for (i in 1:length(celltype_directories)){

FeatureChIPs <- read.table(paste0("FeatureDir/",celltype_directories[i],"_features_integernames/",celltype_directories[i],"_Features_ExtensionFinucane2015.txt"),sep = "\t", header = T, stringsAsFactors = FALSE)
filenum <- FeatureChIPs$new_file_num
        for (j in 1:nrow(FeatureChIPs)){
        count <- count + 1
	print(count)
        ds <- read.table(paste0("FeatureDir/",celltype_directories[i],"_features_integernames/gzip/",filenum[j],".bed.gz"), sep = "\t", header = F, stringsAsFactors = FALSE)
        bed <- unique(cbind(ds[,1:3], 0,0,"+"))
        colnames(bed) <- c('chr','start','end','id','score','strand') #score may be number of reads
        Granges <- with(bed, GRanges(chr, IRanges(start,end), strand, score, id = id))
        assign(paste0("DS",count,"_Granges"),Granges)
        }
}
print(count)
datasets <- paste0("DS",1:count,"_Granges")
dataset_names <- c("label","Conservation",datasets)
save(list = ls(all.names = TRUE), file = "Features.Rdata")

#load("Features.Rdata")

