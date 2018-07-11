#want to find how many peaks there are genome-wide where the score exceeds some threshold and doesn't overlap with
#ChIP data. 
#go back on this decision: Make sure that we are only comparing to motif centered pseudo peaks and not chip. 

library(GenomicRanges)

TFs <- c("Tbet","Gata3","Stat3","Foxp3")
for (tf in 4:length(TFs)){
        TF_to_plot <- TFs[tf]
	dat <- read.table(paste0("/tiffany/IMPACT_manuscript/CandidateLoci/",TF_to_plot,"_CandidateLociBedFile.txt"), sep = " ", header = F, stringsAsFactors = FALSE)
	colnames(dat) <- c("chr", "start", "end", "id", "score", "strand")	
	dat_Granges <- with(dat, GRanges(chr, IRanges(start,end), strand, score, id = id))
	print(sum(width(reduce(dat_Granges))))
}
