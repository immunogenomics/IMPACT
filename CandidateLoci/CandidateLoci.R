#want to find how many peaks there are genome-wide where the score exceeds some threshold and doesn't overlap with
#ChIP data. 
#go back on this decision: Make sure that we are only comparing to motif centered pseudo peaks and not chip. 


TFs <- c("Tbet","Gata3","Stat3","Foxp3")
threshold <- 0.9

for (tf in 1:length(TFs)){
	TF_to_plot <- TFs[tf]

#	TF_chip <- read.table(paste0("../CD4Tcell_MasterTFs/train_test_positive_bed_",TF_to_plot,"only_center.txt"), sep = "\t", 
#			header = F, stringsAsFactors = FALSE) #want all unsampled regions. motif + ChIP. pseudo peaks. 

	TF_chip <- read.table(paste0("/tiffany/IMPACT_manuscript/TargetvNonTargets/",TF_to_plot,"_ChIPseq_Inputs.bed.gz"), sep = "\t", 
                       header = F, stringsAsFactors = FALSE)
	chromosomes <- c(seq(1:22),"X","Y")
	#candidateloci_count <- 0
	count <- 0
	homerbed <- matrix(0,500000,3)
	for (j in 1:length(chromosomes)){
		print(j)
		chr <- paste0("chr",chromosomes[j])
		pred <- read.table(paste0("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/",TF_to_plot,"_IMPACT_scaled_predictions_",chr,".txt.gz"), 
			sep = "\t", header = T, stringsAsFactors = FALSE)

		#zero out training regions so they are not considered.      
		TF_chip_subset <- TF_chip[which(TF_chip[,1]==chr),]	
		starts <- floor(TF_chip_subset[,2]/3)
		ends <- floor(TF_chip_subset[,3]/3)
		for (i in 1:nrow(TF_chip_subset)){
			pred[starts[i]:ends[i],1] <- 0
		}
		#candidateloci_count <- candidateloci_count + length(which(as.numeric(pred[,1]) >= 0.9))	
		
		#want to create bed file of these. taken from FollowUponFP.R local
		w_highscore <- which(as.numeric(pred[,1]) >= threshold)
		if (length(w_highscore)>1){
		    #need to make sure this is one continuous region greater than 0.8
		    newpeak <- 0
		    lastpeakcoord <- 1
		    for (n in 1:(length(w_highscore)-1)){
		      diff <- w_highscore[n+1]-w_highscore[n]
		      if (diff > 1){ #we've found a peak. 
			#starting boundaries of peak, replace 1 with boundary of last peak in j space
		        peakmin <- min(w_highscore[lastpeakcoord:n]) 
		        peakmax <- max(w_highscore[lastpeakcoord:n]) #ending boundaries of peak
		        #translate to coordinates 
		        peakmincoord <- peakmin*3
		        peakmaxcoord <- peakmax*3
		        peakchr <- chr
		        #update lastpeakcoord
                          lastpeakcoord <- n+1
		          newpeak <- 1 + newpeak
		          count <- 1 + count
          	          homerbed[count,1] <- peakchr
	                  homerbed[count,2] <- peakmincoord
		          homerbed[count,3] <- peakmaxcoord
        }
      }
    }
  }

	m <- min(which(homerbed[,1] == 0))
	homerbed2 <- homerbed[1:(m-1),]
	homerbed <- homerbed2		
	homerbed <- cbind(homerbed, paste0("peak_",1:nrow(homerbed)), 0, "+")
 	write.table(homerbed, paste0("/tiffany/IMPACT_manuscript/CandidateLoci/",TF_to_plot,"_CandidateLociBedFile.txt"), quote = F, row.names = FALSE, col.names = FALSE)
	#candidateloci_count <- candidateloci_count*3
	#print(candidateloci_count)
	#write.table(candidateloci_count, paste0("Total_bpCoverage_CandidateLoci_",TF_to_plot,"_IMPACTmodel.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
	#}
}



