args = commandArgs(trailingOnly=TRUE)
TFs <- strsplit(args, split = ",")[[1]]
print(TFs)
#1: tf 


library(GenomicRanges)
library(glmnet)
library(ROCR)



Fill.Classifier.Binary.SNPs <- function(my_Granges, datasets, chromosome){

  classifier_test <- matrix(0,length(my_Granges),length(datasets))

  for (j in 1:length(datasets)){
    print(paste0("working on dataset ",j))
    curr_dataset <- get(datasets[j])
    curr_dataset <- curr_dataset[seqnames(curr_dataset) == chromosome,]
    result <- intersect(my_Granges,curr_dataset)

    if (length(result)>0){
      mastermatrix <- data.frame(result_start=start(result),result_end=end(result),matchstartTPTN=match(start(result),
                 start(my_Granges)),matchendTPTN=match(end(result), end(my_Granges)),TPTN_index=rep(0,length(result)))
      #if there's a match to either the end or the start, that's a 1 
      #if there's no match to either the end of the start of TPTN testrain beds, that's a 0
      indices <- na.omit(c(mastermatrix$matchstartTPTN, mastermatrix$matchendTPTN))
      classifier_test[indices, j] <- 1
    }
  }
  return(classifier_test)
}


load("Features.Rdata")

#TFs <- args[1]
numTFs <- length(TFs) 
for (mastertf in 1:numTFs){
	print(TFs[mastertf])

	      #load training data for each TF
	      ExtraFeatures <- read.table(paste0("../../CD4Tcell_MasterTFs/testGitPipeline/ExtraFeatures_update_",TFs[mastertf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
	      positive_bed <- read.table(paste0("../../CD4Tcell_MasterTFs/testGitPipeline/train_test_positive_bed_",TFs[mastertf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
	      negative_bed <- read.table(paste0("../../CD4Tcell_MasterTFs/testGitPipeline/train_test_negative_bed_",TFs[mastertf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)

	ExtraFeatures_positive <- ExtraFeatures[1:nrow(positive_bed),1]
	ExtraFeatures_negative <- ExtraFeatures[(nrow(positive_bed)+1):nrow(ExtraFeatures),1]
	#remove all the chr_hap6_blah_blah in the negative files. later use 10000 regions to train, not 12000
        #should have been removed in an earlier step.  
	w_p <- which(positive_bed[,1] == "chrM")
        if (length(w_p)>0){positive_bed <- positive_bed[-w_p,]
	ExtraFeatures_positive <- ExtraFeatures_positive[-w_p]}
        w_p <- which(nchar(positive_bed[,1])>5)
        if (length(w_p)>0){positive_bed <- positive_bed[-w_p,]
	ExtraFeatures_positive <- ExtraFeatures_positive[-w_p]}
	w_n <- which(negative_bed[,1] == "chrM")
	if (length(w_n)>0){negative_bed <- negative_bed[-w_n,]
	ExtraFeatures_negative <- ExtraFeatures_negative[-w_n]}
        w_n <- which(nchar(negative_bed[,1])>5)
        if (length(w_n)>0){negative_bed <- negative_bed[-w_n,]
	ExtraFeatures_negative <- ExtraFeatures_negative[-w_n]}

        #make Granges objects
        colnames(positive_bed) <- c('chr','start','end','id','score','strand')
        colnames(negative_bed) <- c('chr','start','end','id','score','strand')
        positive_bed_Granges <- with(positive_bed, GRanges(chr, IRanges(start,end), strand, score, id = id))
        negative_bed_Granges <- with(negative_bed, GRanges(chr, IRanges(start,end), strand, score, id = id))

        numextrafeatures_pluslabel <- 1+1 #label, conservation at center nucleotide

        #build feature matrix for positive and negative sets, then sample from each. 
        #length(datasets)*2 because each feature has a local and regional/distal subfeature 
        classifier_TP_train <- matrix(0,nrow(positive_bed),(length(datasets)*2+numextrafeatures_pluslabel))
        classifier_TP_train[,1] <- 1
        classifier_TP_train[,2] <- ExtraFeatures_positive
        classifier_TN_train <- matrix(0,nrow(negative_bed),(length(datasets)*2+numextrafeatures_pluslabel))
        classifier_TN_train[,1] <- 0
        classifier_TN_train[,2] <- ExtraFeatures_negative

	hg19_sizes <- read.table("../GenomeTracks/hg19.sizes", sep = "\t", header = F, stringsAsFactors = FALSE)
        traintypes <- c("positive","negative")
        regionaltypes <- c("upstream","downstream")
        regionaldists <- c(-1000,1000)
        for (i in 1:2){
                for (j in 1:2){
                dummy <- get(paste0(traintypes[i],"_bed"))
                dummy[,2] <- dummy[,2] + regionaldists[j]
                dummy[,3] <- dummy[,3] + regionaldists[j]
                w <- which(dummy[,2]<=0) #check if upstream start is < 0
                if (length(w)>0){dummy[w,2] <- 1
                dummy[w,3] <- 2}
                m <- match(dummy[,1],hg19_sizes[,1])
		s <- sapply(1:nrow(dummy), function (x) if(dummy[x,3] > hg19_sizes[m[x],2]){dummy[x,3] <- hg19_sizes[m[x],2]
                dummy[x,2] <- dummy[x,3]-1
                }) #if it is past, use end, if not don't do anything
                colnames(dummy) <- c('chr','start','end','id','score','strand')
                dummy_Granges <- with(dummy, GRanges(chr, IRanges(start,end), strand, score, id = id))
                assign(paste0("regional",regionaltypes[j],"_",traintypes[i]), dummy_Granges)
                }
        }
		chromosomes <- paste0("chr",c(1:22,"X","Y"))
	      for (chrom in 1:length(chromosomes)){
                print(chromosomes[chrom])

                #POSITIVE
                #local
                w <- which(seqnames(positive_bed_Granges) == chromosomes[chrom])
		if (length(w)>0){
                dataList <- Fill.Classifier.Binary.SNPs(positive_bed_Granges[w], datasets, chromosome = chromosomes[chrom])
                classifier_TP_train[w,(numextrafeatures_pluslabel+1):(numextrafeatures_pluslabel+length(datasets))] <- dataList
                #distal upstream
                dataList_up <- Fill.Classifier.Binary.SNPs(regionalupstream_positive[w], datasets, chromosome = chromosomes[chrom])
                #distal downstream
                dataList_down <- Fill.Classifier.Binary.SNPs(regionaldownstream_positive[w], datasets, chromosome = chromosomes[chrom])
                #if either of the two. 
                dataList_distal <- dataList_up + dataList_down
                dataList_distal[dataList_distal>=1] <- 1 #change all values > 1 to 1
                classifier_TP_train[w,(numextrafeatures_pluslabel+length(datasets)+1):ncol(classifier_TP_train)] <- dataList_distal
		}
	
                #NEGATIVE
                #local
                w <- which(seqnames(negative_bed_Granges) == chromosomes[chrom])
        	if (length(w)>0){
	        dataList <- Fill.Classifier.Binary.SNPs(negative_bed_Granges[w], datasets, chromosome = chromosomes[chrom])
                classifier_TN_train[w,(numextrafeatures_pluslabel+1):(numextrafeatures_pluslabel+length(datasets))] <- dataList
                #distal upstream
                dataList_up <- Fill.Classifier.Binary.SNPs(regionalupstream_negative[w], datasets, chromosome = chromosomes[chrom])
                #distal downstream
                dataList_down <- Fill.Classifier.Binary.SNPs(regionaldownstream_negative[w], datasets, chromosome = chromosomes[chrom])
                #if either of the two. 
                dataList_distal <- dataList_up + dataList_down
                dataList_distal[dataList_distal>=1] <- 1 #change all values > 1 to 1
                classifier_TN_train[w,(numextrafeatures_pluslabel+length(datasets)+1):ncol(classifier_TN_train)] <- dataList_distal
		}
         }
        
	
      	s_pos <- sample(1:nrow(classifier_TP_train), size = min(c(nrow(classifier_TP_train),1000)), replace = F)
      	classifier_TP_train <- classifier_TP_train[s_pos,]
      	s_neg <- sample(1:nrow(classifier_TN_train), size = 10000, replace = F)
      	classifier_TN_train <- classifier_TN_train[s_neg,]
      
      	classifier <- rbind(classifier_TP_train, classifier_TN_train)	
      	train <- classifier[,-1] 
      	train_labels <- classifier[,1]
      	ENet_fit <- cv.glmnet(x=train[complete.cases(train),], y= train_labels[complete.cases(train)], family = "binomial", type.measure = "auc", alpha = 0.5) #alpha: mixing term, lasso, 1-alpha ridge. 
      
      	assign(paste0(TFs[mastertf],"_IMPACT_fit"), ENet_fit)	

w1 <- which(ls(1) == paste0(TFs[mastertf],"_IMPACT_fit"))
save(list = ls(1)[w1], file = paste0("IMPACT_model_fits_7_",TFs[mastertf],".RData"), envir = .GlobalEnv)


}









