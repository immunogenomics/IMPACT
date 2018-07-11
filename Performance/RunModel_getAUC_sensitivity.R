#List of Functions

library(glmnet)
library(ROCR)
library(GenomicRanges)

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


#### Load_and_Train_Model_AUC_CVglmnet_LabeledData ####
Load_and_Train_Model_AUC_CVglmnet_LabeledData <- function(classifier_all, TP_training_regions, TN_training_regions, specificity_cutoff){
  
  sa_matrix <- matrix(0,2,1)
  train <- classifier_all[1:(TP_training_regions + TN_training_regions),2:ncol(classifier_all)] 
  train_labels <- classifier_all[1:(TP_training_regions + TN_training_regions),1]
  test <- classifier_all[(TP_training_regions + TN_training_regions+1):nrow(classifier_all),2:ncol(classifier_all)] 
  test_labels <- classifier_all[(TP_training_regions + TN_training_regions+1):nrow(classifier_all),1]
  
  ENet_fit <- cv.glmnet(x=train[complete.cases(train),], y= train_labels[complete.cases(train)], family = "binomial", type.measure = "auc", alpha = 0.5)
  ENet_pred_lambdamin <- predict(ENet_fit,test[complete.cases(test),],s="lambda.min", type = "response") #type = response ensures that the scale is from 0 to 1 
  
  pred <- prediction(ENet_pred_lambdamin, test_labels[complete.cases(test)])
  perf <- performance(pred, 'sens', 'spec') #x is spec, y is sens
  ix <- which.min(abs(perf@x.values[[1]] - specificity_cutoff))
  sensitivity <- perf@y.values[[1]][ix]
  auc <- round(slot(performance(pred, measure = "auc"), "y.values")[[1]],4)
  sa_matrix[1,1] <- sensitivity
  sa_matrix[2,1] <- auc 

  test_pos <- test[which(test_labels == 1),]
  test_pos_dim_comp <- nrow(test_pos[complete.cases(test_pos),])

  test_neg <- test[which(test_labels == 0),]
  test_neg_dim_comp <- nrow(test_neg[complete.cases(test_neg),])

  Labels <- c(rep(1,test_pos_dim_comp), rep(0,test_neg_dim_comp))

  newList <- list("betas" = coef(ENet_fit, s = "lambda.min"), "prediction" = ENet_pred_lambdamin, "SensitivityAUC" = sa_matrix, "test_pos_dim_complete" = test_pos_dim_comp, "test_neg_dim_complete" = test_neg_dim_comp, "Labels" = Labels)

  return(newList)
}


load("/tiffany/IMPACT_manuscript/GenomeWide/IMPACT_modelfit/Features.Rdata")

TFs <- c("Tbet","Gata3","Stat3","Foxp3")
numTFs <- length(TFs) #Tbet Th1, Gata3 Th2, Stat3 Th17, Foxp3 Treg
for (tf in 1:numTFs){

	#load training data for each TF
	ExtraFeatures <- read.table(paste0("/tiffany/IMPACT_manuscript/Features/ExtraFeatures_update_",TFs[tf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
	positive_bed <- read.table(paste0("/tiffany/IMPACT_manuscript/Training/train_test_positive_bed_",TFs[tf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)
	negative_bed <- read.table(paste0("/tiffany/IMPACT_manuscript/Training/train_test_negative_bed_",TFs[tf],"only_center.txt"),sep = "\t", header = F, stringsAsFactors = FALSE)

	ExtraFeatures_positive <- ExtraFeatures[1:nrow(positive_bed),1]
	ExtraFeatures_negative <- ExtraFeatures[(nrow(positive_bed)+1):nrow(ExtraFeatures),1]
	
	
	#remove all the chr_hap6_blah_blah in the negative files. later use 10000 regions to train, not 12000
	#should have been removed in an earlier step.  
	w_p <- which(nchar(positive_bed[,1])>5)
        if (length(w_p)>0){positive_bed <- positive_bed[-w_p,]
	ExtraFeatures_positive <- ExtraFeatures_positive[-w_p]}
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

	#split by chromosome, like we do for SNPs, don't use X and Y? 
	hg19_sizes <- read.table("/tiffany/IMPACT_manuscript/Performance/hg19.sizes", sep = "\t", header = F, stringsAsFactors = FALSE)
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
		s <- sapply(1:nrow(dummy), function (x) ifelse(dummy[x,3] > hg19_sizes[m[x],2], hg19_sizes[m[x],2], 0))
		w <- which(s != 0)
		if (length(w)>0){dummy[,3] <- s[w]
		dummy[,2] <- s[w]-1} 
		#w <- which(end(dummy_granges)>max_on_chr) #check if downstream end is > max chromosome length
		#if (length(w)>0){dummy_granges <- dummy_granges[-w]}
		colnames(dummy) <- c('chr','start','end','id','score','strand')
	        dummy_Granges <- with(dummy, GRanges(chr, IRanges(start,end), strand, score, id = id))
		assign(paste0("regional",regionaltypes[j],"_",traintypes[i]), dummy_Granges)		
		}
	}	
	chromosomes <- paste0("chr",c(1:22,"X","Y"))
	for (chrom in 1:length(chromosomes)){
	#for (chrom in 1:2){
		
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
	

	#CV
	numiters <- 10
	Sens_AUC_collection <- matrix(0,2,numiters)
	for (iter in 1:numiters){
	print(paste0("iteration ",iter))
	#training <- 0.8
	#need Load_and_Train_Model_AUC_CVglmnet_LabeledData input matrix to be trainP, trainN, testP, testN
	s_p <- sample(1:nrow(classifier_TP_train), size = 1000, replace = F)
	s_n <- sample(1:nrow(classifier_TN_train), size = 10000, replace = F)
	classifier <- rbind(classifier_TP_train[s_p[1:800],], classifier_TN_train[s_n[1:8000],], classifier_TP_train[s_p[801:1000],],classifier_TN_train[s_n[8001:10000],])
		
	LATM <- Load_and_Train_Model_AUC_CVglmnet_LabeledData(classifier, 800, 8000, 0.99)

	#want AUC and sensitivity at FPR 0.99. 
	Sens_AUC <- LATM$SensitivityAUC
	print(Sens_AUC)
	Sens_AUC_collection[,iter]<- Sens_AUC

	Predictions <- LATM$prediction
	Labels <- LATM$Labels
	predictionmat <- matrix(0,length(Labels),2)
	predictionmat[,1] <- Labels
	predictionmat[,2] <- Predictions
	write.table(predictionmat,paste0("/tiffany/IMPACT_manuscript/Performance/PredictionOutput/Predictions_",TFs[tf],"_trial",iter,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
	Betas <- LATM$betas #want to show that betas don't change much between samplings 
	write.table(Betas[,1],paste0("/tiffany/IMPACT_manuscript/Performance/PredictionOutput/Betas_",TFs[tf],"_trial",iter,".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
	}
	write.table(Sens_AUC_collection, paste0("/tiffany/IMPACT_manuscript/Performance/Sens_AUC_collection_",TFs[tf],".txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
}

  



