library(GenomicRanges)
library(glmnet)
library(ROCR)
remove(list = ls())

args = commandArgs(trailingOnly=TRUE)

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

load("../IMPACT_modelfit/Features.Rdata")
print("loaded features")
TFs <- args[1]
chrom <- args[2]
numTFs <- length(TFs) 
load(paste0("../IMPACT_modelfit/IMPACT_model_fits_7_",TFs,".RData")) #only models here, no other variables. was 2
print("loaded model fits")
chromosomes <- c(1:22,"X","Y") #exclude X and Y chromosomes- not all features have these. why train on these if no h2g done on them? 
for (mtf in 1:numTFs){ #watch out for variables saved in Rdata files overwriting current variables. 
	mymastertf <- TFs[mtf]
	numextrafeatures_pluslabel <- 2 #one for intercept, one for conservation
	#split up into local and regional coefficients 
	myIMPACTfit <- get(paste0(mymastertf,"_IMPACT_fit"))
	mycoef <- coef(myIMPACTfit)
	intercept <- coefficients(myIMPACTfit)[1]
	conservation <- coefficients(myIMPACTfit)[2]

	coef_local <- mycoef[numextrafeatures_pluslabel+1:(length(datasets)+numextrafeatures_pluslabel)] #skipping intercept and conservation
	coef_regional <- mycoef[(length(datasets)+numextrafeatures_pluslabel+1):length(mycoef)]
	
	local_datasets <- datasets[which(coef_local!=0)]
	regional_datasets <- datasets[which(coef_regional!=0)]
	
	local_betas <- coef_local[which(coef_local!=0)]
	regional_betas <- coef_regional[which(coef_regional!=0)]

        #list of coefficients goes like this: intercept, conservation, 397 datasets (local), 397 datasets (regional)
	#now decide which features we actually care about, each model has an intercept, is conservation beta pushed to zero for each model? usually
        #shift w by 2 and remove 1 for intercept 

	#for (chrom in 1:length(chromosomes)){
		mychromosome <- paste0("chr",chromosomes[chrom])
		print(mychromosome)
		load(paste0(mychromosome,"_bedfile.RData"))	
		print("loaded chromosome R data files")
		#load(paste0(mychromosome,"_conservation.RData"))
		#print("loaded conservation R data files")
		test_forannotating <- cbind(genemat,0,0,"+")	
		#ExtraFeatures_allsnps <- Phastcons
		#print(length(ExtraFeatures_allsnps) == nrow(test_forannotating)) 
		colnames(test_forannotating) <- c("chr", "start", "end", "id", "score", "strand")  
		#dataset_names_1 <- c("Conservation",datasets)
		test_forannotating_df <- as.data.frame(test_forannotating)
		w <- which(ls(1) == "test_forannotating") #save space
		remove(list = ls()[w]) 	
		test_forannotating_df$start <- as.numeric(as.character(test_forannotating_df$start))
       		test_forannotating_df$end <- as.numeric(as.character(test_forannotating_df$end))
		test_bed_1_Granges <- makeGRangesFromDataFrame(test_forannotating_df,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field="end",
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)  	
			
		if (length(local_datasets)>0){
		dataList_local <- Fill.Classifier.Binary.SNPs(test_bed_1_Granges, local_datasets, chromosome = mychromosome)
		}else{dataList_local <- NULL
                local_betas <- NULL}
	
		if (length(regional_datasets)>0){
		

		regdist <- 1000
                hg19_sizes <- read.table("hg19.sizes", sep = "\t", header = F, stringsAsFactors = FALSE)
                chrmax <- hg19_sizes[match(mychromosome, hg19_sizes[,1]),2]
                #bug fix: need to make a test_bed_1_Granges upstream and downstream version and test accordingly
                #we know that genome resolution is 3 bp, so we know which of the bedfile rows will be less than 1000, or more than the length 
                #of the chr when adding the reg region
                test_bed_1_Granges_up <- test_bed_1_Granges
                start(test_bed_1_Granges_up)[1:(ceiling(regdist/3))] <- 0
                start(test_bed_1_Granges_up)[(ceiling(regdist/3)):length(test_bed_1_Granges_up)] <- start(test_bed_1_Granges_up)[(ceiling(regdist/3)):length(test_bed_1_Granges_up)]-regdist
                end(test_bed_1_Granges_up) <- start(test_bed_1_Granges_up) + 1

                test_bed_1_Granges_down <- test_bed_1_Granges
                end(test_bed_1_Granges_down)[(length(test_bed_1_Granges_down) - ceiling(regdist/3)):length(test_bed_1_Granges_down)] <- chrmax
                end(test_bed_1_Granges_down)[1:(length(test_bed_1_Granges_down) - floor(regdist/3))] <- end(test_bed_1_Granges_down)[1:(length(test_bed_1_Granges_down) - floor(regdist/3))] + regdist
                start(test_bed_1_Granges_down) <- end(test_bed_1_Granges_down) - 1
                #where 3 is the basepair resolution.            
                #check granges object is good to go: which(end(test_bed_1_Granges_up) - start(test_bed_1_Granges_up) <= 0)
                #which(end(test_bed_1_Granges_down) - start(test_bed_1_Granges_down) <= 0)
		

		dataList_up <- Fill.Classifier.Binary.SNPs(test_bed_1_Granges_up, regional_datasets, chromosome = mychromosome)
		dataList_down <- Fill.Classifier.Binary.SNPs(test_bed_1_Granges_down, regional_datasets, chromosome = mychromosome)
		dataList_distal <- dataList_up + dataList_down
                dataList_distal[dataList_distal>=1] <- 1
		}else{dataList_distal <- NULL
                regional_betas <- NULL}
		
		classifier <- cbind(dataList_local, dataList_distal)
		mybetas <- c(local_betas, regional_betas)
		
		print("done filling out feature matrix")
		ENet_pred_lambdamin <- 1/ ( 1 + exp(-(intercept + mybetas%*%t(classifier))))
		#to add in conservation, simply add conservation_beta * ExtraFeatures_allsnps after intercept term above
		#scale between 0 and 1 after all chromosomes are predicted (to find max and min values) 
		print("done making predictions")

		gz1 <- gzfile(paste0(mymastertf,"_IMPACT_predictions_",mychromosome,".txt.gz"), "w")
	  	write.table(t(ENet_pred_lambdamin), gz1, sep = "\t", quote = F, row.names = FALSE)
  		close(gz1)

		
#}
}


#next thing is to take a chromosome on which we have a target/nontarget gene, and see if predictions in that area agree with shiny app 





