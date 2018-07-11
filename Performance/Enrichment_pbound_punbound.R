##Look for enrichment in label = 1 group vs label = 0 group. 
##is the mean prediction value similar (over positive and negative) within and between TF models? 
TFs <- c("Tbet","Gata3","Stat3","Foxp3")

#run once to calculate maxmin, then rerun to normalize 
#not genomewide annotation yet, need to scale manually here 
maxmin <- matrix(0,4,2)
for (j in 1:4){
  TF_to_plot <- TFs[j]
  iters <- 10
  mean_predictions <- numeric()
  positive_values <- numeric()
  negative_values <- numeric()
  for (i in 1:iters){
    predictions <- read.table(paste0("/tiffany/IMPACT_manuscript/Performance/PredictionOutput/Predictions_",TF_to_plot,"_trial",i,".txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    positive_values <- c(positive_values, predictions[which(predictions[,1] == 1),2])
    negative_values <- c(negative_values, predictions[which(predictions[,1] == 0),2])
    mean_predictions <- c(mean_predictions,mean(predictions[,2]))
  }
    maxmin[j,1] <- min(c(positive_values, negative_values))
    maxmin[j,2] <- max(c(positive_values, negative_values))
}

for (j in 1:4){
  TF_to_plot <- TFs[j]
  iters <- 10
  mean_predictions <- numeric()
  positive_values <- numeric()
  negative_values <- numeric()
  for (i in 1:iters){
    predictions <- read.table(paste0("/tiffany/IMPACT_manuscript/Performance/PredictionOutput/Predictions_",TF_to_plot,"_trial",i,".txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
    predictions[,2] <- (predictions[,2]- maxmin[j,1])/(maxmin[j,2]-maxmin[j,1])
    positive_values <- c(positive_values, predictions[which(predictions[,1] == 1),2])
    negative_values <- c(negative_values, predictions[which(predictions[,1] == 0),2])
    mean_predictions <- c(mean_predictions,mean(predictions[,2]))
  }
  assign(paste0("positive_values_",TF_to_plot),positive_values)
  assign(paste0("negative_values_",TF_to_plot),negative_values)
  print(log(x = t.test(positive_values, negative_values, alternative = "greater")$p.val, base = 10))
  print(mean(mean_predictions)) #NEED TO SCALE want this number to be similar for each TF 
}

pdf("/tiffany/IMPACT_manuscript/Performance/Supp_Fig_1_ProbRegEnrichinClasses.pdf")
boxplot(positive_values_Tbet,negative_values_Tbet, positive_values_Gata3, negative_values_Gata3, positive_values_Stat3, negative_values_Stat3, positive_values_Foxp3, negative_values_Foxp3, 
        col = c("beige", "beige", "plum2", "plum2","palegreen3","palegreen3","slateblue1","slateblue1"), outline = F, 
        names = c("Th1+", "Th1-", "Th2+", "Th2-", "Th17+", "Th17-", "Treg+", "Treg-"), ylab = "IMPACT predictions", main = "Active class enriched for higher predictions")
dev.off()

#this value is very similar
#new: (200 training, 2000 test)
#Tbet: 0.08931407
#Gata3: 0.09115423
#Stat3: 0.08105479
#Foxp3: 0.09217836

#sd AUCs, sensitivity 
AUCs <- numeric()
Sens <- numeric()
for (j in 1:4){
  TF_to_plot <- TFs[j]
  mymat <- read.table(paste0("/tiffany/IMPACT_manuscript/Performance/Sens_AUC_collection_",TF_to_plot,".txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
  AUCs <- c(AUCs, sd(as.numeric(mymat[2,])))
  Sens <- c(Sens, sd(as.numeric(mymat[1,])))
}

#Supplementary Table 2
#Beta variability analysis 10 trial samplings.
TFs <- c("Tbet","Gata3","Stat3","Foxp3")
featureTFmat <- matrix(0,796,4*3) #4 types of measurements
#mean, sd, num non zeros, direction change
for (j in 1:4){ #4TFs
  featuremat <- matrix(0,796,10) #796 = 396*2 + 1 (cons) + 1 (int)
  for (tr in 1:10){ #10 trials
  TF_to_plot <- TFs[j]
  mymat <- read.table(paste0("/tiffany/IMPACT_manuscript/Performance/PredictionOutput/Betas_",TF_to_plot,"_trial",tr,".txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
  featuremat[,tr] <- mymat$V1
  }  
  avgvalue <- sapply(1:nrow(featuremat),function(x) mean(featuremat[x,]))
  stdevs <- sapply(1:nrow(featuremat),function(x) sd(featuremat[x,]))
  num_nonzeros <- sapply(1:nrow(featuremat),function(x) length(which(featuremat[x,]!=0)))
  featureTFmat[,j*3-2] <- avgvalue
  featureTFmat[,j*3-1] <- stdevs
  featureTFmat[,j*3] <- num_nonzeros
}
#Supplementary Table 2
write.table(featureTFmat, "/tiffany/IMPACT_manuscript/Performance/featureTFmat.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)


#for genome-wide annotations later
#only computed betas once. 
#Supplementary Table 1
load("/tiffany/IMPACT_manuscript/GenomeWide/IMPACTmodel/Features.Rdata") 
#make a matrix for all TFs 
mat <- matrix(0,nrow(coef(Tbet_IMPACT_fit)),4)
mat[,1] <- coef(Tbet_IMPACT_fit)[,1]
mat[,2] <- coef(Gata3_IMPACT_fit)[,1]
mat[,3] <- coef(Stat3_IMPACT_fit)[,1]
mat[,4] <- coef(Foxp3_IMPACT_fit)[,1]

#Supplementary Table 1
write.table(mat, "/tiffany/IMPACT_manuscript/Performance/featuremat_genomewidetrack.txt", sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)





