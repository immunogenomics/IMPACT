########################################################################################################
##plotting changing distal parameter
TFplotnames <- c("Tbet Th1","Gata3 Th2","Stat3 Th17","Foxp3 Treg")
TFs <- c("Tbet","Gata3","Stat3","Foxp3")
type <- c(20,100,200,400,1000,2000,20000) #bps
type <- paste0(type," bp")
mean_sd_sens_stuff <- matrix(0,length(type),2)
for (i in 1:length(TFs)){
  model <- read.table(paste0("/tiffany/IMPACT_manuscript/KnockOuts/Sens_collection_",TFs[i],"_changeDistal.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_sens_stuff[,1] <- sapply(1:nrow(model), function(x) mean(as.numeric(model[x,])))
  mean_sd_sens_stuff[,2] <- sapply(1:nrow(model), function(x) sd(as.numeric(model[x,])))
  assign(paste0("mean_sd_sens_stuff_",TFs[i]),mean_sd_sens_stuff) 
}

# TFplotnames <- c("Tbet Th1","Gata3 Th2","Stat3 Th17","Foxp3 Treg")
# TFs <- c("Tbet","Gata3","Stat3","Foxp3")
# type <- c(20,100,200,400,1000,2000,20000) #bps
# type <- paste0(type," bp")
# mean_sd_sens_stuff <- matrix(0,length(type),2)
# for (i in 1:length(TFs)){
#   model <- read.table(paste0("/tiffany/IMPACT_manuscript/KnockOuts/AUC_collection_",TFs[i],"_changeDistal.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
#   mean_sd_sens_stuff[,1] <- sapply(1:nrow(model), function(x) mean(as.numeric(model[x,])))
#   mean_sd_sens_stuff[,2] <- sapply(1:nrow(model), function(x) sd(as.numeric(model[x,])))
#   assign(paste0("mean_sd_sens_stuff_",TFs[i]),mean_sd_sens_stuff) 
# }

TFcolors <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
pdf("/tiffany/IMPACT_manuscript/KnockOuts/ChangingDistalParamPlot.pdf")
par(mfrow=c(2,2))
for(i in 1:length(TFs)){
  info <- get(paste0("mean_sd_sens_stuff_",TFs[i]))
  bplot <- barplot(as.numeric(info[,1]), col = TFcolors[i], ylim = c(0,1), 
                   xpd = F, ylab = "Sensitivity at 0.99 Specificity", main = TFplotnames[i]) #1 if you want auc, 3 if you want sensitivity, ylim = c(0.45,1) if auc
  abline(h = 0, col = "black", lwd = 1)
  text(x = bplot, y = 0-0.02, type, xpd = T, srt = 30, adj = 1, cex = 0.7)
  for (a in 1:nrow(info)){
    segments(x0 = bplot[a], x1 = bplot[a], y0 = info[a,1]-2*info[a,2], y1 = info[a,1]+2*info[a,2],col = "black", lwd = 2) #2 if you want auc sd, 4 if you want sensitivity sd
  }
}
dev.off()

########################################################################################################
#Knock Out feature categories 
TFplotnames <- c("Tbet Th1","Gata3 Th2","Stat3 Th17","Foxp3 Treg")
TFs <- c("Tbet","Gata3","Stat3","Foxp3")
type <- c("no KO","cell-state","H3K4me1","H3K27ac","DNase","ATAC","Promoter","H3K4me3","H3K9ac","H3K27me3") #feature categories
mean_sd_sens_stuff <- matrix(0,length(type),2)
for (i in 1:length(TFs)){
  #original
  model <- read.table(paste0("/tiffany/IMPACT_manuscript/Performance/Sens_AUC_collection_",TFs[i],".txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_sens_stuff[1,1] <- mean(as.numeric(model[1,]))
  mean_sd_sens_stuff[1,2] <- sd(as.numeric(model[1,]))
  #cell-state ko
  model <- read.table(paste0("tiffany/IMPACT_manuscript/KnockOuts/Sens_AUC_collection_",TFs[i],"_KOcellstatefeatures.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_sens_stuff[2,1] <- mean(as.numeric(model[1,]))
  mean_sd_sens_stuff[2,2] <- sd(as.numeric(model[1,]))
  #feature categories ko
  model <- read.table(paste0("tiffany/IMPACT_manuscript/KnockOuts/Sens_collection_",TFs[i],"_kofeatures.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_sens_stuff[3:length(type),1] <- sapply(1:nrow(model), function(x) mean(as.numeric(model[x,])))
  mean_sd_sens_stuff[3:length(type),2] <- sapply(1:nrow(model), function(x) sd(as.numeric(model[x,])))
  assign(paste0("mean_sd_sens_stuff_",TFs[i]),mean_sd_sens_stuff) 
}

# TFplotnames <- c("Tbet Th1","Gata3 Th2","Stat3 Th17","Foxp3 Treg")
# TFs <- c("Tbet","Gata3","Stat3","Foxp3")
# type <- c("no KO","cell-state","H3K4me1","H3K27ac","DNase","ATAC","Promoter","H3K4me3","H3K9ac","H3K27me3") #feature categories
# mean_sd_sens_stuff <- matrix(0,length(type),2)
# for (i in 1:length(TFs)){
#   #original
#   model <- read.table(paste0("tiffany/IMPACT_manuscript/Performance/Sens_AUC_collection_",TFs[i],".txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
#   mean_sd_sens_stuff[1,1] <- mean(as.numeric(model[2,]))
#   mean_sd_sens_stuff[1,2] <- sd(as.numeric(model[2,]))
#   #cell-state ko
#   model <- read.table(paste0("tiffany/IMPACT_manuscript/KnockOuts/Sens_AUC_collection_",TFs[i],"_KOcellstatefeatures.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
#   mean_sd_sens_stuff[2,1] <- mean(as.numeric(model[2,]))
#   mean_sd_sens_stuff[2,2] <- sd(as.numeric(model[2,]))
#   #feature categories ko
#   model <- read.table(paste0("tiffany/IMPACT_manuscript/KnockOuts/AUC_collection_",TFs[i],"_kofeatures.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
#   mean_sd_sens_stuff[3:length(type),1] <- sapply(1:nrow(model), function(x) mean(as.numeric(model[x,])))
#   mean_sd_sens_stuff[3:length(type),2] <- sapply(1:nrow(model), function(x) sd(as.numeric(model[x,])))
#   assign(paste0("mean_sd_sens_stuff_",TFs[i]),mean_sd_sens_stuff) 
# }

TFcolors <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
pdf("tiffany/IMPACT_manuscript/KnockOuts/Knockoutcategories.pdf")
par(mfrow=c(2,2))
for(i in 1:length(TFs)){
  info <- get(paste0("mean_sd_sens_stuff_",TFs[i]))
  bplot <- barplot(as.numeric(info[,1]), col = TFcolors[i], ylim = c(0,1), 
                   xpd = F, ylab = "Sensitivity at 0.99 Specificity", main = TFplotnames[i]) #1 if you want auc, 3 if you want sensitivity, ylim = c(0.45,1) if auc
  abline(h = 0, col = "black", lwd = 1)
  text(x = bplot, y = 0-0.02, type, xpd = T, srt = 30, adj = 1, cex = 0.7)
  for (a in 1:nrow(info)){
    segments(x0 = bplot[a], x1 = bplot[a], y0 = info[a,1]-2*info[a,2], y1 = info[a,1]+2*info[a,2],col = "black", lwd = 2) #2 if you want auc sd, 4 if you want sensitivity sd
  }
}
dev.off()



#compute a p value representing the significance in the difference between the original auc or sensitivity and the knock out version. 
TFplotnames <- c("Tbet Th1","Gata3 Th2","Stat3 Th17","Foxp3 Treg")
TFs <- c("Tbet","Gata3","Stat3","Foxp3")
type <- c("no KO","cell-state","H3K4me1","H3K27ac","DNase","ATAC") 
mean_sd_sens_stuff <- matrix(0,length(type),10)
mean_sd_auc_stuff <- matrix(0,length(type),10)
sens_pvalue <- matrix(0,length(type)-1,1)
auc_pvalue <- matrix(0,length(type)-1,1)
for (i in 1:length(TFs)){
  #original
  model <- read.table(paste0("tiffany/IMPACT_manuscript/Performance/Sens_AUC_collection_",TFs[i],".txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_sens_stuff[1,] <- as.numeric(model[1,])
  mean_sd_auc_stuff[1,] <- as.numeric(model[2,])
  #cell-state ko
  model <- read.table(paste0("tiffany/IMPACT_manuscript/KnockOuts/Sens_AUC_collection_",TFs[i],"_KOcellstatefeatures.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_sens_stuff[2,] <- as.numeric(model[1,])
  mean_sd_auc_stuff[2,] <- as.numeric(model[2,])
  #feature categories ko
  model <- read.table(paste0("tiffany/IMPACT_manuscript/KnockOuts/Sens_collection_",TFs[i],"_kofeatures.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_sens_stuff[3:length(type),] <- sapply(1:4,function(x) as.numeric(model[x,])) #only take "H3K4me1","H3K27ac","DNase","ATAC"
  
  model <- read.table(paste0("tiffany/IMPACT_manuscript/KnockOuts/AUC_collection_",TFs[i],"_kofeatures.txt"), quote = "", header = F, row.names = NULL, stringsAsFactors = FALSE)  
  mean_sd_auc_stuff[3:length(type),] <- sapply(1:4, function(x) as.numeric(model[x,]))
  
  for (ty in 1:(length(type)-1)){
    sens_pvalue[ty,1] <- t.test(mean_sd_sens_stuff[1,], mean_sd_sens_stuff[(ty+1),], alternative = "greater")$p.val
    auc_pvalue[ty,1] <- t.test(mean_sd_auc_stuff[1,], mean_sd_auc_stuff[(ty+1),], alternative = "greater")$p.val
  }

  assign(paste0("sens_pvalue_",TFs[i]),sens_pvalue)
  assign(paste0("auc_pvalue_",TFs[i]),auc_pvalue) 
}




