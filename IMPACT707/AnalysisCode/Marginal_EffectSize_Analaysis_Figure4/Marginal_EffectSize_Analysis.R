# Marginal beta effect size, Fst, heterogeneity concordance 

traitnames <- c("Basophil","BMI","Eosinophil","Hb","Height","Ht","Lymphocyte","MCHC","MCH","MCV","Monocyte","Neutrophil","Platelet","RBC","SBP","WBC","Asthma.ukbb","CAD","RA","PrCa","T2D")  
eas_traitnames <- c("Baso","BMI","Eosino","Hb","Height","Ht","Lym","MCHC","MCH","MCV","Mono","Neutro","Plt","RBC","SBP","WBC","Asthma","CAD","RA","PrCa","T2D")  
eur_traitnames <- c('Basophil','BMI','Eosinophil','Hb','Height','Ht','Lymphocyte','MCHC','MCH','MCV','Monocyte','Neutrophil','Platelet','RBC','SBP','WBC','Asthma.ukbb','CAD','PASS_Rheumatoid_Arthritis','PrCa','T2D.gwascatalog')

for (tr in 1:length(traitnames)){
  ta_bottom <- fread(paste0(traitnames[tr],"_Loci_EUREAS_manyIMPACTs_BETA_bottom95_May13.txt.gz"))
  ta_top <- fread(paste0(traitnames[tr],"_Loci_EUREAS_manyIMPACTs_BETA_top5_May13.txt.gz"))
  ta_ordinary <- fread(paste0(traitnames[tr],"_Loci_EUREAS_manyIMPACTs_BETA_ordinary_May13.txt.gz"))
  
  assign(paste0("ta_ordinary_",eur_traitnames[tr]),ta_ordinary)
  assign(paste0("ta_top_",eur_traitnames[tr]),ta_top)
  assign(paste0("ta_bottom_",eur_traitnames[tr]),ta_bottom)
}

pth <- seq(0,8,0.5) #1:8
pth <- 10^(-pth)
pth_lowers <- c(0,rev(pth)[-length(pth)])
pth_uppers <- rev(pth)

myMAmat_BetaCorr <- matrix(0,length(traitnames),3*length(pth_uppers))
myMAmat_se_BetaCorr <- matrix(0,length(traitnames),3*length(pth_uppers))

myMAmat_2pqCorr <- matrix(0,length(traitnames),3*length(pth_uppers))
myMAmat_se_2pqCorr <- matrix(0,length(traitnames),3*length(pth_uppers))

myMAmat_Fst <- matrix(0,length(traitnames),3*length(pth_uppers))
myMAmat_se_Fst <- matrix(0,length(traitnames),3*length(pth_uppers))

min_val <- 3 #minimum number of significant loci required for correlation

for (tr in 1:length(traitnames)){
  numSNPs <- matrix(0,length(pth_uppers),3)
  ta_ordinary <- get(paste0("ta_ordinary_",eur_traitnames[tr]))
  ta_top <- get(paste0("ta_top_",eur_traitnames[tr]))
  ta_bottom <- get(paste0("ta_bottom_",eur_traitnames[tr]))
  
  #remove all the rows with no beta; because this is how we eliminated the MHC 
  w_remove <- which(is.na(ta_ordinary$Beta_EUR) | is.na(ta_ordinary$Beta_EAS))
  if(length(w_remove) > 0){ta_ordinary <- ta_ordinary[-w_remove,]}
  w_remove <- which(is.na(ta_top$Beta_EUR) | is.na(ta_top$Beta_EAS))
  if(length(w_remove) > 0){ta_top <- ta_top[-w_remove,]}
  w_remove <- which(is.na(ta_bottom$Beta_EUR) | is.na(ta_bottom$Beta_EAS))
  if(length(w_remove) > 0){ta_bottom <- ta_bottom[-w_remove,]}
  
  cormat_BetaCorr <- matrix(0,length(pth_lowers),3)
  cormat_se_BetaCorr <- matrix(0,length(pth_lowers),3)
  
  cormat_2pqCorr <- matrix(0,length(pth_lowers),3)
  cormat_se_2pqCorr <- matrix(0,length(pth_lowers),3)
  
  cormat_Fst <- matrix(0,length(pth_lowers),3)
  cormat_se_Fst <- matrix(0,length(pth_lowers),3)
  
  for (j in 1:length(pth_uppers)){
  
    for (k in 1:3){ #ordinary/IMPACT top 5/IMPACT bottom 95 
      if(k == 1){ta <- ta_ordinary}
      if(k == 2){ta <- ta_top}
      if(k == 3){ta <- ta_bottom}
      ta_pth <- ta[(ta$P_EUR < pth_uppers[j]),]
      
      if(nrow(ta_pth)>min_val){
        
        w_noAlleleInfo <- which(ta_pth$EUR_A1 == "" | ta_pth$EUR_A2 == "" | ta_pth$EAS_A1 == "" | ta_pth$EAS_A2 == "")
        if(length(w_noAlleleInfo) > 0){
          ta_pth <- ta_pth[-w_noAlleleInfo,]
        }
        w_match <- which(paste0(ta_pth$EUR_A1,"_",ta_pth$EUR_A2) == paste0(ta_pth$EAS_A1,"_",ta_pth$EAS_A2))
        w_dontmatch <- which(!paste0(ta_pth$EUR_A1,"_",ta_pth$EUR_A2) == paste0(ta_pth$EAS_A1,"_",ta_pth$EAS_A2))
        w_moreComplicatedThanSwaps <- which(!paste0(ta_pth$EUR_A1[w_dontmatch],"_",ta_pth$EUR_A2[w_dontmatch]) == paste0(ta_pth$EAS_A2[w_dontmatch],"_",ta_pth$EAS_A1[w_dontmatch])) #which don't match after swap
        print(paste0("MoreComplicatedThanSwaps:",length(w_moreComplicatedThanSwaps)))
        if(length(w_moreComplicatedThanSwaps) > 0){
          ta_pth <- ta_pth[-w_dontmatch[w_moreComplicatedThanSwaps],]
          #update
          w_dontmatch <- which(!paste0(ta_pth$EUR_A1,"_",ta_pth$EUR_A2) == paste0(ta_pth$EAS_A1,"_",ta_pth$EAS_A2))
        }
        ta_pth$Beta_EAS[w_dontmatch] <- (-1)*ta_pth$Beta_EAS[w_dontmatch]
        
        if(nrow(ta_pth)>min_val){
          cormat_BetaCorr[j,k] <- cor(ta_pth$Beta_EUR, ta_pth$Beta_EAS)
          cormat_se_BetaCorr[j,k] <- (cor.test(ta_pth$Beta_EUR, ta_pth$Beta_EAS)$conf.int[2] - cormat_BetaCorr[j,k])/1.96
        }else{
          cormat_BetaCorr[j,k] <- NA
          cormat_se_BetaCorr[j,k] <- NA
        }
        
        #heterozygosity
        eur_af <- ta_pth$EUR_AltF
        eas_af <- ta_pth$EAS_AltF
        eur_het <- 2*eur_af*(1-eur_af)
        eas_het <- 2*eas_af*(1-eas_af)
        
        rem <- union(which(is.na(eur_het)), which(is.na(eas_het)))
        if(length(rem) > 0){
          eur_het <- eur_het[-rem]
          eas_het <- eas_het[-rem]
        }
        if(length(eur_het) < 4){
          cormat_2pqCorr[j,k] <- NA
          cormat_se_2pqCorr[j,k] <- NA
        }else{
          cormat_2pqCorr[j,k] <- cor(eur_het, eas_het)
          cormat_se_2pqCorr[j,k] <- (cor.test(eur_het, eas_het)$conf.int[2] - cormat_2pqCorr[j,k])/1.96
        }
        
        #fst
        eur_af <- ta_pth$EUR_AltF
        eas_af <- ta_pth$EAS_AltF
        mean_p <- sapply(1:length(eas_af), function(x) mean(c(eur_af[x],eas_af[x])))
        fst <- (eur_af - eas_af)^2 / (2*mean_p*(1-mean_p))
        cormat_Fst[j,k] <- mean(fst, na.rm = T)
        cormat_se_Fst[j,k] <- sd(fst, na.rm = T)
      }else{ 
        cormat_BetaCorr[j,k] <- NA
        cormat_2pqCorr[j,k] <- NA
        cormat_Fst[j,k] <- NA
        cormat_se_BetaCorr[j,k] <- NA
        cormat_se_2pqCorr[j,k] <- NA
        cormat_se_Fst[j,k] <- NA
      }
    }
  }
  myMAmat_BetaCorr[tr,] <- c(rev(cormat_BetaCorr[,1]),rev(cormat_BetaCorr[,2]),rev(cormat_BetaCorr[,3]))
  myMAmat_se_BetaCorr[tr,] <- c(rev(cormat_se_BetaCorr[,1]),rev(cormat_se_BetaCorr[,2]),rev(cormat_se_BetaCorr[,3]))
  myMAmat_2pqCorr[tr,] <- c(rev(cormat_2pqCorr[,1]),rev(cormat_2pqCorr[,2]),rev(cormat_2pqCorr[,3]))
  myMAmat_se_2pqCorr[tr,] <- c(rev(cormat_se_2pqCorr[,1]),rev(cormat_se_2pqCorr[,2]),rev(cormat_se_2pqCorr[,3]))
  myMAmat_Fst[tr,] <- c(rev(cormat_Fst[,1]),rev(cormat_Fst[,2]),rev(cormat_Fst[,3]))
  myMAmat_se_Fst[tr,] <- c(rev(cormat_se_Fst[,1]),rev(cormat_se_Fst[,2]),rev(cormat_se_Fst[,3]))
}

