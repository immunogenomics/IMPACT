
commontraits_EAS <- c("Baso","BMI","DBP","Eosino","Hb","Height","Ht","Lym","MCHC","MCH","MCV","Mono","Neutro","Plt","RBC","SBP","WBC","Arrhythmia","Asthma","BrCa","CAD","Cataract","CHF","Glaucoma","IS","Osteoporosis","PrCa","RA","T2D")
commontraits_EUR <- c("Basophil","BMI","DBP","Eosinophil","Hb","Height","Ht","Lymphocyte","MCHC","MCH","MCV","Monocyte","Neutrophil","Platelet","UKB2_145K.blood_RED_COUNT","SBP","WBC","AF","Asthma.ukbb","BrCa","CAD","Cataract","CHF","Glaucoma","IS","Osteoporosis","PrCa","PASS_Rheumatoid_Arthritis","T2D")

m_eas <- match(commontraits_EAS, traits_EAS) #get traits_EAS and traits_EUR from 2_Get_tau_star_association_matrix_707by111.R script 
m_eur <- match(commontraits_EUR, traits_EUR)

y <- fread("IMPACT_annotations_PearsonRCorr.txt.gz", header = F)
y <- as.matrix(y)
all(diag(y) == 1)
w <- which(is.na(y), arr.ind = T)
removethese <- w[1:21,1] #these are the same as in other scripts, annotations QCed out- not enough training data. 
kept_these <- c(1:728)[-removethese]
y1 <- y[-removethese,-removethese]
y2 <- y1^2

for (i in 1:length(commontraits_EAS)){
  taustar_disease <- mytau_star_EUR_IMPACT[,m_eur[i]] #relative to EUR 
  taustar_disease <- taustar_disease[-removethese] #707 
  #need 707 by 707 corr structure: go to section titled: pairwise correlation between 728 impact annotations
  
  
  annotsleft <- length(taustar_disease)
  taustar_disease_manipulate <- taustar_disease
  y2_manipulate <- y2 
  Correlation_Clumping <- data.frame()
  count <- 1
  while(annotsleft > 0){
    mymax <- max(taustar_disease_manipulate, na.rm = T)
    print(mymax)
    w_lead_annot <- which(taustar_disease_manipulate == mymax)
    w_corr_annots <- which(y2_manipulate[w_lead_annot,] > 0.5) #includes the lead 
    Correlation_Clumping[count,1] <- kept_these[w_lead_annot]
    Correlation_Clumping[count,2] <- mymax
    Correlation_Clumping[count,3] <- paste0(kept_these[w_corr_annots],collapse = ",")
    count <- count + 1 
    y2_manipulate[w_corr_annots,] <- NA
    y2_manipulate[,w_corr_annots] <- NA
    taustar_disease_manipulate[w_corr_annots] <- NA
    
    annotsleft <- length(which(!is.na(taustar_disease_manipulate)))
  }
  colnames(Correlation_Clumping) <- c("Lead","taustar","CorrAnnots")
  assign(paste0("Correlation_Clumping_",commontraits_EAS[i]),Correlation_Clumping)
}