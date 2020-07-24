pop <- "EUR" #or "EAS"
file_key <- fread(paste0("IMPACTxTRAIT_Top5p_",pop,"_filekey.txt.gz"),header = F)
file_key <- file_key$V1
top5p_h2_file <- fread(paste0("IMPACTxTRAIT_Top5p_",pop,".txt.gz"), header = F)
model <- sapply(1:length(file_key), function(x) strsplit(file_key[x], split = "[.]")[[1]][2])
trait <- sapply(1:length(file_key), function(x) strsplit(file_key[x], split = "[.]")[[1]][3])

if (pop == "EUR"){
  w <- which(trait == "UKB2_145K" | trait == "Asthma" | trait == "RA") #get more definition for these traits
  trait_addition <- sapply(1:length(w), function(x) strsplit(file_key[w[x]], split = "[.]")[[1]][4])
  fulltrait <- paste0(trait[w],".",trait_addition)
  trait[w] <- fulltrait
  unique_traits <- unique(trait)
  
  top5p_EUR_h2 <- matrix(0,728,length(unique_traits))
  top5p_EUR_h2se <- matrix(0,728,length(unique_traits))
  for (i in 1:728){
    w <- which(model == paste0("IMPACT_ENCODE_",i)) 
    if (length(w) > 0){
      for (j in 1:length(unique_traits)){
        top5p_EUR_h2[i,j] <- as.numeric(top5p_h2_file[w[j],3]) #prop h2
        top5p_EUR_h2se[i,j] <- as.numeric(top5p_h2_file[w[j],4]) #prop h2 se
      }
    }else{
      top5p_EUR_h2[i,] <- NA
      top5p_EUR_h2se[i,] <- NA
    }
  }
  toremove <- c(2,28,30,36,40,46,59,51,56,61,54,64) #removed some traits for QC
  top5p_EUR_h2 <- top5p_EUR_h2[,-toremove]
  top5p_EUR_h2se <- top5p_EUR_h2se[,-toremove]
  
}else{
  w <- which(trait == "RA") #get more definition for these traits
  trait_addition <- sapply(1:length(w), function(x) strsplit(file_key[w[x]], split = "[.]")[[1]][4])
  fulltrait <- paste0(trait[w],".",trait_addition)
  trait[w] <- fulltrait
  unique_traits <- unique(trait)
  
  top5p_EAS_h2 <- matrix(0,728,length(unique_traits))
  top5p_EAS_h2se <- matrix(0,728,length(unique_traits))
  for (i in 1:728){
    w <- which(model == paste0("IMPACT_ENCODE_",i)) 
    if (length(w) > 0){
      for (j in 1:length(unique_traits)){
        top5p_EAS_h2[i,j] <- as.numeric(top5p_h2_file[w[j],3]) #prop h2
        top5p_EAS_h2se[i,j] <- as.numeric(top5p_h2_file[w[j],4]) #prop h2 se 
      }
    }else{
      top5p_EAS_h2[i,] <- NA
      top5p_EAS_h2se[i,] <- NA
    }
  }
  toremove <- c(2,29,38,41,44) #removed some traits for QC
  top5p_EAS_h2 <- top5p_EAS_h2[,-toremove] 
  top5p_EAS_h2se <- top5p_EAS_h2se[,-toremove] 
}

