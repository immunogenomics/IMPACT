#Acquire matrix of 707 IMPACT annotaitons by 111 summary statistics 

library(gplots)
library(scales)
impact_sds <- read.table("IMPACT_Kawakami_ENCODE_sds.txt", stringsAsFactors = F, header = F, sep = "\t")
TF_cell_type_pairs_Kawa <- read.table("TF_cell_type_pairs_Kawa.txt", stringsAsFactors = F, header = F, sep = " ") 
mat_Kawa <- read.table("df_Kawa.txt", stringsAsFactors = F, header = T, sep = " ") 
IMPACT_Kawakami_nicenames <- read.table("IMPACT_Kawakami_nicenames.txt", stringsAsFactors = F, header = F, sep = "\t")

pop <- "EUR" #either run with EUR or EAS 
modelname <- "IMPACT_ENCODE_" 
filemodelname <- "IMPACT" 
key <- read.table(paste0("IMPACTxTRAIT_Kawakami",filemodelname,"_",pop,"_filekey.txt.gz"), stringsAsFactors = F, header = F, sep = "\t") #IMPACT
mydat <- read.table(paste0("IMPACTxTRAIT_Kawakami",filemodelname,"_",pop,".txt.gz"), stringsAsFactors = F, header = F, sep = "\t")

removethese <- c(3,4,5,6,31,48,75,96,120,124,132,134,139,181,257,258,290,301,719,727,728) #annotations that were QCed out 
removethese_eas <- c(3,4,5,6,31,48,75,96,132,134,139,181,257,258,290,301,719,727,728)

length(key) == nrow(mydat) #should be TRUE 

#load s-ldsc results
model <- sapply(1:length(key), function(x) strsplit(key[x], split = "[.]")[[1]][1])
trait <- sapply(1:length(key), function(x) strsplit(key[x], split = "[.]")[[1]][2])
#add more info about traits, when . separates from useful info 
if (pop == "EUR"){
  w <- which(trait == "UKB2_145K" | trait == "Asthma" | trait == "RA") #get more definition for these traits
  trait_addition <- sapply(1:length(w), function(x) strsplit(key[w[x]], split = "[.]")[[1]][3])
  fulltrait <- paste0(trait[w],".",trait_addition)
  trait[w] <- fulltrait
}else{
  w <- which(trait == "RA") #get more definition for these traits
  trait_addition <- sapply(1:length(w), function(x) strsplit(key[w[x]], split = "[.]")[[1]][3])
  fulltrait <- paste0(trait[w],".",trait_addition)
  trait[w] <- fulltrait
}
length(unique(model)) #707 impact
length(unique(trait)) #traits
unique_models <- unique(model)
unique_traits <- unique(trait)

#trait heritability 
mydat_h2 <- read.table(paste0("IMPACTxTRAIT_Kawakami",filemodelname,"_",pop,"_h2.txt.gz"), stringsAsFactors = F, header = F, sep = "\t")
myh2key_h2 <- read.table(paste0("IMPACTxTRAIT_Kawakami",filemodelname,"_",pop,"_h2_key.txt.gz"), stringsAsFactors = F, header = F, sep = "\t")
nrow(mydat_h2) == length(key)
if (modelname == "IMPACT_ENCODE_" & pop == "EUR"){
  length(removethese)*length(unique_traits) == nrow(myh2key_h2) - nrow(mydat_h2)
  myh2key_h2_annotindex <- sapply(1:nrow(myh2key_h2), function(x) strsplit(myh2key_h2$V1[x], split = ",")[[1]][1])
  m <- match(as.numeric(myh2key_h2_annotindex), removethese)
  myh2key_h2 <- myh2key_h2[which(is.na(m)),1]
  mm <- match("2,Asthma.ukbb",myh2key_h2)
  myh2key_h2 <- myh2key_h2[-mm]
  length(myh2key_h2) == nrow(mydat_h2)
  length(myh2key_h2) == length(key)
  length(myh2key_h2) == nrow(mydat)
}
if (modelname == "IMPACT_ENCODE_" & pop == "EAS"){
  
  s <- sapply(1:length(unique_models), function(x) strsplit(unique_models[x],"_")[[1]][3])
  m <- match(1:728,as.numeric(s))
  removethese_eas <- which(is.na(m))
  
  length(removethese_eas)*length(unique_traits) == nrow(myh2key_h2) - nrow(mydat_h2)
  
  myh2key_h2_annotindex <- sapply(1:nrow(myh2key_h2), function(x) strsplit(myh2key_h2$V1[x], split = ",")[[1]][1])
  m <- match(as.numeric(myh2key_h2_annotindex), removethese_eas)
  
  myh2key_h2 <- myh2key_h2[which(is.na(m)),1]
  length(myh2key_h2) == nrow(mydat_h2)
  length(myh2key_h2) == length(key)
  length(myh2key_h2) == nrow(mydat)
}

#tau_star
proportions <- mydat[match(unique_models, model),2] #find one occurrence of each of the models, take the prop from the output 
if(modelname == "IMPACT_ENCODE_"){
  hyphenparse <- 3
}else{hyphenparse <- 2}
unique_model_values <- as.numeric(sapply(1:length(unique_models), function(x) strsplit(unique_models[x], split = "_")[[1]][hyphenparse]) )
#a <- cbind(proportions,unique_model_values)
#a1 <- a[order(a[,2],decreasing = F),]
proportions_ordered <- proportions[match(1:728, unique_model_values)] #want this to be 728, because NAs are ok to have 
head(proportions_ordered)
tail(proportions_ordered)

master_scaled_IMPACT_Kawakami <- as.data.frame(matrix(0, ncol=3+length(unique_traits)*8, nrow=728))
master_scaled_IMPACT_Kawakami[,1] <- TF_cell_type_pairs_Kawa$V1
master_scaled_IMPACT_Kawakami[,2] <- 1:728
master_scaled_IMPACT_Kawakami[,3] <- proportions_ordered

#annotation results; here's where we reorder from 1:728 
for (i in 1:728){
  w <- which(model == paste0(modelname,i)) 
  if (length(w) > 0){
    counter <- 4 #how many columns need to skip before putting enrichment / tau output
    for (j in 1:length(unique_traits)){
      master_scaled_IMPACT_Kawakami[i,counter:(counter+7)] <- mydat[w[j],c(3:10)] #how many values - 1
      counter <- counter + 8 #how many values adding from ldsc output 
    }
  }else{
    master_scaled_IMPACT_Kawakami[i,4:ncol(master_scaled_IMPACT_Kawakami)] <- NA
  }
}

traitcolnames <- c()
for (i in 1:length(unique_traits)){
  traitcolnames <- c(traitcolnames, paste0(c("h2Explained:","h2se:","Enrichment:","EnrichmentSE:","EnrichmentP:","TauCoef:","TauCoefSE:","TauZ:"), unique(trait)[i]))
}
colnames(master_scaled_IMPACT_Kawakami) <- c("Model","Index","ChipAnnotSize", traitcolnames)

#caution: this is tau, not tau star
g <- grep("TauCoef:", names(master_scaled_IMPACT_Kawakami), value = TRUE)
mytau <- master_scaled_IMPACT_Kawakami[,g]
print(dim(mytau)) #728 x length(unique_traits), need to transform each of these values

g <- grep("TauCoefSE:", names(master_scaled_IMPACT_Kawakami), value = TRUE)
mytause <- master_scaled_IMPACT_Kawakami[,g]
print(dim(mytause)) #728 x length(unique_traits), need to transform each of these values

g <- grep("Enrichment:", names(master_scaled_IMPACT_Kawakami), value = TRUE)
myenr <- master_scaled_IMPACT_Kawakami[,g]
g <- grep("EnrichmentSE:", names(master_scaled_IMPACT_Kawakami), value = TRUE)
myenrse <- master_scaled_IMPACT_Kawakami[,g]
g <- grep("EnrichmentP:", names(master_scaled_IMPACT_Kawakami), value = TRUE)
myenrp <- master_scaled_IMPACT_Kawakami[,g]

s <- sapply(1:nrow(mydat_h2), function(x) strsplit(mydat_h2$V1[x], split = ": "))
h2 <- sapply(1:length(s), function(x) as.numeric(strsplit(s[[x]][2], split = " ")[[1]][1]))
h2_se <- c()
for (i in 1:length(s)){
  p1 <- strsplit(s[[i]][2], split = "[(]")[[1]][2]
  p2 <- strsplit(p1, split = "[)]")[[1]][1]
  h2_se[i] <- p2
}

h2_matrix <- matrix(0,728,length(unique_traits))
count <- 1
for (i in 1:728){
  if(modelname == "IMPACT_ENCODE_" & pop == "EUR"){rm1 <- removethese}
  if(modelname == "IMPACT_ENCODE_" & pop == "EAS"){rm1 <- removethese_eas}
  
  if (i %in% rm1){
    h2_matrix[i,] <- NA #and we don't advance the count 
  }else{
    h2_matrix[i,] <- h2[count:(length(unique_traits)+count-1)] #i = 1, 1:728*1
    count <- count + length(unique_traits)
  }
}

h2_se_matrix <- matrix(0,728,length(unique_traits))
count <- 1
for (i in 1:728){
  if(modelname == "IMPACT_ENCODE_" & pop == "EUR"){rm1 <- removethese}
  if(modelname == "IMPACT_ENCODE_" & pop == "EAS"){rm1 <- removethese_eas}
  
  if (i %in% rm1){
    h2_se_matrix[i,] <- NA #and we don't advance the count 
  }else{
    h2_se_matrix[i,] <- h2_se[count:(length(unique_traits)+count-1)] #i = 1, 1:728*1
    count <- count + length(unique_traits)
  }
}

if (pop == "EAS"){M <- 5469053}else{M <- 5961159}
mytau_star <- mytau
for (i in 1:nrow(mytau)){
  for (j in 1:ncol(mytau)){
    if (filemodelname == "ChIP"){
      mysd <- sqrt(proportions_ordered[i]*(1 - proportions_ordered[i]))
    }else{
      mysd <- impact_sds$V1[match(i, impact_sds$V1) + 1]
    }
    mytau_star[i,j] <- (M * mysd) / h2_matrix[1,j] * mytau[i,j]
  }
}

mytau_star_se <- mytause
for (i in 1:nrow(mytause)){
  for (j in 1:ncol(mytause)){
    if (filemodelname == "ChIP"){
      mysd <- sqrt(proportions_ordered[i]*(1 - proportions_ordered[i]))
    }else{
      mysd <- impact_sds$V1[match(i, impact_sds$V1) + 1]
    }
    mytau_star_se[i,j] <- (M * mysd) / h2_matrix[1,j] * mytause[i,j]
  }
}

mytau_star_logp <- mytau_star
for (i in 1:nrow(mytause)){
  for (j in 1:ncol(mytause)){
    if (is.na(mytau_star[i,j])){
      mytau_star_logp[i,j] <- NA
    }else{
      if (mytau_star[i,j] < 0){
        mytau_star_logp[i,j] <- log(x = pnorm(q = 0,mean = mytau_star[i,j], sd = mytau_star_se[i,j], lower.tail = F), base = 10) 
      }else{
        mytau_star_logp[i,j] <- -log(x = pnorm(q = 0,mean = mytau_star[i,j], sd = mytau_star_se[i,j], lower.tail = T), base = 10)
      }
    }
  }
}

#there are more sumstats in the IMPACTxTRAIT_KawakamiIMPACT* files than used in manuscript, e.g. dups of the same phenotype; remove these (those with less power)
if (pop == "EUR"){
    toremove <- c(2,28,30,36,40,46,59,51,56,61,54,64)
    unique_traits <- unique_traits[-toremove] 
    assign("traits_EUR",unique_traits)
    mytau_star_logp <- mytau_star_logp[,-toremove]
    h2_matrix <- h2_matrix[,-toremove]
    h2_se_matrix <- h2_se_matrix[,-toremove]
    assign("h2_matrix_EUR",h2_matrix)
    assign("h2_se_matrix_EUR",h2_se_matrix)
    mytau_star_se <- mytau_star_se[,-toremove]
    mytau_star <- mytau_star[,-toremove]
    myenr <- myenr[,-toremove]
    myenrse <- myenrse[,-toremove]
    myenrp <- myenrp[,-toremove]
}else{
  toremove <- c(2,13,30,39,42,45)
  unique_traits <- unique_traits[-toremove] 
  unique_traits[which(unique_traits == "RA.Kanai")] <- "RA"
  assign("traits_EAS",unique_traits)
  mytau_star_logp <- mytau_star_logp[,-toremove]
  h2_matrix <- h2_matrix[,-toremove]
  h2_se_matrix <- h2_se_matrix[,-toremove]
  assign("h2_matrix_EAS",h2_matrix)
  assign("h2_se_matrix_EAS",h2_se_matrix)
  mytau_star_se <- mytau_star_se[,-toremove]
  mytau_star <- mytau_star[,-toremove]
  myenr <- myenr[,-toremove]
  myenrse <- myenrse[,-toremove]
  myenrp <- myenrp[,-toremove]
}

#apply 5% FDR
mytau_star_logp_na <- mytau_star_logp
mytau_star_logp_na[is.na(mytau_star_logp_na)] <- 0 #all tau*
assign(paste0("mytau_star_logp_na_",pop,"_",filemodelname),mytau_star_logp_na) #store this information for later 

mytau_star_se[is.na(mytau_star_se)] <- 0
mytau_star[is.na(mytau_star)] <- 0
assign(paste0("mytau_star_se_",pop,"_",filemodelname),mytau_star_se) #store this information for later
assign(paste0("mytau_star_",pop,"_",filemodelname),mytau_star) #store this information for later 

assign(paste0("myenrse_",pop,"_",filemodelname),myenrse)
assign(paste0("myenr_",pop,"_",filemodelname),myenr)
assign(paste0("myenrp_",pop,"_",filemodelname),myenrp)

a <- cbind(mytau_star_EUR_IMPACT, mytau_star_EAS_IMPACT)
b <- cbind(mytau_star_se_EUR_IMPACT, mytau_star_se_EAS_IMPACT)

mytau_star_p <- a
for (i in 1:nrow(a)){
  for (j in 1:ncol(a)){
    if (is.na(a[i,j])){
      mytau_star_p[i,j] <- NA
    }else{
      if (a[i,j] < 0){
        mytau_star_p[i,j] <- pnorm(q = 0,mean = a[i,j], sd = b[i,j], lower.tail = F) #want this to be signed so remove the - part 
      }else{
        mytau_star_p[i,j] <- pnorm(q = 0,mean = a[i,j], sd = b[i,j], lower.tail = T)
      }
    }
  }
}

mytau_star_p_num <- as.numeric(unlist(mytau_star_p))
mytau_star_p_num_fdr <- p.adjust(mytau_star_p_num, method = "fdr", n = length(mytau_star_p_num))
mytau_star_p_fdr <- matrix(mytau_star_p_num_fdr,nrow = nrow(a), ncol=ncol(a))
length(which(mytau_star_p_fdr < 0.05)) #7843 associations 
mytau_star_p_EUR_IMPACT <- mytau_star_p_fdr[,1:length(traits_EUR)]
mytau_star_p_EAS_IMPACT <- mytau_star_p_fdr[,((length(traits_EUR)+1):ncol(mytau_star_p_fdr))]


