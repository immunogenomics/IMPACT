#Compute Mean and St Dev Genome-wide, from per chromosome values 
sds <- read.table("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/StDevChrWidePredictionsperTF.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
means <- read.table("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/MeanChrWidePredictionsperTF.txt", sep = "\t", header = F, stringsAsFactors = FALSE)
sizes <- read.table("/tiffany/IMPACT_manuscript/GenomeWide/GenomeTracks/hg19.sizes", sep = "\t", header = F, stringsAsFactors = FALSE)

#pooled s = sq rt of ([sum from 1 to k of (nk - 1)*sk^2 + nk(ybark - ybarall)^2 ]/ [(sum from 1 to k of nk) - 1]
#https://stats.stackexchange.com/questions/55999/is-it-possible-to-find-the-combined-standard-deviation

#do this for T-bet. sds col 1 means col 1.
for (tfi in 1:4){
  
tf_index <- tfi
#pool means first
meansum <- 0
for (i in 1:24){
  meansum <- meansum + sizes$V2[i]*means[i,tf_index]
}
new_mean <- meansum/sum(as.numeric(sizes$V2))

#step 1. 
#sum from 1 to k (chromosomes) of (nk - 1)*sk^2 + nk(ybark - ybarall)^2
numer <- 0
denom <- 0
for (i in 1:24){
  numer <- numer + (sizes$V2[i]-1)*(sds[i,tf_index])^2 + sizes$V2[i]*(means[i,tf_index] - new_mean)^2
  denom <- denom + sizes$V2[i]
}
new_sd <- sqrt(numer/(denom-1))

print(new_mean)
print(new_sd)

}
