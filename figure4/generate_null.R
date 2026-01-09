#!/usr/bin/env Rscript 

setwd("/scratch/midway2/rozennpineau/drought/generate_null/")
daily_af <- read.table("daily_af_43loci.txt")
fq_full <- read.table("/scratch/midway2/rozennpineau/drought/two_pulse_flexible_prop_2/anc_call_frequency_site1-786261.txt")
fq <- fq_full[-c(1,2),]

#how many sites per bin for the data?
max(daily_af[1,]) - min(daily_af[1,])
bins <- hist(as.numeric(daily_af[1,]), breaks = 100, xlim = c(min(daily_af[1,]),max(daily_af[1,])))
iteration <- which(bins$counts > 0) #31 different starting frequencies
n_sites <- bins$counts[bins$counts > 0]

#Keep sites that have the same starting frequency range
fq_fil <- fq[ , which(fq[1,] >= min(daily_af[1,]) & fq[1,] <= max(daily_af[1,]))] 

af_change_list <- vector("list", length(n_sites))

#make pools of same range bins in a list
for (i in 1:length(iteration)) {
  
  idx <- which(bins$breaks[iteration[i]] <= fq_fil[1,] & fq_fil[1,] <= bins$breaks[iteration[i]+1])
  print(length(idx))
  if (length(idx)>0) { af_change_list[[i]] <- t(fq_fil[20,idx] - fq_fil[1,idx])}

  }
#export list
sink('af_change_list_43.txt')
print(af_change_list)
sink()


#
#L <- 1000 # number of permutations
#mean_af_change <- NA
#write.table(mean_af_change, "mean_af_change_null.txt", sep = "\t", row.names = F, col.names = F)
#for (i in 1:L) {
  
  #sample X site from each bin
 # mean_af_change <- mean(af_change_list[[i]][sample(length(af_change_list[[i]]), n_sites[i])]) #see if we have a bug here with not enough samples to sample randomly 
 # write.table(mean_af_change, "mean_af_change_null.txt", sep = "\t", row.names = F, col.names = F, append = T)
#}
