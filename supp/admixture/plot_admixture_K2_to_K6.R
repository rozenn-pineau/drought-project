# Plot PCA res from drought selection experiment data. 
rm(list= ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(data.table)
library(extrafont)
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/2.admixture/admixture_out")

# ---------------------------------------------------------------------------- #
# Load data
# ---------------------------------------------------------------------------- #

# upload sample info
samplelist <- read.table("../samps.txt", header=F)
sample_order <- samplelist$V1 # keep only the first column

# read in admixture results
k1 <- read.table("merged_numericChr_LDpruned.1.Q", sep = " ")
k2 <- read.table("merged_numericChr_LDpruned_mod.2.Q", sep = " ")
k3 <- read.table("merged_numericChr_LDpruned_mod.3.Q", sep = " ")
k4 <- read.table("merged_numericChr_LDpruned_mod.4.Q", sep = " ")
k5 <- read.table("merged_numericChr_LDpruned_mod.5.Q", sep = " ")
k6 <- read.table("merged_numericChr_maf05_LDpruned.6.Q", sep = " ")
k6_samp_order <- read.table("merged_numericChr_maf05_LDpruned.fam", sep = " ")

#remove outliers from k1 and k6
k6 <- k6[-c(which(k6_samp_order$V1  %in%   c("P16_Nat_1_T",  "P12_Nat_14_T") )),]
dt <- data.frame(samp = sample_order,
                 k2 = k2,
                 k3 = k3,
                 k4 = k4,
                 k5 = k5, 
                 k6 = k6)

summary(dt)
#colnames(sample_order) <- "samp"
#write.table(samplelist, "sample_info.csv", sep = ",")
sample_info <- read.table("../sample_info.csv", header=T, sep = ",")

# ---------------------------------------------------------------------------- #
# Prep data
# ---------------------------------------------------------------------------- #

# order sample_info based on longitude, then order sample_order (and admixture data) based on sample_info
sample_info_ordered <- sample_info[order(sample_info$long),]

# join datasets together
fulldt <- inner_join(sample_info_ordered,dt, by=c("samp"))

summary(fulldt)


tub_col <- "#76528BFF"
rud_col <- "#CBCE91FF"




# ---------------------------------------------------------------------------- #
# Admixture plot
# ---------------------------------------------------------------------------- #

# full admixture plots (k=1 to k=10)
pdf(file="/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/2.admixture/figure2_admx.pdf", 
    bg = "transparent", width=12, height=4, family = "Times New Roman")

par(family = "Times New Roman")

barplot(t(cbind(fulldt$k2.V1, fulldt$k2.V2)), col = c(tub_col, rud_col), space = 0, border = c("#929292"),
        ylab = "% ancestry", main = "")


dev.off()


#k=2 to k=6
pal <- colour_values(1:6, palette = "viridis")
pdf(file="/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/2.admixture/k1-6_admx.pdf", 
    bg = "transparent", width=8, height=8, family = "Times New Roman")

par(mfrow = c(4,1), family ="Times New Roman", mar = c(2, 4, 1, 1), cex.lab = 2)

barplot(t(cbind(fulldt$k3.V1, fulldt$k3.V2, fulldt$k3.V3)), col = pal[1:3], space = 0, border = NA,
        ylab = "% ancestry", main = "K=3")

barplot(t(cbind(fulldt$k4.V1, fulldt$k4.V2, fulldt$k4.V3, fulldt$k4.V4)),  col = pal[1:4], space = 0, border =NA,
        ylab = "% ancestry", main = "K=4")

barplot(t(cbind(fulldt$k5.V1, fulldt$k5.V2, fulldt$k5.V3, fulldt$k5.V4, fulldt$k5.V5)),  col = pal[1:5], space = 0, border =NA,
        ylab = "% ancestry", main = "K=5")

barplot(t(cbind(fulldt$k6.V1, fulldt$k6.V2, fulldt$k6.V3, fulldt$k6.V4, fulldt$k6.V5, fulldt$k6.V6)),  col =pal, space = 0, border = NA,
        ylab = "% ancestry", main = "K=6")

dev.off()

# ---------------------------------------------------------------------------- #
# CV plots
# ---------------------------------------------------------------------------- #

# Plot cross validation error (find which model is "best")
pdf(file="/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/2.admixture/CV.pdf", 
    bg = "transparent", width=4, height=4, family = "Times New Roman")

cv <- fread("cv.txt")
plot(c(2:6), cv$V4[c(1:4,6)], type = "b", pch = 16, 
     xlim = c(0,7), ylim = c(0, 0.4), 
     cex = 1.5, col = "#525252",
     ylab = "Cross Validation error", 
     xlab = "K")

dev.off()

# ---------------------------------------------------------------------------- #
#Extract individuals that are pure var rudis, pure var tuberculatus, and hybrid
# ---------------------------------------------------------------------------- #
#dataset is ordered by longitude: lower longitudes to higher longitudes, so West to East
#so we expect var rudis to be to the West (lower longitudes)
(mean(fulldt$k2.V2)) #%of var.rudis 0.67
(mean(fulldt$k2.V1)) #%of var. tuberculatus 0.33

keep_rudis <- which(fulldt$k2.V2 > .99) #66
keep_tub <- which(fulldt$k2.V1 > .99) #32
keep_hyb <- which(fulldt$k2.V2 > .4 & fulldt$k2.V1 > 0.4) #19
#not the same number of individuals for each, so randomly select 19 from var rudis and var tub 
keep_rudis_19 <- sample(keep_rudis, 19)
keep_tub_19 <- sample(keep_tub, 19)

#check
fulldt$k2.V1[keep_tub] #ok
fulldt$k2.V2[keep_rudis] #ok
fulldt$k2.V1[keep_hyb] #ok

#export individual IDs for each
info_export <- data.frame( samp_name  = fulldt$samp[c(keep_rudis_19, keep_tub_19, keep_hyb)], 
                           admxtr  = fulldt$k2.V1[c(keep_rudis_19, keep_tub_19, keep_hyb)], 
                           ancestry = c(rep("var_rudis", 19), rep("var_tuberculatus", 19), rep("hybrid", 19)) )

write.table(info_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/6.trajectories/samp_names_to_compare_ancestry.txt",
      sep = "\t", col.names = T, row.names = F) 










