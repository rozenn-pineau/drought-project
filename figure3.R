# Plot results from the gwas analysis done on ancestry mapped sites and day to full wilt with correction for population structure
rm(list= ls())

library(data.table)
library(dplyr)
library(plyr)
library(tidyverse)
library(ggplot2)
library(extrafont)
library(reshape2)
library(cowplot)
library(colourvalues)

# ----------------------------------------------------------------------------- #
# F3A
# ----------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------- #
# read in GWAS probabilities
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/8.ancestry_hmm/gwas/gemma_cutoff09_miss02/")
lmm_gemma_anc <- read.table("ancestry_corrected_inflated_gemma_gwas.assoc.txt", header= T, sep = "\t") 
#modify rs col to get chrom position information
lmm_gemma_anc[c("chrom", "pos")] <- as.numeric(str_split_fixed(lmm_gemma_anc$rs, ':', 2))

#length(which(lmm_gemma_anc$FDR < 0.05)) 893 --> correct file

#read in clumped lead SNPs
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/8.ancestry_hmm/1.results/trajectories/")
genotype <- read.table("drought_adapted_43clumps_GT.txt")
lead_snps <- genotype[,c(1,2)]



# ----------------------------------------------------------------------------- #
# Manhattan plot
# ----------------------------------------------------------------------------- #
idx2 <- which(lmm_gemma_anc$FDR == max(lmm_gemma_anc$FDR[which(lmm_gemma_anc$FDR < 0.05)]))

#colors
lmm_gemma_anc$cols <- alpha("gray", 0.5)
lmm_gemma_anc$cols_stroke <- "gray"
lmm_gemma_anc$cols_additional <- NA
#lmm_gemma_anc$cols[lmm_gemma_anc$FDR < 0.05] <- alpha("#fc9272",0.3) #893

#assign color to each + make outer circle black for white points
for (i in 1:length(lead_snps$V1)) {
  lmm_gemma_anc$cols[lmm_gemma_anc$chrom == lead_snps$V1[i] & lmm_gemma_anc$pos==lead_snps$V2[i]] <- "red" #colvec[i]
  lmm_gemma_anc$cols_additional[lmm_gemma_anc$chrom == lead_snps$V1[i] & lmm_gemma_anc$pos==lead_snps$V2[i]] <- "red"
  lmm_gemma_anc$cols_stroke[lmm_gemma_anc$chrom == lead_snps$V1[i] & lmm_gemma_anc$pos==lead_snps$V2[i]]  <- "black"
}




p1 <- lmm_gemma_anc %>% #filter(chrom==1) %>% 
  ggplot(aes(pos/1000000,-log10(p_wald))) +
  geom_point(size = 2, shape = 21, colour = lmm_gemma_anc$cols_stroke, fill = lmm_gemma_anc$cols, stroke = 0.01) +
  geom_point(size = 2, shape = 16, colour = lmm_gemma_anc$cols_additional) +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Genomic position (Mb)") + ylab("-log10(p-value)") +
  ylim(0,10) +
  theme(axis.text.x = element_text(size = 8, angle = 90),
        text=element_text(size=16,  family="Times New Roman")) +
  #geom_hline(yintercept = -log10(lmm_gemma_anc$p_wald[idx]),lty="dashed",lwd=.5, col = "#929292") +
  geom_hline(yintercept = -log10(6.5830432740603e-05),lty="dashed",lwd=.5, col = "red") 

# ----------------------------------------------------------------------------- #
# F3B
# ----------------------------------------------------------------------------- #

#goal : to plot the ancestry per site per chromosome per individual

#plot colors
tub_col <- "#76528BFF"
rud_col <- "#CBCE91FF"
het_col <- "#44A3BB"

# values is a matrix with the ancestry information for each site (rows) and each individual (columns)
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/8.ancestry_hmm/1.results/")
values <- read.table("two_pulse_flexible_prop_2_values.txt", sep = "\t", header = T)
#rm outliers
values <- values[,-c(which(colnames(values) == "P16_Nat_1_T" | colnames(values) == "P12_Nat_14_T"))] #282 samples

#transform file from large to long format
df_long <- melt(values, id.vars = c( "chrom", "position"))

#upload df
#df_long <- melt(values[1:100,], id.vars = c( "chrom", "position")) #when testing
colnames(df_long) <- c("chrom", "pos", "variable", "GT")

#upload longitude info
long <- read.table("k2_structure_res.csv", sep = ",", header = T)
long_info <- data.frame(variable = long$samp, long = long$long)
#remove the two outliers
long_info <- long_info[-c(which(long_info$variable == "P16_Nat_1_T" | long_info$variable == "P12_Nat_14_T")),] #280 samples

#subsample (more than 196 000 000 lines, the figure stops rendering at chromosome 3 if I do not subsample)
df_long <- df_long[sort(sample(dim(df_long)[1], 500000, replace = F )),]

#merge 
dataset <- merge(df_long, long_info, by = "variable" )
#sort samples in values based on long
dataset <- as.data.frame(dataset[order(dataset$long),]) # 280 77011
dataset$GT <- as.character(as.numeric(dataset$GT))
dataset$pos <- as.numeric(dataset$pos)


p2 <- dataset %>% 
  #ggplot(aes(pos, variable, color=GT, fill=GT)) + x = reorder(category, -value), y = value
  ggplot(aes(x = pos, y = reorder(variable, -long), color=GT, fill=GT)) +
  geom_tile() +
  facet_grid(~chrom, scales = "free_x",space = "free_x") +
  scale_colour_manual(values = c(tub_col, het_col, rud_col)) +
  scale_fill_manual(values = c(tub_col, het_col, rud_col)) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),  # Adjust size to your preference
    axis.text.x = element_blank(), #element_text(size = 6,angle = 90),
    text=element_text(size=16,  family="Times New Roman")) + #,
    #panel.spacing = unit(0, "lines"),           # No space between facets
    #panel.border = element_rect(color = "black", fill = NA, size = .5)) + # Black outline
  labs(fill="Ancestry", y="Individual", x="") +
  guides(color="none")


# ----------------------------------------------------------------------------- #
# plot
# ----------------------------------------------------------------------------- #

pdf(file="figure3.pdf", family = "Times New Roman",
    bg = "transparent", width=10, height=5)


plot_grid(p2,p1, ncol = 1, nrow = 2, align = "hv", rel_heights = c(1,2))


dev.off()