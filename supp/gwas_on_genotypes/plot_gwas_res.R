# Plot results from the gwas analysis
library(data.table)
library(ggplot2)
library(dplyr)
library(qqman)
library(tidyverse)
library(extrafont)

# Theoretical quantile function for -log10(p)
qlog10 <- function(p.values) {
  theoretical <- rank(p.values)/length(p.values)
  return(-log10(theoretical))
}

 
# ----------------------------------------------------------------------------- #
# read in GEMMA output
# ----------------------------------------------------------------------------- #

# gwas with kinship matrix (corrected for pop structure)
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/3.gwas/output")
list.files()
res <- fread("merged_numericChr_nooutliers_maf01_gwaswanc.assoc.txt",header=T) # one phenotype was tested - day to full wilt
res$FDR_corr <- p.adjust(res$p_lrt, method = "fdr") # Given a set of p-values, returns p-values adjusted using one of several methods
res$bon <- p.adjust(res$p_lrt,method = "bonferroni")
res$qlogp <- qlog10(res$p_lrt)
res <- res[order(-qlogp),]
summary(res$chr)

# uncorrected gwas
res_vanilla <- fread("merged_numericChr_nooutliers_maf01_gwasnoanc.assoc.txt",header=T) #no covariate GEMMA GWA results
res_vanilla$FDR_corr <- p.adjust(res$p_wald,method = "fdr")
res_vanilla$bon<-p.adjust(res$p_wald,method = "bonferroni")

res_vanilla$qlogp<-qlog10(res_vanilla$p_wald)
res_vanilla<-res_vanilla[order(-qlogp),]


# ----------------------------------------------------------------------------- #
# QQPlots
# ----------------------------------------------------------------------------- #
pdf(file=paste("qqplot_corrected.pdf", sep="" ), bg = "transparent", width=4, height=4, family = "Times New Roman")

# QQplot with base 10
ggplot(res, aes(x = qlogp, y = -log10(p_lrt))) + 
  geom_path() +
  geom_abline(intercept = 0, slope = 1)
# to test if there is enrichment in low p-values compared to the null, which are uniformly distributed under the null

dev.off()

# QQplot with base 10 - pwald test for uncorrected gwas
pdf(file=paste("qqplot_uncorrected.pdf", sep="" ), bg = "transparent", width=4, height=4, family = "Times New Roman")

ggplot(res_vanilla, aes(x = qlogp, y = -log10(p_wald))) + 
  geom_path() +
  geom_abline(intercept = 0, slope = 1)

dev.off()

## calculates lambda
 # The genomic inflation factor \(\lambda\) is defined as the ratio of the median 
# of the empirically observed distribution of the test statistic to the expected median,
# thus quantifying the extent of the bulk inflation and the excess false positive rate.

chisq <-qchisq(1-res$p_lrt,1)
lambda <- median(chisq) / qchisq(0.5, 1)
lambda

SE.median <- qchisq(0.975, 1) * (1.253 * ( sd(chisq) / sqrt( length(chisq) ) ) )
lambda + (SE.median / qchisq(0.5, 1))
lambda - (SE.median / qchisq(0.5, 1))

#plot

pdf(file=paste("gwas_uncorrected_subset.pdf", sep="" ), bg = "transparent", width=10, height=3, family = "Times New Roman")

#subsample
subres <- res_vanilla[sample(1:dim(res_vanilla)[1], 500000),]

p1<- subres %>%
  ggplot(aes(ps/1000000,-log10(p_wald), color=af)) +
  geom_point(alpha=.5) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,10) +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none") +
  geom_hline(yintercept = -log(0.05),lty="dashed",lwd=1.25,color="grey") 
  #geom_hline(yintercept = 6.39597,lty="dashed",lwd=.5) +
  #geom_hline(yintercept = 8.260492,lty="dashed",lwd=.25,color="grey")  #plot GW manhattan 

p1 

dev.off()


# res_vanilla %>% filter(chr==5) %>% 
#   ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
#   geom_rect(xmin=23.5,xmax=30,ymin=0,ymax=15, alpha=.5, fill="lightgrey",color="lightgrey") +
#   geom_point(alpha=.5) +
#   #facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
#   theme_classic() +
#   xlab("Position in Mb") +
#   ylim(0,15) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   geom_hline(yintercept = 6.39597,lty="dashed",lwd=.5) +
#   geom_hline(yintercept = 8.260492,lty="dashed",lwd=.25,color="grey")  #just chromosome 5, outlining EPSPS amplification


# plot
summary(res$FDR_corr)
gwa_FDR <- res[res$FDR_corr < 0.05,]
gwa_FDR <- res[res$bon < 0.05,]
gwa_lrt <- res[res$p_lrt < 0.05,]
max(gwa_FDR$p_wald)

pdf(file=paste("gwas_corrected_subset.pdf", sep="" ), bg = "transparent", width=10, height=3, family = "Times New Roman")
#subsample
subres <- res[sample(1:dim(res)[1], 500000),]
p2 <- subres %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_point(alpha=.5) +
  facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
  theme_classic() + 
  ylim(0,10) +
  xlab("Position in Mb") +
  theme(axis.text.x = element_text(angle = 90),legend.position = "none") +
  geom_hline(yintercept = -log(0.05),lty="dashed",lwd=1.25,color="grey") 
  #geom_hline(yintercept = ,lty="dashed",lwd=.5) +
  
p2
dev.off()





res %>% filter(chr==1) %>% 
  ggplot(aes(ps/1000000,-log10(p_lrt), color=af)) +
  geom_rect(xmin=23.5,xmax=30,ymin=0,ymax=15, alpha=.5, fill="lightgrey",color="lightgrey") +
  geom_point(alpha=.5) +
  #facet_grid(. ~ chr, scales = "free_x", space = "free_x",) +
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,15) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = 6.39597,lty="dashed",lwd=.5) +
  geom_hline(yintercept = 8.260492,lty="dashed",lwd=.25,color="grey")  #just chromosome 1, where high prop of sig SNPs map

#to plot no covariate and TSR manhattans together
library(lemon)
grid_arrange_shared_legend(p1,p2, nrow=2, ncol=1, position="right")
