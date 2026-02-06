#!/usr/bin/env Rscript

# qqplot and Manhattan plots from CMH scan
rm(list= ls())

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
# read in regression results files
# ----------------------------------------------------------------------------- #

#new runs with -miss 0.2 and cutoff at 0.9
cmh_res <- read.table("FDRdrought", header= T, sep = " ")
#make new column for chrom values
cmh_res[c("scaffold","chrom")] <- as.numeric(str_split_fixed(cmh_res$CHR, '_', 2))


df_vec <- c("cmh_res") #if more than one dataset
for (d in 1:length(df_vec)) {

  df <- get(df_vec[d])
  N <- dim(df)[1]
  if (length(df$P) >0 ){
  # df$FDR <- p.adjust(df$pval_cmh, method = "fdr", n = N)
  # df$bon <- p.adjust(df$pval, method = "bonferroni", n = N)
  df$qlogp <- qlog10(df$P)
  df <- df[order(-df$qlogp),]
  }

  #update dataset
  assign(df_vec[d], df)

}
# ----------------------------------------------------------------------------- #
# QQPlots
# ----------------------------------------------------------------------------- #

pdf(file = "cmh_qqplot_pvalue.pdf",
    bg = "transparent", width=5, height=5, family = "Times New Roman")

# QQplot with base 10 - lrt test for ancestry corrected gemma gwas
ggplot(cmh_res, aes(x = qlogp, y = -log10(P))) +
  geom_path() + ylim(0,10) + xlim(0,10) +
  geom_abline(intercept = 0, slope = 1) +
  geom_hline(yintercept = -log10(0.01), lty="dashed",lwd=.5, col = "red")# + 
  #geom_vline(xintercept = -log10(0.015), lty="dashed",lwd=.5, col = "red")

dev.off()

# ----------------------------------------------------------------------------- #
# Manhattan plot
# ----------------------------------------------------------------------------- #
#subsample randomly to plot 
cmh_res_sub <- cmh_res[sample(1:dim(cmh_res)[1],500000, replace = F),]

#CMH scan
pdf(file="manhattan_cmh_pvalue.pdf",
    bg = "transparent", width=10, height=3, family = "Times New Roman")

p1<- cmh_res_sub %>%
  ggplot(aes(BP/1000000,-log10(FDR_p))) +
  geom_point(alpha=.5, col = "gray") +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  #facet_wrap(~ chrom, scales = "free_x")+
  theme_classic() +
  xlab("Position in Mb") +
  ylim(0,7.5) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_hline(yintercept = -log10(0.05),lty="dashed",lwd=.5, col = "red")
  # geom_hline(yintercept = -log10(lmm_gemma_anc$p_wald[idx2]),lty="dashed",lwd=.5, col = "#929292")
p1

dev.off()
