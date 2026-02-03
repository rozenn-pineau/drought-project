# Plot LD decay rate
rm(list= ls())

# libraries
library(extrafont)
library(minpack.lm)

#upload data
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/12.LD_decay/")
ld_decay <- read.table("genotype_LDdecay.stat.gz",
                       header = F, sep = "\t")
summary(ld_decay)
colnames(ld_decay) <- c("#Dist",   "Mean_r^2", "Mean_D'", "Sum_r^2", "Sum_D'",  "NumberPairs")


#Fit exponential decay function
x <- ld_decay$`#Dist`
y <- ld_decay$`Mean_r^2`
c0 <- min(y)
a0 <- max(y) - c0
b0 <- 1 / (max(x) - min(x))

fit <- nlsLM(
  y ~ a * exp(-b * x) + c,
  start = list(a = a0, b = b0, c = c0),
  lower = c(a = 0, b = 0, c = 0)
)


#plot
pdf("ld_decay.pdf", bg = "transparent", width=5, height=4, family = "Times New Roman")

plot(ld_decay$`#Dist`, ld_decay$`Mean_r^2`, pch = 16, ylab = expression('Mean r'^2), xlab = "Distance (bp)",
     xlim = c(0,50000))
lines(x, predict(fit), col = "red", lwd = 2)
text(30000,0.1,"decay rate b = 0.00011")

dev.off()

#Half-decay distance
log(2)/0.000106 #6539.124

summary(ld_decay$`Mean_r^2`)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.3452  0.3795  0.4132  0.4638  0.4880  1.0000
summary(ld_decay$NumberPairs) # 1556 


#How many loci per clump?
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/12.LD_decay/")
clump_info <- read.table("clump_summary.txt",
                         sep = "\t", header = T)

pdf(file="gene_per_clump_info.pdf", bg = "transparent", width=6, height=4, family = "Times New Roman")

par(family="Times New Roman", mfrow = c(1,2))
plot(clump_info$n_variants, clump_info$length_bp, pch = 16, cex = 1.2,
     xlab = "Number of SNP per clump", 
     ylab = "Length of clump (bp)")

summary(clump_info$n_variants)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9.0    37.0    76.0   117.5   141.0   776.0

summary(clump_info$length_bp)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6360   27889   62751   77004  109111  196411 

#How many genes per clump?
gene_info <- read.table("number_of_overlapping_genes_in_clumps.txt",
                        sep = "\t", header = F)
hist(gene_info[,6], breaks = 15, xlab = "Number of genes per clump", main = "")

dev.off()

summary(gene_info[,6])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1.000   2.500   5.000   6.442   9.000  25.000 
