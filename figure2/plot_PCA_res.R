# 07.24.24
# Plot PCA res from drought selection experiment data. 
# Data from the drought experiment: all samples were subjected to the drought (no control was sequenced).

rm(list= ls())

# libraries
library(colourvalues)
library(extrafont)
library(ggplot2)
#library(tidyft)
library(magrittr)
#library(tidyr)
library(dplyr)
library(lsmeans)
library(car)

# ----------------------------------------------------------------------------- #
# initial PCA with outliers
# ----------------------------------------------------------------------------- #

# Upload data
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/1.pca/1.pca_v1")
pca_output <- read.table("drought_eigenvec.txt", sep = "\t", header = T)
summary(pca_output)

# plot pca
pca_output$env_cols <- color_values(pca_output$env, palette = "viridis")
plot(pca_output$PC1, pca_output$PC2, pch = 16, col = pca_output$env_cols) # I wonder what is going on with PC1 - 

# we have two samples that seem to be outliers - maybe other species. I will remove them and run the PCA with Plink again. 
# identify the outliers
idx <- which(pca_output$PC2 < -0.1) # samples 44 and 173
pca_output$samp[idx] # "P16_Nat_1_T"  "P12_Nat_14_T"
# the two outliers are natural samples!

# modify the .fam file to re-run the analyses (in bash)

# ----------------------------------------------------------------------------- #
# PCA without outliers
# ----------------------------------------------------------------------------- #

# Upload data
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/1.pca/2.rm_outliers")
pca_output <- read.table("merged_numericChr_LDpruned_pca.eigenvec.txt", header = T)
summary(pca_output)
variance <- read.table("merged_numericChr_LDpruned_pca.eigenval.txt")
variance_pct <- variance/sum(variance) * 100 # turn to percentage

# sample info
sample_info <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/2.admixture/sample_info.csv", header=T, sep = ",")
# remove outliers
idx <- c(which(sample_info$samp== c("P16_Nat_1_T")), which(sample_info$samp== c("P12_Nat_14_T")) )
sample_info <- sample_info[-idx,]
# order sample_info based on longitude, then order pca samples (and admixture data) based on sample_info
sample_info_ordered <- sample_info[order(sample_info$long),]
sample_info_ordered <- sample_info_ordered[,c(1,6,7)]
# join datasets together
pca_output <- inner_join(sample_info_ordered,pca_output, by=c("samp"))



# ----------------------------------------------------------------------------- #
# plots
# ----------------------------------------------------------------------------- #


# plot pca by environmemt
par(family = "Times New Roman", mfrow = c(1,2))
pca_output$env_cols <- color_values(pca_output$env, palette = "spectral")
plot(pca_output$PC1, pca_output$PC2, pch = 16, col = pca_output$env_cols, 
     xlab = paste("PC1, ",round(variance_pct$V1[1],2),"%"), 
     ylab = paste("PC2, ",round(variance_pct$V1[2],2),"%")) 

legend("topright", c("Ag", "Nat"), pch = 16, col = c(unique(pca_output$env_cols)), box.lwd = 0)

plot(pca_output$PC3, pca_output$PC4, pch = 16, col = pca_output$env_cols,
     xlab = paste("PC3, ",round(variance_pct$V1[3],2),"%"), 
     ylab = paste("PC4, ",round(variance_pct$V1[4],2),"%"))


# color pca by pair
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/")
pdf("pca_colored_by_pair_supp.pdf", bg = "transparent", width=8, height=5, family = "Times New Roman")

par(family = "Times New Roman", mfrow = c(1,2))
pca_output$pair_cols <- color_values(pca_output$pair, palette = "purples")
plot(pca_output$PC1, pca_output$PC2, pch = 21,
     col ="#252525", bg = pca_output$pair_cols, lwd=0.2,
     xlab = paste("PC1, ",round(variance_pct$V1[1],2),"%"), 
     ylab = paste("PC2, ",round(variance_pct$V1[2],2),"%")) 

plot(pca_output$PC3, pca_output$PC4, pch = 21, 
     col ="#252525", bg = pca_output$pair_cols, lwd=0.2,
     xlab = paste("PC3, ",round(variance_pct$V1[3],2),"%"), 
     ylab = paste("PC4, ",round(variance_pct$V1[4],2),"%"))


#legend("top", inset=c(-0.2,0), xpd = T, horiz=TRUE,
#       c(unique(pca_output$pair)), pch = 16, col = c(unique(pca_output$pair_cols)), box.lwd = 0)

dev.off()


setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/")
pdf("pca_colored_by_pair_supp.pdf", bg = "transparent", width=8, height=5, family = "Times New Roman")

par(family = "Times New Roman", mfrow = c(1,2))



# color pca by longitude
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/")

pdf("pca_long_pc1-2.pdf", bg = "transparent", width=6, height=5, family = "Times New Roman")
par(family="Times New Roman", cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1)

pca_output$long <- as.factor(as.character(pca_output$long))
pca_output$long_cols <- color_values(pca_output$long, palette = "purples")
pchvec <- rep(21,length(pca_output$PC1))
pchvec[pca_output$env == "Nat"] <- 22
plot(pca_output$PC1, pca_output$PC2, pch = pchvec, cex = 1.5,
     col ="#252525", bg = pca_output$long_cols, lwd=0.2,
     xlab = paste("PC1, ",round(variance_pct$V1[1],2),"%"), 
     ylab = paste("PC2, ",round(variance_pct$V1[2],2),"%")) 

# plot(pca_output$PC3, pca_output$PC4, pch = 21,cex = 1.5,
#      col ="#252525", bg = pca_output$long_cols, lwd=0.2,
#      xlab = paste("PC3, ",round(variance_pct$V1[3],2),"%"),
#      ylab = paste("PC4, ",round(variance_pct$V1[4],2),"%"))

dev.off()

# ----------------------------------------------------------------------------- #
# PC1 vs PC3
# ----------------------------------------------------------------------------- #
pdf("pca_long_pc1-4_opposed.pdf", bg = "transparent", width=8, height=6, family = "Times New Roman")
par(family="Times New Roman", cex.lab=1.5, cex.axis=1.5, cex.main=1, cex.sub=1, mfrow = c(1,2))

plot(pca_output$PC1, pca_output$PC3, pch = 21, cex = 1.5,
     col ="#252525", bg = pca_output$long_cols, lwd=0.2,
     xlab = paste("PC1, ",round(variance_pct$V1[1],2),"%"), 
     ylab = paste("PC3, ",round(variance_pct$V1[3],2),"%"))

plot(pca_output$PC2, pca_output$PC3, pch = 21, cex = 1.5,
     col ="#252525", bg = pca_output$long_cols, lwd=0.2,
     xlab = paste("PC2, ",round(variance_pct$V1[2],2),"%"), 
     ylab = paste("PC3, ",round(variance_pct$V1[3],2),"%"))

dev.off()


# ----------------------------------------------------------------------------- #
# export PC values
# ----------------------------------------------------------------------------- #
df_to_export <- pca_output[,-c(28:30)]
write.table(df_to_export, "/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/1.pca/2.rm_outliers/PCA_output_noout.txt",
            sep = "\t", col.names = T, row.names = F)



# ----------------------------------------------------------------------------- #
# get data ready for gwas
# ----------------------------------------------------------------------------- #
library(gsheet)
phenos <- gsheet2tbl('docs.google.com/spreadsheets/d/1cQCLJCZzQtytUC_PaClJhNBud_S9-x6XCSAUuNy8LsQ/edit#gid=1855741146')
phenos_mod <- data.frame(phenos$Label, phenos$`Day of Full Wilt`)
# filter phenos based on sample_info
phenos_mod_fil <- matrix(NA, dim(sample_info_ordered)[1], 2)
for (i in 1:dim(sample_info_ordered)[1]) {
  
  idx <- which(sample_info_ordered$samp[i] == phenos_mod$phenos.Label)
  phenos_mod_fil[i,1] <- phenos_mod[idx,1]
  phenos_mod_fil[i,2] <- phenos_mod[idx,2]
  
}

# export phenos
#write.table(phenos_mod_fil, "/Users/rozenn/GaTech Dropbox/CoS/BioSci/BioSci-Ratcliff/Rozenn/36.Chicago/2.DroughtProject/1.analyses/data/3.gwas/drought_phenos.txt",
#            sep = "\t", col.names = F, row.names = F, quote = FALSE)


# ----------------------------------------------------------------------------- #
# ld50 as a function of PC1
# ----------------------------------------------------------------------------- #

# find mean PC1 value per population
N <- length(levels(as.factor(pca_output$pair))) * length(levels(as.factor(pca_output$env))) 
pc1_stats <- pca_output %>% group_by(pair, env, long_cols) %>% summarise(meanPC1 = mean(PC1), stdPC1 = sd(PC1)/sqrt(N))
#how many samples per group?
num_samp_per_group <- pca_output %>% count(pair, env)
num_samp_per_group[order(match( num_samp_per_group$pair, unique(pca_output$pair))), ]
pc1_stats$pairP <- pc1_stats$pair
pc1_stats$pair <- gsub("P","", pc1_stats$pair)
pc1_stats$pair <- as.numeric(pc1_stats$pair)

# correlate with LD50
ld50 <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/4.phenotypes/ld50.txt",
                   sep = "\t", header = T)


summary(ld50$ld50)

fulldt <- inner_join(ld50, pc1_stats, by=c("pair", "env"))

pdf("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/pc1_vs_ld50.pdf",
    bg = "transparent", width=5, height=4, family = "Times New Roman")

par(family="Times New Roman", cex.lab=1.5, cex.axis=1.5, cex.main=1, cex.sub=1, mfrow = c(1,1))
plot("n", xlim = c(-.15,.1), ylim = c(5,15), 
     xlab = "Mean PC1 score (per population)", 
     ylab = "LD50")

# simple linear regression
lm.res <- lm(fulldt$ld50~fulldt$meanPC1)
summary(lm.res)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2.35501 -0.74622  0.00732  0.81123  1.89271 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     10.9271     0.2481   44.04  < 2e-16 ***
#   fulldt$meanPC1  16.0810     4.1770    3.85 0.000999 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.155 on 20 degrees of freedom
# Multiple R-squared:  0.4256,	Adjusted R-squared:  0.3969 
# F-statistic: 14.82 on 1 and 20 DF,  p-value: 0.0009992

# multiple linear regression with longitude as covariate
lm.res <- lm(ld50 ~ meanPC1 + long + env, data = fulldt)

#with PC1 in place of ancestry
#multiple linear regression : Drought ~ PC1 + Long + Lat  + Env + Env:PC1 + Long:PC1
lm.res <- lm(ld50 ~ meanPC1 + long + lat + env + long:meanPC1 , data = fulldt) # + env:meanPC1
Anova(lm.res, type = "III") # R2 = 0.6697
summary(lm.res)
#interaction + env:meanPC1 term removed because not significant

#drop ancestry to calculate partial R2
lm.res <- lm(ld50 ~long + lat + env , data = fulldt) # + env:meanPC1
summary(lm.res) # 0.3567

0.6697-0.3567 #=0.313

lines(c(-.15, .1), 10.9271 + 16.0810 * c(-.15, .1), lwd = 2, col = "gray", lty = 2)

points(fulldt$meanPC1, fulldt$ld50, pch = 21,  cex = 1.5,
       bg = fulldt$long_cols, col = "black", lwd = 0.3)

arrows(fulldt$meanPC1 + fulldt$stdPC1, fulldt$ld50,  fulldt$meanPC1 - fulldt$stdPC1, fulldt$ld50, code = 3, angle = 90, length = 0, col = "gray")
arrows(fulldt$meanPC1, fulldt$ld50_minus_stderr,  fulldt$meanPC1, fulldt$ld50_plus_stderr, code = 3, angle = 90, length = 0, col = "gray")


dev.off()

#note: why are the standard errors so small for higher PC values?
#from pc1 versus pc2 plot, it seems like there is more variaiton in the southeatstern populations


# ----------------------------------------------------------------------------- #
# ld50 as a function of admixture
# ----------------------------------------------------------------------------- #

pca_output$long <- as.factor(as.character(pca_output$long))
pca_output$long_cols <- color_values(pca_output$long, palette = "purples")
library(magrittr)
library(dplyr)
k2 <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/2.admixture/k2_structure_res.csv", 
                 sep = ",", header = T)
summary(k2)
k2$long_cols <- color_values(k2$long, palette = "purples")
phenos$samp <- phenos$Label
k2_wilt <- merge(k2, phenos, by = "samp")

# ----------------------------------------------------------------------------- #
# individual-level plot
# ----------------------------------------------------------------------------- #

pdf("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/admixture_wilt.pdf",
    bg = "transparent", width=6, height=5, family = "Times New Roman")


par(family="Times New Roman", cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1, mfrow = c(1,1))
plot("n", xlim = c(0,1), ylim = c(0,23), 
     xlab = expression(paste(italic("var. rudis"), " proportion per sample")), 
     ylab = "Day to full wilt")

#multiple linear regression : Drought ~ Long + Lat + Ancestry + Env + Env:Ancestry + Long:Ancestry
lm_res <- lm(`Day of Full Wilt` ~ long + lat + k2.V1 + env + long:k2.V1, data = k2_wilt) # + meank2:env
summary(lm_res) #multiple R2 = 0.6617
Anova(lm_res, type="III")

#for plot
lm1 <- lm(`Day of Full Wilt`~ k2.V1, data = k2_wilt)
summary(lm1)

lines(c(0, 1),  lm1$coefficients[1] + lm1$coefficients[2] * c(0, 1), lwd = 2, col = "#353535", lty = 2)

pchvec <- rep(21, length(k2_wilt$samp))
pchvec[k2_wilt$env == "Nat"] <- 22
points(k2_wilt$k2.V1, k2_wilt$`Day of Full Wilt`, pch = pchvec,  cex = 1.5,
       bg = k2_wilt$long_cols, col = "#353535", lwd = 0.3)

text(0.3,2, expression(paste("r"^2 , " = 0.065, p-value = 1.41e-05")))

dev.off()

# ----------------------------------------------------------------------------- #
# pop-level plot
# ----------------------------------------------------------------------------- #

#summarize admixture by pair and environment
k2_stats <- k2 %>% group_by(pair, env, long_cols) %>% summarise(meank2 = mean(k2.V1), stdk2 = sd(k2.V1)/sqrt(N)) #proportion var.rudis
k2_stats$pairP <- k2_stats$pair
k2_stats$pair <- as.numeric(gsub("P", "", k2_stats$pair))

fulldt <- inner_join(ld50, k2_stats, by=c("pair", "env"))

pdf("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/admixture_vs_ld50.pdf",
bg = "transparent", width=6, height=5, family = "Times New Roman")

par(family="Times New Roman", cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1, mfrow = c(1,1))
plot("n", xlim = c(0,1), ylim = c(5,15), 
     xlab = expression(paste("Mean ", italic("var. rudis"), " proportion")), 
     ylab = expression(paste("Drought LD"[50], " (days)")))


#multiple linear regression : Drought ~ Long + Lat + Ancestry + Env + Env:Ancestry + Long:Ancestry
lm_res <- lm(ld50 ~ long + lat + meank2 + env + long:meank2, data = fulldt) # + meank2:env
summary(lm_res) #multiple R2 = 0.6617
#for partial R2 for mean ancestry
summary(lm(ld50 ~ long + lat + env, data = fulldt)) #multiple R2 = 0.3567
0.6617 - 0.3567 # 0.305 partial R2 for ancestry
Anova(lm_res, type="III")
#dropped interaction term that was not significant

#pull out main effects at different longitudes
lsm <- lsmip(lm_res, meank2 ~ long, plot = F, at = list(meank2 = 0.5, long = c(-96,-93,-90,-86,-82)))
lsm

#drop ancestry to calculate R2 without ancestry
# lm_res <- lm(ld50 ~ long + lat + env, data = fulldt) # + meank2:env
# summary(lm_res) #multiple R2 = 0.3567
#partial R2
0.6617-0.3567# = 0.305

#long varies from
-42.54637 - 0.54116*min(fulldt$long) #8.98 (southwest)
-42.54637 - 0.54116*max(fulldt$long) #2.79 (northeast)



#least square means 
lsm <- lsmeans(lm_res, ~  meank2 + env)
print(lsm)

#for plot
lm1 <- lm(ld50 ~ meank2, data = fulldt)
summary(lm1)
Anova(lm1)
lsmeans(lm1, ~ meank2)

lines(c(0, 1),  lm1$coefficients[1] + lm1$coefficients[2] * c(0, 1), lwd = 2, col = "#353535", lty = 2)


N <- length(fulldt$stdk2)
arrows(fulldt$meank2 + fulldt$stdk2/sqrt(N), fulldt$ld50,  fulldt$meank2 - fulldt$stdk2/sqrt(N), fulldt$ld50, code = 3, angle = 90, length = 0, col = "gray")
arrows(fulldt$meank2, fulldt$ld50_minus_stderr,  fulldt$meank2, fulldt$ld50_plus_stderr, code = 3, angle = 90, length = 0, col = "#353535")

pchvec <- rep(21, length(fulldt$ld50))
pchvec[fulldt$env == "Nat"] <- 22
points(fulldt$meank2, fulldt$ld50, pch = pchvec,  cex = 1.5,
       bg = fulldt$long_cols, col = "#353535", lwd = 0.3)


dev.off()



#plot projections
#pull out main effects at different longitudes
lm_res <- lm(ld50 ~ long + lat + meank2 + env + long:meank2, data = fulldt) # + meank2:env
long_groups <- c(-94.5,-91.5,-88,-84)
lsm <- lsmip(lm_res, meank2 ~ long, plot = F, at = list(meank2 = c(0,1), long = long_groups))
lsm

pdf(file=paste("Supp_amdx_drought_long.pdf", sep="" ), bg = "transparent", width=5, height=4, family = "Times New Roman")

par(family= "Times New Roman")
long_col <- color_values(c(-100,long_groups), palette = "purples") #take the middle of each longitude group

#Plot for linear predictions of temperature
plot(lsm$meank2[1:2], lsm$yvar[1:2], type = "l", col = long_col[2], lwd =2,
     xlab = "Mean ancestry", ylab = "Linear prediction", ylim = c(0,15), lty = 1, xlim = c(0,1))
lines(lsm$meank2[3:4], lsm$yvar[3:4], type = "l", col = long_col[3], lwd =2)
lines(lsm$meank2[5:6], lsm$yvar[5:6], type = "l", col = long_col[4], lwd =2)
lines(lsm$meank2[7:8], lsm$yvar[7:8], type = "l", col = long_col[5], lwd =2)

legend("bottomright", legend = long_groups, col = long_col[2:5],
  bty = "n", lty = c(1,1,1,1), lwd = 2, title = "Longitude", cex = 0.8, border=F, ncol=1)


dev.off()

#make legend
fulldt_ordered <- fulldt[order(fulldt$long), ]
legend_image <- as.raster(matrix(unique(fulldt_ordered$long_cols), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'legend title')
text(x=1.5, y = seq(0,1), labels = c(round(fulldt_ordered$long[1]),round(fulldt_ordered$long[22])))
rasterImage(legend_image, 0, 0, 1,1)


# ----------------------------------------------------------------------------- #
# multiple regression models
# ----------------------------------------------------------------------------- #

#join admixture and PC values
fulldt <- inner_join(pca_output, k2, by=c("samp"))

# what variables explain the most how the data cluster for the first principal component?
# variables: Geographic (lat, long), Origin (nat versus ag, in pairs), admixture
lmres1 <- lm(PC1 ~ lat.y + long.y + env.y + k2.V2 + env.y:k2.V2, data = fulldt) 
Anova(lmres1, type = "III") #0.9893
summary(lmres1) #Multiple R-squared:  0.9893

#correlation coef value 
cor.test(fulldt$PC1, fulldt$k2.V2)
#partial R2
lmres1 <- lm(PC1 ~ lat.y + long.y + env.y, data = fulldt)
summary(lmres1) #Multiple R-squared:  0.7132  


# Individual-level analyses
lmres2 <- lm(PC1 ~ lat.y + long.y + env.y + PC1 + long.y:PC1, data = fulldt) 
Anova(lmres2, type = "III")

lmres2 <- lm(k2.V2 ~ lat.y + long.y + env.y + long.y:k2.V2, data = fulldt) 
Anova(lmres2, type = "III")


#without admixture
lmres1 <- lm(PC1 ~ lat.y + long.y + env.y , data = fulldt) 
summary(lmres1) #multiple R2 = 0.7132

# Call:
#   lm(formula = PC1 ~ lat.y + long.y + env.y, data = fulldt)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.077285 -0.020216 -0.002158  0.028624  0.106078 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.2436326  0.1655072   1.472    0.142    
# lat.y       -0.0236474  0.0026457  -8.938  < 2e-16 ***
#   long.y      -0.0079076  0.0008045  -9.829  < 2e-16 ***
#   env.yNat    -0.0168571  0.0038571  -4.370 1.76e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.03224 on 276 degrees of freedom
# Multiple R-squared:  0.7132,	Adjusted R-squared:   0.71 
# F-statistic: 228.7 on 3 and 276 DF,  p-value: < 2.2e-16

#Partial R2 for latitude
lmres1 <- lm(PC1 ~ long.y + env.y, data = fulldt) 
summary(lmres1)
0.7132 - 0.6301 #= 0.0831
#Partial R2 for longitude
lmres1 <- lm(PC1 ~ lat.y, data = fulldt) 
summary(lmres1)
0.7132 - 0.5955 #= 0.1177
#Partial R2 forenv
lmres1 <- lm(PC1 ~ long.y + lat.y, data = fulldt) 
summary(lmres1)
0.7132 - 0.6933 #= 0.0199



lmres2 <- lm(PC2 ~ lat.y + long.y + env.y + k2.V2 + env.y:k2.V2, data = fulldt) # lat + 
summary(lmres2)
Anova(lmres2)

# Call:
#   lm(formula = PC2 ~ lat.y + long.y + env.y + k2.V2 + env.y:k2.V2, 
#      data = fulldt)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.18482 -0.01290 -0.00073  0.01122  0.32365 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    -0.156534   0.273018  -0.573   0.5669    
# lat.y          -0.005079   0.005060  -1.004   0.3164    
# long.y         -0.003804   0.001476  -2.576   0.0105 *  
#   env.yNat        0.014500   0.008679   1.671   0.0959 .  
# k2.V2           0.120833   0.020597   5.866 1.28e-08 ***
#   env.yNat:k2.V2 -0.146768   0.018344  -8.001 3.52e-14 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05242 on 274 degrees of freedom
# Multiple R-squared:  0.247,	Adjusted R-squared:  0.2332 
# F-statistic: 17.97 on 5 and 274 DF,  p-value: 2.054e-15


#without admixture
lmres2 <- lm(PC2 ~ lat.y + long.y + env.y, data = fulldt) # lat + 
summary(lmres2)

# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.116927 -0.021718 -0.008751  0.016451  0.286964 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.0811175  0.2986152  -0.272    0.786    
# lat.y       -0.0007813  0.0047735  -0.164    0.870    
# long.y      -0.0014120  0.0014515  -0.973    0.332    
# env.yNat    -0.0286747  0.0069592  -4.120    5e-05 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.05816 on 276 degrees of freedom
# Multiple R-squared:  0.06625,	Adjusted R-squared:  0.0561 
# F-statistic: 6.528 on 3 and 276 DF,  p-value: 0.0002795


#Partial R2 for env
lmres1 <- lm(PC2 ~ long.y + lat.y, data = fulldt) 
summary(lmres1)
0.06625 - 0.008816 #0.057434


# ----------------------------------------------------------------------------- #
# Least square means on LD50 with admixture and geography
# ----------------------------------------------------------------------------- #

#linear model
lm_res <- lm(ld50 ~ env + long_group , data = ld50)
summary(lm_res)

#least square means 
lsm <- lsmeans(lm_res, ~  env + long_group)
print(lsm)

# ----------------------------------------------------------------------------- #
# multiple regression models
# ----------------------------------------------------------------------------- #
colnames(phenos_mod) <- c("samp", "wilt")
fulldt_wilt <- merge(fulldt, phenos_mod, by = "samp" )

# Individual-level analyses
lmres <- lm( wilt ~ lat.y + long.y + env.y + PC1 + long.y:PC1, data = fulldt_wilt) 
Anova(lmres, type = "II")

lmres <- lm(wilt ~ lat.y + long.y + env.y + k2.V2 + long.y:k2.V2, data = fulldt_wilt) 
Anova(lmres, type = "II")

