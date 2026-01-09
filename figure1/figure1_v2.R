#this script : plots raw phenotypic data for survival through time,
#calculates LD50 from non linear fit,
#maps LD50 value in space.

library(gsheet)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(data.table)
library(car)
library(lsmeans)
library(extrafont)
library(colourvalues)
library(lsmeans)

rm(list= ls())

# ------------------------------------------------------------------------------------- #
# prep data
# ------------------------------------------------------------------------------------- #

setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/4.phenotypes")
phenos <- fread("phenos_drought.csv")
phenos$Pair[phenos$Pair == "P12"] <- "P10" #fix pair swap during labeling
phenos$Pair<-as.numeric(gsub("P","",phenos$Pair)) #change to numeric
meta <-fread("meta.csv") #metadata

#parse and process
head(meta)
names(meta)[c(9,10)]<-c("Env","Pair") # modify header
meta<-meta[c(1:25),c(1:11)]

#merge with metadata
phenos_wmeta<-inner_join(phenos,meta, by=c("Pair","Env")) # join metadata and phenotypic data in a single dataset based on Pair and Env columns
phenos_wmeta$Lat<-as.numeric(phenos_wmeta$Lat)
(check<-phenos_wmeta %>% group_by(Treatment, Env) %>% drop_na(`Stem Width (mm)`) %>% summarise(n=n())) #sample size within each treatment, by env
phenos_wmeta<-phenos_wmeta %>% filter(is.na(`Day of First Wilt`)|`Day of First Wilt`!="drop") %>% filter(is.na(`Day of First Wilt`)|`Day of First Wilt`!="?") #remove samples that should be dropped from experiment. we had a bunch of cheaters that grew roots under the plastic and got access to water without permission!
sum(check$n)
phenos_wmeta<-phenos_wmeta %>% filter(is.na(`Day of First Wilt`)|`Day of First Wilt`!="?") %>% filter(is.na(`Day of First Wilt`)|`Day of First Wilt`!="drop") # drop nas

(check<-phenos_wmeta %>% group_by(Treatment, Env) %>% drop_na(`Stem Width (mm)`) %>% summarise(n=n()))
sum(check$n) #final sample size of full exp

phenos_wmeta$`Day of First Wilt`<-as.numeric(phenos_wmeta$`Day of First Wilt`)
phenos_wmeta$`Day of First Wilt`<-as.numeric(phenos_wmeta$`Day of First Wilt`)
sort(unique(phenos_wmeta$`Day of First Wilt`)) #days at which plants were recorded to die
phenos_wmeta$`Day of Full Wilt` <- as.numeric(phenos_wmeta$`Day of Full Wilt`)

# ------------------------------------------------------------------------------------- #
# plot raw data grouped by longitude
# ------------------------------------------------------------------------------------- #
#save figure
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/4.phenotypes/")
pdf("raw_survival_proportion.pdf", 
    bg = "white", width=8, height=4.5, family = "Times New Roman")

#day to full wilt by longitude group
par(family = "Times New Roman",  cex.lab=1.5, cex.axis=1)
plot(1, pch = "", xlab = "Time (days)", ylab = "Survival proportion", xlim = c(0,20), ylim = c(0,1))
#add grid
#abline(v = c(seq(0,20,by=2)), lty = 1, col = "grey", lwd = 0.5)
#abline(h = c(0,0.2,0.4,0.6,0.8,1),  lty = 1, col = "grey", lwd = 0.5)

#make longitude + env groups
long_loop <- c(-96,-93,-90,-86,-82)
long_cols <- color_values(long_loop, palette = "purples")

env_loop <- c("Nat", "Ag")
env_types <- c(1,2)
pch_vec <- c(22, 21)

for (e in 1:2) {
  
  pch_val <- pch_vec[e]
  
  for (a in 1:4) { #loop through longitudes
    
    idx <- which( phenos_wmeta$Long > long_loop[a] & phenos_wmeta$Long  <= long_loop[a+1] & phenos_wmeta$Env == env_loop[e] & !is.na(phenos_wmeta$`Day of Full Wilt`) ) # find indices for each sub group (divided by environment and longitude here)
    sub <- phenos_wmeta[idx,]
    n_days <- 20
    n_ind <- length(idx)
    surv_dt <- matrix(0, n_ind, n_days)
    
    for (n in 1:n_ind) {
      surv_dt[n, 1:sub$`Day of Full Wilt`[n]] <- 1
    }
    #calculate proportion surviving per group
    surv_prop <- colSums(surv_dt)/n_ind
    
    lines(0:n_days, c(1,surv_prop), col = long_cols[a+1], lty = e)
    points(0:n_days, c(1,surv_prop), pch = pch_val, col = "#353535", bg = long_cols[a+1] )
    #add se : formula sqrt[p(1-p)/n)
    se <- sqrt(surv_prop*(1-surv_prop)/n_ind)
    arrows(0:n_days, c(1,surv_prop) + se,
           0:n_days, c(1,surv_prop) - se, code = 3, 
           angle = 90, length= 0.01, col =  long_cols[a+1])
    
  }
}

# legend("bottomleft", title="Longitude         Env", 
#        legend=c("(-96,-93] ","(-93,-90] ","(-90,-86]", "(-86,-82]","Ag","Nat"), 
#        col=c(long_cols[2:5], "black", "black"), pch=c(16,16,16,16,16,15), lty = c(NA,NA,NA,NA,2,1),
#        bty="n", border=F, ncol=2)

dev.off()

#save legend for longitude
pdf("legend_longitude.pdf", 
    bg = "white", width=6, height=3, family = "Times New Roman")

par(family = "Times New Roman", cex = 1.3)
scale_data <- data.frame(a = 4, x = min(phenos_wmeta$Long):max(phenos_wmeta$Long))
scale_data$cols <- colour_values(scale_data$x, palette = "purples")
barplot(height = scale_data$a, col = scale_data$cols, space = 0, border = NA, names.arg = round(scale_data$x), yaxt='n', ann=FALSE) # 

dev.off()

# ------------------------------------------------------------------------------------- #
# compare day of full wilt between pairs of nat versus env samples
# ------------------------------------------------------------------------------------- #
# https://web.archive.org/web/20160527230917/http://socserv.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Nonlinear-Regression.pdf


pair_vec <- as.numeric(levels(as.factor(phenos_wmeta$Pair)))
store <- matrix(NA, length(pair_vec))
ld50_ag <- matrix(NA, length(pair_vec), 3) #record LD50, L50-standard error and LD50+standard error
ld50_nat <- matrix(NA, length(pair_vec), 3)
lats <- matrix(NA, length(pair_vec))
longs <- matrix(NA, length(pair_vec))
ag_surv <- c()
nat_surv <- c()

ag_col <- rgb(241/255, 80/255, 103/255)
nat_col <- rgb(76/255, 175/255, 202/255)

pdf("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/ld50_ag_nat.pdf", bg = "transparent", width=5, height=4, family = "Times New Roman")

par(family="Times New Roman", cex.lab=1.3, cex.axis=1.3, cex.main=1, cex.sub=1)
plot(1, type="n", xlab="Days", ylab="Survival proportion", xlim=c(0, 25), ylim=c(0, 1)) # plot to be filled
col_pair <- colour_values(unique(phenos_wmeta$Pair), palette = "blues")
lines(c(1,25), rep(0.5,2), col = "gray", lwd = 2, lty = 2) #col_pair[p]

for (p in 1 : length(pair_vec)) {
  
  # group data per pair and environment
  ag_vec <- phenos_wmeta$`Day of Full Wilt`[which(phenos_wmeta$Pair == pair_vec[p] & phenos_wmeta$Env == "Ag" & phenos_wmeta$Treatment == "T")]
  nat_vec <- phenos_wmeta$`Day of Full Wilt`[which(phenos_wmeta$Pair == pair_vec[p] & phenos_wmeta$Env == "Nat" & phenos_wmeta$Treatment == "T")]
  
  # store survival curves info for AG
  cur <- as.data.frame(table(ag_vec))
  cur$ag_vec <- as.numeric(as.character(cur$ag_vec))
  cur <- rbind(c(1,0), cur) # Constrain the curve to start at 100% survival at day 1
  cur$cum <- cumsum(cur$Freq)
  cur$surv <- 1-cur$cum/max(cur$cum)
  cur$pair <- pair_vec[p]
  ag_surv <- rbind(ag_surv, cur)
  # fit logistic
  x <- cur$ag_vec
  y <- cur$surv
  log.ss <- nls( y ~ SSlogis(x, phi1, phi2, phi3))
  phi1 <- summary(log.ss)$coef[1]
  phi2 <- summary(log.ss)$coef[2]
  phi3 <- summary(log.ss)$coef[3]
  #std error
  phi1_stderr <- summary(log.ss)$coef[4]
  phi2_stderr <- summary(log.ss)$coef[5]
  phi3_stderr <- summary(log.ss)$coef[6]
  
  # calculate LD50 based on parameters (equation is y = phi1/(1+exp((phi2-x)/phi3)), see above link)
  ld50_ag[p, 1] <- phi2 - phi3*log((phi1/0.5 - 1)) 
  ld50_ag[p, 2] <- (phi2-phi2_stderr) - (phi3-phi3_stderr)*log(( (phi1-phi1_stderr)/0.5 - 1))
  ld50_ag[p, 3] <- (phi2+phi2_stderr) - (phi3+phi3_stderr)*log(( (phi1+phi1_stderr)/0.5 - 1)) 
  
  # plot
  points(x,y, pch = 16, col = ag_col, cex = 1) #col_pair[p]
  lines(1:20, phi1/(1+exp((phi2-c(1:20))/phi3)), col = ag_col , lwd = 1.8) #col_pair[p]
  
  # store survival curves info for NAT
  cur <- as.data.frame(table(nat_vec))
  cur$nat_vec <- as.numeric(as.character(cur$nat_vec))
  cur <- rbind(c(1,0), cur) # Constrain the curve to start at 100% survival at day 1
  cur$cum <- cumsum(cur$Freq)
  cur$surv <- 1-cur$cum/max(cur$cum)
  cur$pair <- pair_vec[p]
  nat_surv <- rbind(nat_surv, cur)
  # fit logistic
  x <- cur$nat_vec
  y <- cur$surv
  log.ss <- nls( y ~ SSlogis(x, phi1, phi2, phi3))
  phi1 <- 1
  phi2 <- summary(log.ss)$coef[2]
  phi3 <- summary(log.ss)$coef[3]
  #std error
  phi1_stderr <- summary(log.ss)$coef[4]
  phi2_stderr <- summary(log.ss)$coef[5]
  phi3_stderr <- summary(log.ss)$coef[6]
  
  # calculate LD50 based on parameters (equation is y = phi1/(1+exp((phi2-x)/phi3)), see above link)
  ld50_nat[p, 1] <- phi2 - phi3*log((phi1/0.5 - 1)) 
  ld50_nat[p, 2] <- (phi2-phi2_stderr) - (phi3-phi3_stderr)*log(( (phi1-phi1_stderr)/0.5 - 1))
  ld50_nat[p, 3] <- (phi2+phi2_stderr) - (phi3+phi3_stderr)*log(( (phi1+phi1_stderr)/0.5 - 1)) 
  
  # calculate mean survival difference b/w environments per pair
  store[p] <- mean(ag_vec, na.rm = T) - mean(nat_vec, na.rm = T)
  
  # record coordinate
  lats[p] <- phenos_wmeta$Lat[which(phenos_wmeta$Pair == pair_vec[p] & phenos_wmeta$Env == "Ag" & phenos_wmeta$Treatment == "T")][1]
  longs[p] <- phenos_wmeta$Long[which(phenos_wmeta$Pair == pair_vec[p] & phenos_wmeta$Env == "Ag" & phenos_wmeta$Treatment == "T")][1]
  
  # plot
  points(x,y, pch = 16, col = nat_col, cex = 1) #col_pair[p]
  lines(1:20, phi1/(1+exp((phi2-c(1:20))/phi3)), col = nat_col, lwd = 1.8) #col_pair[p]
   
  
}

legend("topright", title ="Environment", legend = c("Ag", "Nat"), col = c(ag_col, nat_col), pch = 16, bty = "n")

dev.off()


# ------------------------------------------------------------------------------------- #
# plot example for supp figure
# ------------------------------------------------------------------------------------- #

pdf("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/ld50example.pdf", 
    bg = "transparent", width=5, height=4, family = "Times New Roman")


par(family="Times New Roman", cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1)
plot(1, type="n", xlab="Days", ylab="Survival proportion", xlim=c(0, 25), ylim=c(0, 1)) # plot to be filled
lines(1:20, phi1/(1+exp((phi2-c(1:20))/phi3)), lwd = 1.8) 


val <- phi2 - phi3*log((phi1/0.5 - 1))
lines(c(1,val), rep(0.5,2), col = "gray", lwd = 2, lty = 2) #horizontal line
arrows(val, 0.5, val, -0.03, col = "gray", lwd = 2, lty = 1, length = 0.1, angle = 45) #vertical line at ld50 value
points(x,y, pch = 16, cex = 1) 

#add equation text
text(20,0.8,expression(y == frac(phi["1"],1+e^(frac(phi["2"]-x,phi["3"])))), cex = 1.5)
# /(1+exp((phi2-x)/phi3)))
dev.off()


# ------------------------------------------------------------------------------------- #
# export ld50 table
# ------------------------------------------------------------------------------------- #

# ld50 <- data.frame(ld50 = rbind(ld50_ag, ld50_nat), 
#                    pair = rep(pair_vec, 2), 
#                    lat = rep(lats, 2),
#                    long = rep(longs,2),
#                    env = c(rep("Ag", length(pair_vec)), rep("Nat", length(pair_vec))))
# colnames(ld50) <- c("ld50", "ld50_minus_stderr", "ld50_plus_stderr", "pair",   "lat",    "long",   "env")
#write.table(ld50, "ld50.txt", sep = "\t", 
#            col.names = T, row.names = F)



# ------------------------------------------------------------------------------------- #
# upload ld50 info that was created with above code
# ------------------------------------------------------------------------------------- #
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/4.phenotypes")
ld50 <- fread("ld50.txt") 

# ------------------------------------------------------------------------------------- #
# linear regression b/w LD50 and longitude
# ------------------------------------------------------------------------------------- #

#multiple linear regression model
lm_res <- lm(ld50 ~ env + long + lat, data = ld50)
summary(lm_res)

# Call:
#   lm(formula = ld50 ~ env + lat + long, data = ld50)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -2.9031 -0.7686  0.1204  0.8127  2.1297 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)   7.6730    23.2270   0.330   0.7446  
# envNat       -0.9251     0.5283  -1.751   0.0952 .
# lat          -0.2501     0.3681  -0.679   0.5047  
# long         -0.1534     0.1130  -1.357   0.1899  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.294 on 20 degrees of freedom
# Multiple R-squared:  0.3658,	Adjusted R-squared:  0.2707 
# F-statistic: 3.846 on 3 and 20 DF,  p-value: 0.02529

lm_res <- lm(ld50 ~ env + long , data = ld50)
summary(lm_res) #Multiple R-squared:  0.3512
#removing lat improves R^2, p value

#partial R2 for long
lm_res <- lm(ld50 ~ env  , data = ld50)
summary(lm_res) #multiple R2 = 0.09725

#partial R2 for env
lm_res <- lm(ld50 ~ long  , data = ld50)
summary(lm_res) #Multiple R-squared:  0.254


lm_res <- lm(ld50 ~ long , data = ld50)
summary(lm_res )

# > lm_res <- lm(ld50 ~ long , data = ld50)
# > summary(lm_res )
# 
# Call:
#   lm(formula = ld50 ~ long, data = ld50)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.5054 -0.5578  0.1470  0.8863  2.3460 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -7.90315    6.90045  -1.145    0.264  
# long        -0.21104    0.07712  -2.737    0.012 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 1.338 on 22 degrees of freedom
# Multiple R-squared:  0.254,	Adjusted R-squared:   0.22 
# F-statistic: 7.489 on 1 and 22 DF,  p-value: 0.01205

#Given that the model with both long and env has a better adjusted R-squared and a lower residual standard error, 
#I include env (despite its marginal significance) helps improve the model's explanatory power.
lm_res <- lm(ld50 ~ env + long, data = ld50)
# Generate predicted values bases on model
ld50$predicted_ld50 <- predict(lm_res, newdata = ld50)
ld50$long_cols <- color_values(ld50$long, palette = "purples")
pch_vec <- matrix(21,dim(ld50)[1],1)
pch_vec[ld50$env == "Nat"] <- 22

# ------------------------------------------------------------------------------------- #
# Figure 1 - panels C&D
# ------------------------------------------------------------------------------------- #


pdf("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/1.ld50_pca/F1C&D.pdf", 
    bg = "white", width=8, height=4, family = "Times New Roman")

par(family="Times New Roman", cex.lab=2, cex.axis=1.5)

layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))

plot(ld50$long, ld50$ld50, pch = pch_vec, xlab = "Longitude", ylab = expression(paste("Drought LD"[50], " (days)")), col = "#353535", bg =  ld50$long_cols, cex = 1.8, ylim = c(7,14))

#add grid
#abline(v = seq(-96,-84,by = 2) , col = "#353535", lwd = 0.5)
#abline(h = seq(7,14,by = 1) , col = "#353535", lwd = 0.5)
points(ld50$long[ld50$env == "Ag"][c(1,12)], ld50$predicted_ld50[ld50$env == "Ag"][c(1,12)], lty = 2, col = "#353535",  lwd = 2, type = "l") #Ag
lines(ld50$long[ld50$env == "Nat"], ld50$predicted_ld50[ld50$env == "Nat"], col = "#353535",  lwd = 2) #Nat
arrows(ld50$long, ld50$ld50_minus_stderr,  ld50$long, ld50$ld50_plus_stderr, code = 3, angle = 90, length = 0, col = "#353535")
points(ld50$long, ld50$ld50, pch = pch_vec, col = "#353535", bg = ld50$long_cols, cex = 1.8, ylim = c(7,14))

legend("bottomleft", legend = c("Ag", "Nat"), lty = c(2,1), col = "#353535", bty = "n", lwd = 2, cex = 1.5) #pt.bg = c("#e7298a", "#02818a")


# ------------------------------------------------------------------------------------- #
# LD50 as a function of Env (Ag/Nat) and colored by longitude
# ------------------------------------------------------------------------------------- #

par(family="Times New Roman", cex.lab=2, cex.axis=1.5)
ld50$long_cols <- color_values(ld50$long, palette = "purples")

jitter <- 0
store <- c()
#empty plot
plot(1, 
     xlim = c(0.5,2.6),  ylim=c(7.3, 14), 
     xlab = "Habitat", ylab = "LD50", xaxt = "n", yaxt = "n")
axis(1, at = c(1,2), labels = c("Ag", "Nat") )
#add pairs
for (p in levels(as.factor(ld50$pair)) ) {
  
  idx <- which(ld50$pair == p )
  lines(c(1 + jitter, 2+ jitter), ld50$ld50[idx], type = "b", pch = c(21, 22), cex = 1.8,
        bg = ld50$long_cols[idx], col = "#353535")
  jitter <- jitter + 0.01
  store <- c(store, ld50$ld50[idx][1] - ld50$ld50[idx][2])
  
}
mean_data <- mean(store)
# add box plots (summaries)
ag_idx <- which(ld50$env == "Ag")
nat_idx <- which(ld50$env == "Nat")
boxplot(ld50$ld50[ag_idx], add=TRUE, at = 0.85, boxwex = 0.15, border = "#505050", col = "#D3D3D3", yaxt = "n")
boxplot(ld50$ld50[nat_idx], add=TRUE, at = 2.25, boxwex = 0.15, border = "#505050", col = "#D3D3D3", yaxt = "n")


dev.off()

#Wilcoxon test for paired samples (nat/ag)
wilcox.test(ld50$ld50[ag_idx],ld50$ld50[nat_idx], paired = TRUE, alternative = "greater")
#if not paired
wilcox.test(ld50$ld50[ag_idx],ld50$ld50[nat_idx], alternative = "greater")



# Compare to random - randomize ENV
P <- 10000
means_rand <- matrix(NA, P)
for (i in 1 : P) {
  
  phenos_wmeta$Env_rand <- sample(phenos_wmeta$Env)
  store_rand <- matrix(NA, length(pair_vec))
  
  for (p in 1 : length(pair_vec)) {
    
    ag_vec <- phenos_wmeta$`Day of Full Wilt`[which(phenos_wmeta$Pair == pair_vec[p] & phenos_wmeta$Env_rand == "Ag" & phenos_wmeta$Treatment == "T")]
    nat_vec <- phenos_wmeta$`Day of Full Wilt`[which(phenos_wmeta$Pair == pair_vec[p] & phenos_wmeta$Env_rand == "Nat" & phenos_wmeta$Treatment == "T")]
    store_rand[p] <- mean(ag_vec, na.rm = T) - mean(nat_vec, na.rm = T)
    
  }
  
  means_rand[i] <- mean(store_rand)
  
}

pdf("permutations_env_ld50.pdf", bg = "transparent", width=4, height=3, family = "Times New Roman")

par(family="Times New Roman", cex.lab=1.3, cex.axis=1.3)
hist(means_rand, 50, xlab = "LD50 difference Ag/Nat", main = "")
abline(v=mean_data, col = "red", lwd = 2)
legend("topleft", legend = c("data", "random"), pch = 15, col = c("red", "gray"), bty = "n")

# main = "Difference of day to wilt is explained by environment"
box()

dev.off()


# ------------------------------------------------------------------------------------- #
# multiple linear regression and least square means
# ------------------------------------------------------------------------------------- #
# to be able to calculate least square means on the longitude, I need to create a factor that groups longitude together (three bins)

#linear model with longitude and env
lm_res <- lm(ld50 ~ env + long , data = ld50)
summary(lm_res) #multiple R2 = 0.3512

#linear model with longitude only
lm_res <- lm(ld50 ~ long , data = ld50)
summary(lm_res) #multiple R2 = 0.3512

#linear model
lm_res <- lm(ld50 ~ env + long + lat , data = ld50)
summary(lm_res) #multiple R2 = 0.3658

#partial R2 for longitude
0.3658 - 0.3074 = 0.0584


#least square means # computes the adjusted mean (controlled for other factors in the linear model)
#lsmin is to treat a quantitative factor as a qualitative one by making categories.
#considers continuous predictors as categorical predictors.
lsm <- lsmip(lm_res, env ~ long, at=list(long=c(-96, -92, -88, -84)), plot=F)
print(lsm)

#Calculate least square means for env together
lsm <- lsmip(lm_res, ~ long, at=list(long=c(-96, -84)), plot=F)
print(lsm)
12.0 + 1.96* 0.790 #for west
12.0 - 1.96* 0.790 #for west
10.1 + 1.96*0.666 #for east
10.1 - 1.96*0.666 #for east
#Calculate least square means for long together
lsm <- lsmip(lm_res, ~ env, at=list(env=c("Ag", "Nat")), plot=F)
print(lsm)
11.4 + 1.96*0.374 #for ag
11.4 - 1.96*0.374 #for ag
10.5 + 1.96*0.374 #for nat
10.5 - 1.96*0.374 #for nat



#least square mean in Ag
mean(lsm$yvar[c(1,3,5,7)]) #11.52 (mean)
mean(lsm$SE[c(1,3,5,7)]) #0.6075838 (sd)
#Nat 
mean(lsm$yvar[c(1,3,5,7)+1]) #10.59 (mean)

#calculating the confidence interval
confint(lm_res, level = 0.95)


#calculate stats
#zscore " (mean obs - mean exp)/sd exp
(mean_data - mean(means_rand)) / sd(means_rand)
#pvalue
length( which(means_rand >mean_data) ) / P # 0.003 - yes
# ------------------------------------------------------------------------------------- #
# Map LD50 values in space!! - interactive map
# ------------------------------------------------------------------------------------- #
library(leaflet)
#list of providers for maps here : https://leaflet-extras.github.io/leaflet-providers/preview/
ag_idx <- which(ld50$env == "Ag")
nat_idx <- which(ld50$env == "Nat")
ld50$long_cols <- color_values(ld50$long, palette = "purples")

# Upload herbarium samples
herb <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/5.herbarium/data/herbarium_samps.txt", sep = "\t", header = T)
herb_fil <- herb[herb$Year>=1880,]

herb_fil$cols <- colour_values(herb_fil$Year, palette = "magma")
min(herb_fil$Long, na.rm=T)

m <- leaflet() %>%
  addProviderTiles(providers$CartoDB.Positron)  %>% 
  setView(lng = min(herb_fil$Long, na.rm=T), lat = min(herb_fil$Lat, na.rm=T), zoom = 6)

# Jitter lat lon for Ag samples so that they are easier to see on the map
ld50$long_mod <- ld50$long
ld50$long_mod[ag_idx] <- ld50$long_mod[ag_idx] - 0.1
ld50$lat_mod <- ld50$lat
ld50$lat_mod[ag_idx] <- ld50$lat_mod[ag_idx] - 0.1

m <- addCircleMarkers(m, lng = ld50$long_mod[ag_idx], lat = ld50$lat_mod[ag_idx], weight = 3 , 
                      color = "gray", fill = TRUE, fillColor = ld50$long_cols[ag_idx], #fillColor = ld50$colors[ag_idx]
                      radius = 11, opacity = 1,   fillOpacity = 1) 

m <- addRectangles(m, lng1 = ld50$long[nat_idx], lat1 = ld50$lat[nat_idx], weight = 3 , 
                   lng2 = ld50$long[nat_idx]+0.5, lat2 = ld50$lat[nat_idx]+0.4, 
                      color = "gray", fill = TRUE, fillColor = ld50$long_cols[nat_idx] , opacity = 1,   fillOpacity = 1) #ld50$colors[nat_idx]
m

m <- addScaleBar(m, position = c( "bottomleft"), options = scaleBarOptions(imperial = FALSE))


m <- addCircleMarkers(m, lng = herb_fil$Long[herb_fil$Nat.Ag.Dist == "Ag"], lat = herb_fil$Lat[herb_fil$Nat.Ag.Dist == "Ag"], 
                      color = herb_fil$cols[herb_fil$Nat.Ag.Dist == "Ag"], radius = 5, opacity = 1, weight = 0,fill = TRUE, fillOpacity = 1)
m <- addRectangles(m, lng1 = herb_fil$Long[herb_fil$Nat.Ag.Dist == "Nat"], lat1 = herb_fil$Lat[herb_fil$Nat.Ag.Dist == "Nat"], weight = 3 , 
                   lng2 = herb_fil$Long[herb_fil$Nat.Ag.Dist == "Nat"]+0.2, lat2 = herb_fil$Lat[herb_fil$Nat.Ag.Dist == "Nat"]+0.15, 
                   color = herb_fil$cols[herb_fil$Nat.Ag.Dist == "Nat"], fill = TRUE, fillColor = herb_fil$cols[herb_fil$Nat.Ag.Dist == "Nat"], opacity = 1,   fillOpacity = 1) #ld50$colors[nat_idx]


m

# opacity = outer circle color strentgh, higher is stronger

# add legend - scales
#LD50 
pdf(file=paste("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/4.phenotypes/scale_herb.pdf", sep="" ), bg = "transparent", width=6, height=3, family = "Times New Roman")

par(family = "Times New Roman", cex = 1.3)
scale_data <- data.frame(a = 4, x = min(herb_fil$Year):max(herb_fil$Year))
scale_data$cols <- colour_values(scale_data$x, palette = "magma")
barplot(height = scale_data$a, col = scale_data$cols, space = 0, border = NA, names.arg = round(scale_data$x), yaxt='n', ann=FALSE) # 

dev.off()



# ------------------------------------------------------------------------------------- #
# calculate growth rates
# ------------------------------------------------------------------------------------- #

#how does stem width and plant height correlate? these were measured prior to drought
plot(phenos_wmeta$`Stem Width (mm)`, phenos_wmeta$`Plant Height (cm)`)

#subtract height/width after drought to before drought and compare between treatments
##stem width
phenos_wmeta$growth.rate.stem<-phenos_wmeta$`Stem Width 1W-Post-Drought` - phenos_wmeta$`Stem Width (mm)`
hist(phenos_wmeta$growth.rate.stem,breaks=20) #1 week growth rate

phenos_wmeta<-phenos_wmeta %>% mutate(stemrate_tilldeath=phenos_wmeta$`Stem Width at ""Death""` - phenos_wmeta$`Stem Width (mm)`) %>%
  mutate(stemrate_combined = coalesce(stemrate_tilldeath,growth.rate.stem))
hist(phenos_wmeta$stemrate_combined,breaks=20,)

##plant height
phenos_wmeta$growth_rate_height<-(phenos_wmeta$`Plant Height 1W-Post-Drought` - phenos_wmeta$`Plant Height (cm)` )/7
hist(phenos_wmeta$growth_rate_height,breaks=100) #1 week growth rate, early life history

phenos_wmeta$growth_rate_height_2<-(phenos_wmeta$`Plant Height 2W` - phenos_wmeta$`Plant Height 1W-Post-Drought`)/7
hist(phenos_wmeta$growth_rate_height_2,breaks=100) #2 week growth rate, mid life history 

phenos_wmeta <- phenos_wmeta %>% mutate(heightrate_tilldeath=(phenos_wmeta$`Plant Height at Death`-phenos_wmeta$`Plant Height (cm)`)/7) %>%
  mutate(heightrate_combined = coalesce(heightrate_tilldeath,growth_rate_height)) #combine to get more complete information (note measured across diff times, perhaps not best approach)
hist(phenos_wmeta$heightrate_combined,breaks=20)


# ------------------------------------------------------------------------------------- #
# multiple linear regression w/ growth rate
# ------------------------------------------------------------------------------------- #

# linear regression
summary(lm(growth_rate_height ~ Env + Long, data = phenos_wmeta)) #not significant
summary(lm(growth_rate_height ~ Long , data = phenos_wmeta)) #marginally significant


summary(lm(growth_rate_height_2 ~ Env +  Long, data = phenos_wmeta)) #marginally significant
summary(lm(growth_rate_height_2 ~ Long , data = phenos_wmeta)) #marginally significant


lm_res <- lm(growth_rate_height_2 ~ Env +  Long, data = phenos_wmeta)
# Generate predicted values bases on model
phenos_wmeta$predicted_growth <- predict(lm_res, newdata = phenos_wmeta)

pdf("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/4.phenotypes/growth_rate_longitude.pdf", 
    bg = "white", width=6, height=5, family = "Times New Roman")

par(family="Times New Roman", cex.lab=1.5, cex.axis=1.5)

phenos_wmeta$env_cols <- "#e7298a"
phenos_wmeta$env_cols[phenos_wmeta$Env == "Nat"] <- "#4eb3d3"

plot(phenos_wmeta$Long, phenos_wmeta$growth_rate_height_2, pch =21, xlab = "Longitude", ylab = "Growth rate", ylim = c(0,8))
#add grid
abline(v = seq(-96,-84,by = 2) , col = "gray", lwd = 0.5)
abline(h = seq(0,8,by = 1) , col = "gray", lwd = 0.5)
lines(phenos_wmeta$Long[phenos_wmeta$Env == "Ag"], phenos_wmeta$predicted_growth[phenos_wmeta$Env == "Ag"], col = "#ce1256",  lwd = 2) #Ag
lines(phenos_wmeta$Long[phenos_wmeta$Env == "Nat"], phenos_wmeta$predicted_growth[phenos_wmeta$Env == "Nat"], col = "#02818a",  lwd = 2) #Nat
points(phenos_wmeta$Long, phenos_wmeta$growth_rate_height_2, pch =21, col = "gray", bg = phenos_wmeta$env_cols, cex = 1.2)

legend("topleft", legend = c("Ag", "Nat"), pch = 21, pt.bg = c("#e7298a", "#02818a"), col = "gray", bty = "n")

dev.off()


# ------------------------------------------------------------------------------------- #
# Trade-off survival to drought (LD50) and growth
# ------------------------------------------------------------------------------------- #









