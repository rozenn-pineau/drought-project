#!/usr/bin/env Rscript

library(dplyr)
library(extrafont)
library(geosphere)
library(data.table)
library(car)
library(scales)
library(paletteer)
library(lsmeans)
library(ggplot2)

rm(list=ls())

setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/5.herbarium/data/54samples/")
tmax <- read.table("tmax_NOAA.txt", sep = "\t", header = T)
n_snp <- dim(tmax[,-c(1:6)])[2]

# Define parameters
max_distance_km <- 100      # Maximum distance between pairs
genomic_cols <- 7:39        # Columns containing genomic site data
n_perm <- 100
n_ind <- dim(tmax)[1]

# ================================
# FUNCTION: SNP-SPECIFIC PAIR SELECTION
# ================================

select_pairs_for_snp <- function(clean_data, snp_col, min_pairs = 3) {
  # Calculate all pairwise combinations
  pairs <- expand.grid(i = 1:nrow(clean_data), j = 1:nrow(clean_data)) %>%
    filter(i < j) %>%
    mutate(
      sample1 = clean_data$sample[i],
      sample2 = clean_data$sample[j],
      year1 = clean_data$year[i],
      year2 = clean_data$year[j],
      lon_sample1 = clean_data$samp_lon[i],
      lon_sample2 = clean_data$samp_lon[j],
      lat_sample1 = clean_data$samp_lat[i],
      lat_sample2 = clean_data$samp_lat[j],
      tmax_sample1 = clean_data$tmax[i],
      tmax_sample2 = clean_data$tmax[j],
      GT_samp_1 = clean_data[i,snp_col],
      GT_samp_2 = clean_data[j,snp_col],
      year_diff = abs(clean_data$year[i] - clean_data$year[j]),
      distance_km = distVincentyEllipsoid(
        cbind(clean_data$samp_lon[i], clean_data$samp_lat[i]),
        cbind(clean_data$samp_lon[j], clean_data$samp_lat[j])
      ) / 1000
    )
  
  # Filter based on criteria AND SNP-specific data availability
  valid_pairs <- pairs %>%
    filter(
      year_diff > 0,                                    # Different years
      distance_km <= max_distance_km,                   # Within distance threshold
      !is.na(clean_data[i, snp_col]),                  # Sample i has data for this SNP
      !is.na(clean_data[j, snp_col])                   # Sample j has data for this SNP
    ) %>%
    arrange(distance_km)  # Prioritize closer pairs
  
  # If we don't have enough pairs, return empty
  if(nrow(valid_pairs) < min_pairs) {
    return(data.frame())
  }
  
  # Greedy selection for non-overlapping pairs by snp
  selected <- data.frame()
  used_samples <- c()
  
  # Loop through each valid pair
  for(k in 1:nrow(valid_pairs)) {
    # Check if both samples in this pair are unused
    if(!valid_pairs$i[k] %in% used_samples & !valid_pairs$j[k] %in% used_samples) {
      # Select this pair
      selected <- rbind(selected, valid_pairs[k,])
      # Mark both samples as used
      used_samples <- c(used_samples, valid_pairs$i[k], valid_pairs$j[k])
    }
  }
  
  return(selected)
}

# ================================
# FUNCTION: SNP analysis based on optimal pairs calculated site-by-site
# ================================

analyze_all_snps_comprehensive <- function(clean_data, use_randomized_temp = FALSE) {
  snp_correlations <- c()
  snp_pvalues <- c()
  snp_sample_sizes <- c()
  year_diff_vec <- c()
  
  if(!use_randomized_temp) {
    cat("Analyzing", n_snp, "SNPs with observed temperature data...\n")
  }
  
  for(snp_idx in 1:n_snp) {
    snp_col <- genomic_cols[snp_idx]
    
    # Get optimal pairs for this SNP
    pairs_for_snp <- select_pairs_for_snp(clean_data, snp_col)
    
    if(nrow(pairs_for_snp) == 0) {
      if(!use_randomized_temp) {
        cat("  SNP", snp_idx, ": No valid pairs found\n")
      }
      snp_correlations <- c(snp_correlations, NA)
      snp_pvalues <- c(snp_pvalues, NA)
      snp_sample_sizes <- c(snp_sample_sizes, 0)
      next
    }
    
    if(!use_randomized_temp) {
      cat("  SNP", snp_idx, ": Found", nrow(pairs_for_snp), "pairs")
    }
    
    snp_sample_sizes <- c(snp_sample_sizes, nrow(pairs_for_snp))
    
    # Calculate temperature and allele frequency changes
    temp_changes <- c()
    af_changes <- c()
    
    for(p in 1:nrow(pairs_for_snp)) {
      i <- pairs_for_snp$i[p]
      j <- pairs_for_snp$j[p]
      
      year_diff <- clean_data$year[i] - clean_data$year[j]
      year_diff_vec <- c(year_diff_vec, abs(year_diff))
      
      if(year_diff < 0) {
        if(use_randomized_temp) {
          temp_change <- clean_data$tmax_rdn[j] - clean_data$tmax_rdn[i]
        } else {
          temp_change <- clean_data$tmax[j] - clean_data$tmax[i]
        }
        af_change <- clean_data[j, snp_col] - clean_data[i, snp_col]
      } else {
        if(use_randomized_temp) {
          temp_change <- clean_data$tmax_rdn[i] - clean_data$tmax_rdn[j]
        } else {
          temp_change <- clean_data$tmax[i] - clean_data$tmax[j]
        }
        af_change <- clean_data[i, snp_col] - clean_data[j, snp_col]
      }
      
      temp_changes <- c(temp_changes, temp_change)
      af_changes <- c(af_changes, af_change)
    }
    
    # Calculate correlation between temperature and allele frequency changes
    tryCatch({
      # Check if both variables have variation (not all identical values)
      if(var(temp_changes, na.rm = TRUE) > 0 && var(af_changes, na.rm = TRUE) > 0) {
        # Calculate Pearson correlation and significance test
        cor_result <- cor.test(temp_changes, af_changes, method = "spearman",exact=F)
        correlation <- cor_result$estimate
        p_value <- cor_result$p.value
        
        # Store results in vectors
        snp_correlations <- c(snp_correlations, correlation)
        snp_pvalues <- c(snp_pvalues, p_value)
        
        # Print results (unless using randomized temperature data)
        if(!use_randomized_temp) {
          if(p_value < 0.05 && correlation > 0) {
            cat(" → SIGNIFICANT positive correlation:", round(correlation, 3), "(p =", round(p_value, 4), ")\n")
          } else if(correlation > 0) {
            cat(" → Positive correlation:", round(correlation, 3), "(p =", round(p_value, 3), ")\n")
          } else {
            cat(" → Negative correlation:", round(correlation, 3), "(p =", round(p_value, 3), ")\n")
          }
        }
      } else {
        # No variation in data - correlation is undefined
        snp_correlations <- c(snp_correlations, NA)
        snp_pvalues <- c(snp_pvalues, NA)
        if(!use_randomized_temp) {
          cat(" → No variation in data\n")
        }
      }
    }, error = function(e) {
      # Handle analysis failures
      snp_correlations <- c(snp_correlations, NA)
      snp_pvalues <- c(snp_pvalues, NA)
      if(!use_randomized_temp) {
        cat(" → Analysis failed\n")
      }
    })
  }
  
  # Remove NAs for analysis (moved outside the for loop)
  valid_indices <- !is.na(snp_correlations)
  correlations_clean <- snp_correlations[valid_indices]
  pvalues_clean <- snp_pvalues[valid_indices]
  
  if(!use_randomized_temp) {
    cat("\nSummary:\n")
    cat("  Total SNPs with data:", sum(valid_indices), "\n")
    cat("  Significant positive correlations:", sum(pvalues_clean < 0.05 & correlations_clean > 0, na.rm = TRUE), "\n")
    cat("  Total positive correlations:", sum(correlations_clean > 0, na.rm = TRUE), "\n")
    cat("  Mean correlation:", round(mean(correlations_clean, na.rm = TRUE), 3), "\n")
    cat("  Pairs per SNP - Range:", min(snp_sample_sizes[snp_sample_sizes > 0]), "-", max(snp_sample_sizes),
        ", Mean:", round(mean(snp_sample_sizes[snp_sample_sizes > 0]), 1), "\n")
    cat("  Mean year difference:", round(mean(year_diff_vec, na.rm = TRUE), 3), "\n")
  }
  
  return(list(
    correlations = snp_correlations,
    correlations_clean = correlations_clean,
    pvalues = snp_pvalues,
    pvalues_clean = pvalues_clean,
    sample_sizes = snp_sample_sizes,
    n_analyzed = sum(valid_indices),
    pairs_for_snp = pairs_for_snp,
    
    # Test Statistics
    stat_n_significant = sum(pvalues_clean < 0.05 & correlations_clean > 0, na.rm = TRUE),
    stat_n_positive = sum(correlations_clean > 0, na.rm = TRUE),  # ADD THIS LINE
    stat_sum_positive_cors = sum(correlations_clean[correlations_clean > 0], na.rm = TRUE),
    stat_median_correlation = median(correlations_clean, na.rm = TRUE),
    stat_mean_correlation = mean(correlations_clean, na.rm = TRUE),
    stat_mean_year_diff = round(mean(year_diff_vec, na.rm = TRUE), 3),
    stat_mean_year_diff_vec = year_diff_vec
  ))
    
}

# ================================
# STEP 1: ANALYZE OBSERVED DATA
# ================================

cat("Analyzing observed data...\n")

clean_data <- tmax %>%
  filter(!is.na(sample), !is.na(samp_lon))

observed_results <- analyze_all_snps_comprehensive(clean_data, use_randomized_temp = FALSE)

cat("Observed results:\n")
cat("  SNPs analyzed:", observed_results$n_analyzed, "\n")
cat("  Significant positive (p<0.05):", observed_results$stat_n_significant, "\n")
cat("  Sum of positive correlations:", round(observed_results$stat_sum_positive_cors, 3), "\n")
cat("  Median correlation:", round(observed_results$stat_median_correlation, 3), "\n")
cat("  Mean correlation:", round(observed_results$stat_mean_correlation, 3), "\n")
cat("  Mean temporal diff:", round(observed_results$stat_mean_year_diff, 3), "\n")

mean(observed_results$stat_mean_year_diff_vec[1:54])
sd(observed_results$stat_mean_year_diff_vec)/sqrt(54)
# ================================
# STEP 2: PERMUTATION TESTS
# ================================

cat("\nRunning", n_perm, "permutations...\n")

# Storage for different test statistics
null_sum_positive_cors <- c()
null_median_correlations <- c()
null_mean_correlations <- c()
null_all_correlations <- list()
null_n_significant <- c()        
null_n_positive <- c()          

for(perm in 1:n_perm) {
  clean_data$tmax_rdn <- clean_data$tmax[sample(1:nrow(clean_data))]
  
  null_result <- analyze_all_snps_comprehensive(clean_data, use_randomized_temp = TRUE)
  
  null_sum_positive_cors <- c(null_sum_positive_cors, null_result$stat_sum_positive_cors)
  null_median_correlations <- c(null_median_correlations, null_result$stat_median_correlation)
  null_mean_correlations <- c(null_mean_correlations, null_result$stat_mean_correlation)
  null_n_significant <- c(null_n_significant, null_result$stat_n_significant)     
  null_n_positive <- c(null_n_positive, sum(null_result$correlations_clean > 0, na.rm = TRUE)) 
  
  null_all_correlations[[perm]] <- null_result$correlations_clean
  
  if(perm %% 100 == 0) cat("Completed", perm, "permutations\n")
}

# ================================
# STEP 3: STATISTICAL TESTS
# ================================
cat("\n=== COMPREHENSIVE STATISTICAL TESTS ===\n")

# Extract observed correlations for testing
observed_cors <- observed_results$correlations_clean
pooled_null_cors <- unlist(null_all_correlations)

# One-sample test: are observed correlations significantly > 0?
if(length(observed_cors) > 0) {
  wilcox_observed <- wilcox.test(observed_cors, mu = 0, alternative = "greater")
  cat("1. Wilcoxon test - observed correlations vs 0:\n")
  cat("   P-value:", round(wilcox_observed$p.value, 4), "\n")
  
  # Compare to null expectation
  null_wilcox_pvals <- c()
  for(perm in 1:n_perm) {
    if(length(null_all_correlations[[perm]]) > 0) {
      null_wilcox <- wilcox.test(null_all_correlations[[perm]], mu = 0, 
                                 alternative = "greater", exact = FALSE)
      null_wilcox_pvals <- c(null_wilcox_pvals, null_wilcox$p.value)
    }
  }
  # How often do null data give p-values as extreme as observed?
  p7 <- mean(null_wilcox_pvals <= wilcox_observed$p.value, na.rm = TRUE)
  cat("   Null-adjusted p-value:", round(p7, 4), "\n")
}

# Test 2: Sum of positive correlations
p3 <- mean(null_sum_positive_cors >= observed_results$stat_sum_positive_cors)
cat("2. Sum of positive correlations:\n")
cat("   Observed:", round(observed_results$stat_sum_positive_cors, 3), 
    "| Null mean:", round(mean(null_sum_positive_cors), 3), 
    "| P-value:", round(p3, 4), "\n")

# Test 3: Median correlation
p6 <- mean(null_median_correlations >= observed_results$stat_median_correlation)
cat("3. Median correlation:\n")
cat("   Observed:", round(observed_results$stat_median_correlation, 3), 
    "| Null mean:", round(mean(null_median_correlations), 3), 
    "| P-value:", round(p6, 4), "\n")

# Test 3: Median correlation
px <- mean(null_mean_correlations >= observed_results$stat_mean_correlation)
cat("3. Mean correlation:\n")
cat("   Observed:", round(observed_results$stat_mean_correlation, 3), 
    "| Null mean:", round(mean(null_mean_correlations), 3), 
    "| P-value:", round(px, 4), "\n")

# Test 5: Number of significant correlations
p_sig <- mean(null_n_significant >= observed_results$stat_n_significant)
cat("5. Number of significant correlations:\n")
cat("   Observed:", observed_results$stat_n_significant, 
    "| Null mean:", round(mean(null_n_significant), 1), 
    "| P-value:", round(p_sig, 4), "\n")

# Test 6: Number of positive correlations (all, not just significant)
observed_n_positive <- sum(observed_results$correlations_clean > 0, na.rm = TRUE)
p_pos <- mean(null_n_positive >= observed_n_positive)
cat("6. Number of positive correlations (all):\n")
cat("   Observed:", observed_n_positive, 
    "| Null mean:", round(mean(null_n_positive), 1), 
    "| P-value:", round(p_pos, 4), "\n")


# ================================
# STEP 4: VISUALIZATION
# ================================


# Prepare data with permutation IDs
cors_plot_data_detailed <- data.frame(
  Correlation = c(observed_results$correlations_clean, unlist(null_all_correlations)),
  Type = c(rep("Observed", length(observed_results$correlations_clean)),
           rep(paste0("Null_", rep(1:n_perm, sapply(null_all_correlations, length))), 
               lengths = sapply(null_all_correlations, length))),
  Permutation = c(rep(0, length(observed_results$correlations_clean)),
                  rep(1:n_perm, sapply(null_all_correlations, length)))
)

p2_plot2 <- ggplot(cors_plot_data_detailed, aes(x = Correlation)) +
  # Show each null permutation as a faint line
  geom_density(data = subset(cors_plot_data_detailed, Permutation > 0),
               aes(group = Permutation), alpha = 0.05, color = "gray", size = 0.3) +
  # Show observed as bold line
  geom_density(data = subset(cors_plot_data_detailed, Permutation == 0),
               color = "#E31A1C", size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", alpha = 0.7) +
  labs(
    subtitle = paste("Wilcoxon p-value:", round(wilcox_observed$p.value, 4), 
                     "\nNull-adjusted p-value:", round(p7, 4)),
    x = "Correlation coefficient",
    y = "Density"
  ) +
  theme_minimal() + 
  theme(text=element_text(size=12,  family="Times New Roman"),    
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size=8))

p2_plot2

# Prepare data for sum of positive correlations
sum_positive_data <- data.frame(
  Value = c(observed_results$stat_sum_positive_cors, null_sum_positive_cors),
  Type = c("Observed", rep("Null", length(null_sum_positive_cors)))
)

p3_plot <- ggplot(sum_positive_data[sum_positive_data$Type == "Null",], aes(x = Value, fill = Type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 30) +
  geom_density(aes(color = Type), alpha = 0.8, size = 1.2) +
  scale_fill_manual(values = c("Observed" = "#E31A1C", "Null" = "gray")) +
  scale_color_manual(values = c("Observed" = "#E31A1C", "Null" = "gray")) +
  geom_vline(aes(xintercept = observed_results$stat_sum_positive_cors), 
             color = "#E31A1C", size = 1) +
  labs(
   # title = "Distribution of sum of positive correlations",
    subtitle = paste("\nPermuted p-value:", round(p3, 4), "| Observed:", round(observed_results$stat_sum_positive_cors, 3)),
    x = "Sum of positive correlations",
    y = "Density"
  ) +
  theme_minimal() +
  theme(text=element_text(size=12,  family="Times New Roman"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size=8)
  ) 
p3_plot


# Prepare data for median correlation
median_data <- data.frame(
  Value = c(observed_results$stat_median_correlation, null_median_correlations),
  Type = c("Observed", rep("Null", length(null_median_correlations)))
)

p4_plot <- ggplot(median_data[median_data$Type == "Null",], aes(x = Value, fill = Type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 30) +
  geom_density(aes(color = Type), alpha = 0.8, size = 1.2) +
  scale_fill_manual(values = c("Observed" = "#E31A1C", "Null" = "gray")) +
  scale_color_manual(values = c("Observed" = "#E31A1C", "Null" = "gray")) +
  geom_vline(aes(xintercept = observed_results$stat_median_correlation), 
             color = "#E31A1C", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", alpha = 0.7) +
  labs(
   # title = "Distribution of median correlations",
    subtitle = paste("\nPermuted p-value:", round(p6, 4), "| Observed:", round(observed_results$stat_median_correlation, 3)),
    x = "Median correlation",
    y = "Density"
  ) +
  theme_minimal() +
  theme(text=element_text(size=12,  family="Times New Roman"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size=8),
    legend.position = "none"
  )

p4_plot

#Prep Bar plot for visulalizing number of sig positive coresamp_lations
sig_counts <- table(null_n_significant)
sig_props <- prop.table(sig_counts)

bar_data <- data.frame(
  Count = as.numeric(names(sig_counts)),
  Proportion = as.numeric(sig_props),
  Frequency = as.numeric(sig_counts)
)

p_sig_bar <- ggplot(bar_data, aes(x = Count, y = Proportion)) +
  geom_col(fill = "gray", alpha = 0.7, width = 0.8) +
  geom_vline(xintercept = observed_results$stat_n_significant, 
             color = "#E31A1C", size = 1, linetype = "solid") +
  scale_x_continuous(breaks = 0:max(null_n_significant), 
                     limits = c(-0.5, max(c(null_n_significant, observed_results$stat_n_significant)) + 0.5)) +
  labs(
    subtitle = paste("\nPermuted p-value:", round(p_sig, 4), "| Observed:", observed_results$stat_n_significant),
    x = "Number of significant correlations",
    y = "Proportion of permutations"
  ) +
  theme_minimal() +
  theme(text=element_text(size=12,  family="Times New Roman"),
    plot.subtitle = element_text(size=8),
    panel.grid.minor.x = element_blank()
  )


#lemon::grid_arrange_shared_legend(p2_plot2,p4_plot,p_sig_bar, nrow=1)
library(cowplot)

pdf(file="/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/5.herbarium/tracking_permutations.pdf", 
    bg = "transparent", width=7, height=2.5, family = "Times New Roman")


plot_grid(p2_plot2,p4_plot,p_sig_bar, nrow=1)

dev.off()


#locus at position 10 in the dataset has the highest correlation
sort(observed_results$correlations_clean)
which(observed_results$correlations_clean >= 0.5 & observed_results$correlations_clean <= 0.6)

pairs <- select_pairs_for_snp(clean_data, snp_col = 15)


# ================================
# Composite plot
# ================================

# ------------------------------------------------- map
#for sample 1 is earlier
pairs$pch_samp1 <- 21
pairs$pch_samp2 <- 23
#if sample 2 is earlier than sample 1
pairs$pch_samp1[pairs$year1 >= pairs$year2] <- 21 #sample 2 is earlier
pairs$pch_samp2[pairs$year1 > pairs$year2] <- 23 #sample 2 is earlier

#color by temperature
n <- length(pairs$tmax_sample1)*2
cols <- colour_values(c(pairs$tmax_sample1,pairs$tmax_sample2), palette = "magma")
pairs$col_samp1 <- cols[1:length(pairs$tmax_sample1)]
pairs$col_samp2 <- cols[(length(pairs$tmax_sample1)+1):n]

pdf(file=paste("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/5.herbarium/Figure5A.pdf", sep="" ), 
    bg = "transparent", width=4, height=8, family = "Times New Roman")

par(family="Times New Roman", cex.axis = 1, cex.lab = 1.5, mfrow = c(2,1))
plot(0, xlim = c(min(pairs$lon_sample1), max(pairs$lon_sample1)), ylim = c(min(pairs$lat_sample1), max(pairs$lat_sample1)), 
     xlab = "Longitude (\u00B0)", ylab = "Latitude (\u00B0)", pch = "")
grid()
for (p in 1:length(pairs$lon_sample1)){
  
  lines(c(pairs$lon_sample1[p], pairs$lon_sample2[p]),
       c(pairs$lat_sample1[p], pairs$lat_sample2[p]), lty = 1)
  
  points(c(pairs$lon_sample1[p], pairs$lon_sample2[p]),
        c(pairs$lat_sample1[p], pairs$lat_sample2[p]), 
        pch = c(pairs$pch_samp1[p], pairs$pch_samp2[p]),
        bg = c(pairs$col_samp1[p], pairs$col_samp2[p]),
        col = "#353535",
        cex = 1.5)

}

# ------------------------------------------------- genotype

#color by distance
pairs$col_dist <- colour_values(pairs$distance_km, palette = "viridis")
pairs$tmax_change <- NA
pairs$af_change <- NA
pairs$pch_size <- NA

#calculate AF change and temperature change
for (p in 1:length(pairs$lon_sample1)){
  
  #if year sample 1 is after year sample 2
  if (pairs$year1[p] > pairs$year2[p]) {#year 1 - year 2
    #tmax difference
    pairs$tmax_change <- pairs$tmax_sample1 - pairs$tmax_sample2
    #AF difference
    pairs$af_change <- pairs$GT_samp_1 - pairs$GT_samp_2

  }
  else {
    #tmax difference
    pairs$tmax_change <- pairs$tmax_sample2 - pairs$tmax_sample1
    #AF difference
    pairs$af_change <- pairs$GT_samp_2 - pairs$GT_samp_1
  }
  
  #how many years apart? coded by pch size
  if (pairs$year_diff[p] < 20) { pairs$pch_size[p] <- 1 }
  else if (pairs$year_diff[p] >= 20 & pairs$year_diff[p] < 40) {pairs$pch_size[p] <- 1.3}
  else if (pairs$year_diff[p] >= 40 & pairs$year_diff[p] < 60) {pairs$pch_size[p] <- 1.6}
  else {pairs$pch_size[p] <- 2}
  
}


par(family="Times New Roman", cex.axis = 1, cex.lab = 1.5)
plot(0, ylim = c(-.5,.5), 
     xlim = c(-4,4.5), 
     xlab = "Temperature change (\u00B0C)", ylab = "Allele frequency change", pch = "")
grid()

#plot
for (p in 1:length(pairs$lon_sample1)){
  
  points(pairs$tmax_change, pairs$af_change,
         pch = 21,
         bg = pairs$col_dist,
         col = "#353535",
         cex = pairs$pch_size)

}

#add curve
lm_res <- lm(pairs$af_change ~ pairs$tmax_change)
lm.res <- summary(lm_res)

polygon(c(c(-4,4.5),rev(c(-4,4.5))),
        c(lm.res$coefficients[1]-lm.res$coefficients[3] + (lm.res$coefficients[2]-lm.res$coefficients[4])*c(-4,4.5),
          rev(lm.res$coefficients[1]+lm.res$coefficients[3] + (lm.res$coefficients[2]+lm.res$coefficients[4])*c(-4,4.5))),
        col = scales::alpha("#959595", 0.15), border = FALSE)
lines(c(-4,4.5), lm.res$coefficients[1] + lm.res$coefficients[2]*c(-4,4.5), lty = 2, lwd = 2, col = "red")

dev.off()



#Export the legends
pdf(file=paste("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/figures/5.herbarium/Legend_dist.pdf", sep="" ), 
    bg = "transparent", width=3, height=4, family = "Times New Roman")

par(family= "Times New Roman", cex.axis = 2, cex.lab = 2.5)

#Export legend for precipitation
lg <- color_values(c(min(pairs$dist):max(pairs$dist)), palette = "viridis")
legend_image <- as.raster(matrix(lg , ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'Max. temperature')
text(x=1.5, y = seq(0,1), labels = c("0 km","99 km"))
rasterImage(legend_image, 0, 0, 1,1)

dev.off()

