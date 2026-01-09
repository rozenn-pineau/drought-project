#!/usr/bin/env Rscript
rm(list=ls())

library(dplyr)
library(extrafont)
library(geosphere)
library(data.table)
library(car)
library(scales)
library(paletteer)
library(lsmeans)
library(ggplot2)

#Upload precipitation data
prcp <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/1.DroughtProject/1.analyses/data/5.herbarium/merged/prcp/prcp_NOAA.txt", sep = "\t", header = T)
prcp$prcp <- prcp$prcp/10 #change from tenth of mm to mm
n_snp <- dim(prcp[,-c(1:6)])[2]

# Define parameters
max_distance_km <- 100      # Maximum distance between pairs
genomic_cols <- 7:39        # Columns containing genomic site data
n_perm <- 100
n_ind <- dim(prcp)[1]

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

analyze_all_snps_comprehensive <- function(clean_data, use_randomized_prcp = FALSE) {
  snp_correlations <- c()
  snp_pvalues <- c()
  snp_sample_sizes <- c()
  
  if(!use_randomized_prcp) {
    cat("Analyzing", n_snp, "SNPs with observed precip data...\n")
  }
  
  for(snp_idx in 1:n_snp) {
    snp_col <- genomic_cols[snp_idx]
    
    # Get optimal pairs for this SNP
    pairs_for_snp <- select_pairs_for_snp(clean_data, snp_col)
    
    if(nrow(pairs_for_snp) == 0) {
      if(!use_randomized_prcp) {
        cat("  SNP", snp_idx, ": No valid pairs found\n")
      }
      snp_correlations <- c(snp_correlations, NA)
      snp_pvalues <- c(snp_pvalues, NA)
      snp_sample_sizes <- c(snp_sample_sizes, 0)
      next
    }
    
    if(!use_randomized_prcp) {
      cat("  SNP", snp_idx, ": Found", nrow(pairs_for_snp), "pairs")
    }
    
    snp_sample_sizes <- c(snp_sample_sizes, nrow(pairs_for_snp))
    
    # Calculate precipitation and allele frequency changes
    prcp_changes <- c()
    af_changes <- c()
    
    for(p in 1:nrow(pairs_for_snp)) {
      i <- pairs_for_snp$i[p]
      j <- pairs_for_snp$j[p]
      
      year_diff <- clean_data$year[i] - clean_data$year[j]
      
      if(year_diff < 0) {
        if(use_randomized_prcp) {
          prcp_change <- clean_data$prcp_rdn[j] - clean_data$prcp_rdn[i]
        } else {
          prcp_change <- clean_data$prcp[j] - clean_data$prcp[i]
        }
        af_change <- clean_data[j, snp_col] - clean_data[i, snp_col]
      } else {
        if(use_randomized_prcp) {
          prcp_change <- clean_data$prcp_rdn[i] - clean_data$prcp_rdn[j]
        } else {
          prcp_change <- clean_data$prcp[i] - clean_data$prcp[j]
        }
        af_change <- clean_data[i, snp_col] - clean_data[j, snp_col]
      }
      
      prcp_changes <- c(prcp_changes, prcp_change)
      af_changes <- c(af_changes, af_change)
    }
    
    # Calculate correlation between precipitation and allele frequency changes
    tryCatch({
      # Check if both variables have variation (not all identical values)
      if(var(prcp_changes, na.rm = TRUE) > 0 && var(af_changes, na.rm = TRUE) > 0) {
        # Calculate Pearson correlation and significance test
        cor_result <- cor.test(prcp_changes, af_changes, method = "spearman",exact=F)
        correlation <- cor_result$estimate
        p_value <- cor_result$p.value
        
        # Store results in vectors
        snp_correlations <- c(snp_correlations, correlation)
        snp_pvalues <- c(snp_pvalues, p_value)
        
        # Print results (unless using randomized precipitation data)
        if(!use_randomized_prcp) {
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
        if(!use_randomized_prcp) {
          cat(" → No variation in data\n")
        }
      }
    }, error = function(e) {
      # Handle analysis failures
      snp_correlations <- c(snp_correlations, NA)
      snp_pvalues <- c(snp_pvalues, NA)
      if(!use_randomized_prcp) {
        cat(" → Analysis failed\n")
      }
    })
  }
  
  # Remove NAs for analysis (moved outside the for loop)
  valid_indices <- !is.na(snp_correlations)
  correlations_clean <- snp_correlations[valid_indices]
  pvalues_clean <- snp_pvalues[valid_indices]
  
  if(!use_randomized_prcp) {
    cat("\nSummary:\n")
    cat("  Total SNPs with data:", sum(valid_indices), "\n")
    cat("  Significant positive correlations:", sum(pvalues_clean < 0.05 & correlations_clean > 0, na.rm = TRUE), "\n")
    cat("  Total positive correlations:", sum(correlations_clean > 0, na.rm = TRUE), "\n")
    cat("  Mean correlation:", round(mean(correlations_clean, na.rm = TRUE), 3), "\n")
    cat("  Pairs per SNP - Range:", min(snp_sample_sizes[snp_sample_sizes > 0]), "-", max(snp_sample_sizes),
        ", Mean:", round(mean(snp_sample_sizes[snp_sample_sizes > 0]), 1), "\n")
  }
  
  return(list(
    correlations = snp_correlations,
    correlations_clean = correlations_clean,
    pvalues = snp_pvalues,
    pvalues_clean = pvalues_clean,
    sample_sizes = snp_sample_sizes,
    n_analyzed = sum(valid_indices),
    
    # Test Statistics
    stat_n_significant = sum(pvalues_clean < 0.05 & correlations_clean > 0, na.rm = TRUE),
    stat_n_positive = sum(correlations_clean > 0, na.rm = TRUE),  # ADD THIS LINE
    stat_sum_positive_cors = sum(correlations_clean[correlations_clean > 0], na.rm = TRUE),
    stat_median_correlation = median(correlations_clean, na.rm = TRUE),
    stat_mean_correlation = mean(correlations_clean, na.rm = TRUE)          
  ))
    
}

# ================================
# STEP 1: ANALYZE OBSERVED DATA
# ================================

cat("Analyzing observed data...\n")

clean_data <- prcp %>%
  filter(!is.na(samp_lat), !is.na(samp_lon))

observed_results <- analyze_all_snps_comprehensive(clean_data, use_randomized_prcp = FALSE)

cat("Observed results:\n")
cat("  SNPs analyzed:", observed_results$n_analyzed, "\n")
cat("  Significant positive (p<0.05):", observed_results$stat_n_significant, "\n")
cat("  Sum of positive correlations:", round(observed_results$stat_sum_positive_cors, 3), "\n")
cat("  Median correlation:", round(observed_results$stat_median_correlation, 3), "\n")
cat("  Mean correlation:", round(observed_results$stat_mean_correlation, 3), "\n")

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
  clean_data$prcp_rdn <- clean_data$prcp[sample(1:nrow(clean_data))]
  
  null_result <- analyze_all_snps_comprehensive(clean_data, use_randomized_prcp = TRUE)
  
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
               aes(group = Permutation), alpha = 0.05, color = "#1F78B4", size = 0.3) +
  # Show observed as bold line
  geom_density(data = subset(cors_plot_data_detailed, Permutation == 0),
               color = "#E31A1C", size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", alpha = 0.7) +
  labs(
    subtitle = paste("Wilcoxon p-value:", round(wilcox_observed$p.value, 4), 
                     "\nNull-adjusted p-value:", round(p7, 4)),
    x = "Correlation Coefficient",
    y = "Density"
  ) +
  theme_minimal()

p2_plot2

# Prepare data for sum of positive correlations
sum_positive_data <- data.frame(
  Value = c(observed_results$stat_sum_positive_cors, null_sum_positive_cors),
  Type = c("Observed", rep("Null", length(null_sum_positive_cors)))
)

p3_plot <- ggplot(sum_positive_data[sum_positive_data$Type == "Null",], aes(x = Value, fill = Type)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 30) +
  geom_density(aes(color = Type), alpha = 0.8, size = 1.2) +
  scale_fill_manual(values = c("Observed" = "#E31A1C", "Null" = "#1F78B4")) +
  scale_color_manual(values = c("Observed" = "#E31A1C", "Null" = "#1F78B4")) +
  geom_vline(aes(xintercept = observed_results$stat_sum_positive_cors), 
             color = "#E31A1C", size = 1) +
  labs(
   # title = "Distribution of Sum of Positive Correlations",
    subtitle = paste("\nPermuted p-value:", round(p3, 4), "| Observed:", round(observed_results$stat_sum_positive_cors, 3)),
    x = "Sum of Positive Correlations",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
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
  scale_fill_manual(values = c("Observed" = "#E31A1C", "Null" = "#1F78B4")) +
  scale_color_manual(values = c("Observed" = "#E31A1C", "Null" = "#1F78B4")) +
  geom_vline(aes(xintercept = observed_results$stat_median_correlation), 
             color = "#E31A1C", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", alpha = 0.7) +
  labs(
   # title = "Distribution of Median Correlations",
    subtitle = paste("\nPermuted p-value:", round(p6, 4), "| Observed:", round(observed_results$stat_median_correlation, 3)),
    x = "Median Correlation",
    y = "Density"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    legend.position = "none"
  )

p4_plot

#Prep Bar plot for visulalizing number of sig positive corelations
sig_counts <- table(null_n_significant)
sig_props <- prop.table(sig_counts)

bar_data <- data.frame(
  Count = as.numeric(names(sig_counts)),
  Proportion = as.numeric(sig_props),
  Frequency = as.numeric(sig_counts)
)

p_sig_bar <- ggplot(bar_data, aes(x = Count, y = Proportion)) +
  geom_col(fill = "#1F78B4", alpha = 0.7, width = 0.8) +
  geom_vline(xintercept = observed_results$stat_n_significant, 
             color = "#E31A1C", size = 1, linetype = "solid") +
  scale_x_continuous(breaks = 0:max(null_n_significant), 
                     limits = c(-0.5, max(c(null_n_significant, observed_results$stat_n_significant)) + 0.5)) +
  labs(
    subtitle = paste("\nPermuted p-value:", round(p_sig, 4), "| Observed:", observed_results$stat_n_significant),
    x = "Number of Significant Correlations",
    y = "Proportion of Permutations"
  ) +
  theme_minimal() +
  theme(
    plot.subtitle = element_text(size = 12),
    panel.grid.minor.x = element_blank()
  )


#lemon::grid_arrange_shared_legend(p2_plot2,p4_plot,p_sig_bar, nrow=1)
library(cowplot)
plot_grid(p2_plot2,p4_plot,p_sig_bar, nrow=1)


