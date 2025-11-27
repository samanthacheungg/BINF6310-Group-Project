# ==============================================================================
# Figure 3b: Correlation between genes detected and generation at 36h
# ==============================================================================
# Author: Anaita Lashley
# Date: November 26, 2024
# 
# This script creates a scatter plot showing the relationship between 
# generation number and number of genes detected in 36-hour yeast cells.
# Expected result: Positive correlation (R ~ 0.57 in paper)
# ==============================================================================

# Load required libraries
library(ggplot2)  # For creating plots
library(dplyr)    # For data manipulation (filtering, merging)

# ==============================================================================
# STEP 1: LOAD DATA FILES
# ==============================================================================

# Load the metadata file containing sample information
# This file has: sample_id, subgroup (2h, 16h/S, 16h/F, etc.), generation
metadata <- read.csv("processed_data/metadata.csv")

# Load the QC metrics file containing quality control measurements
# This file has: cell (SRR ID), genes_detected, mapped_reads, etc.
qc_metrics <- read.delim("processed_data/filtered/qc_filtered_cell_qc_metrics.txt")

# Load the mapping file that links sample IDs to SRR IDs
# This file has: srr_id, sample_id, correlation
mapping <- read.csv("processed_data/srr_to_sample_correlation_mapping.csv")

# ==============================================================================
# STEP 2: MERGE DATA FILES
# ==============================================================================

# First merge: Add SRR IDs to metadata using the mapping file
# by.x = "sample_id" means use this column from metadata
# by.y = "sample_id" means use this column from mapping
# all.x = TRUE keeps all rows from metadata even if no match
metadata_with_srr <- merge(metadata, mapping, by = "sample_id", all.x = TRUE)

# Second merge: Add QC metrics to metadata using SRR IDs
# by.x = "srr_id" means use this column from metadata_with_srr
# by.y = "cell" means use this column from qc_metrics
# Result: One dataframe with sample_id, subgroup, generation, genes_detected
full_metadata <- merge(metadata_with_srr, qc_metrics, 
                       by.x = "srr_id", 
                       by.y = "cell")

# Print summary to verify data loaded correctly
cat("Total cells loaded:", nrow(full_metadata), "\n")
cat("Cells per subgroup:\n")
print(table(full_metadata$subgroup))

# ==============================================================================
# STEP 3: FILTER DATA FOR 36-HOUR CELLS
# ==============================================================================

# Filter to keep only 36h cells (both 36h/S and 36h/F)
# grepl("36h", subgroup) finds any subgroup containing "36h"
# %>% is the pipe operator - it passes data to the next function
metadata_36h <- full_metadata %>% filter(grepl("36h", subgroup))

# Print how many 36h cells we have
cat("\n36h cells:", nrow(metadata_36h), "\n")

# ==============================================================================
# STEP 4: CALCULATE CORRELATION
# ==============================================================================

# Perform Pearson correlation test between generation and genes_detected
# This calculates:
#   - R value (correlation coefficient): measures strength of relationship
#   - p-value: measures statistical significance
cor_36h <- cor.test(metadata_36h$generation, metadata_36h$genes_detected)

# Print the correlation results
cat(sprintf("36h correlation: R = %.2f, p = %.2e\n", 
            cor_36h$estimate, cor_36h$p.value))
cat("Paper reported: R = 0.57, p = 2.6e-06\n")

# ==============================================================================
# STEP 5: CREATE THE SCATTER PLOT
# ==============================================================================

# Create the plot using ggplot2
fig3b <- ggplot(metadata_36h, aes(x = generation, y = genes_detected)) +
    
    # Add red points for each cell
    # size = 3 controls point size
    # alpha = 0.7 makes points slightly transparent (0=invisible, 1=solid)
    geom_point(color = "red", size = 3, alpha = 0.7) +
    
    # Add blue linear regression line with grey confidence interval
    # method = "lm" means linear model (straight line)
    # se = TRUE shows the shaded confidence interval
    # fill = "grey80" colors the confidence interval
    # linewidth = 1.2 controls line thickness
    geom_smooth(method = "lm", se = TRUE, fill = "grey80", 
                color = "blue", linewidth = 1.2) +
    
    # Add vertical dashed line at mean generation (divides slow/fast)
    # mean(metadata_36h$generation) calculates the average generation
    # linetype = "dashed" makes it a dashed line
    geom_vline(xintercept = mean(metadata_36h$generation), 
               linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # Add text annotation showing R and p-value in top-right corner
    # x = Inf, y = Inf positions text at the plot edges
    # hjust/vjust fine-tune the exact position
    # sprintf formats numbers: %.2f = 2 decimal places, %.2e = scientific notation
    # \n creates a line break between R and p
    annotate("text", x = Inf, y = Inf, 
             label = sprintf("R = %.2f\np = %.2e", 
                           cor_36h$estimate, cor_36h$p.value),
             hjust = 1.2, vjust = 1.5, size = 5) +
    
    # Add axis labels and plot title
    labs(x = "Generation at 36 hr", 
         y = "Number of Genes Detected",
         title = "36h: Genes Detected vs Generation") +
    
    # Use a clean minimal theme with 14-point font
    theme_minimal(base_size = 14)

# ==============================================================================
# STEP 6: SAVE THE PLOT
# ==============================================================================

# Save as PDF (vector format - good for publications, scales without pixelation)
# width and height are in inches
ggsave("figures/figure3b_correlation_36h.pdf", fig3b, width = 6, height = 5)

# Save as PNG (raster format - good for presentations and documents)
# dpi = 300 means 300 dots per inch (high quality)
ggsave("figures/figure3b_correlation_36h.png", fig3b, width = 6, height = 5, dpi = 300)

# Print success message
cat("\n✓ Figure 3b saved successfully!\n")
cat("  - figures/figure3b_correlation_36h.pdf\n")
cat("  - figures/figure3b_correlation_36h.png\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
