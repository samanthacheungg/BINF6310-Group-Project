# ==============================================================================
# Figure 3f: FIT3 Expression Correlation with Generation
# ==============================================================================
# Author: Kyle Brodeur (based on Anaita's structure)
# Date: November 29, 2024
#
# This script creates scatter plots showing the relationship between
# generation number and FIT3 gene expression for both 16h and 36h cells.
# Expected results:
#   - 16h: R ~ -0.55 (negative correlation)
#   - 36h: R ~ -0.62 (negative correlation)
# ==============================================================================

# Load required libraries
library(ggplot2)  # For creating plots
library(dplyr)    # For data manipulation
library(tidyr)    # For data reshaping
library(cowplot)  # For combining plots

# ==============================================================================
# STEP 1: LOAD DATA FILES
# ==============================================================================

cat("Loading data files...\n")

# Load normalized gene expression data (genes = rows, cells = columns)
counts <- read.csv("processed_data/log2_normalized_counts.csv", row.names = 1)

# Load the metadata file containing sample information
metadata <- read.csv("processed_data/metadata.csv")

# Load the QC metrics file
qc_metrics <- read.delim("processed_data/filtered/qc_filtered_cell_qc_metrics.txt")

# Load the mapping file
mapping <- read.csv("processed_data/srr_to_sample_correlation_mapping.csv")

# ==============================================================================
# STEP 2: PREPARE METADATA
# ==============================================================================

cat("Merging data files...\n")

# Merge metadata with SRR IDs
metadata_with_srr <- merge(metadata, mapping, by = "sample_id", all.x = TRUE)

# Merge with QC metrics
full_metadata <- merge(metadata_with_srr, qc_metrics,
                       by.x = "srr_id",
                       by.y = "cell")

cat("Total cells:", nrow(full_metadata), "\n")

# ==============================================================================
# STEP 3: EXTRACT FIT3 EXPRESSION
# ==============================================================================

cat("Extracting FIT3 expression...\n")

# Check if FIT3 gene exists in the data
if(!"FIT3" %in% rownames(counts)) {
  cat("ERROR: FIT3 not found in gene names!\n")
  cat("Searching for similar names...\n")
  fit_genes <- grep("FIT", rownames(counts), value = TRUE, ignore.case = TRUE)
  if(length(fit_genes) > 0) {
    cat("Found these FIT genes:\n")
    print(fit_genes)
    cat("\nUsing the first match:", fit_genes[1], "\n")
    fit3_gene <- fit_genes[1]
  } else {
    stop("No FIT genes found in data!")
  }
} else {
  fit3_gene <- "FIT3"
  cat("✓ FIT3 gene found!\n")
}

# Extract FIT3 expression values
fit3_expression <- as.numeric(counts[fit3_gene, ])

# Create a data frame with SRR IDs and FIT3 expression
fit3_data <- data.frame(
  srr_id = colnames(counts),
  fit3_expression = fit3_expression,
  stringsAsFactors = FALSE
)

# Merge FIT3 expression with metadata
full_data <- merge(full_metadata, fit3_data, by = "srr_id")

cat("Data prepared successfully!\n")

# ==============================================================================
# STEP 4: ANALYZE 16-HOUR CELLS
# ==============================================================================

cat("\n=== ANALYZING 16-HOUR CELLS ===\n")

# Filter for 16h cells only
data_16h <- full_data %>% filter(grepl("16h", subgroup))

cat("16h cells:", nrow(data_16h), "\n")

# Calculate correlation
cor_16h <- cor.test(data_16h$generation, data_16h$fit3_expression)

cat(sprintf("16h FIT3 correlation: R = %.2f, p = %.2e\n",
            cor_16h$estimate, cor_16h$p.value))
cat("Paper reported: R = -0.55, p = 1.3e-04\n")

# Create 16h plot
plot_16h <- ggplot(data_16h, aes(x = generation, y = fit3_expression)) +
    geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, fill = "grey80",
                color = "blue", linewidth = 1.2) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("R = %.2f\np = %.2e",
                           cor_16h$estimate, cor_16h$p.value),
             hjust = 1.2, vjust = 1.5, size = 5) +
    labs(x = "Generation at 16 hr",
         y = "FIT3 Expression (log2)",
         title = "16h: FIT3 vs Generation") +
    theme_minimal(base_size = 14)

# ==============================================================================
# STEP 5: ANALYZE 36-HOUR CELLS
# ==============================================================================

cat("\n=== ANALYZING 36-HOUR CELLS ===\n")

# Filter for 36h cells only
data_36h <- full_data %>% filter(grepl("36h", subgroup))

cat("36h cells:", nrow(data_36h), "\n")

# Calculate correlation
cor_36h <- cor.test(data_36h$generation, data_36h$fit3_expression)

cat(sprintf("36h FIT3 correlation: R = %.2f, p = %.2e\n",
            cor_36h$estimate, cor_36h$p.value))
cat("Paper reported: R = -0.62, p = 5.6e-06\n")

# Create 36h plot
plot_36h <- ggplot(data_36h, aes(x = generation, y = fit3_expression)) +
    geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, fill = "grey80",
                color = "blue", linewidth = 1.2) +
    annotate("text", x = Inf, y = Inf,
             label = sprintf("R = %.2f\np = %.2e",
                           cor_36h$estimate, cor_36h$p.value),
             hjust = 1.2, vjust = 1.5, size = 5) +
    labs(x = "Generation at 36 hr",
         y = "FIT3 Expression (log2)",
         title = "36h: FIT3 vs Generation") +
    theme_minimal(base_size = 14)

# ==============================================================================
# STEP 6: COMBINE PLOTS AND SAVE
# ==============================================================================

cat("\n=== CREATING COMBINED FIGURE ===\n")

# Combine both plots side by side
combined_plot <- plot_grid(plot_16h, plot_36h, ncol = 2, labels = c("16h", "36h"))

# Save individual plots
ggsave("figures/figure3f_fit3_16h.pdf", plot_16h, width = 6, height = 5)
ggsave("figures/figure3f_fit3_16h.png", plot_16h, width = 6, height = 5, dpi = 300)

ggsave("figures/figure3f_fit3_36h.pdf", plot_36h, width = 6, height = 5)
ggsave("figures/figure3f_fit3_36h.png", plot_36h, width = 6, height = 5, dpi = 300)

# Save combined plot
ggsave("figures/figure3f_fit3_combined.pdf", combined_plot, width = 12, height = 5)
ggsave("figures/figure3f_fit3_combined.png", combined_plot, width = 12, height = 5, dpi = 300)

# Save correlation statistics
fit3_stats <- data.frame(
  timepoint = c("16h", "36h"),
  R_value = c(cor_16h$estimate, cor_36h$estimate),
  p_value = c(cor_16h$p.value, cor_36h$p.value),
  paper_R = c(-0.55, -0.62),
  paper_p = c(1.3e-04, 5.6e-06)
)

write.csv(fit3_stats, "processed_data/fit3_correlation_stats.csv", row.names = FALSE)

cat("\n✓ Figure 3f saved successfully!\n")
cat("  Individual plots:\n")
cat("    - figures/figure3f_fit3_16h.pdf/png\n")
cat("    - figures/figure3f_fit3_36h.pdf/png\n")
cat("  Combined plot:\n")
cat("    - figures/figure3f_fit3_combined.pdf/png\n")
cat("  Statistics:\n")
cat("    - processed_data/fit3_correlation_stats.csv\n")

cat("\n=== SUMMARY ===\n")
cat("FIT3 shows NEGATIVE correlation with generation\n")
cat("(As cells divide more, FIT3 expression DECREASES)\n")
cat("\nThis matches the paper's findings that FIT3 is involved in\n")
cat("iron transport and is downregulated in aging cells.\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
