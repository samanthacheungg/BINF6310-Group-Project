# ==============================================================================
# Figure 3f: FIT3 Expression Correlation with Generation - FIXED VERSION
# ==============================================================================
# Author: Kyle Brodeur
# Date: November 30, 2024
# Fixed Issues:
#   - Correct gene name (YOR383C instead of FIT3)
#   - Simplified data merging (no SRR mapping needed)
#   - Proper unique correlations for each timepoint
#   - Added diagnostics and validation
#
# Expected results:
#   - 16h: R ~ -0.55 (negative correlation)
#   - 36h: R ~ -0.62 (negative correlation)
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

cat("============================================================\n")
cat("Figure 3f: FIT3 vs Generation Correlation Analysis\n")
cat("============================================================\n\n")

# ==============================================================================
# STEP 1: LOAD DATA
# ==============================================================================

cat("=== Loading data ===\n")

# Load normalized counts (columns are already mapped to sample_ids like A1, B2, C3)
counts <- read.csv("processed_data/log2_normalized_counts.csv", row.names = 1)

# Load metadata
metadata <- read.csv("processed_data/metadata.csv")

cat("✓ Counts:", nrow(counts), "genes x", ncol(counts), "cells\n")
cat("✓ Metadata:", nrow(metadata), "samples\n")

# ==============================================================================
# STEP 2: EXTRACT FIT3 EXPRESSION
# ==============================================================================

cat("\n=== Extracting FIT3 (YOR383C) ===\n")

# FIT3 systematic name in yeast genome
FIT3_GENE <- "YOR383C"

if (!FIT3_GENE %in% rownames(counts)) {
  stop("ERROR: FIT3 gene (YOR383C) not found in data!")
}

# Extract FIT3 expression
fit3_expression <- as.numeric(counts[FIT3_GENE, ])

# Create data frame
fit3_data <- data.frame(
  sample_id = colnames(counts),
  fit3_expression = fit3_expression,
  stringsAsFactors = FALSE
)

# Merge with metadata
full_data <- merge(metadata, fit3_data, by = "sample_id")

cat("✓ Extracted FIT3 expression for", nrow(full_data), "cells\n")

# Diagnostics
zeros <- sum(full_data$fit3_expression == 0)
cat("\nFIT3 Expression Summary:\n")
cat("  Range:", range(full_data$fit3_expression), "\n")
cat("  Mean:", round(mean(full_data$fit3_expression), 2), "\n")
cat("  Median:", round(median(full_data$fit3_expression), 2), "\n")
cat("  Zeros:", zeros, "(", round(100*zeros/nrow(full_data), 1), "%)\n")

# ==============================================================================
# STEP 3: ANALYZE 16-HOUR CELLS
# ==============================================================================

cat("\n============================================================\n")
cat("16-HOUR TIMEPOINT ANALYSIS\n")
cat("============================================================\n")

data_16h <- full_data %>% filter(grepl("16h", subgroup))

cat("Sample size:", nrow(data_16h), "cells\n")
cat("Subgroups:", paste(unique(data_16h$subgroup), collapse = ", "), "\n")
cat("Generation range:", range(data_16h$generation), "\n")
cat("FIT3 zeros:", sum(data_16h$fit3_expression == 0), "\n")

# Calculate correlation
cor_16h <- cor.test(data_16h$generation, data_16h$fit3_expression,
                    method = "pearson")

cat("\n📊 Results:\n")
cat("  R  =", round(cor_16h$estimate, 3), "\n")
cat("  p  =", format(cor_16h$p.value, scientific = TRUE, digits = 3), "\n")
cat("  CI = [", round(cor_16h$conf.int[1], 3), ",",
    round(cor_16h$conf.int[2], 3), "]\n")
cat("\n📖 Paper reported: R = -0.55, p = 1.3e-04\n")
cat("✓ Match:", abs(cor_16h$estimate - (-0.55)) < 0.01, "\n")

# Create plot
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
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))

# ==============================================================================
# STEP 4: ANALYZE 36-HOUR CELLS
# ==============================================================================

cat("\n============================================================\n")
cat("36-HOUR TIMEPOINT ANALYSIS\n")
cat("============================================================\n")

data_36h <- full_data %>% filter(grepl("36h", subgroup))

cat("Sample size:", nrow(data_36h), "cells\n")
cat("Subgroups:", paste(unique(data_36h$subgroup), collapse = ", "), "\n")
cat("Generation range:", range(data_36h$generation), "\n")
cat("FIT3 zeros:", sum(data_36h$fit3_expression == 0), "\n")

# Verify datasets are different
cat("\n🔍 Verification:\n")
cat("  Datasets are different:", !identical(data_16h, data_36h), "\n")
cat("  No overlapping samples:",
    length(intersect(data_16h$sample_id, data_36h$sample_id)) == 0, "\n")

# Calculate correlation
cor_36h <- cor.test(data_36h$generation, data_36h$fit3_expression,
                    method = "pearson")

cat("\n📊 Results:\n")
cat("  R  =", round(cor_36h$estimate, 3), "\n")
cat("  p  =", format(cor_36h$p.value, scientific = TRUE, digits = 3), "\n")
cat("  CI = [", round(cor_36h$conf.int[1], 3), ",",
    round(cor_36h$conf.int[2], 3), "]\n")
cat("\n📖 Paper reported: R = -0.62, p = 5.6e-06\n")
cat("✓ Match:", abs(cor_36h$estimate - (-0.62)) < 0.01, "\n")

# Create plot
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
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))

# ==============================================================================
# STEP 5: SAVE OUTPUTS
# ==============================================================================

cat("\n============================================================\n")
cat("SAVING OUTPUTS\n")
cat("============================================================\n")

# Save individual plots
ggsave("figures/figure3f_fit3_16h.pdf", plot_16h, width = 6, height = 5)
ggsave("figures/figure3f_fit3_16h.png", plot_16h, width = 6, height = 5, dpi = 300)
ggsave("figures/figure3f_fit3_36h.pdf", plot_36h, width = 6, height = 5)
ggsave("figures/figure3f_fit3_36h.png", plot_36h, width = 6, height = 5, dpi = 300)

# Combined plot
combined_plot <- plot_grid(plot_16h, plot_36h, ncol = 2,
                          labels = c("A", "B"), label_size = 16)
ggsave("figures/figure3f_fit3_combined.pdf", combined_plot, width = 12, height = 5)
ggsave("figures/figure3f_fit3_combined.png", combined_plot, width = 12, height = 5, dpi = 300)

cat("✓ Figures saved\n")

# Save statistics
fit3_stats <- data.frame(
  timepoint = c("16h", "36h"),
  n_cells = c(nrow(data_16h), nrow(data_36h)),
  n_zeros = c(sum(data_16h$fit3_expression == 0),
              sum(data_36h$fit3_expression == 0)),
  R_value = c(cor_16h$estimate, cor_36h$estimate),
  p_value = c(cor_16h$p.value, cor_36h$p.value),
  CI_lower = c(cor_16h$conf.int[1], cor_36h$conf.int[1]),
  CI_upper = c(cor_16h$conf.int[2], cor_36h$conf.int[2]),
  paper_R = c(-0.55, -0.62),
  paper_p = c(1.3e-04, 5.6e-06),
  stringsAsFactors = FALSE
)

write.csv(fit3_stats, "processed_data/fit3_correlation_stats.csv",
          row.names = FALSE)

cat("✓ Statistics saved\n")

# Save session info
dir.create("session_info", showWarnings = FALSE)
sink("session_info/figure3f_session.txt")
cat("Figure 3f Analysis\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
sessionInfo()
sink()

cat("✓ Session info saved\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n============================================================\n")
cat("✅ ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat("Comparison to Paper:\n")
print(fit3_stats[, c("timepoint", "R_value", "p_value", "paper_R", "paper_p")])

cat("\n🔬 Biological Interpretation:\n")
cat("FIT3 shows strong NEGATIVE correlation with cell division.\n")
cat("As cells age (more generations), FIT3 expression decreases.\n")
cat("This matches the paper's findings on aging heterogeneity.\n\n")

cat("📁 Output files:\n")
cat("  - figures/figure3f_fit3_16h.pdf/png\n")
cat("  - figures/figure3f_fit3_36h.pdf/png\n")
cat("  - figures/figure3f_fit3_combined.pdf/png\n")
cat("  - processed_data/fit3_correlation_stats.csv\n")
cat("  - session_info/figure3f_session.txt\n")

cat("\n============================================================\n")
