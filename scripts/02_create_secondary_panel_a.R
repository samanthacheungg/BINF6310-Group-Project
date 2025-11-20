#!/usr/bin/env Rscript

# SECONDARY PANEL A: QC Metrics Visualization
# Author: Anaita Lashley
# Date: 2025-11-19
#
# NOTE: This script must be run interactively in R due to renv/conda conflicts
# Do NOT run with: Rscript scripts/02_create_secondary_panel_a.R
# Instead run with: R
# Then inside R: source("scripts/02_create_secondary_panel_a.R")

cat("=== SECONDARY PANEL A: QC Metrics Visualization ===\n\n")

# Load required packages
cat("Loading packages...\n")
suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
})
cat("✓ Packages loaded\n\n")

# Load QC metadata (all 136 cells, including those that failed)
cat("Loading QC data...\n")
metadata <- read.table("qc/qc_filtered_cell_qc_metrics.txt", 
                      header = TRUE, sep = "\t")

cat(sprintf("Loaded %d cells\n", nrow(metadata)))
cat(sprintf("  - %d passing QC\n", sum(metadata$pass_all)))
cat(sprintf("  - %d failing QC\n\n", sum(!metadata$pass_all)))

# Add QC status labels for plotting
metadata$QC_Status <- ifelse(metadata$pass_all, "Pass QC", "Fail QC")
metadata$QC_Status <- factor(metadata$QC_Status, levels = c("Pass QC", "Fail QC"))

# === PANEL 1: Genes Detected ===
cat("Creating Panel 1: Genes Detected...\n")

p1 <- ggplot(metadata, aes(x = "All Cells", y = genes_detected, fill = QC_Status)) +
    geom_violin(alpha = 0.7, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(0.9)) +
    geom_hline(yintercept = 1000, linetype = "dashed", color = "red", linewidth = 1) +
    scale_fill_manual(values = c("Pass QC" = "#4ECDC4", "Fail QC" = "#FF6B6B")) +
    labs(title = "Genes Detected per Cell", 
         y = "Genes Detected", 
         x = "", 
         fill = "QC Status") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

cat("✓ Panel 1 complete\n")

# === PANEL 2: Total Mapped Reads ===
cat("Creating Panel 2: Total Reads...\n")

p2 <- ggplot(metadata, aes(x = "All Cells", y = mapped_reads, fill = QC_Status)) +
    geom_violin(alpha = 0.7, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(0.9)) +
    geom_hline(yintercept = 40000, linetype = "dashed", color = "red", linewidth = 1) +
    scale_y_log10(labels = scales::comma) +
    scale_fill_manual(values = c("Pass QC" = "#4ECDC4", "Fail QC" = "#FF6B6B")) +
    labs(title = "Total Mapped Reads", 
         y = "Reads (log10)", 
         x = "", 
         fill = "QC Status") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

cat("✓ Panel 2 complete\n")

# === PANEL 3: ERCC Spike-in Proportion ===
cat("Creating Panel 3: ERCC Proportion...\n")

p3 <- ggplot(metadata, aes(x = "All Cells", y = ercc_fraction, fill = QC_Status)) +
    geom_violin(alpha = 0.7, position = position_dodge(0.9)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(0.9)) +
    geom_hline(yintercept = 0.74, linetype = "dashed", color = "red", linewidth = 1) +
    scale_fill_manual(values = c("Pass QC" = "#4ECDC4", "Fail QC" = "#FF6B6B")) +
    labs(title = "ERCC Spike-in Proportion", 
         y = "ERCC Fraction", 
         x = "", 
         fill = "QC Status") +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          legend.position = "bottom",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())

cat("✓ Panel 3 complete\n")

# === COMBINE PANELS ===
cat("\nCombining panels...\n")

combined <- grid.arrange(
    p1, p2, p3, 
    ncol = 3,
    top = grid::textGrob("SECONDARY PANEL A: Quality Control Metrics",
                        gp = grid::gpar(fontsize = 20, fontface = "bold"))
)

# === SAVE FIGURES ===
cat("Saving figures...\n")

ggsave("figures/SECONDARY_PANEL_A_qc_metrics.pdf", 
       combined, 
       width = 18, height = 6, dpi = 300)
cat("✓ Saved: figures/SECONDARY_PANEL_A_qc_metrics.pdf\n")

ggsave("figures/SECONDARY_PANEL_A_qc_metrics.png", 
       combined, 
       width = 18, height = 6, dpi = 300)
cat("✓ Saved: figures/SECONDARY_PANEL_A_qc_metrics.png\n")

cat("\n✅ SECONDARY PANEL A COMPLETE!\n\n")
cat("Summary:\n")
cat(sprintf("  - Total cells plotted: %d\n", nrow(metadata)))
cat(sprintf("  - Pass QC: %d (%.1f%%)\n", 
            sum(metadata$pass_all), 
            100*sum(metadata$pass_all)/nrow(metadata)))
cat(sprintf("  - Fail QC: %d (%.1f%%)\n", 
            sum(!metadata$pass_all), 
            100*sum(!metadata$pass_all)/nrow(metadata)))
cat("\nFigures saved in figures/ folder\n")
