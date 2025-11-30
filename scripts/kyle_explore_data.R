
# ==============================================================================
# Exploring the Yeast Aging Data
# By: Kyle Brodeur
# Date: November 29, 2024
# Purpose: Understand what data we have before running correlation analyses
# ==============================================================================

# Clear the workspace
rm(list = ls())

# Load required libraries
library(dplyr)

# Set working directory to the project folder
setwd("/Users/brodeur.k/Desktop/BINF6310-Group-Project")

# ==============================================================================
# STEP 1: Load the gene expression data
# ==============================================================================
cat("Loading gene expression data...\n")

# Read the normalized counts (genes = rows, cells = columns)
counts <- read.csv("processed_data/log2_normalized_counts.csv", row.names = 1)

cat("Data loaded successfully!\n")
cat("Number of genes:", nrow(counts), "\n")
cat("Number of cells:", ncol(counts), "\n\n")

# Look at first few rows
cat("First 5 genes in first 3 cells:\n")
print(counts[1:5, 1:3])

# ==============================================================================
# STEP 2: Load the metadata (info about each cell)
# ==============================================================================
cat("\n\nLoading cell metadata...\n")

metadata <- read.csv("processed_data/metadata_filtered.csv", row.names = 1)

cat("Metadata loaded!\n")
cat("Number of cells:", nrow(metadata), "\n\n")

cat("Columns in metadata:\n")
print(colnames(metadata))

cat("\n\nFirst 5 cells:\n")
print(head(metadata, 5))

# ==============================================================================
# STEP 3: Count how many cells we have per timepoint
# ==============================================================================
cat("\n\n=== HOW MANY CELLS PER TIMEPOINT? ===\n")
print(table(metadata$timepoint))

# ==============================================================================
# STEP 4: Calculate "genes detected" for each cell
# ==============================================================================
cat("\n\n=== CALCULATING GENES DETECTED ===\n")
cat("(This counts how many genes are expressed in each cell)\n\n")

# Count genes with expression > 0 in each cell
genes_detected <- colSums(counts > 0)

cat("Range of genes detected:\n")
print(summary(genes_detected))

# Add to metadata
metadata$genes_detected <- genes_detected[rownames(metadata)]

# ==============================================================================
# STEP 5: Look for the FIT3 gene
# ==============================================================================
cat("\n\n=== CHECKING FOR FIT3 GENE ===\n")

if("FIT3" %in% rownames(counts)) {
  cat("✓ FIT3 gene found in the data!\n")
  fit3_expression <- as.numeric(counts["FIT3", ])
  cat("FIT3 expression range:\n")
  print(summary(fit3_expression))
} else {
  cat("✗ FIT3 not found with exact name.\n")
  cat("Searching for similar names...\n")
  fit3_matches <- grep("FIT", rownames(counts), value = TRUE, ignore.case = TRUE)
  print(fit3_matches)
}

# ==============================================================================
# STEP 6: Save the updated metadata
# ==============================================================================
cat("\n\n=== SAVING UPDATED METADATA ===\n")

write.csv(metadata, "processed_data/metadata_with_genes_detected.csv")

cat("✓ Saved to: processed_data/metadata_with_genes_detected.csv\n")
cat("\n✓ Data exploration complete! You're ready for the next step.\n")
