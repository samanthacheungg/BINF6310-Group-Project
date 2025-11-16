library(DESeq2)

counts_file <- "../qc/qc_filtered_gene_counts.txt"  # adjust if needed

cat("Loading counts from:", counts_file, "\n")

df <- read.delim(
  counts_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Gene IDs as rownames
rownames(df) <- df$Geneid

# featureCounts: columns 7+ are count columns
count_mat <- as.matrix(df[, 7:ncol(df)])
mode(count_mat) <- "integer"

cat("Count matrix dimensions:",
    nrow(count_mat), "genes ×", ncol(count_mat), "cells\n")

# 2. Build minimal colData -------------------------------

samples <- colnames(count_mat)

# Dummy condition (DESeq2 just needs something for the design)
coldata <- data.frame(
  row.names = samples,
  condition = factor(rep("A", length(samples)))
)

# 3. Create DESeqDataSet and estimate size factors -------

cat("Creating DESeqDataSet...\n")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = coldata,
  design    = ~ 1   # no groups; just size factors
)

# Drop genes with zero total counts
dds <- dds[rowSums(DESeq2::counts(dds)) > 0, ]

cat("Estimating size factors...\n")
dds <- DESeq2::estimateSizeFactors(dds)

# 4. Extract normalized counts ---------------------------

norm_counts <- DESeq2::counts(dds, normalized = TRUE)

cat("Writing DESeq2 normalized counts to normalized_counts.txt\n")

write.table(
  norm_counts,
  file = "normalized_counts.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# 5. Compute log2NormCounts (log2(norm + 1)) ------------- 

log2_norm <- log2(norm_counts + 1)

cat("Writing log2-normalized counts to log2_normalized_counts.txt\n")

write.table(
  log2_norm,
  file = "log2_normalized_counts.txt",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

cat("Done.\n")
