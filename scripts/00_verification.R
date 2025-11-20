#!/usr/bin/env Rscript

cat("=== QC VERIFICATION REPORT ===\n\n")
cat(sprintf("Date: %s\n", Sys.Date()))
cat(sprintf("Analyst: Anaita Lashley\n\n"))

cat("=== FILE CHECK ===\n")

# Check critical files exist
critical_files <- c(
    "counts/gene_counts.txt",
    "qc/qc_filtered_gene_counts.txt",
    "qc/qc_filtered_cell_qc_metrics.txt",
    "normalize/log2_normalized_counts.txt"
)

file_descriptions <- c(
    "Raw counts",
    "Filtered counts",
    "Cell QC metrics",
    "Log2 normalized counts"
)

for (i in 1:length(critical_files)) {
    file <- critical_files[i]
    desc <- file_descriptions[i]
    
    if (file.exists(file)) {
        size <- file.info(file)$size / 1024^2
        cat(sprintf("✓ %s: %s (%.1f MB)\n", desc, file, size))
    } else {
        cat(sprintf("✗ %s: %s MISSING\n", desc, file))
    }
}

cat("\n=== QC VALIDATION ===\n")

# Load the QC metadata
metadata <- read.table("qc/qc_filtered_cell_qc_metrics.txt", 
                      header = TRUE, sep = "\t")

cat(sprintf("Total cells in file: %d\n", nrow(metadata)))

# Filter to only passing cells
if ("pass_all" %in% colnames(metadata)) {
    passing_cells <- sum(metadata$pass_all)
    failing_cells <- nrow(metadata) - passing_cells
    
    cat(sprintf("  - Cells PASSING all QC: %d\n", passing_cells))
    cat(sprintf("  - Cells FAILING QC: %d\n", failing_cells))
    
    # Filter to only passing cells
    metadata_passed <- metadata[metadata$pass_all == TRUE, ]
    
    cat(sprintf("\nAfter filtering: %d cells (Expected: 125)\n", nrow(metadata_passed)))
    
    # Calculate stats on passing cells only
    mean_genes_passed <- mean(metadata_passed$genes_detected)
    cat(sprintf("Mean genes detected: %.0f (Paper: 2202)\n", mean_genes_passed))
    
    pct_diff <- abs(mean_genes_passed - 2202) / 2202 * 100
    
    if (pct_diff < 10) {
        cat(sprintf("✓ Within 10%% of paper (%.1f%% difference)\n", pct_diff))
    } else {
        cat(sprintf("⚠ Differs from paper by %.1f%%\n", pct_diff))
    }
    
    # QC breakdown
    cat("\n=== QC Criteria Breakdown ===\n")
    cat(sprintf("Pass mapped reads (>40,000): %d / %d\n", 
                sum(metadata$pass_mapped), nrow(metadata)))
    cat(sprintf("Pass genes detected (>1,000): %d / %d\n", 
                sum(metadata$pass_genes), nrow(metadata)))
    cat(sprintf("Pass ERCC fraction (<0.74): %d / %d\n", 
                sum(metadata$pass_ercc), nrow(metadata)))
    
    # Final validation
    if (nrow(metadata_passed) == 125) {
        cat("\n✅ QC VALIDATION PASSED - 125 cells ready for analysis!\n")
    } else {
        cat(sprintf("\n⚠️ Cell count mismatch: %d cells (expected 125)\n", 
                   nrow(metadata_passed)))
    }
    
} else {
    cat("✗ ERROR: No 'pass_all' column found in metadata\n")
}

cat("\n=== END VERIFICATION ===\n")
