##############################
## qc_filter_cells.R
## QC to reproduce yeast aging scRNA-seq filtering:
## - remove cells with:
##   < 1000 yeast genes detected
##   < 40,000 total mapped reads
##   ERCC / total mapped > 0.74
##############################

#### 0. Paths ####
counts_file   <- "../counts/gene_counts.txt"
bam_dir       <- "../alignment/bam"
output_prefix <- "qc_filtered"

cat("\n=== Loading count matrix ===\n")

#### 1. Read gene-level counts from featureCounts ####
counts_raw <- read.delim(
  counts_file,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

count_col_idx    <- 7:ncol(counts_raw)
count_col_paths  <- colnames(counts_raw)[count_col_idx]

# rownames = Gene IDs
rownames(counts_raw) <- counts_raw$Geneid

count_mat <- as.matrix(counts_raw[, count_col_idx])

# clean sample names
sample_names <- basename(count_col_paths)
sample_names <- sub("_sorted\\.bam$", "", sample_names)
colnames(count_mat) <- sample_names

total_cells <- ncol(count_mat)
cat("Total cells in dataset:", total_cells, "\n\n")

#### 2. Per-cell metrics from counts ####
cat("=== Computing per-cell metrics from counts ===\n")

# ERCC detection
ercc_rows <- grepl("^ERCC", rownames(count_mat))
ercc_reads <- colSums(count_mat[ercc_rows, , drop = FALSE])

# Yeast gene detection
genes_detected <- colSums(count_mat[!ercc_rows, , drop = FALSE] > 0)

cat("Summary of yeast genes detected per cell:\n")
print(summary(genes_detected))

#### 3. Total mapped reads per cell from BAMs ####
cat("\n=== Reading mapped reads from BAM files ===\n")

bam_files <- file.path(bam_dir, paste0(sample_names, "_sorted.bam"))

get_mapped_reads <- function(bam) {
  cmd <- paste("samtools flagstat", shQuote(bam))
  x   <- tryCatch(system(cmd, intern = TRUE), error = function(e) character(0))
  line <- x[grepl(" mapped \\(", x)]
  if (length(line) == 0L) return(NA_integer_)
  as.integer(sub(" .*", "", line))
}

mapped_reads <- vapply(bam_files, get_mapped_reads, integer(1))
cat("Summary of total mapped reads per cell:\n")
print(summary(mapped_reads))

#### 4. Apply QC criteria ####
cat("\n=== Applying QC Filters ===\n")

ercc_fraction <- ercc_reads / mapped_reads

# thresholds from paper
pass_mapped <- mapped_reads   >= 40000
pass_genes  <- genes_detected >= 1000
pass_ercc   <- ercc_fraction  <= 0.74

cat("Cells failing mapped_reads < 40,000 :", sum(!pass_mapped), "\n")
cat("Cells failing genes_detected < 1000:", sum(!pass_genes),  "\n")
cat("Cells failing ERCC fraction > 0.74 :", sum(!pass_ercc),   "\n\n")

qc_df <- data.frame(
  cell           = sample_names,
  mapped_reads   = mapped_reads,
  genes_detected = genes_detected,
  ercc_reads     = ercc_reads,
  ercc_fraction  = ercc_fraction,
  pass_mapped    = pass_mapped,
  pass_genes     = pass_genes,
  pass_ercc      = pass_ercc,
  stringsAsFactors = FALSE
)

qc_df$pass_all <- with(qc_df, pass_mapped & pass_genes & pass_ercc)

cat("Cells passing mapped_read filter   :", sum(pass_mapped), "\n")
cat("Cells passing gene_detected filter :", sum(pass_genes),  "\n")
cat("Cells passing ERCC filter          :", sum(pass_ercc),   "\n")
cat("Cells passing ALL filters          :", sum(qc_df$pass_all), "\n\n")

#### 5. Filter counts matrix
cells_keep <- qc_df$cell[qc_df$pass_all]
count_mat_filt <- count_mat[, cells_keep, drop = FALSE]

#### 6. Write outputs ####
cat("=== Writing output files ===\n")

write.table(
  qc_df,
  file = paste0(output_prefix, "_cell_qc_metrics.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# use first 6 annotation cols + selected cell columns
keep_idx_in_samples <- match(cells_keep, sample_names)
keep_idx_in_counts  <- count_col_idx[keep_idx_in_samples]
counts_filt <- counts_raw[, c(1:6, keep_idx_in_counts)]

write.table(
  counts_filt,
  file = paste0(output_prefix, "_gene_counts.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nQC complete.\n")
cat("Cells before filtering:", total_cells, "\n")
cat("Cells after filtering :", length(cells_keep), "\n\n")

