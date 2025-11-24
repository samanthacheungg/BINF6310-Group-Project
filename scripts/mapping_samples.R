#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ---- file paths ----
log_file <- "../processed_data/log2_normalized_counts.txt"
map_file <- "../processed_data/srr_to_sample_correlation_mapping.csv"
gse_file <- "../processed_data/GEO data/GSE210032_Log2NormCounts_1W-A37B43C45_H2_NoERCC_CSV.csv"
out_file <- "../processed_data/log2norm_mapped_GSEgenes_alphabetical.csv"

# ---- 1. Read log2-normalized counts ----
log2dat <- read_tsv(log_file, show_col_types = FALSE)
colnames(log2dat)[1] <- "Geneid"

# ---- 2. Extract SRR IDs from BAM path column names ----
bam_cols <- colnames(log2dat)[-1]
srr_ids <- gsub(".*(SRR[0-9]+).*", "\\1", bam_cols)

# ---- 3. Map SRR → sample names ----
map <- read_csv(map_file, show_col_types = FALSE)  # must contain srr_id + sample_id

sample_names <- map$sample_id[match(srr_ids, map$srr_id)]

if (any(is.na(sample_names))) {
  warning("Some SRR IDs could not be mapped:\n",
          paste(srr_ids[is.na(sample_names)], collapse=", "))
}

# Apply sample names
colnames(log2dat) <- c("Geneid", sample_names)

# ---- 4. Read paper GSE file ----
gse <- read_csv(gse_file, show_col_types = FALSE)
gse_genes <- gse$Geneid

# ---- 5. Restrict log2 data to only GSE genes (same order as paper) ----
log2_sub <- log2dat %>%
  filter(Geneid %in% gse_genes) %>%
  arrange(match(Geneid, gse_genes))

# ---- 6. Reorder sample columns alphabetically ----
sample_cols <- sort(colnames(log2_sub)[-1])  # alphabetical
log2_alph <- log2_sub[, c("Geneid", sample_cols)]

# ---- 7. Write output ----
write_csv(log2_alph, out_file)
cat("Wrote:", out_file, "\n")

