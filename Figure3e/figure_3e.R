# -----------------------
# 0. Load libraries
# -----------------------
library(tidyr)
library(dplyr)
library(tibble)   # for rownames_to_column
library(ggplot2)
library(patchwork) # for combining plots

# -----------------------
# 1. Read counts
# -----------------------
counts <- read.table(
  "/Users/kirensadiq/Desktop/BINF6310-Group-Project/Figure 3e/log2_normalized_counts.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# Clean counts column names
colnames(counts) <- gsub("alignment/bam/|_sorted.bam", "", colnames(counts))

# -----------------------
# 2. Read full metadata mapping
# -----------------------
meta <- read.csv(
  "/Users/kirensadiq/Desktop/BINF6310-Group-Project/Figure 3e/metadata.csv",
  stringsAsFactors = FALSE
)

# -----------------------
# 3. Pivot counts to long format
# -----------------------
counts_long <- counts %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "SRR_ID",
    values_to = "log2_count"
  ) %>%
  left_join(meta, by = "SRR_ID")

# Check data
head(counts_long)
table(counts_long$subgroup)  # should show counts per subgroup

# -----------------------
# 4. Define gene sets
# -----------------------
FIT3_genes <- c("FIT3")
HAC1_genes <- c("HAC1")
mito_genes_11 <- c("ATG27", "COR1", "MRPL28", "COX4", "ISC1", "UFD4",
                    "GDB1", "SAP185", "LYS4", "CAB5", "NCS6")

# -----------------------
# 5. Colors for subgroups
# -----------------------
fill_colors <- c(
  "2h" = "#B0B0B0",
  "16h/S" = "#66C2A5",
  "16h/F" = "#66C2A5",
  "36h/S" = "#3288BD",
  "36h/F" = "#3288BD"
)

# -----------------------
# 6. Plotting function
# -----------------------
plot_gene_set <- function(genes, title){
  counts_long %>%
    filter(gene %in% genes) %>%
    ggplot(aes(x = subgroup, y = log2_count, fill = subgroup)) +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = fill_colors) +
    labs(title = title, x = "", y = "log2(NormCounts)")
}

# -----------------------
# 7. Generate plots
# -----------------------
p1 <- plot_gene_set(FIT3_genes, "FIT3")
p2 <- plot_gene_set(HAC1_genes, "HAC1")
p3 <- plot_gene_set(mito_genes_11, "Mitochondrion (11 genes)")

# -----------------------
# 8. Display plots
# -----------------------

print(p1 + p2 + p3)
