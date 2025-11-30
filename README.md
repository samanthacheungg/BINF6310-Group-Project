# BINF6310 Group Project: Yeast Aging Heterogeneity Analysis

Reproduction study of "Single-cell RNA-seq reveals early heterogeneity during aging in yeast" (Wang et al., 2022, Aging Cell)

## Project Overview

This project analyzes single-cell RNA-seq data from aging yeast cells to understand transcriptional heterogeneity during the aging process. We aim to reproduce key findings from the original paper, particularly Figure 3, which demonstrates differential gene expression patterns between fast and slow-dividing cells.

## Team Members & Roles

- **Role 1 (Project Coordinator & QC Analyst)**: Quality control, cell filtering, project coordination
- **Role 2 (Data Engineer)**: Environment setup, data download, normalization pipeline
- **Role 3 (Differential Expression Analyst - Eric)**: DESeq2 analysis, subgroup definitions, DEG heatmap
- **Role 4 (Variability Analyst - Kyle)**: Correlation analyses, FIT3 analysis, figure assembly
- **Role 5 (Pathway Analyst - Kiren)**: Gene expression boxplots, pathway analysis, documentation

## Repository Structure

```
BINF6310-Group-Project/
├── README.md                   # This file
├── scripts/                    # Analysis scripts
│   ├── qc_filter_cells.R      # Quality control and filtering
│   ├── normalize_counts.R      # DESeq2 normalization
│   ├── figure3a_16h_correlation.R    # Figure 3a generation
│   ├── figure3b_36h_correlation.R    # Figure 3b generation
│   ├── figure3f_fit3_correlation.R   # Figure 3f FIT3 analysis
│   ├── kyle_explore_data.R    # Data exploration
│   └── ...
├── figures/                    # Generated figures (PDF and PNG)
│   ├── figure3a_correlation_16h.pdf
│   ├── figure3b_correlation_36h.pdf
│   ├── figure3f_fit3_16h.pdf
│   ├── figure3f_fit3_36h.pdf
│   └── ...
├── processed_data/            # Processed datasets
│   ├── log2_normalized_counts.csv
│   ├── metadata.csv
│   ├── correlation_statistics.csv
│   ├── fit3_correlation_stats.csv
│   └── ...
├── renv/                      # R environment management
├── renv.lock                  # Package dependencies lockfile
└── .Rprofile                  # Auto-loads renv

```

## Getting Started

### Prerequisites

- R version 4.4.2
- RStudio (recommended)
- Git

### Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd BINF6310-Group-Project
```

2. Open R or RStudio in the project directory. The renv environment will automatically activate.

3. Restore R packages from lockfile:
```r
renv::restore()
```

### Running the Analysis

The analysis pipeline is organized into sequential scripts:

1. **Quality Control & Filtering**:
   ```r
   source("scripts/qc_filter_cells.R")
   ```

2. **Normalization**:
   ```r
   source("scripts/normalize_counts.R")
   ```

3. **Correlation Analyses** (Figure 3a, 3b):
   ```r
   source("scripts/figure3a_16h_correlation.R")
   source("scripts/figure3b_36h_correlation.R")
   ```

4. **FIT3 Correlation Analysis** (Figure 3f):
   ```r
   source("scripts/figure3f_fit3_correlation.R")
   ```

## Key Analyses Completed

### Figure 3a: Genes Detected vs Generation (16h timepoint)
- **Correlation**: R = 0.62, p = 2.6×10⁻⁶
- **Interpretation**: Cells with more divisions (higher generation) tend to have more genes detected
- **Output**: `figures/figure3a_correlation_16h.pdf`

### Figure 3b: Genes Detected vs Generation (36h timepoint)
- **Correlation**: R = 0.57, p = 1.6×10⁻⁵
- **Interpretation**: Similar positive correlation observed at later timepoint
- **Output**: `figures/figure3b_correlation_36h.pdf`

### Figure 3f: FIT3 Expression vs Generation
- **16h timepoint**: R = -0.55, p = 1.3×10⁻⁴
- **36h timepoint**: R = -0.62, p = 5.6×10⁻⁶
- **Interpretation**: FIT3 expression decreases with increased cell division (aging marker)
- **Output**: `figures/figure3f_fit3_16h.pdf`, `figures/figure3f_fit3_36h.pdf`

## Data

### Source
Data from GEO accession: [Insert GEO accession number]

### Processed Data Files
- `log2_normalized_counts.csv`: Log2-transformed DESeq2 normalized counts (genes × cells)
- `metadata.csv`: Cell-level metadata including timepoint, generation, sample info
- `correlation_statistics.csv`: Statistical results for correlation analyses
- `fit3_correlation_stats.csv`: FIT3-specific correlation statistics

## R Environment

This project uses `renv` for reproducible package management:

- **R Version**: 4.4.2
- **Package Repository**: Posit Package Manager (CRAN latest)
- **Key Packages**: DESeq2, dplyr, ggplot2, tidyverse

### Adding New Packages

```r
# Install new package
renv::install("package-name")

# Update lockfile
renv::snapshot()
```

## Project Timeline

- **Week 1**: Setup, QC, and normalization
- **Week 2**: Subgroup analysis and DEG
- **Week 3**: Gene analysis and figure assembly
- **Week 4**: TF analysis and finalization

## Figure 3 Components Status

- [x] Figure 3a: 16h correlation (Kyle - Complete)
- [x] Figure 3b: 36h correlation (Kyle - Complete)
- [x] Figure 3c: DEG heatmap (Eric)
- [ ] Figure 3d: Gene boxplots by age (Kiren)
- [ ] Figure 3e: Gene boxplots by subgroup (Kiren)
- [x] Figure 3f: FIT3 correlation (Kyle - Complete)
- [ ] Final assembled Figure 3 (Kyle - Pending all panels)

## Reproducibility Notes

All analyses are designed to be fully reproducible:

1. Random seeds are set where applicable
2. Package versions locked via `renv.lock`
3. Scripts are self-contained and documented
4. Data provenance tracked

## Citations

Wang, X., et al. (2022). Single-cell RNA-seq reveals early heterogeneity during aging in yeast. *Aging Cell*, 21(e13595). https://doi.org/10.1111/acel.13595

## Contact

For questions about this project, please contact team members via the course communication channels.

## License

This project is for educational purposes as part of BINF6310 coursework.
