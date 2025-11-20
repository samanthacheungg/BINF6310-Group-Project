# How to Create SECONDARY PANEL A

## Quick Start

The figure creation script cannot be run with `Rscript` due to renv/conda package conflicts.

### Method 1: Interactive R (Recommended)

bash
# 1. Activate environment
conda activate R_env

# 2. Navigate to project
cd /courses/BINF6310.202610/students/lashley.an/BINF6310-Group-Project

# 3. Open R
R

# 4. Run the script
source("scripts/02_create_secondary_panel_a.R")

# 5. Exit R
quit()


### Method 2: Run Commands Directly

See `scripts/02_create_secondary_panel_a.R` for the complete code.

## Required Packages

- ggplot2
- dplyr  
- gridExtra

Install with:
r
install.packages(c("ggplot2", "dplyr", "gridExtra"), repos="https://cloud.r-project.org")


## Output

- `figures/SECONDARY_PANEL_A_qc_metrics.pdf` (high-res, for reports)
- `figures/SECONDARY_PANEL_A_qc_metrics.png` (high-res, for slides)

## What the Figure Shows

3-panel visualization of QC metrics:
1. **Genes Detected:** Distribution with 1,000 gene threshold
2. **Total Reads:** Distribution with 40,000 read threshold  
3. **ERCC Proportion:** Distribution with 0.74 threshold

Cells are color-coded:
- Blue: Pass all QC criteria (125 cells)
- Red: Fail QC (11 cells)
