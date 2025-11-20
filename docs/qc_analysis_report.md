# QC Analysis Report

**Analyst:** Anaita Lashley  
**Date:** November 19, 2025  
**Project:** BINF6310 Group Project - Yeast Aging Reproducibility Study

---

## Executive Summary

Quality control filtering was performed on 136 single yeast cells following the methodology described in Wang et al. (2022). After applying three filtering criteria, **125 cells (91.9%) passed QC**, exactly matching the paper's reported cell count.

---

## Methods

### QC Filtering Criteria

Three quality control filters were applied to remove low-quality cells:

1. **Genes Detected > 1,000**
   - Rationale: Cells with <1,000 detected genes likely represent failed library preparations or damaged cells
   
2. **Total Mapped Reads > 40,000**
   - Rationale: Low read depth results in noisy gene expression measurements
   
3. **ERCC Spike-in Proportion < 0.74**
   - Rationale: High ERCC proportions indicate cell lysis or poor RNA capture

Cells passing all three criteria were retained for downstream analysis.

---

## Results

### Overall QC Statistics

| Metric | Value |
|--------|-------|
| **Total cells sequenced** | 136 |
| **Cells passing QC** | 125 (91.9%) |
| **Cells failing QC** | 11 (8.1%) |
| **Mean genes detected (passing cells)** | 2,278 |
| **Paper reported (mean genes)** | 2,202 |
| **Difference from paper** | +3.4% |

### QC Criteria Breakdown

| Criterion | Cells Passing | Cells Failing |
|-----------|---------------|---------------|
| Mapped reads > 40,000 | 136 (100%) | 0 (0%) |
| Genes detected > 1,000 | 125 (91.9%) | 11 (8.1%) |
| ERCC proportion < 0.74 | 134 (98.5%) | 2 (1.5%) |
| **All criteria (AND)** | **125 (91.9%)** | **11 (8.1%)** |

### Primary Failure Mode

The primary reason for QC failure was **insufficient genes detected** (11 cells failed this criterion). All cells had sufficient read depth, and only 2 cells had excessive ERCC spike-in proportions.

---

## Validation Against Paper

### Cell Count Validation

| Metric | Our Analysis | Wang et al. 2022 | Match? |
|--------|--------------|------------------|--------|
| Total cells after QC | 125 | 125 | ✅ Yes |
| Mean genes detected | 2,278 | 2,202 | ✅ Within 10% |
| Percent difference | +3.4% | baseline | ✅ Acceptable |

**Conclusion:** Our QC filtering successfully reproduced the paper's cell count (125 cells) with mean gene detection within 3.4% of the reported value.

---

## Data Quality Assessment

### Strengths
- ✅ High proportion of cells passed QC (91.9%)
- ✅ All cells had sufficient sequencing depth (>40k reads)
- ✅ Low ERCC contamination (98.5% of cells <0.74 threshold)
- ✅ Mean genes detected close to paper (3.4% difference)

### Observations
- 11 cells failed due to low gene detection (<1,000 genes)
- 2 cells had elevated ERCC proportions (likely lysed cells)
- No cells failed due to insufficient read depth (good sequencing quality)

---

## Visualization

Quality control metrics are visualized in **SECONDARY PANEL A** (`figures/SECONDARY_PANEL_A_qc_metrics.pdf`), which displays:

1. **Panel 1:** Distribution of genes detected per cell
2. **Panel 2:** Distribution of total mapped reads per cell (log scale)
3. **Panel 3:** Distribution of ERCC spike-in proportions

Cells are color-coded:
- **Blue (Pass QC):** 125 cells meeting all criteria
- **Red (Fail QC):** 11 cells failing one or more criteria

Red dashed lines indicate QC thresholds.

---

## Files Generated

### Data Files
- `qc/qc_filtered_gene_counts.txt` - Gene expression matrix (125 cells × ~6,000 genes)
- `qc/qc_filtered_cell_qc_metrics.txt` - Cell metadata with QC flags
- `normalize/log2_normalized_counts.txt` - Log2-normalized expression values

### Analysis Scripts
- `scripts/00_verification.R` - QC validation script
- `scripts/02_create_secondary_panel_a.R` - Figure generation script

### Figures
- `figures/SECONDARY_PANEL_A_qc_metrics.pdf` - High-resolution QC visualization (300 dpi)
- `figures/SECONDARY_PANEL_A_qc_metrics.png` - Web-friendly version (300 dpi)

### Documentation
- `docs/how_to_run_figure_creation.md` - Instructions for reproducing figures
- `docs/qc_analysis_report.md` - This report

---

## Reproducibility Notes

### Environment
- **R version:** 4.4.3
- **Key packages:** ggplot2, dplyr, gridExtra, DESeq2
- **Conda environment:** R_env
- **Compute cluster:** Northeastern Discovery

### Challenges Encountered
- renv/conda package conflicts prevented running scripts with `Rscript`
- **Solution:** Scripts must be run interactively in R using `source()`
- See `docs/how_to_run_figure_creation.md` for detailed instructions

---

## Conclusions

1. ✅ **QC filtering successfully reproduced:** 125 cells identified, matching Wang et al. 2022 exactly
2. ✅ **Data quality is high:** 91.9% cell pass rate, mean genes within 3.4% of paper
3. ✅ **Primary QC failures:** Low gene detection (11 cells), not sequencing depth issues
4. ✅ **Ready for downstream analysis:** Filtered dataset available for Figure 3 reproduction

The QC validation confirms that the processed data is suitable for downstream analyses and reproduces the paper's filtering approach with high fidelity.

---

## References

Wang, J., Sang, Y., Jin, S., et al. (2022). Single-cell RNA-seq reveals early heterogeneity during aging in yeast. *Aging Cell*, 21(11), e13712. https://doi.org/10.1111/acel.13712
