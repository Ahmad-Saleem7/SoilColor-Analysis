# Workflow: Soil Data — Pearson Correlation Matrix
## Stage 6

## Objective
Compute publication-ready Pearson correlation matrices for the z-score standardized soil dataset. Produce a full variable matrix and a color-variable sub-matrix.

## Input
| File | Location |
|------|----------|
| `soil_standardized.csv` | `.tmp/05_soil_homogeneity/` |

## Tool
```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/06_correlation_soil.R")
```

## Method
- **Pearson correlation** on z-score standardized variables
- **p-values**: `Hmisc::rcorr()` with pairwise complete observations
- **Multiple testing**: Benjamini-Hochberg FDR correction at α=0.05
- Missing data handled via pairwise complete observations

## Outputs

### A — Full Matrix
All numeric soil variables (soil physical + chemical + color).
Expected ~18–22 variables depending on data completeness after Stage 4 filtering.

### B — Color Sub-Matrix
Variables matching color-related patterns: RGB (R, G, B), xyY (x, y, Y), Munsell (Chroma, Value).
Focuses the analysis on how different color measurement systems relate to each other.

## Figure Specifications
- Color palette: RdBu diverging
- Lower triangle; hierarchical clustering
- Cell labels: r (2 dp); blank = p ≥ 0.05 (FDR)
- Full matrix: ~10×9.5 in at 300 DPI
- Color sub-matrix: 7×6.5 in at 300 DPI

## Outputs (`.tmp/06_correlation_soil/`)
| File | Content |
|------|---------|
| `soil_correlation_matrix_full.png` | Full soil correlation heatmap |
| `soil_correlation_matrix_color.png` | Color variables sub-matrix |
| `soil_corr_full_r.csv` | Full Pearson r matrix |
| `soil_corr_full_pval_fdr.csv` | Full FDR-adjusted p-values |
| `soil_corr_color_r.csv` | Color sub-matrix r |
| `soil_corr_color_pval_fdr.csv` | Color sub-matrix FDR p-values |
