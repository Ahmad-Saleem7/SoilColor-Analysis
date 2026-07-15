# Workflow: Combined Correlation Matrix (Soil + Soil)
## Stage 7

## Objective
Merge the normalized soil color dataset and the standardized soil dataset on shared RGB color space variables. Compute a combined Pearson correlation matrix and a scatter matrix to compare how color channels behave across soil and soil measurement contexts.

## Inputs
| File | Location |
|------|----------|
| `normalized_data.csv` | Project root |
| `soil_standardized.csv` | `.tmp/05_soil_homogeneity/` |

## Tool
```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/07_combined_correlation.R")
```

## Combination Strategy
Both datasets share RGB color space. Combined by row-binding:

| Variable | Soil Source | Soil Source |
|----------|-------------|-------------|
| R | `R_norm` (z-scaled) | `rgb_r` (already z-scored) |
| G | `G` (z-scaled) | `rgb_g` |
| B | `B` (z-scaled) | `rgb_b` |
| Lightness | `L` (z-scaled, LAB) | `xyY.Y` or Munsell Value |

A `Source` indicator column (Soil/Soil) is added. Combined n = 76 + 97 = 173.

All soil variables are z-scored in this tool to match the z-score scale of the soil standardized dataset.

## Method
- **Pearson correlation** on shared z-scored color variables
- **p-values**: `Hmisc::rcorr()` with pairwise complete observations
- **Multiple testing**: Benjamini-Hochberg FDR at α=0.05

## Scientific Rationale
The combined correlation reveals whether RGB channel co-variation is consistent across soil and soil color measurement contexts, and whether Lightness/luminance follows similar patterns in both domains.

## Outputs (`.tmp/07_combined_correlation/`)
| File | Content |
|------|---------|
| `combined_correlation_matrix.png` | Publication-ready Pearson heatmap |
| `combined_scatter_matrix.png` | GGally scatter matrix colored by Source |
| `combined_dataset.csv` | Row-bound Soil + Soil data (with Source) |
| `combined_corr_r.csv` | Pearson r matrix |
| `combined_corr_pval_fdr.csv` | FDR-adjusted p-values |
