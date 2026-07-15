# Workflow: Soil Color Correlation Matrix
## Stage 3 — Parametric Pearson Correlation

## Objective
Compute a publication-ready Pearson correlation matrix for all 9 soil color variables using the normalized dataset from Stages 1–2. Identify significant pairwise correlations with FDR correction.

## Input
| File | Location |
|------|----------|
| `normalized_data.csv` | Project root |

## Tool
```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/03_correlation_plant.R")
```

## Variables (9 color channels)
| Variable | Color Space | Source Column |
|----------|-------------|---------------|
| R | RGB | `R_norm` (OrderNorm-transformed) |
| G | RGB | `G` (original — normal) |
| B | RGB | `B` (original — normal) |
| H | HSB | `H_norm` |
| S | HSB | `S_norm` |
| Br | HSB | `Br_norm` |
| L | LAB | `L` (original — normal) |
| a | LAB | `a_norm` |
| b | LAB | `b_norm` |

## Method
- **Pearson correlation** (justified: all variables normalized in Stages 1–2)
- **p-values**: `Hmisc::rcorr()` — pairwise complete observations (n=76, no missing data)
- **Multiple testing**: Benjamini-Hochberg FDR correction at α=0.05

## Figure Specifications
- Color palette: RdBu diverging (#2166AC → white → #D6604D)
- Lower triangle only; hierarchical clustering order
- Cell labels: correlation coefficient (2 dp)
- Non-significant cells (FDR p ≥ 0.05): blanked
- Resolution: 300 DPI | Size: 8 × 7.5 in

## Outputs (`.tmp/03_correlation_plant/`)
| File | Content |
|------|---------|
| `plant_correlation_matrix.png` | Publication-ready heatmap |
| `plant_corr_r.csv` | 9×9 Pearson r matrix |
| `plant_corr_pval_raw.csv` | 9×9 raw p-values |
| `plant_corr_pval_fdr.csv` | 9×9 FDR-adjusted p-values |
