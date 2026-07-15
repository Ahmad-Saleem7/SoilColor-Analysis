# Workflow: Soil Data — Distribution Analysis & Normalization
## Stage 4

## Objective
Characterize the statistical distribution of all numeric variables in `soilcolor_database.xlsx`. Normalize non-normal variables using the best available transformation. Produce a normalized soil dataset for downstream homogeneity and correlation analyses.

## Input
| File | Location | Sheet |
|------|----------|-------|
| `soilcolor_database.xlsx` | Project root | "all data" |

## Tool
```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/04_soil_dist_norm.R")
```

## Data Cleaning
1. Column names cleaned to snake_case using `janitor::clean_names()`
2. Non-numeric and identifier columns excluded (ID, coordinates, text categories)
3. Columns with >50% missing values dropped
4. Variables with <15 valid observations skipped

## Tests Applied
- Shapiro-Wilk + Anderson-Darling normality tests
- Skewness and excess kurtosis
- AIC-based distribution fitting: normal, log-normal, gamma, Weibull

## Normalization Transformations (applied to NON-NORMAL variables)
log1p_shift → sqrt_shift → Box-Cox → Yeo-Johnson → OrderNorm
Best selected by max(min(SW_p, AD_p)).

## Missing Data
- Pairwise complete observations used throughout
- Missing values preserved as NA in output dataset

## Outputs (`.tmp/04_soil_dist_norm/`)
| File | Content |
|------|---------|
| `soil_normalized.csv` | Original + `*_norm` columns for transformed variables |
| `soil_distribution_summary.csv` | Per-variable normality tests, skewness, AIC, verdict |
| `soil_transformation_summary.csv` | Best transformation and post-transformation p-values |
| `<var>_dist.png` | Per-variable histogram + Q-Q plot |
| `00_distributions_page*.png` | Composite panels (8 variables per page) |
