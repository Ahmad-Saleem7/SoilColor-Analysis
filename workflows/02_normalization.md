# Workflow: Data Normalization
## Stage 2 of 5 — Soil Pathology Color Analysis (MSU)

## Objective
For each variable flagged `NON-NORMAL` in Stage 1, systematically apply a suite of mathematical transformations and identify which (if any) achieves normality. Produce a normalized dataset and before/after visualizations. Assign each variable to either the parametric or non-parametric analysis path for Stages 4 and 5.

## Prerequisite
Stage 1 must have been run first:
```r
source("tools/01_distribution_analysis.R")
```
The tool reads `.tmp/01_distributions/distribution_summary.csv` to identify which variables need normalization.

## Tool
```
tools/02_normalization.R
```

Run from the project root:
```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/02_normalization.R")
```

Packages are auto-installed: `readxl`, `bestNormalize`, `nortest`, `ggplot2`, `patchwork`, `dplyr`, `moments`, `MASS`.

## Transformations Applied (in order)

| Name | Method | Handles |
|------|--------|---------|
| `log1p_shift` | log(x + shift + 1) | Right-skewed, positive data; shift ensures strict positivity |
| `sqrt_shift` | sqrt(x + shift) | Mild right skew; less aggressive than log |
| `boxcox` | Box-Cox optimal lambda (MASS) | Flexible power transformation for strictly positive data |
| `yeo_johnson` | Yeo-Johnson transform (bestNormalize) | Handles negative values; more general than Box-Cox |
| `order_norm` | Rank-based inverse normal (bestNormalize) | Last resort; forces normality but loses original scale |

Selection criterion: transformation that maximizes `min(SW_p, AD_p)`. If the winner achieves both p > 0.05, the variable is classified `NORMAL`; otherwise `REQUIRES-NONPARAMETRIC`.

## Outputs
All saved to `.tmp/02_normalization/`:

| File | Description |
|------|-------------|
| `normalized_data.csv` | Full dataset with original columns plus `<var>_norm` columns for transformed values |
| `transformation_summary.csv` | Best transformation per variable, original and post-transformation p-values, final status |
| `<var>_normalization.png` | 2×2 panel: original histogram + Q-Q (top) vs. transformed histogram + Q-Q (bottom) |

## Plot Color Legend
| Color | Meaning |
|-------|---------|
| Red | Original (non-normal) |
| Green | Transformed — achieved normality |
| Orange | Transformed — still non-parametric |

## Decision Gate

| Final_Status | Downstream Tests |
|--------------|-----------------|
| `NORMAL` | Pearson correlation (Stage 4), ANOVA + Tukey HSD (Stage 5) |
| `REQUIRES-NONPARAMETRIC` | Spearman correlation (Stage 4), Kruskal-Wallis + Dunn's test (Stage 5) |

## How to Use the Normalized Dataset
In Stages 3–5, use `<var>_norm` columns (from `normalized_data.csv`) for variables that were transformed to normality, and the original `<var>` columns for variables that were already normal in Stage 1. Non-parametric tests always use the original scale.

## Edge Cases
- **All variables already normal**: Tool detects empty `non_normal` list and exits cleanly. Skip to Stage 3.
- **Transformation errors**: Each transform is wrapped in `tryCatch`; failures are logged and that method is skipped.
- **order_norm always works** but ranks data, losing the original measurement scale. Only use it if interpretability of the scale is less important than satisfying downstream test assumptions.
- **If even order_norm fails to normalize**: Flag the variable as `REQUIRES-NONPARAMETRIC` and use rank-based tests throughout.

## Common Issues
| Problem | Fix |
|---------|-----|
| `Stage 1 output not found` | Run `source("tools/01_distribution_analysis.R")` first |
| `bestNormalize` install fails | Install manually: `install.packages("bestNormalize")` |
| Transformed Q-Q still looks curved | Accept `REQUIRES-NONPARAMETRIC` status; non-parametric tests are valid and robust |
