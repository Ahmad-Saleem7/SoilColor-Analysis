# Workflow: Soil Data — Homogeneity Assessment & Standardization
## Stage 5

## Objective
Test whether soil variables have equal variance across Carbonate Class groups (NC/SC/EC) using Levene's test. Apply z-score standardization to all numeric variables to place them on a common scale for correlation analysis.

## Inputs
| File | Location |
|------|----------|
| `soil_normalized.csv` | `.tmp/04_soil_dist_norm/` |
| `soilcolor_database.xlsx` | Project root (for Carbonate Class grouping) |

## Tool
```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/05_soil_homogeneity.R")
```

## Grouping Variable
**Carbonate Class** (NC = Non-Calcareous, SC = Slightly Calcareous, EC = Extremely Calcareous).
Scientifically meaningful: carbonate content strongly influences soil color, organic matter, and pH-related properties.

## Levene's Test
- `car::leveneTest(variable ~ Carbonate_Class, center = median)`
- Median-based (robust to non-normality)
- Minimum group size: 5 observations
- Verdict: HOMOGENEOUS (p ≥ 0.05) or HETEROGENEOUS (p < 0.05)

## Z-Score Standardization
Applied to **all** numeric variables regardless of Levene result. Rationale: the dataset mixes bulk density (kg/l), carbon percentages (0–15%), and RGB values (0–255) — vastly different units and scales. Z-scoring ensures each variable contributes equally to the correlation analysis.

Z-score: `z = (x − mean(x)) / sd(x)` computed on non-missing values.

## Outputs (`.tmp/05_soil_homogeneity/`)
| File | Content |
|------|---------|
| `soil_standardized.csv` | Z-scored numeric variables + Carbonate_Class label |
| `homogeneity_summary.csv` | Per-variable Levene F-stat, p-value, verdict |
| `00_homogeneity_page*.png` | Boxplots by Carbonate Class group (4 per row) |

## Note
`soil_standardized.csv` is the input for both Stage 6 (soil correlation) and Stage 7 (combined correlation).
