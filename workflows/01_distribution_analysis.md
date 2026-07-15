# Workflow: Distribution Analysis
## Stage 1 of 5 — Soil Pathology Color Analysis (MSU)

## Objective
Characterize the statistical distribution of each color measurement variable. Determine whether variables are normally distributed and identify the best-fitting parametric distribution. This decision gates all downstream analysis choices (parametric vs. non-parametric).

## Inputs
| Item | Path | Required |
|------|------|----------|
| Color data | `color_analysis All samples.xlsx` (project root) | Yes |

## Tool
```
tools/01_distribution_analysis.R
```

Run from the project root in R or RStudio:
```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/01_distribution_analysis.R")
```

Packages are auto-installed on first run: `readxl`, `fitdistrplus`, `nortest`, `ggplot2`, `patchwork`, `dplyr`, `moments`.

## Variables Analyzed
The 9 mean color measurements:

| Variable | Color Space | Description |
|----------|-------------|-------------|
| R | RGB | Red channel mean |
| G | RGB | Green channel mean |
| B | RGB | Blue channel mean |
| H | HSB | Hue mean |
| S | HSB | Saturation mean |
| Br | HSB | Brightness mean |
| L | LAB | Lightness mean |
| a | LAB | a* chroma mean |
| b | LAB | b* chroma mean |

## Tests Applied
| Test | When | Interpretation |
|------|------|----------------|
| Shapiro-Wilk (W) | All variables; n=76 | p > 0.05 = consistent with normality |
| Anderson-Darling | All variables; more powerful for n>50 | p > 0.05 = consistent with normality |
| Skewness | All variables | Near 0 = symmetric; >1 or <-1 = strongly skewed |
| Excess kurtosis | All variables | Near 0 = mesokurtic; >1 = heavy tails |

**Normality decision rule**: A variable is classified `NORMAL` only if **both** SW p > 0.05 **and** AD p > 0.05. Otherwise it is `NON-NORMAL`.

## Distribution Fitting
Candidate distributions fitted via maximum likelihood (`fitdistrplus::fitdist`):
- Normal (`norm`) — fitted on full data
- Log-normal (`lnorm`) — fitted on positive values
- Gamma (`gamma`) — fitted on positive values
- Weibull (`weibull`) — fitted on positive values

Model selection is by lowest AIC. The best-fit distribution is overlaid on histograms.

## Outputs
All saved to `.tmp/01_distributions/`:

| File | Description |
|------|-------------|
| `distribution_summary.csv` | One row per variable: descriptives, SW/AD test statistics and p-values, skewness, excess kurtosis, AIC for each candidate distribution, best-fit label, normality verdict |
| `<var>_distribution.png` | Histogram with fitted density + Normal Q-Q plot for each variable |
| `00_all_distributions_composite.png` | 9-panel composite figure (publication-ready) |

## Decision Gate
The `Normality` column in `distribution_summary.csv` controls downstream tools:

| Verdict | Next Action |
|---------|-------------|
| `NORMAL` | Variable uses **parametric** tests in Stages 4 & 5 (Pearson, ANOVA, Tukey HSD) |
| `NON-NORMAL` | Variable passes to Stage 2 (normalization); if un-normalizable → non-parametric tests |

## Edge Cases
- **All variables normal**: Skip Stage 2; proceed directly to Stage 3.
- **Fitting failure** (e.g., negative values for log-normal): Tool catches errors; that distribution is excluded from AIC comparison.
- **Small p but mild skewness**: Trust the tests over visual impression — both SW and AD must pass.

## Common Issues
| Problem | Fix |
|---------|-----|
| `Cannot read Excel file` | Confirm `color_analysis All samples.xlsx` is in the project root and not open in Excel |
| Package install fails | Run `install.packages(c("readxl","fitdistrplus","nortest","ggplot2","patchwork","dplyr","moments"))` manually |
| Column name mismatch | Check actual column names in the Excel file; update `MEAN_VARS` in the script if different |
