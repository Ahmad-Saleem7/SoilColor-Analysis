# Master Pipeline — Soil Pathology Color Analysis (MSU)

## Objective
Complete statistical analysis of color space measurements (RGB, HSB, LAB) from soil pathology samples received from Michigan State University. Determine variable distributions, normalize non-normal data, assess composition and homogeneity, quantify relationships, and perform group comparisons.

## Data
- **Input file**: `color_analysis All samples.xlsx` (project root)
- **Samples**: 76 rows
- **Mean variables (9)**: R, G, B, H, S, Br, L, a, b
- **StdDev variables (9)**: Rstdv, Gstdv, Bstdv, Hstdev, Sstdev, Brstdev, Lstdev, astdev, bastdev
- **Grouping variable**: To be added when full MSU dataset arrives (treatment / disease level / variety)

## How to Run
Open R or RStudio, set the working directory to the project root, then source tools in order:

```r
setwd("c:/Users/PMLS/Desktop/pp_msu")
source("tools/01_distribution_analysis.R")
source("tools/02_normalization.R")
# Stages 3-5 added after MSU grouped data arrives
```

---

## Pipeline Stages

### Stage 1 — Distribution Analysis
- **Tool**: `tools/01_distribution_analysis.R`
- **Workflow**: `workflows/01_distribution_analysis.md`
- **Input**: `color_analysis All samples.xlsx`
- **Outputs** (in `.tmp/01_distributions/`):
  - `distribution_summary.csv` — normality test results, skewness, kurtosis, AIC values, best-fit distribution
  - `<var>_distribution.png` — per-variable histogram + Q-Q plot
  - `00_all_distributions_composite.png` — 9-panel composite figure
- **Decision gate**: Variables tagged `NORMAL` take the parametric path. Variables tagged `NON-NORMAL` go to Stage 2.

### Stage 2 — Normalization
- **Tool**: `tools/02_normalization.R`
- **Workflow**: `workflows/02_normalization.md`
- **Input**: `distribution_summary.csv` + raw Excel file
- **Outputs** (in `.tmp/02_normalization/`):
  - `normalized_data.csv` — original dataset plus `*_norm` transformed columns
  - `transformation_summary.csv` — best transformation per variable, before/after p-values
  - `<var>_normalization.png` — before/after histogram + Q-Q comparison
- **Decision gate**: Variables achieving normality after transformation join the parametric path. Variables flagged `REQUIRES-NONPARAMETRIC` use non-parametric equivalents throughout.

### Stage 3 — Composition Analysis *(awaiting MSU grouped data)*
- **Tool**: `tools/03_composition_analysis.R`
- Levene's test (robust) and Bartlett's test for homogeneity of variance
- Z-score or min-max standardization across color spaces
- PCA for dimensionality reduction and variance composition visualization

### Stage 4 — Relationship Analysis *(awaiting MSU grouped data)*
- **Tool**: `tools/04_correlation_analysis.R`
- Pearson correlation for parametric variables
- Spearman correlation for non-parametric variables
- Corrplot heatmap; cross-color-space correlation analysis

### Stage 5 — Comparison Tests *(awaiting MSU grouped data)*
- **Tool**: `tools/05_comparison_tests.R`
- Parametric path: one-way ANOVA + Tukey HSD / Bonferroni / Dunnett post-hoc
- Non-parametric path: Kruskal-Wallis + Dunn's test with Bonferroni correction
- Effect size reporting (eta-squared for ANOVA, epsilon-squared for Kruskal-Wallis)

---

## Decision Tree

```
Raw Data (color_analysis All samples.xlsx)
          |
          v
    Stage 1: Distribution Analysis
          |
          +-- NORMAL ──────────────────────────────┐
          |                                        |
          +-- NON-NORMAL ──> Stage 2: Normalize    |
                                  |                |
                            NORMAL after transform─┤ PARAMETRIC PATH
                            REQUIRES-NONPARAMETRIC─┤ NON-PARAMETRIC PATH
                                                   |
                                                   v
                                         Stage 3: Composition
                                                   |
                                                   v
                                         Stage 4: Relationships
                                         (Pearson | Spearman)
                                                   |
                                                   v
                                         Stage 5: Comparisons
                                         (ANOVA+Tukey | KW+Dunn)
```

## Output Convention
- All intermediate and final outputs go to `.tmp/<stage>/`
- `.tmp/` is not committed to version control (see `.gitignore`)
- Re-running any tool is safe — it overwrites previous outputs
- Final publication figures should be exported from `.tmp/` to a cloud destination (Google Drive, Slides) for sharing
