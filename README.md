# Soil Pathology Color Analysis — Michigan State University

A reproducible statistical pipeline for multi-color-space analysis of soil pathology samples. All computation is in R, producing publication-ready figures and structured CSV outputs at each stage.

---

## Overview

This repository analyzes color measurements (RGB, HSB, LAB, xyY, Munsell) collected from two sources:

- **Soil samples** — 76 specimens measured across 9 mean color variables (R, G, B, H, S, Br, L, a, b), sourced from Michigan State University
- **Soil samples** — 95 specimens with soil physical/chemical properties alongside color measurements (RGB, xyY, Munsell Chroma/Value), sourced from a multi-institution dataset (SchoolLondon, NRI/UoG, Earthwatch)

The pipeline tests for normality, applies the best-fit normalization transformation where needed, assesses compositional homogeneity, and produces Pearson correlation matrices — all at publication quality (300 DPI, RdBu palette, FDR-corrected).

---

## Tech Stack

| Layer | Technology |
|---|---|
| Language | R 4.4.x |
| Data I/O | `readxl`, `janitor` |
| Distribution fitting | `fitdistrplus` (AIC comparison: normal, log-normal, gamma, Weibull) |
| Normality testing | `nortest` (Anderson-Darling), base R `shapiro.test` |
| Normalization transforms | `bestNormalize` (OrderNorm, Yeo-Johnson), `MASS` (Box-Cox) |
| Homogeneity testing | `car::leveneTest` (center = median) |
| Correlation | `Hmisc::rcorr` (Pearson, pairwise complete) + `p.adjust` (BH/FDR) |
| Visualization | `ggplot2`, `ggcorrplot`, `GGally`, `patchwork` |
| Pipeline scripting | R scripts with automatic package bootstrapping |

---

## Repository Structure

```
SoilColor-Analysis/
├── data/                           # Ignored in git
│   ├── raw/                        # Original Excel databases
│   │   ├── plant_color_data.xlsx
│   │   └── soil_color_database.xlsx
│   └── processed/                  # Intermediate normalized data
│       └── normalized_data.csv
│
├── output/                         # Generated plots and CSVs (Ignored in git)
│   ├── 01_distributions/
│   ├── 02_normalization/
│   ├── 03_correlation_plant/
│   ├── 04_soil_dist_norm/
│   ├── 05_soil_homogeneity/
│   ├── 06_correlation_soil/
│   └── 07_combined_correlation/
│
├── scripts/
│   ├── pipeline/                   # R scripts — one per pipeline stage
│   │   ├── 01_plant_distribution.R
│   │   ├── 02_plant_normalization.R
│   │   ├── 03_plant_correlation.R
│   │   ├── 04_soil_distribution_normalization.R
│   │   ├── 05_soil_homogeneity.R
│   │   ├── 06_soil_correlation.R
│   │   └── 07_combined_correlation.R
│   └── setup/                      # Environment setup scripts
│       └── install_packages.R
│
├── workflows/                      # Markdown SOPs — stage objectives, inputs, edge cases
│   ├── 00_master_pipeline.md
│   ├── 01_distribution_analysis.md
│   ├── 02_normalization.md
│   ├── 03_correlation_plant.md
│   ├── 04_soil_dist_norm.md
│   ├── 05_soil_homogeneity.md
│   ├── 06_correlation_soil.md
│   └── 07_combined_correlation.md
│
├── .gitignore
└── README.md
```

---

## Pipeline

```
Soil data (76 samples)            Soil data (95 samples)
        │                                   │
        ▼                                   ▼
Stage 1: Distribution analysis     Stage 4: Distribution analysis
  Shapiro-Wilk + Anderson-Darling    Shapiro-Wilk + Anderson-Darling
  AIC: normal/lnorm/gamma/Weibull    AIC: normal/lnorm/gamma/Weibull
        │                                   │
        ▼ (non-normal vars)                 ▼ (non-normal vars)
Stage 2: Normalization             Stage 4: Normalization (same tool)
  log1p · sqrt · Box-Cox             OrderNorm selected for all 24/25
  Yeo-Johnson · OrderNorm            non-normal variables
        │                                   │
        ▼                                   ▼
Stage 3: Pearson correlation       Stage 5: Levene's homogeneity test
  9×9 matrix (soil color vars)      by Carbonate Class (5 groups)
  FDR-corrected (BH)                 → z-score standardization (all vars)
  34/36 pairs significant                    │
                                             ▼
                                   Stage 6: Pearson correlation
                                     25×25 full matrix (119/300 sig.)
                                     8×8 color sub-matrix
                                             │
                           ┌─────────────────┘
                           ▼
                  Stage 7: Combined correlation
                    Row-bind soil + soil on R, G, B, Lightness
                    n = 170 | 4×4 matrix | all 6 pairs significant
```

---

## Running the Pipeline

Ensure that you have your raw data placed in `data/raw/`:
- `data/raw/plant_color_data.xlsx`
- `data/raw/soil_color_database.xlsx`

Run scripts in order from the project root in R or from the terminal:

```r
setwd("C:/path/to/SoilColor-Analysis")

# Setup (if packages are not installed yet)
source("scripts/setup/install_packages.R")

# Soil Analysis
source("scripts/pipeline/01_plant_distribution.R")
source("scripts/pipeline/02_plant_normalization.R")
source("scripts/pipeline/03_plant_correlation.R")

# Soil Analysis
source("scripts/pipeline/04_soil_distribution_normalization.R")
source("scripts/pipeline/05_soil_homogeneity.R")
source("scripts/pipeline/06_soil_correlation.R")

# Combined Analysis
source("scripts/pipeline/07_combined_correlation.R")
```

Each tool is self-contained — it reads its inputs from `data/raw/` or `data/processed/`, bootstraps required packages if missing, and writes all results to `output/<stage>/`. Re-running any tool safely overwrites previous outputs.

---

## Key Results

| Stage | Output | Result |
|---|---|---|
| Soil distribution | 9 variables tested | 3 normal (G, B, L) · 6 non-normal |
| Soil normalization | 6 non-normal vars | All normalized via OrderNorm |
| Soil correlation | 9×9 Pearson matrix | 34/36 pairs significant (FDR p < 0.05) |
| Soil distribution | 25 variables tested | 1 normal · 24 non-normal |
| Soil normalization | 24 non-normal vars | 20 normalized · 4 non-parametric |
| Soil homogeneity | Levene's test, 5 groups | 16 homogeneous · 9 heterogeneous |
| Soil correlation | 25×25 Pearson matrix | 119/300 pairs significant |
| Combined | 4×4 on R, G, B, Lightness (n=170) | All 6 pairs significant |
