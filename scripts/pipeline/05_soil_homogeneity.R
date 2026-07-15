# ============================================================
# Tool: Soil Data — Homogeneity Assessment + Z-score Standardization
# Project: Soil Pathology Color Analysis — Michigan State University
# Stage: 5
# Input:  output/04_soil_dist_norm/soil_normalized.csv
#         soilcolor_database.xlsx (for grouping variable)
# Output: output/05_soil_homogeneity/
# ============================================================

# ---------- Package Bootstrap ----------
# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))
pkgs <- c("car", "readxl", "janitor", "ggplot2", "patchwork", "dplyr", "tidyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(car); library(readxl); library(janitor)
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
})

# ---------- Config ----------
NORM_FILE  <- "output/04_soil_dist_norm/soil_normalized.csv"
EXCEL_FILE <- "data/raw/soil_color_database.xlsx"
SHEET      <- "all data"
OUT_DIR    <- "output/05_soil_homogeneity"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_VAR <- "carbonate_class"   # after janitor::clean_names()
ALPHA     <- 0.05
MIN_GRP   <- 5                   # minimum group size for Levene's test

# ---------- Load Normalized Data ----------
if (!file.exists(NORM_FILE))
  stop("Stage 4 output not found. Run scripts/pipeline/04_soil_distribution_normalization.R first.")

cat("Loading normalized data:", NORM_FILE, "\n")
df_norm <- read.csv(NORM_FILE, stringsAsFactors = FALSE)

# ---------- Attach Grouping Variable ----------
# Same 3-row header structure as in Tool 04
cat("Loading grouping variable from:", EXCEL_FILE, "\n")
hdr_raw  <- suppressMessages(
  read_excel(EXCEL_FILE, sheet = SHEET, skip = 1, n_max = 1, col_names = FALSE)
)
headers  <- as.character(unlist(hdr_raw[1, ]))
na_pos   <- is.na(headers) | trimws(headers) == ""
headers[na_pos] <- paste0("col_", which(na_pos))
raw_data <- suppressMessages(
  read_excel(EXCEL_FILE, sheet = SHEET, skip = 3, col_names = FALSE)
)
if (ncol(raw_data) != length(headers)) {
  n <- min(ncol(raw_data), length(headers))
  raw_data <- raw_data[, seq_len(n), drop = FALSE]
  headers  <- headers[seq_len(n)]
}
names(raw_data) <- headers
df_raw   <- janitor::clean_names(raw_data)

# Find the carbonate class column
grp_col  <- names(df_raw)[grep(GROUP_VAR, names(df_raw), ignore.case = TRUE)][1]
if (is.na(grp_col)) stop("Carbonate Class column not found. Check column names.")
cat("Grouping column found:", grp_col, "\n")

# Trim df_raw to match df_norm rows (same order assumed)
n_rows   <- min(nrow(df_norm), nrow(df_raw))
df_norm  <- df_norm[seq_len(n_rows), ]
groups   <- df_raw[[grp_col]][seq_len(n_rows)]
groups   <- trimws(as.character(groups))

# Remove rows where group is NA
valid_idx <- !is.na(groups) & groups != "" & groups != "NA"
df_norm   <- df_norm[valid_idx, ]
groups    <- groups[valid_idx]
groups    <- factor(groups)

cat("Groups:", paste(levels(groups), collapse = ", "),
    "| n =", nrow(df_norm), "\n")
print(table(groups))
cat("\n")

# Check minimum group sizes
grp_sizes <- table(groups)
if (any(grp_sizes < MIN_GRP)) {
  small_grps <- names(grp_sizes)[grp_sizes < MIN_GRP]
  cat("WARNING: Groups with < ", MIN_GRP, " obs:", paste(small_grps, collapse = ", "),
      "— Levene's test may be unreliable for these.\n\n")
}

# ---------- Select Numeric Variables for Analysis ----------
# Use *_norm columns where they exist; fall back to original numeric columns
all_cols   <- names(df_norm)
orig_num   <- all_cols[sapply(df_norm, is.numeric) & !grepl("_norm$", all_cols)]
norm_suffix <- paste0(orig_num, "_norm")

# For each original numeric variable: use _norm version if it exists
use_cols <- sapply(orig_num, function(v) {
  nm <- paste0(v, "_norm")
  if (nm %in% all_cols) nm else v
})
use_cols   <- unname(use_cols)
label_cols <- sub("_norm$", "", use_cols)   # clean display labels

analysis_df <- df_norm[, use_cols, drop = FALSE]
names(analysis_df) <- label_cols

# Remove columns with zero variance or all-NA
ok <- sapply(analysis_df, function(col) {
  x <- col[is.finite(col)]
  length(x) >= MIN_GRP & sd(x, na.rm = TRUE) > 0
})
analysis_df <- analysis_df[, ok, drop = FALSE]
cat("Variables for homogeneity test:", ncol(analysis_df), "\n\n")

# ---------- Levene's Test Per Variable ----------
cat("Running Levene's test (center = median) per variable...\n\n")

levene_records <- list()

for (v in names(analysis_df)) {
  x   <- analysis_df[[v]]
  idx <- is.finite(x) & !is.na(groups)
  xi  <- x[idx]; gi <- groups[idx]

  if (length(unique(gi)) < 2 || length(xi) < MIN_GRP * 2) {
    cat(sprintf("  [SKIP] %-30s — insufficient data\n", v))
    next
  }

  lt <- tryCatch(
    car::leveneTest(xi ~ gi, center = median),
    error = function(e) NULL
  )
  if (is.null(lt)) {
    cat(sprintf("  [SKIP] %-30s — Levene error\n", v)); next
  }

  f_val <- lt[["F value"]][1]
  p_val <- lt[["Pr(>F)"]][1]
  verd  <- if (!is.na(p_val) && p_val < ALPHA) "HETEROGENEOUS" else "HOMOGENEOUS"

  cat(sprintf("  %-30s  F=%.3f  p=%.4f  -> %s\n", v, f_val, p_val, verd))

  levene_records[[v]] <- data.frame(
    Variable  = v,
    F_stat    = round(f_val, 4),
    p_value   = round(p_val, 4),
    Verdict   = verd,
    stringsAsFactors = FALSE
  )
}

homog_summ <- bind_rows(levene_records)
out_homog  <- file.path(OUT_DIR, "homogeneity_summary.csv")
write.csv(homog_summ, out_homog, row.names = FALSE)
cat("\nHomogeneity summary saved:", out_homog, "\n\n")

n_het <- sum(homog_summ$Verdict == "HETEROGENEOUS")
n_hom <- sum(homog_summ$Verdict == "HOMOGENEOUS")
cat("Homogeneous    :", n_hom, "variables\n")
cat("Heterogeneous  :", n_het, "variables\n\n")

# ---------- Z-Score Standardization (all variables) ----------
# Z-score is applied regardless of Levene result so that all variables
# (spanning different measurement scales and units) are comparable
# for the subsequent correlation analysis.

cat("Applying z-score standardization to all numeric variables...\n")

std_df <- as.data.frame(
  scale(analysis_df),   # scale() handles NAs column-wise
  stringsAsFactors = FALSE
)
std_df[["carbonate_class"]] <- as.character(groups)

out_std <- file.path(OUT_DIR, "soil_standardized.csv")
write.csv(std_df, out_std, row.names = FALSE)
cat("Standardized dataset saved:", out_std, "\n\n")

# ---------- Visualization: Boxplots by Carbonate Class ----------
cat("Generating group boxplots...\n")

VERD_COLORS <- c(HOMOGENEOUS = "#2980B9", HETEROGENEOUS = "#E74C3C")
GROUP_FILLS <- c(NC = "#F1C40F", SC = "#E67E22", EC = "#8E44AD",
                 # fallback
                 `Non-Calcareous` = "#F1C40F",
                 `Slightly Calcareous` = "#E67E22",
                 `Extremely Calcareous` = "#8E44AD")

make_boxplot <- function(var_nm, verd) {
  x_df <- data.frame(
    val = analysis_df[[var_nm]],
    grp = groups
  ) %>% filter(is.finite(val))

  bord <- VERD_COLORS[[verd]]

  ggplot(x_df, aes(x = grp, y = val, fill = grp)) +
    geom_boxplot(outlier.size = 1.2, outlier.alpha = 0.5,
                 color = bord, linewidth = 0.6) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(title    = var_nm,
         subtitle = paste("Levene:", verd),
         x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(plot.title    = element_text(face = "bold", size = 9,
                                       color = bord),
          plot.subtitle = element_text(size = 7, color = "grey45"),
          axis.text.x   = element_text(angle = 30, hjust = 1))
}

bp_vars <- homog_summ$Variable
chunks  <- split(bp_vars, ceiling(seq_along(bp_vars) / 8))

for (i in seq_along(chunks)) {
  panels <- mapply(
    function(v, verd) make_boxplot(v, verd),
    v    = chunks[[i]],
    verd = homog_summ$Verdict[homog_summ$Variable %in% chunks[[i]]],
    SIMPLIFY = FALSE
  )
  composite <- wrap_plots(panels, ncol = 4) +
    plot_annotation(
      title    = paste0("Soil Homogeneity — Levene's Test by Carbonate Class (Page ", i, ")"),
      subtitle = "Blue border = HOMOGENEOUS  |  Red border = HETEROGENEOUS",
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 9, color = "grey40"))
    )
  comp_png <- file.path(OUT_DIR, sprintf("00_homogeneity_page%02d.png", i))
  ggsave(comp_png, composite,
         width  = 14,
         height = 4 * ceiling(length(chunks[[i]]) / 4),
         dpi    = 150)
  cat("Boxplot page", i, "saved:", comp_png, "\n")
}

# ---------- Final Report ----------
cat("\n========================================\n")
cat("STAGE 5 COMPLETE — Soil Homogeneity & Standardization\n")
cat("========================================\n")
cat("Homogeneous   :", n_hom, "| Heterogeneous:", n_het, "\n")
cat("All variables z-score standardized for downstream correlation.\n")
cat("Outputs in:", OUT_DIR, "\n")
cat("  soil_standardized.csv    - ready for Stage 6 correlation\n")
cat("  homogeneity_summary.csv  - Levene test results\n")
cat("\nNext: source('scripts/pipeline/06_soil_correlation.R')\n")
