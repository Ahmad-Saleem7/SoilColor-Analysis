# ============================================================
# Tool: Correlation Matrix — Soil Data (Parametric)
# Project: Soil Pathology Color Analysis — Michigan State University
# Stage: 6
# Input:  output/05_soil_homogeneity/soil_standardized.csv
# Output: output/06_correlation_soil/
# ============================================================

# ---------- Package Bootstrap ----------
# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))
pkgs <- c("Hmisc", "ggcorrplot", "ggplot2", "dplyr", "RColorBrewer")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(Hmisc); library(ggcorrplot); library(ggplot2)
  library(dplyr); library(RColorBrewer)
})

# ---------- Config ----------
STD_FILE <- "output/05_soil_homogeneity/soil_standardized.csv"
OUT_DIR  <- "output/06_correlation_soil"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05

# Munsell & color variable patterns (for sub-matrix)
COLOR_PATTERNS <- c("rgb", "xy_y", "xyy", "munsell", "chroma", "value",
                    "colour", "color", "hue")

# ---------- Load Standardized Data ----------
if (!file.exists(STD_FILE))
  stop("Stage 5 output not found. Run scripts/pipeline/05_soil_homogeneity.R first.")

cat("Loading:", STD_FILE, "\n")
df <- read.csv(STD_FILE, stringsAsFactors = FALSE)

# Select only numeric columns (drop carbonate_class group column)
num_df  <- df %>% select(where(is.numeric))
num_vars <- names(num_df)
cat("Numeric variables:", length(num_vars), "\n")
cat(paste(num_vars, collapse = ", "), "\n\n")

# ---------- Helper: Readable Variable Labels ----------
# Shorten common soil variable names for cleaner plot axis labels
LABEL_MAP <- c(
  bulk_density_kg_l               = "BD (kg/l)",
  total_carbon_percent            = "Total C",
  total_nitrogen_percent          = "Total N",
  c_n_ratio                       = "C:N",
  organic_matter_percent          = "OM (%)",
  soil_organic_carbon_soc_percent = "SOC (%)",
  organic_carbon_stock_t_ha       = "C Stock",
  organic_matter_loi_percent_w_w  = "OM LOI",
  sand_2_00_0_063mm_percent_w_w   = "Sand",
  silt_0_063_0_002mm_percent_w_w  = "Silt",
  clay_0_002mm_percent_w_w        = "Clay",
  loi_percent_som                 = "LOI-SOM",
  water_content                   = "H2O",
  bulk_density_g_cm3              = "BD (g/cm3)",
  om_percent                      = "OM-alt",
  ca_co3                          = "CaCO3",
  col_22                          = "Col.22",
  what_is_the_soil_chroma         = "Chroma",
  what_is_the_soil_value          = "Value",
  xy_y_x                          = "xyY.x",
  xy_y_y                          = "xyY.y",
  xy_y_y_2                        = "xyY.Y",
  rgb_r                           = "RGB.R",
  rgb_g                           = "RGB.G",
  rgb_b                           = "RGB.B"
)

prettify <- function(nms) {
  sapply(nms, function(nm) {
    if (nm %in% names(LABEL_MAP)) LABEL_MAP[[nm]]
    else gsub("_", " ", nm)
  }, USE.NAMES = FALSE)
}

# ---------- Pearson Correlation Helper ----------
compute_corr <- function(mat) {
  res <- Hmisc::rcorr(as.matrix(mat), type = "pearson")
  r   <- res$r
  p   <- res$P
  # FDR correction
  p_flat <- as.vector(p)
  non_na <- !is.na(p_flat)
  p_adj_flat <- rep(NA_real_, length(p_flat))
  p_adj_flat[non_na] <- p.adjust(p_flat[non_na], method = "BH")
  p_adj <- matrix(p_adj_flat, nrow = nrow(p), dimnames = dimnames(p))
  list(r = r, p_raw = p, p_adj = p_adj, n = res$n)
}

# ---------- Plot Helper ----------
rdbu <- c("#2166AC", "white", "#D6604D")

make_corrplot <- function(r, p_adj, title, subtitle, out_file,
                          width = 10, height = 9.5, dpi = 300) {
  p <- ggcorrplot(
    r,
    p.mat         = p_adj,
    hc.order      = TRUE,
    type          = "lower",
    method        = "square",
    lab           = TRUE,
    lab_size      = 3.6,
    digits        = 2,
    colors        = rdbu,
    outline.color = "white",
    insig         = "blank",
    sig.level     = ALPHA,
    tl.cex        = 9,
    tl.col        = "black",
    ggtheme       = theme_bw(base_size = 11)
  ) +
    labs(title = title, subtitle = subtitle, fill = "r") +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, size = 8.5, color = "grey45"),
      legend.title  = element_text(face = "bold", size = 10),
      axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold",
                                    color = "black", size = 9),
      axis.text.y   = element_text(face = "bold", color = "black", size = 9),
      panel.grid    = element_blank()
    )
  ggsave(out_file, p, width = width, height = height, dpi = dpi)
  cat("Saved:", out_file, "\n")
  p
}

# ===== A. FULL SOIL CORRELATION MATRIX =====
cat("Computing full Pearson correlation matrix...\n")
corr_full <- compute_corr(num_df)

# Apply readable labels
colnames(corr_full$r)     <- prettify(colnames(corr_full$r))
rownames(corr_full$r)     <- prettify(rownames(corr_full$r))
colnames(corr_full$p_adj) <- prettify(colnames(corr_full$p_adj))
rownames(corr_full$p_adj) <- prettify(rownames(corr_full$p_adj))

n_sig_full <- sum(corr_full$p_adj[lower.tri(corr_full$p_adj)] < ALPHA, na.rm = TRUE)
n_tot_full <- sum(lower.tri(corr_full$p_adj))
cat(sprintf("Full matrix: %d variables | %d significant pairs / %d\n\n",
            ncol(num_df), n_sig_full, n_tot_full))

write.csv(round(corr_full$r,     4), file.path(OUT_DIR, "soil_corr_full_r.csv"))
write.csv(round(corr_full$p_adj, 4), file.path(OUT_DIR, "soil_corr_full_pval_fdr.csv"))

n_full <- ncol(num_df)
fig_dim <- max(9, n_full * 0.5)

make_corrplot(
  corr_full$r, corr_full$p_adj,
  title    = "Soil Data — Pearson Correlation Matrix (All Variables)",
  subtitle = paste0("n = 97 (pairwise complete)  |  z-score standardized  |",
                    "  FDR-corrected (BH)  |  Blanked: p ≥ ", ALPHA),
  out_file = file.path(OUT_DIR, "soil_correlation_matrix_full.png"),
  width    = fig_dim, height = fig_dim - 0.5
)

# ===== B. COLOUR SUB-MATRIX =====
cat("\nComputing color sub-matrix...\n")
color_cols <- num_vars[sapply(num_vars, function(nm)
  any(sapply(COLOR_PATTERNS, function(p) grepl(p, nm, ignore.case = TRUE))))]

cat("Color variables:", paste(color_cols, collapse = ", "), "\n")

if (length(color_cols) >= 3) {
  color_df <- num_df[, color_cols, drop = FALSE]
  corr_col <- compute_corr(color_df)

  colnames(corr_col$r)     <- prettify(colnames(corr_col$r))
  rownames(corr_col$r)     <- prettify(rownames(corr_col$r))
  colnames(corr_col$p_adj) <- prettify(colnames(corr_col$p_adj))
  rownames(corr_col$p_adj) <- prettify(rownames(corr_col$p_adj))

  write.csv(round(corr_col$r,     4), file.path(OUT_DIR, "soil_corr_color_r.csv"))
  write.csv(round(corr_col$p_adj, 4), file.path(OUT_DIR, "soil_corr_color_pval_fdr.csv"))

  make_corrplot(
    corr_col$r, corr_col$p_adj,
    title    = "Soil Color Variables — Pearson Correlation Matrix",
    subtitle = paste0("xyY  •  RGB  •  Munsell  |  n = 97 (pairwise complete)",
                      "  |  FDR-corrected (BH)  |  Blanked: p ≥ ", ALPHA),
    out_file = file.path(OUT_DIR, "soil_correlation_matrix_color.png"),
    width    = 7, height = 6.5
  )
} else {
  cat("Not enough color variables for sub-matrix (found:", length(color_cols), ")\n")
}

# ---------- Final Report ----------
cat("\n========================================\n")
cat("STAGE 6 COMPLETE — Soil Correlation\n")
cat("========================================\n")
cat("Significant pairs (full matrix, FDR):", n_sig_full, "/", n_tot_full, "\n")
cat("Outputs in:", OUT_DIR, "\n")
cat("\nNext: source('scripts/pipeline/07_combined_correlation.R')\n")
