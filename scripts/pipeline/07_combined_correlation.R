# ============================================================
# Tool: Combined Correlation Matrix (Soil + Soil Datasets)
# Project: Soil Pathology Color Analysis — Michigan State University
# Stage: 7
# Inputs: data/processed/normalized_data.csv
#         output/05_soil_homogeneity/soil_standardized.csv
# Output: output/07_combined_correlation/
# ============================================================

# ---------- Package Bootstrap ----------
# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))
pkgs <- c("Hmisc", "ggcorrplot", "ggplot2", "dplyr", "RColorBrewer", "GGally")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(Hmisc); library(ggcorrplot); library(ggplot2)
  library(dplyr); library(RColorBrewer); library(GGally)
})

# ---------- Config ----------
PLANT_FILE <- "data/processed/normalized_data.csv"
SOIL_FILE  <- "output/05_soil_homogeneity/soil_standardized.csv"
OUT_DIR    <- "output/07_combined_correlation"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05

# ---------- Load Datasets ----------
if (!file.exists(SOIL_FILE))
  stop("Stage 5 output not found. Run scripts/pipeline/05_soil_homogeneity.R first.")

cat("Loading soil color data:", PLANT_FILE, "\n")
soil_color <- read.csv(PLANT_FILE, stringsAsFactors = FALSE)
cat("  Rows:", nrow(soil_color), "| Cols:", ncol(soil_color), "\n")

cat("Loading soil standardized data:", SOIL_FILE, "\n")
soil_db  <- read.csv(SOIL_FILE, stringsAsFactors = FALSE)
cat("  Rows:", nrow(soil_db), "| Cols:", ncol(soil_db), "\n\n")

# ---------- Extract Common Color Variables ----------
# Both datasets share the RGB color space.
# Soil Color: R_norm (normalized R), G (original, normal), B (original, normal)
#             L (LAB Lightness)
# Soil Database:  rgb_r, rgb_g, rgb_b (z-score standardized)
#                 xyY.Y luminance (closest to Lightness)

# --- From Soil Color Data ---
cat("Extracting from soil color data...\n")
if (!all(c("R_norm", "G", "B", "L") %in% names(soil_color)))
  stop("Expected columns R_norm, G, B, L not found in soil color data.")

# Z-score the soil color RGB variables (G, B, L) for comparability with soil database z-scores
# R is already OrderNorm-transformed; z-score all for consistent scale
soil_color_rgb <- data.frame(
  R         = scale(soil_color$R_norm)[, 1],
  G         = scale(soil_color$G)[, 1],
  B         = scale(soil_color$B)[, 1],
  Lightness = scale(soil_color$L)[, 1],
  Source    = "Soil Color",
  stringsAsFactors = FALSE
)
cat("  Soil Color: R, G, B, Lightness (L-LAB) | n =", nrow(soil_color_rgb), "\n")

# --- From Soil Database Data ---
cat("Extracting from soil database data...\n")

# Locate RGB columns (janitor clean_names produces rgb_r / rgb_g / rgb_b)
find_col <- function(df, patterns) {
  for (pat in patterns) {
    hits <- names(df)[grepl(pat, names(df), ignore.case = TRUE)]
    if (length(hits) > 0) return(hits[1])
  }
  return(NA_character_)
}

r_col <- find_col(soil_db, c("^rgb_r$", "rgb\\.r", "rgb_r"))
g_col <- find_col(soil_db, c("^rgb_g$", "rgb\\.g", "rgb_g"))
b_col <- find_col(soil_db, c("^rgb_b$", "rgb\\.b", "rgb_b"))

# Lightness proxy: xyY.Y (luminance) or Munsell Value
y_col <- find_col(soil_db, c("xy_y_y_2", "xyy_y_2", "xy_y_y_3",
                           "xy_y_y$", "xyy_y$", "y_3"))
if (is.na(y_col)) y_col <- find_col(soil_db, c("what_is_the_soil_value",
                                              "munsell_value", "value"))

cat("  Soil Database RGB columns: R=", r_col, " G=", g_col, " B=", b_col, "\n")
cat("  Soil Database Lightness:  ", y_col, "\n")

missing_cols <- c(r_col, g_col, b_col)[is.na(c(r_col, g_col, b_col))]
if (length(missing_cols) > 0)
  stop("Could not find soil RGB columns. Check column names in soil_standardized.csv.")

soil_db_rgb <- data.frame(
  R         = if (!is.na(r_col)) as.numeric(soil_db[[r_col]]) else NA_real_,
  G         = if (!is.na(g_col)) as.numeric(soil_db[[g_col]]) else NA_real_,
  B         = if (!is.na(b_col)) as.numeric(soil_db[[b_col]]) else NA_real_,
  Lightness = if (!is.na(y_col)) as.numeric(soil_db[[y_col]]) else NA_real_,
  Source    = "Soil Database",
  stringsAsFactors = FALSE
)

# Retain only rows with at least R/G/B complete
complete_rgb <- complete.cases(soil_db_rgb[, c("R", "G", "B")])
soil_db_rgb     <- soil_db_rgb[complete_rgb, ]
cat("  Soil Database: R, G, B",
    if (!is.na(y_col)) ", Lightness (xyY.Y/Munsell Value)" else " (no Lightness column found)",
    "| n =", nrow(soil_db_rgb), "\n\n")

# ---------- Combine ----------
combined <- rbind(soil_color_rgb, soil_db_rgb)
cat("Combined dataset: n =", nrow(combined),
    "(Soil Color:", sum(combined$Source == "Soil Color"),
    "| Soil Database:", sum(combined$Source == "Soil Database"), ")\n\n")

# Drop Lightness if entirely NA
if (all(is.na(combined$Lightness))) {
  combined$Lightness <- NULL
  cat("NOTE: Lightness dropped (all-NA in combined dataset)\n\n")
}

write.csv(combined, file.path(OUT_DIR, "combined_dataset.csv"), row.names = FALSE)

# ---------- Pearson Correlation on Shared Variables ----------
shared_vars <- names(combined)[names(combined) != "Source"]
cor_mat     <- as.matrix(combined[, shared_vars])

cat("Computing Pearson correlation on:", paste(shared_vars, collapse = ", "), "\n")

cor_res <- Hmisc::rcorr(cor_mat, type = "pearson")
r_mat   <- cor_res$r
p_mat   <- cor_res$P

# FDR correction
p_flat     <- as.vector(p_mat)
non_na     <- !is.na(p_flat)
p_adj_flat <- rep(NA_real_, length(p_flat))
p_adj_flat[non_na] <- p.adjust(p_flat[non_na], method = "BH")
p_adj <- matrix(p_adj_flat, nrow = nrow(p_mat), dimnames = dimnames(p_mat))

cat("\nPearson correlation matrix (r):\n")
print(round(r_mat, 3))
cat("\nFDR-adjusted p-values (BH):\n")
print(round(p_adj, 4))

n_sig <- sum(p_adj[lower.tri(p_adj)] < ALPHA, na.rm = TRUE)
n_tot <- sum(lower.tri(p_adj))
cat(sprintf("\nSignificant pairs (FDR p < %.2f): %d / %d\n\n", ALPHA, n_sig, n_tot))

write.csv(round(r_mat, 4), file.path(OUT_DIR, "combined_corr_r.csv"))
write.csv(round(p_mat, 4), file.path(OUT_DIR, "combined_corr_pval_raw.csv"))
write.csv(round(p_adj, 4), file.path(OUT_DIR, "combined_corr_pval_fdr.csv"))

# ---------- Publication-Ready Correlation Matrix ----------
rdbu <- c("#2166AC", "white", "#D6604D")

p_corr <- ggcorrplot(
  r_mat,
  p.mat         = p_adj,
  hc.order      = TRUE,
  type          = "lower",
  method        = "square",
  lab           = TRUE,
  lab_size      = 5,
  digits        = 2,
  colors        = rdbu,
  outline.color = "white",
  insig         = "blank",
  sig.level     = ALPHA,
  tl.cex        = 13,
  tl.col        = "black",
  ggtheme       = theme_bw(base_size = 13)
) +
  labs(
    title    = "Combined Soil Color & Soil Database — Pearson Correlation Matrix",
    subtitle = paste0(
      "n = ", nrow(combined),
      " (Soil Color: ", sum(combined$Source == "Soil Color"),
      " | Soil Database: ", sum(combined$Source == "Soil Database"), ")",
      "  |  z-score standardized  |  FDR-corrected (BH)",
      "  |  Blanked: p ≥ ", ALPHA
    ),
    fill = "r"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey45"),
    legend.title  = element_text(face = "bold", size = 10),
    axis.text     = element_text(face = "bold", color = "black", size = 13),
    panel.grid    = element_blank()
  )

out_corr <- file.path(OUT_DIR, "combined_correlation_matrix.png")
ggsave(out_corr, p_corr, width = 7, height = 6.5, dpi = 300)
cat("Correlation matrix saved:", out_corr, "\n")

# ---------- Scatter Matrix by Source ----------
cat("Generating scatter matrix by Source...\n")

scatter_df        <- combined
scatter_df$Source <- factor(scatter_df$Source, levels = c("Soil Color", "Soil Database"))

source_colors <- c(`Soil Color` = "#2980B9", `Soil Database` = "#E67E22")

p_scatter <- GGally::ggpairs(
  scatter_df,
  columns  = seq_along(shared_vars),
  mapping  = aes(color = Source, alpha = 0.5),
  upper    = list(continuous = GGally::wrap("cor", size = 4, method = "pearson")),
  lower    = list(continuous = GGally::wrap("points", size = 1.2)),
  diag     = list(continuous = GGally::wrap("densityDiag", alpha = 0.5)),
  title    = "Soil Color vs. Soil Database Color — Scatter Matrix"
) +
  scale_color_manual(values = source_colors) +
  scale_fill_manual(values  = source_colors) +
  theme_bw(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 13),
    strip.text  = element_text(face = "bold")
  )

out_scatter <- file.path(OUT_DIR, "combined_scatter_matrix.png")
ggsave(out_scatter, p_scatter,
       width  = max(7, length(shared_vars) * 2.5),
       height = max(6, length(shared_vars) * 2.5),
       dpi    = 300)
cat("Scatter matrix saved:", out_scatter, "\n")

# ---------- Final Report ----------
cat("\n========================================\n")
cat("STAGE 7 COMPLETE — Combined Correlation\n")
cat("========================================\n")
cat("Variables:", paste(shared_vars, collapse = ", "), "\n")
cat("n combined:", nrow(combined), "\n")
cat("Significant pairs:", n_sig, "/", n_tot, "\n")
cat("Outputs in:", OUT_DIR, "\n")
cat("  combined_correlation_matrix.png  - publication-ready heatmap\n")
cat("  combined_scatter_matrix.png      - scatter matrix by Source\n")
cat("  combined_corr_r.csv              - correlation coefficients\n")
cat("  combined_corr_pval_fdr.csv       - FDR-adjusted p-values\n")
