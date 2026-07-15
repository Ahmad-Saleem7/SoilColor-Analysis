# ============================================================
# Tool: Data Normalization
# Project: Soil Pathology Color Analysis — Michigan State University
# Stage: 2 of 5
# Prerequisite: scripts/pipeline/01_soil_distribution.R must have run first
# Output: output/02_normalization/
# ============================================================

# ---------- Package Bootstrap ----------
# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))
pkgs <- c("readxl", "bestNormalize", "nortest", "ggplot2",
          "patchwork", "dplyr", "moments", "MASS")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(bestNormalize)
  library(nortest)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(moments)
  library(MASS)
})

# ---------- Config ----------
DATA_FILE   <- "data/raw/soil_color_data.xlsx"
PHASE1_CSV  <- "output/01_distributions/distribution_summary.csv"
OUT_DIR     <- "output/02_normalization"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MEAN_VARS <- c("R", "G", "B", "H", "S", "Br", "L", "a", "b")
ALPHA     <- 0.05

# ---------- Load Inputs ----------
if (!file.exists(PHASE1_CSV))
  stop("Stage 1 output not found. Run scripts/pipeline/01_soil_distribution.R first.")

df     <- tryCatch(
  read_excel(DATA_FILE, sheet = 1),
  error = function(e) stop("Cannot read Excel file: ", conditionMessage(e))
)

phase1     <- read.csv(PHASE1_CSV, stringsAsFactors = FALSE)
non_normal <- phase1$Variable[phase1$Normality == "NON-NORMAL"]

cat("Non-normal variables identified in Stage 1:",
    if (length(non_normal) == 0) "none" else paste(non_normal, collapse = ", "), "\n\n")

if (length(non_normal) == 0) {
  cat("All variables are normally distributed. No normalization required.\n")
  cat("Proceed to Stage 3 (composition analysis).\n")
  quit(save = "no")
}

# ---------- Normality Check Helper ----------
check_norm <- function(x) {
  x <- as.numeric(na.omit(x))
  if (length(x) < 3) return(list(sw_p = NA, ad_p = NA, is_normal = FALSE, skewness = NA))
  sw <- shapiro.test(x)
  ad <- nortest::ad.test(x)
  list(
    sw_p      = sw$p.value,
    ad_p      = ad$p.value,
    is_normal = (sw$p.value > ALPHA) & (ad$p.value > ALPHA),
    skewness  = moments::skewness(x)
  )
}

# ---------- Transformation Suite ----------
# Each function returns the transformed vector.
# Offset ensures strict positivity for log/sqrt.

apply_log1p_shift <- function(x) {
  shift <- if (min(x) <= 0) abs(min(x)) + 1e-6 else 0
  log1p(x + shift)
}

apply_sqrt_shift <- function(x) {
  shift <- if (min(x) < 0) abs(min(x)) + 1e-6 else 0
  sqrt(x + shift)
}

apply_boxcox <- function(x) {
  shift <- if (min(x) <= 0) abs(min(x)) + 1e-6 else 0
  x_s   <- x + shift
  bc    <- MASS::boxcox(x_s ~ 1, plotit = FALSE)
  lam   <- bc$x[which.max(bc$y)]
  if (abs(lam) < 1e-4) log(x_s) else (x_s^lam - 1) / lam
}

apply_yeo_johnson <- function(x) {
  obj <- bestNormalize::yeojohnson(x)
  predict(obj)
}

apply_order_norm <- function(x) {
  obj <- bestNormalize::orderNorm(x)
  predict(obj)
}

TRANSFORMS <- list(
  log1p_shift = apply_log1p_shift,
  sqrt_shift  = apply_sqrt_shift,
  boxcox      = apply_boxcox,
  yeo_johnson = apply_yeo_johnson,
  order_norm  = apply_order_norm
)

# ---------- Process Each Non-Normal Variable ----------
transform_records <- list()
normalized_cols   <- list()

for (v in non_normal) {
  x_orig <- as.numeric(df[[v]])
  orig   <- check_norm(x_orig)

  cat("---", v, "---\n")
  cat(sprintf("  Original: SW p=%.4f | AD p=%.4f | Skew=%.3f\n",
              orig$sw_p, orig$ad_p, orig$skewness))

  best_combined <- -1
  best_name     <- "none"
  best_x        <- x_orig
  best_sw       <- orig$sw_p
  best_ad       <- orig$ad_p
  best_is_norm  <- FALSE

  for (t_name in names(TRANSFORMS)) {
    x_t    <- tryCatch(TRANSFORMS[[t_name]](x_orig), error = function(e) NULL)
    if (is.null(x_t) || any(!is.finite(x_t))) {
      cat(sprintf("  %-14s -> ERROR or non-finite values\n", t_name))
      next
    }
    res    <- check_norm(x_t)
    combined <- min(res$sw_p, res$ad_p, na.rm = TRUE)

    cat(sprintf("  %-14s -> SW p=%.4f | AD p=%.4f%s\n",
                t_name, res$sw_p, res$ad_p,
                if (isTRUE(res$is_normal)) "  [NORMAL]" else ""))

    if (combined > best_combined) {
      best_combined <- combined
      best_name     <- t_name
      best_x        <- x_t
      best_sw       <- res$sw_p
      best_ad       <- res$ad_p
      best_is_norm  <- isTRUE(res$is_normal)
    }
  }

  final_status <- if (best_is_norm) "NORMAL" else "REQUIRES-NONPARAMETRIC"
  cat(sprintf("  => Best: %s | Status: %s\n\n", best_name, final_status))

  transform_records[[v]] <- data.frame(
    Variable             = v,
    Best_Transform       = best_name,
    SW_p_original        = round(orig$sw_p, 4),
    AD_p_original        = round(orig$ad_p, 4),
    Skewness_original    = round(orig$skewness, 3),
    SW_p_transformed     = round(best_sw, 4),
    AD_p_transformed     = round(best_ad, 4),
    Final_Status         = final_status,
    stringsAsFactors     = FALSE
  )

  normalized_cols[[v]] <- best_x
}

# ---------- Save Normalized Dataset ----------
norm_df <- df
for (v in names(normalized_cols)) {
  norm_df[[paste0(v, "_norm")]] <- normalized_cols[[v]]
}

out_data <- "data/processed/normalized_data.csv"
write.csv(norm_df, out_data, row.names = FALSE)
cat("Normalized dataset saved:", out_data, "\n")

# ---------- Save Transformation Summary ----------
transform_summary <- bind_rows(transform_records)
out_summary <- file.path(OUT_DIR, "transformation_summary.csv")
write.csv(transform_summary, out_summary, row.names = FALSE)
cat("Transformation summary saved:", out_summary, "\n\n")

cat("=== Transformation Results ===\n")
print(transform_summary[, c("Variable", "Best_Transform",
                             "SW_p_original", "SW_p_transformed",
                             "AD_p_transformed", "Final_Status")])

# ---------- Before / After Visualization ----------
cat("\nGenerating before/after plots...\n")

FILL_BEFORE <- "#E74C3C"   # red = non-normal
FILL_AFTER_NORM   <- "#27AE60"   # green = achieved normality
FILL_AFTER_NONPAR <- "#F39C12"   # orange = still non-parametric

for (v in non_normal) {
  x_orig   <- as.numeric(df[[v]])
  x_trans  <- normalized_cols[[v]]
  rec      <- transform_records[[v]]
  fill_a   <- if (rec$Final_Status == "NORMAL") FILL_AFTER_NORM else FILL_AFTER_NONPAR

  df_o <- data.frame(x = x_orig)
  df_t <- data.frame(x = x_trans)

  # Original histogram
  p_hist_o <- ggplot(df_o, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 15,
                   fill = FILL_BEFORE, color = "white", alpha = 0.75) +
    geom_density(linewidth = 0.8) +
    labs(title    = paste(v, "— Original"),
         subtitle = paste0("SW p=", rec$SW_p_original,
                           "  |  AD p=", rec$AD_p_original),
         x = v, y = "Density") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", color = FILL_BEFORE))

  # Original Q-Q
  p_qq_o <- ggplot(df_o, aes(sample = x)) +
    stat_qq(color = FILL_BEFORE, size = 1.8, alpha = 0.75) +
    stat_qq_line(color = "#2C3E50", linewidth = 0.9) +
    labs(title = "Q-Q (Original)") +
    theme_classic(base_size = 10)

  # Transformed histogram
  p_hist_t <- ggplot(df_t, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 15,
                   fill = fill_a, color = "white", alpha = 0.75) +
    geom_density(linewidth = 0.8) +
    labs(title    = paste0(v, " — ", rec$Best_Transform),
         subtitle = paste0("SW p=", rec$SW_p_transformed,
                           "  |  AD p=", rec$AD_p_transformed,
                           "  |  ", rec$Final_Status),
         x = paste0(v, " (transformed)"), y = "Density") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold", color = fill_a))

  # Transformed Q-Q
  p_qq_t <- ggplot(df_t, aes(sample = x)) +
    stat_qq(color = fill_a, size = 1.8, alpha = 0.75) +
    stat_qq_line(color = "#2C3E50", linewidth = 0.9) +
    labs(title = paste0("Q-Q (", rec$Best_Transform, ")")) +
    theme_classic(base_size = 10)

  combined_plot <- (p_hist_o | p_qq_o) / (p_hist_t | p_qq_t) +
    plot_annotation(
      title    = paste("Normalization:", v),
      subtitle = paste("Best transformation:", rec$Best_Transform,
                       "| Final status:", rec$Final_Status),
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 9, color = "grey40"))
    )

  out_png <- file.path(OUT_DIR, paste0(v, "_normalization.png"))
  ggsave(out_png, plot = combined_plot, width = 10, height = 8, dpi = 150)
  cat("  Saved:", out_png, "\n")
}

# ---------- Final Report ----------
achieved_norm <- transform_summary$Variable[transform_summary$Final_Status == "NORMAL"]
needs_nonpar  <- transform_summary$Variable[transform_summary$Final_Status == "REQUIRES-NONPARAMETRIC"]

cat("\n========================================\n")
cat("STAGE 2 COMPLETE — Normalization\n")
cat("========================================\n")
cat("Normalized successfully :", paste(achieved_norm, collapse = ", "),
    if (length(achieved_norm) == 0) "none" else "", "\n")
cat("Requires non-parametric :", paste(needs_nonpar, collapse = ", "),
    if (length(needs_nonpar) == 0) "none" else "", "\n")
cat("\nNote: *_norm columns in normalized_data.csv hold transformed values.\n")
cat("Non-parametric variables use Spearman correlation (Stage 4)\n")
cat("and Kruskal-Wallis + Dunn's test (Stage 5).\n")
cat("Outputs in:", OUT_DIR, "\n")
