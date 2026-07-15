# ============================================================
# Tool: Soil Data — Distribution Analysis + Normalization
# Project: Soil Pathology Color Analysis — Michigan State University
# Stage: 4
# Input:  data/raw/soil_color_database.xlsx (sheet "all data")
# Output: output/04_soil_dist_norm/
# ============================================================

# ---------- Package Bootstrap ----------
# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))
pkgs <- c("readxl", "janitor", "fitdistrplus", "nortest", "bestNormalize",
          "MASS", "moments", "ggplot2", "patchwork", "dplyr")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(readxl); library(janitor); library(fitdistrplus); library(nortest)
  library(bestNormalize); library(MASS); library(moments)
  library(ggplot2); library(patchwork); library(dplyr)
})

# ---------- Config ----------
DATA_FILE <- "data/raw/soil_color_database.xlsx"
SHEET     <- "all data"
OUT_DIR   <- "output/04_soil_dist_norm"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

MIN_VALID  <- 15     # minimum non-missing observations required for analysis
MISS_THRESH <- 0.50  # drop columns with > 50% missing
ALPHA       <- 0.05

DIST_CANDS  <- c("norm", "lnorm", "gamma", "weibull")

# ---------- Load & Clean ----------
# File structure: Row 1 = institution labels, Row 2 = column names,
#                Row 3 = alternative NRM names, Rows 4+ = actual data
cat("Loading:", DATA_FILE, "| Sheet:", SHEET, "\n")

# Step 1: Extract column names from row 2
hdr_raw  <- suppressMessages(
  read_excel(DATA_FILE, sheet = SHEET, skip = 1, n_max = 1, col_names = FALSE)
)
headers <- as.character(unlist(hdr_raw[1, ]))
na_pos  <- is.na(headers) | trimws(headers) == ""
headers[na_pos] <- paste0("col_", which(na_pos))

# Step 2: Read actual data from row 4 onwards (skip 3 header rows)
raw_data <- suppressMessages(
  read_excel(DATA_FILE, sheet = SHEET, skip = 3, col_names = FALSE)
)

# Align column count
if (ncol(raw_data) < length(headers)) {
  headers <- headers[seq_len(ncol(raw_data))]
} else if (ncol(raw_data) > length(headers)) {
  headers <- c(headers, paste0("extra_", seq_len(ncol(raw_data) - length(headers))))
}
names(raw_data) <- headers

df  <- janitor::clean_names(raw_data)

cat("Raw dimensions:", nrow(df), "rows x", ncol(df), "cols\n")
cat("Column names after clean_names():\n")
print(names(df))
cat("\n")

# ---------- Select Numeric Columns ----------
# ID/name columns to always exclude
EXCLUDE_PATTERNS <- c("^id$", "^sample_name$", "^sample$",
                       "^x$", "^y$",   # spatial coords — not soil properties
                       "carbonate_class", "textural_classification",
                       "what_is_the_colour", "select_munsell_chart",
                       "soil_texture", "vegetation_cover",
                       "field_practices", "field_assessment")

is_excl <- function(nm) {
  any(sapply(EXCLUDE_PATTERNS, function(pat) grepl(pat, nm, ignore.case = TRUE)))
}

num_df <- df %>%
  select(where(is.numeric)) %>%
  select(-any_of(Filter(is_excl, names(df))))

# Drop columns with > MISS_THRESH missing
miss_frac  <- colMeans(is.na(num_df))
num_df     <- num_df[, miss_frac <= MISS_THRESH, drop = FALSE]
num_vars   <- names(num_df)

cat("Numeric variables selected (", length(num_vars), "):\n",
    paste(num_vars, collapse = ", "), "\n\n")

# ---------- Per-Variable Analysis Function ----------
analyze_var <- function(var_nm, x_raw) {
  x <- as.numeric(x_raw)
  x <- x[is.finite(x)]

  if (length(x) < MIN_VALID) {
    cat("  [SKIP]", var_nm, "— only", length(x), "valid obs (min:", MIN_VALID, ")\n")
    return(NULL)
  }

  sw  <- tryCatch(shapiro.test(x), error = function(e) list(statistic = NA, p.value = NA))
  ad  <- tryCatch(nortest::ad.test(x), error = function(e) list(statistic = NA, p.value = NA))
  sk  <- moments::skewness(x)
  kt  <- moments::kurtosis(x) - 3

  # Distribution fitting
  x_pos     <- x[x > 0]
  aic_vals  <- setNames(rep(NA_real_, length(DIST_CANDS)), DIST_CANDS)
  fit_list  <- list()

  for (dist in DIST_CANDS) {
    fd <- if (dist == "norm") x else x_pos
    if (length(fd) < 4) next
    tryCatch({
      f <- fitdist(fd, dist, method = "mle")
      fit_list[[dist]] <- f
      aic_vals[[dist]] <- f$aic
    }, error = function(e) NULL)
  }

  best_dist <- if (any(!is.na(aic_vals))) names(which.min(aic_vals)) else "norm"
  is_normal <- !is.na(sw$p.value) && !is.na(ad$p.value) &&
               sw$p.value > ALPHA && ad$p.value > ALPHA

  list(
    variable  = var_nm,
    n         = length(x),
    mean      = mean(x),
    sd        = sd(x),
    min       = min(x), max = max(x),
    skewness  = sk,   kurtosis = kt,
    sw_stat   = unname(sw$statistic),
    sw_p      = sw$p.value,
    ad_stat   = unname(ad$statistic),
    ad_p      = ad$p.value,
    normality = ifelse(is_normal, "NORMAL", "NON-NORMAL"),
    aic_norm    = aic_vals["norm"],
    aic_lnorm   = aic_vals["lnorm"],
    aic_gamma   = aic_vals["gamma"],
    aic_weibull = aic_vals["weibull"],
    best_dist   = best_dist,
    data        = x,
    fits        = fit_list
  )
}

# ---------- Run Analysis ----------
cat("Running distribution analysis on", length(num_vars), "variables...\n")
results     <- lapply(num_vars, function(v) analyze_var(v, num_df[[v]]))
names(results) <- num_vars
results     <- Filter(Negate(is.null), results)
valid_vars  <- names(results)

# ---------- Summary Table ----------
summary_df <- bind_rows(lapply(results, function(r) data.frame(
  Variable    = r$variable,
  N_valid     = r$n,
  Mean        = round(r$mean, 4),
  SD          = round(r$sd,   4),
  Skewness    = round(r$skewness, 3),
  ExKurtosis  = round(r$kurtosis, 3),
  SW_W        = round(r$sw_stat, 4),
  SW_p        = round(r$sw_p,    4),
  AD_stat     = round(r$ad_stat, 4),
  AD_p        = round(r$ad_p,    4),
  Normality   = r$normality,
  BestFit     = r$best_dist,
  stringsAsFactors = FALSE
)))

out_summ <- file.path(OUT_DIR, "soil_distribution_summary.csv")
write.csv(summary_df, out_summ, row.names = FALSE)
cat("\nDistribution Summary:\n")
print(summary_df[, c("Variable","N_valid","Skewness","SW_p","AD_p","Normality","BestFit")])
cat("\nSaved:", out_summ, "\n\n")

# ---------- Normalization of Non-Normal Variables ----------
TRANSFORMS <- list(
  log1p_shift = function(x) { s <- max(0, -min(x)) + 1e-6; log1p(x + s) },
  sqrt_shift  = function(x) { s <- max(0, -min(x)) + 1e-6; sqrt(x + s) },
  boxcox      = function(x) {
    s  <- max(0, -min(x)) + 1e-6
    xs <- x + s
    bc <- MASS::boxcox(xs ~ 1, plotit = FALSE)
    lam <- bc$x[which.max(bc$y)]
    if (abs(lam) < 1e-4) log(xs) else (xs^lam - 1) / lam
  },
  yeo_johnson = function(x) predict(bestNormalize::yeojohnson(x)),
  order_norm  = function(x) predict(bestNormalize::orderNorm(x))
)

check_norm <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(list(sw_p = NA, ad_p = NA, is_normal = FALSE))
  sw <- tryCatch(shapiro.test(x), error = function(e) list(p.value = 0))
  ad <- tryCatch(nortest::ad.test(x), error = function(e) list(p.value = 0))
  list(sw_p = sw$p.value, ad_p = ad$p.value,
       is_normal = sw$p.value > ALPHA & ad$p.value > ALPHA)
}

non_normal   <- summary_df$Variable[summary_df$Normality == "NON-NORMAL"]
norm_records <- list()
norm_cols    <- list()

cat("Normalizing", length(non_normal), "non-normal variables...\n\n")

for (v in non_normal) {
  x_orig <- as.numeric(num_df[[v]])
  x_orig <- ifelse(is.finite(x_orig), x_orig, NA_real_)
  cat("---", v, "---\n")

  best_comb <- -1; best_nm <- "none"; best_x <- x_orig
  best_sw   <- 0;  best_ad <- 0;      best_norm <- FALSE

  for (t_nm in names(TRANSFORMS)) {
    x_t <- tryCatch(TRANSFORMS[[t_nm]](x_orig[!is.na(x_orig)]), error = function(e) NULL)
    if (is.null(x_t) || !all(is.finite(x_t))) {
      cat(sprintf("  %-14s -> ERROR/non-finite\n", t_nm)); next
    }
    # Re-expand to full length (NAs where original was NA)
    x_full <- rep(NA_real_, length(x_orig))
    x_full[!is.na(x_orig)] <- x_t
    res <- check_norm(x_t)
    cmb <- min(res$sw_p, res$ad_p, na.rm = TRUE)
    cat(sprintf("  %-14s -> SW p=%.4f | AD p=%.4f%s\n",
                t_nm, res$sw_p, res$ad_p,
                if (isTRUE(res$is_normal)) "  [NORMAL]" else ""))
    if (cmb > best_comb) {
      best_comb <- cmb; best_nm <- t_nm; best_sw <- res$sw_p
      best_ad   <- res$ad_p; best_norm <- isTRUE(res$is_normal)
      best_x    <- x_full
    }
  }

  status <- if (best_norm) "NORMAL" else "REQUIRES-NONPARAMETRIC"
  cat(sprintf("  => Best: %s | Status: %s\n\n", best_nm, status))

  norm_records[[v]] <- data.frame(
    Variable         = v,
    Best_Transform   = best_nm,
    SW_p_original    = round(results[[v]]$sw_p,  4),
    SW_p_transformed = round(best_sw, 4),
    AD_p_transformed = round(best_ad, 4),
    Final_Status     = status,
    stringsAsFactors = FALSE
  )
  norm_cols[[v]] <- best_x
}

# ---------- Build & Save Normalized Dataset ----------
norm_df <- num_df
for (v in names(norm_cols)) {
  norm_df[[paste0(v, "_norm")]] <- norm_cols[[v]]
}

out_norm <- file.path(OUT_DIR, "soil_normalized.csv")
write.csv(norm_df, out_norm, row.names = FALSE)
cat("Normalized dataset saved:", out_norm, "\n")

if (length(norm_records) > 0) {
  trans_summ <- bind_rows(norm_records)
  out_trans  <- file.path(OUT_DIR, "soil_transformation_summary.csv")
  write.csv(trans_summ, out_trans, row.names = FALSE)
  cat("Transformation summary saved:", out_trans, "\n\n")
  cat("=== Transformation Results ===\n")
  print(trans_summ[, c("Variable","Best_Transform","SW_p_original",
                        "SW_p_transformed","Final_Status")])
}

# ---------- Visualization: Distribution Panels ----------
cat("\nGenerating distribution plots...\n")

STATUS_COLOR <- c(NORMAL = "#1a7a1a", "NON-NORMAL" = "#c0392b")
HIST_FILL    <- c(NORMAL = "#2ecc71", "NON-NORMAL" = "#e74c3c")

make_panel <- function(r) {
  status <- r$normality
  x_df   <- data.frame(x = r$data)
  xseq   <- seq(min(r$data), max(r$data), length.out = 400)

  p_hist <- ggplot(x_df, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 15,
                   fill = HIST_FILL[[status]], color = "white", alpha = 0.75) +
    geom_density(linewidth = 0.7) +
    labs(title    = r$variable,
         subtitle = paste0(status, "  SW p=", round(r$sw_p, 3),
                           "  Skew=", round(r$skewness, 2)),
         x = NULL, y = "Density") +
    theme_classic(base_size = 9) +
    theme(plot.title = element_text(face = "bold", size = 9,
                                    color = STATUS_COLOR[[status]]),
          plot.subtitle = element_text(size = 7, color = "grey40"))

  p_qq <- ggplot(x_df, aes(sample = x)) +
    stat_qq(color = HIST_FILL[[status]], size = 1.5, alpha = 0.7) +
    stat_qq_line(color = "#2C3E50", linewidth = 0.8) +
    labs(title = NULL, x = "Theoretical", y = "Sample") +
    theme_classic(base_size = 9)

  p_hist | p_qq
}

# Individual plots
for (v in valid_vars) {
  p   <- make_panel(results[[v]])
  png <- file.path(OUT_DIR, paste0(v, "_dist.png"))
  ggsave(png, p, width = 8, height = 3.5, dpi = 150)
}

# Composite (max 16 variables per page to keep readable)
chunks <- split(valid_vars, ceiling(seq_along(valid_vars) / 8))
for (i in seq_along(chunks)) {
  panels    <- lapply(results[chunks[[i]]], make_panel)
  composite <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title    = paste0("Soil Data — Distribution Analysis (Page ", i, ")"),
      subtitle = "Green = NORMAL | Red = NON-NORMAL",
      theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                       plot.subtitle = element_text(size = 9, color = "grey40"))
    )
  comp_png <- file.path(OUT_DIR, sprintf("00_distributions_page%02d.png", i))
  ggsave(comp_png, composite, width = 14, height = 4.5 * ceiling(length(chunks[[i]]) / 2),
         dpi = 150)
  cat("Composite page", i, "saved:", comp_png, "\n")
}

# ---------- Final Report ----------
cat("\n========================================\n")
cat("STAGE 4 COMPLETE — Soil Distribution & Normalization\n")
cat("========================================\n")
norm_v <- summary_df$Variable[summary_df$Normality == "NORMAL"]
nonnm_v <- summary_df$Variable[summary_df$Normality == "NON-NORMAL"]
cat("Normal    :", paste(norm_v,  collapse = ", "), "\n")
cat("Non-normal:", paste(nonnm_v, collapse = ", "), "\n")
cat("Outputs in:", OUT_DIR, "\n")
cat("\nNext: source('scripts/pipeline/05_soil_homogeneity.R')\n")
