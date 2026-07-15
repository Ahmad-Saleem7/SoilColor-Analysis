# ============================================================
# Tool: Distribution Analysis
# Project: Soil Pathology Color Analysis — Michigan State University
# Stage: 1 of 5
# Output: output/01_distributions/
# ============================================================

# ---------- Package Bootstrap ----------
# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))
pkgs <- c("readxl", "fitdistrplus", "nortest", "ggplot2",
          "patchwork", "dplyr", "moments")

for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cran.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(readxl)
  library(fitdistrplus)
  library(nortest)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(moments)
})

# ---------- Config ----------
DATA_FILE <- "data/raw/soil_color_data.xlsx"
OUT_DIR   <- "output/01_distributions"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 9 mean variables across RGB / HSB / LAB color spaces
MEAN_VARS <- c("R", "G", "B", "H", "S", "Br", "L", "a", "b")

# Color-space grouping (for labelling)
COLOR_SPACE <- c(R = "RGB", G = "RGB", B = "RGB",
                 H = "HSB", S = "HSB", Br = "HSB",
                 L = "LAB", a = "LAB", b = "LAB")

DIST_CANDIDATES <- c("norm", "lnorm", "gamma", "weibull")

ALPHA <- 0.05   # significance level for normality decision

# ---------- Load Data ----------
cat("Loading data from:", DATA_FILE, "\n")
df <- tryCatch(
  read_excel(DATA_FILE, sheet = 1),
  error = function(e) stop("Cannot read Excel file: ", conditionMessage(e))
)

missing_vars <- setdiff(MEAN_VARS, names(df))
if (length(missing_vars) > 0)
  stop("Expected columns not found in data: ", paste(missing_vars, collapse = ", "))

cat("Data loaded:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Analysing variables:", paste(MEAN_VARS, collapse = ", "), "\n\n")

# ---------- Per-Variable Analysis ----------
analyze_variable <- function(var_name, raw_vec) {

  x <- na.omit(as.numeric(raw_vec))

  # Normality tests
  sw_res  <- shapiro.test(x)
  ad_res  <- nortest::ad.test(x)

  # Descriptive stats
  sk <- moments::skewness(x)
  kt <- moments::kurtosis(x)    # excess kurtosis = kt - 3

  # Distribution fitting (requires strictly positive data)
  x_pos     <- x[x > 0]
  fit_list  <- list()
  aic_vals  <- setNames(rep(NA_real_, length(DIST_CANDIDATES)), DIST_CANDIDATES)

  for (dist in DIST_CANDIDATES) {
    fit_data <- if (dist == "norm") x else x_pos
    if (length(fit_data) < 4) next
    tryCatch({
      f <- fitdist(fit_data, dist, method = "mle")
      fit_list[[dist]]  <- f
      aic_vals[[dist]]  <- f$aic
    }, error = function(e) NULL)
  }

  best_dist <- if (any(!is.na(aic_vals))) names(which.min(aic_vals)) else "norm"

  # Normality decision: both SW and AD must pass
  is_normal <- (sw_res$p.value > ALPHA) & (ad_res$p.value > ALPHA)

  list(
    variable    = var_name,
    color_space = COLOR_SPACE[[var_name]],
    n           = length(x),
    mean        = mean(x),
    sd          = sd(x),
    median      = median(x),
    min         = min(x),
    max         = max(x),
    skewness    = sk,
    kurtosis    = kt - 3,   # excess kurtosis
    sw_stat     = unname(sw_res$statistic),
    sw_p        = sw_res$p.value,
    ad_stat     = unname(ad_res$statistic),
    ad_p        = ad_res$p.value,
    normality   = ifelse(is_normal, "NORMAL", "NON-NORMAL"),
    aic_norm    = aic_vals["norm"],
    aic_lnorm   = aic_vals["lnorm"],
    aic_gamma   = aic_vals["gamma"],
    aic_weibull = aic_vals["weibull"],
    best_dist   = best_dist,
    data        = x,
    fits        = fit_list
  )
}

results <- setNames(
  lapply(MEAN_VARS, function(v) analyze_variable(v, df[[v]])),
  MEAN_VARS
)

# ---------- Summary Table ----------
summary_df <- bind_rows(lapply(results, function(r) {
  data.frame(
    Variable    = r$variable,
    ColorSpace  = r$color_space,
    N           = r$n,
    Mean        = round(r$mean, 3),
    SD          = round(r$sd, 3),
    Median      = round(r$median, 3),
    Min         = round(r$min, 3),
    Max         = round(r$max, 3),
    Skewness    = round(r$skewness, 3),
    ExKurtosis  = round(r$kurtosis, 3),
    SW_W        = round(r$sw_stat, 4),
    SW_p        = round(r$sw_p, 4),
    AD_stat     = round(r$ad_stat, 4),
    AD_p        = round(r$ad_p, 4),
    Normality   = r$normality,
    AIC_Normal  = round(r$aic_norm, 2),
    AIC_LogNorm = round(r$aic_lnorm, 2),
    AIC_Gamma   = round(r$aic_gamma, 2),
    AIC_Weibull = round(r$aic_weibull, 2),
    BestFit     = r$best_dist,
    stringsAsFactors = FALSE
  )
}))

out_csv <- file.path(OUT_DIR, "distribution_summary.csv")
write.csv(summary_df, out_csv, row.names = FALSE)

cat("=== Distribution Summary ===\n")
print(summary_df[, c("Variable", "ColorSpace", "Skewness", "ExKurtosis",
                      "SW_p", "AD_p", "Normality", "BestFit")])
cat("\nSummary saved to:", out_csv, "\n\n")

# ---------- Visualization Helpers ----------

# Color palette
STATUS_COLOR <- c(NORMAL = "#1a7a1a", "NON-NORMAL" = "#c0392b")
HIST_FILL    <- c(NORMAL = "#2ecc71", "NON-NORMAL" = "#e74c3c")
SPACE_COLOR  <- c(RGB = "#3498DB", HSB = "#9B59B6", LAB = "#E67E22")

make_density_curve <- function(r, xseq) {
  f <- r$fits[[r$best_dist]]
  if (is.null(f)) return(NULL)
  params <- as.list(f$estimate)

  # For lnorm/gamma/weibull the density is defined on positive support;
  # clip xseq to positive values for those distributions
  if (r$best_dist != "norm") xseq <- xseq[xseq > 0]

  y <- tryCatch(
    do.call(paste0("d", r$best_dist), c(list(xseq), params)),
    error = function(e) NULL
  )
  if (is.null(y)) return(NULL)
  data.frame(x = xseq, y = y)
}

make_panel <- function(r) {
  status <- r$normality
  fill_c <- HIST_FILL[[status]]
  xseq   <- seq(min(r$data), max(r$data), length.out = 500)
  x_df   <- data.frame(x = r$data)

  # --- Histogram + fitted density ---
  p_hist <- ggplot(x_df, aes(x = x)) +
    geom_histogram(aes(y = after_stat(density)), bins = 15,
                   fill = fill_c, color = "white", alpha = 0.75) +
    geom_density(color = "black", linewidth = 0.7, linetype = "solid")

  curve_df <- make_density_curve(r, xseq)
  if (!is.null(curve_df)) {
    p_hist <- p_hist +
      geom_line(data = curve_df, aes(x = x, y = y),
                color = "#2C3E50", linewidth = 1.1, linetype = "dashed")
  }

  p_hist <- p_hist +
    labs(
      title    = paste0(r$variable, "  [", r$color_space, "]"),
      subtitle = paste0(status, "  |  SW p=", round(r$sw_p, 3),
                        "  |  Skew=", round(r$skewness, 2),
                        "  |  Best fit: ", r$best_dist),
      x = r$variable, y = "Density"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold",
                                   color = STATUS_COLOR[[status]], size = 11),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )

  # --- Q-Q Plot ---
  p_qq <- ggplot(x_df, aes(sample = x)) +
    stat_qq(color = fill_c, size = 1.8, alpha = 0.75) +
    stat_qq_line(color = "#2C3E50", linewidth = 0.9) +
    labs(
      title    = "Normal Q-Q",
      subtitle = paste0("AD p=", round(r$ad_p, 3)),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_classic(base_size = 10) +
    theme(plot.subtitle = element_text(size = 8, color = "grey40"))

  p_hist | p_qq
}

# ---------- Generate Individual Plots ----------
cat("Generating per-variable plots...\n")
for (v in MEAN_VARS) {
  p <- make_panel(results[[v]])
  out_png <- file.path(OUT_DIR, paste0(v, "_distribution.png"))
  ggsave(out_png, plot = p, width = 10, height = 4, dpi = 150)
  cat("  Saved:", out_png, "\n")
}

# ---------- Composite Panel (all 9 variables) ----------
cat("\nBuilding composite panel...\n")
panels    <- lapply(results, make_panel)
composite <- wrap_plots(panels, ncol = 2) +
  plot_annotation(
    title    = "Distribution Analysis — Soil Pathology Color Data (MSU)",
    subtitle = paste0(
      "Green title = NORMAL  |  Red title = NON-NORMAL  |  Dashed curve = best-fit distribution\n",
      "Normality criterion: Shapiro-Wilk AND Anderson-Darling p > ", ALPHA
    ),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )
  )

out_composite <- file.path(OUT_DIR, "00_all_distributions_composite.png")
ggsave(out_composite, plot = composite, width = 16, height = 20, dpi = 150)
cat("Composite saved:", out_composite, "\n\n")

# ---------- Final Report ----------
normal_vars   <- summary_df$Variable[summary_df$Normality == "NORMAL"]
nonnormal_vars <- summary_df$Variable[summary_df$Normality == "NON-NORMAL"]

cat("========================================\n")
cat("STAGE 1 COMPLETE — Distribution Analysis\n")
cat("========================================\n")
cat("Normal variables    :", paste(normal_vars,    collapse = ", "),
    if (length(normal_vars) == 0) "none" else "", "\n")
cat("Non-normal variables:", paste(nonnormal_vars, collapse = ", "),
    if (length(nonnormal_vars) == 0) "none" else "", "\n")
cat("\nNext step: source('scripts/pipeline/02_soil_normalization.R')\n")
cat("Outputs in:", OUT_DIR, "\n")
