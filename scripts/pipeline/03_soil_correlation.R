# ============================================================
# Tool: Correlation Matrix — Soil Color Data (Parametric)
# Project: Soil Pathology Color Analysis — Michigan State University
# Stage: 3
# Input:  data/processed/normalized_data.csv
# Output: output/03_correlation_soil/
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
DATA_FILE <- "data/processed/normalized_data.csv"
OUT_DIR   <- "output/03_correlation_soil"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ALPHA <- 0.05

# ---------- Load & Assemble 9-Variable Dataset ----------
cat("Loading:", DATA_FILE, "\n")
df <- read.csv(DATA_FILE, stringsAsFactors = FALSE)

# Variable selection:
#   G, B, L  -> original values (confirmed NORMAL in Stage 1)
#   R, H, S, Br, a, b -> OrderNorm-transformed *_norm columns (Stage 2)
cor_data <- data.frame(
  R  = df$R_norm,
  G  = df$G,
  B  = df$B,
  H  = df$H_norm,
  S  = df$S_norm,
  Br = df$Br_norm,
  L  = df$L,
  a  = df$a_norm,
  b  = df$b_norm,
  check.names = FALSE
)

cat("Variables:", paste(names(cor_data), collapse = ", "),
    "| n =", nrow(cor_data), "\n\n")

# ---------- Pearson Correlation + p-values ----------
cor_res <- Hmisc::rcorr(as.matrix(cor_data), type = "pearson")
r_mat   <- cor_res$r
p_mat   <- cor_res$P

# FDR correction (Benjamini-Hochberg) — applied to off-diagonal elements
p_flat     <- as.vector(p_mat)
non_na     <- !is.na(p_flat)
p_adj_flat <- rep(NA_real_, length(p_flat))
p_adj_flat[non_na] <- p.adjust(p_flat[non_na], method = "BH")
p_adj <- matrix(p_adj_flat, nrow = nrow(p_mat), dimnames = dimnames(p_mat))

# ---------- Console Summary ----------
cat("Pearson correlation matrix (r):\n")
print(round(r_mat, 3))
cat("\nFDR-adjusted p-values (BH):\n")
print(round(p_adj, 4))

n_sig <- sum(p_adj[lower.tri(p_adj)] < ALPHA, na.rm = TRUE)
n_tot <- sum(lower.tri(p_adj))
cat(sprintf("\nSignificant pairs (FDR p < %.2f): %d / %d\n\n", ALPHA, n_sig, n_tot))

# ---------- Save CSV Outputs ----------
write.csv(round(r_mat, 4), file.path(OUT_DIR, "soil_corr_r.csv"))
write.csv(round(p_mat,  4), file.path(OUT_DIR, "soil_corr_pval_raw.csv"))
write.csv(round(p_adj,  4), file.path(OUT_DIR, "soil_corr_pval_fdr.csv"))

# ---------- Publication-Ready Correlation Plot ----------
rdbu_colors <- c("#2166AC", "white", "#D6604D")   # standard RdBu

p_corrplot <- ggcorrplot(
  r_mat,
  p.mat         = p_adj,
  hc.order      = TRUE,          # hierarchical clustering order
  type          = "lower",
  method        = "square",
  lab           = TRUE,
  lab_size      = 4.2,
  digits        = 2,
  colors        = rdbu_colors,
  outline.color = "white",
  insig         = "blank",       # blank non-significant cells
  sig.level     = ALPHA,
  tl.cex        = 11,
  tl.col        = "black",
  ggtheme       = theme_bw(base_size = 12)
) +
  labs(
    title    = "Soil Color Data — Pearson Correlation Matrix",
    subtitle = paste0(
      "RGB  •  HSB  •  LAB Color Spaces  |  n = ", nrow(cor_data),
      "  |  FDR-corrected (BH)  |  Blanked: p ≥ ", ALPHA
    ),
    fill = "r"
  ) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey45"),
    legend.title  = element_text(face = "bold", size = 10),
    axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold",
                                  color = "black", size = 11),
    axis.text.y   = element_text(face = "bold", color = "black", size = 11),
    panel.grid    = element_blank()
  )

out_png <- file.path(OUT_DIR, "soil_correlation_matrix.png")
ggsave(out_png, plot = p_corrplot, width = 8, height = 7.5, dpi = 300)
cat("Plot saved:", out_png, "\n")

# ---------- Final Report ----------
cat("\n========================================\n")
cat("STAGE 3 COMPLETE — Soil Correlation\n")
cat("========================================\n")
cat("Outputs in:", OUT_DIR, "\n")
cat("  soil_corr_r.csv          - correlation coefficients\n")
cat("  soil_corr_pval_raw.csv   - raw p-values\n")
cat("  soil_corr_pval_fdr.csv   - FDR-adjusted p-values\n")
cat("  soil_correlation_matrix.png - publication-ready figure\n")
cat("\nNext: source('scripts/pipeline/04_soil_distribution_normalization.R')\n")
