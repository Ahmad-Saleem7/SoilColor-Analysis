# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

pkgs <- c("Hmisc", "ggcorrplot", "janitor", "car", "GGally", "RColorBrewer")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing:", p, "\n")
    install.packages(p, repos = "https://cran.r-project.org", quiet = FALSE)
  } else {
    cat("Already installed:", p, "\n")
  }
}
cat("All packages ready.\n")
