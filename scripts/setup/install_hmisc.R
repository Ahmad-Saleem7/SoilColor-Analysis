# Setup user library path to avoid permission errors
local_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_lib, .libPaths()))

options(timeout = 300)   # 5 minute download timeout
pkgs <- c("Hmisc", "gridExtra")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    cat("Installing:", p, "\n")
    install.packages(p, repos = "https://cran.rstudio.com/", quiet = FALSE)
  } else {
    cat("Already installed:", p, "\n")
  }
}
cat("Done.\n")

