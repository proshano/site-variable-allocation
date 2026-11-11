# Load the packages and functions used by every script. Run the scripts from
# this folder so that paths such as R/ and data/ resolve.
required_packages <- c(
  "digest", "dplyr", "jsonlite", "lme4", "readr", "sandwich", "survival"
)
is_installed <- vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)
if (!all(is_installed)) {
  stop(
    "Install these packages first: ",
    paste(required_packages[!is_installed], collapse = ", "),
    ". See README.md."
  )
}

# The model formulas in R/models.R use survival's functions directly.
suppressPackageStartupMessages(library(survival))

for (source_file in c(
  "config", "scenarios", "simulate", "fit_helpers", "models", "summaries",
  "run", "planning_tables", "simulation_tables", "diagnostic_tables"
)) {
  source(file.path("R", paste0(source_file, ".R")))
}
