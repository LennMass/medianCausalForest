#!/usr/bin/env Rscript
# Run ONCE, from the project root, in a fresh session:  Rscript setup.R
# Do not source this from analysis scripts. Installing packages while they
# are attached is what corrupts a lazy-load database.

# parTreat is a public repo. Clear any stale token so an expired credential
# cannot trigger an HTTP 401 during resolve.
# Sys.unsetenv("GITHUB_PAT"); Sys.unsetenv("GITHUB_TOKEN"); Sys.unsetenv("GH_TOKEN")

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

cran <- c(
  "data.tree", "DescTools", "devtools", "dplyr", "ggplot2", "grf",
  "here", "Matching", "partykit", "randomForest", "ranger", "rpart",
  "sandwich", "SimDesign", "speff2trial", "kableExtra"
)

pak::pkg_install(c(cran, "github::michaelpollmann/parTreat"),
                 upgrade = FALSE, dependencies = NA, ask = FALSE)

pak::local_install(here::here("01_code", "01_ct_dev_folders", "robustCausalTree"),
                   upgrade = FALSE, dependencies = NA, ask = FALSE)