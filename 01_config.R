#!/usr/bin/env Rscript
# Install all dependencies with pak. Run once from the project root:
#   Rscript setup.R
# upgrade = FALSE keeps installed packages in place. dependencies = NA
# installs hard dependencies only (Depends, Imports, LinkingTo), no Suggests.

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

cran <- c(
  "data.tree", "DescTools", "devtools", "dplyr", "ggplot2", "grf",
  "here", "Matching", "partykit", "randomForest", "ranger", "rpart",
  "sandwich", "SimDesign", "speff2trial", "xgboost"
)

# CRAN packages plus parTreat from GitHub, in one resolve.
pak::pkg_install(c(cran, "github::michaelpollmann/parTreat"),
                 upgrade = FALSE, dependencies = NA, ask = FALSE)

# robustCausalTree from local source ('here' is installed above).
rct <- here::here("core", "robustCausalTree")
pak::local_install(rct, upgrade = FALSE, dependencies = NA, ask = FALSE)



