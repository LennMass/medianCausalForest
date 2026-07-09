#!/usr/bin/env Rscript
# Install all dependencies with pak. Run once from the project root:
#   Rscript setup.R
# upgrade = FALSE keeps installed packages in place. dependencies = NA
# installs hard dependencies only (Depends, Imports, LinkingTo), no Suggests.


# declare root
here::i_am("01_config.R")


# Load necessary packages that we want to have in the search path
loadpacks <- c(
  "tidyverse", "here", "randomForest", "SimDesign", "sandwich",
  "DescTools", "parTreat", "rpart", "robustCausalTree"
)

suppressPackageStartupMessages(
  invisible(lapply(loadpacks, library, character.only = TRUE))
)

# Project functions
source(here::here("simulation", "helper.R"))
source(here::here("simulation", "dgp.R"))
source(here::here("simulation", "estimators.R"))


