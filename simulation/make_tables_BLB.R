# Create table for BLB analysis

# ------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------

source(here::here("01_config.R"))

results_dir <-here::here("simulation/raw_results_BLB")
path_out <- "output/sim_tables/"


# ------------------------------------------------------------------
# Find and load all S1--S4 .rds files
# ------------------------------------------------------------------

files <- list.files(
  path = results_dir,
  pattern = "^results_S[1-4]_.*\\.rds$",
  full.names = TRUE
)

# Load files and add scenario based on the filename
results_all <- map_dfr(files, function(file) {
  
  # Extract scenario (S1, S2, S3, S4) from filename
  scenario <- str_extract(basename(file), "S[1-4]")
  
  # Load data and add scenario column
  readRDS(file) %>%
    mutate(scenario = scenario)
})

# ------------------------------------------------------------------
# Summarise results
# ------------------------------------------------------------------

results_summary <- results_all %>%
  group_by(scenario, method, blb_scaling) %>%
  summarise(
    across(
      c(
        rmse,
        absbias,
        coverage,
        width,
        bias,
        empirical_sd,
        miss_low,
        miss_high
      ),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  arrange(scenario, method, blb_scaling) %>%
  mutate(scenario = factor(scenario, levels = paste0("S", 1:4)),
         blb_scaling = if_else(blb_scaling, "Yes", "No"), 
         method = case_when(method == "CT" ~ "CF (MSE)",
                            method == "MAD" ~ "CF (MAD)",
                            method == "LMS" ~ "CF (LMS)",
                            method == "MSD" ~ "CF (MSD)"))



tab <- results_summary %>%
  arrange(scenario, method, blb_scaling) %>%
  transmute(
    Scenario   = scenario,
    Method     = method,
    Scaling    = blb_scaling,
    Bias       = sprintf("%.3f", bias),
    Coverage   = sprintf("%.3f", coverage),
    Width      = sprintf("%.2f", width),
    `Miss low` = sprintf("%.3f", miss_low),
    `Miss high`= sprintf("%.3f", miss_high)
  )

# save table
saveRDS(tab, here::here(paste0(path_out, "summary_BLB_scaling.rds")))

# create latex code
tab %>%
  kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "",
      align    = c("l", "l", "l", rep("r", 5)),
      caption  = "",
      label    = "blb_scaling_ci") %>%
  #collapse_rows(columns = 1:2, valign = "top", latex_hline = "major") %>%
  kable_styling(latex_options = c("hold_position"), font_size = 9)
