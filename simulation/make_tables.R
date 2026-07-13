###############################################################################
# Aggregation of results: collapse per-replication rows into per-scenario summaries
###############################################################################

source(here::here("01_config.R"))

path_in <- "simulation/raw_results"
path_out <- "output/sim_tables/"


# Read every results file
files <- list.files(here::here(path_in), pattern = "^results_S[0-9]+_seed.*\\.rds$",
                    full.names = TRUE)
all_results <- do.call(rbind, lapply(files, readRDS))

# order methods
method_order <- c(
  "CF (MSE)",
  "CF (LMS)",
  "CF (MAD)",
  "CF (MSD)",
  "Diff.Mean",
  "R.Adj",
  "Lin.OLS",
  "EIF",
  "WAQ"
)

# renaming
all_results <- all_results %>%
  mutate(method = case_when(
    method == "CT" ~ "CF (MSE)",
    method == "MSD" ~ "CF (MSD)",
    method == "MAD" ~ "CF (MAD)",
    method == "LMS" ~ "CF (LMS)",
    method == "DM" ~ "Diff.Mean",
    method == "Lin" ~ "Lin.OLS",
    method == "RosenbaumAdj" ~ "R.Adj",
    method == "EIF" ~ "EIF",
    method == "WAQ" ~ "WAQ"
    )
  ) %>%
  dplyr::filter(method %in% method_order)
  




# ---- CATE results (full, unordered) -------------------------------------------------
cate_table <- all_results %>%
  dplyr::filter(estimand == "CATE") %>%
  group_by(scenario, method) %>%
  summarise(
    n_reps   = n(),    
    RMSE     = mean(rmse,     na.rm = TRUE),
    RMSE_se  = sd(rmse,       na.rm = TRUE) / sqrt(sum(!is.na(rmse))),
    
    AbsBias  = mean(absbias,  na.rm = TRUE),
    Bias_se  = sd(absbias,    na.rm = TRUE) / sqrt(sum(!is.na(absbias))),
    
    Coverage = mean(coverage, na.rm = TRUE),
    Cov_se   = sd(coverage,   na.rm = TRUE) / sqrt(sum(!is.na(coverage))),
    
    Width    = mean(width,    na.rm = TRUE),
    Wid_se   = sd(width,      na.rm = TRUE) / sqrt(sum(!is.na(width))),
    
    .groups  = "drop"
  ) %>%
  arrange(scenario, method) %>%
  mutate(
    across(c(RMSE, AbsBias, Coverage, Width), ~ round(.x, 2)),
    across(ends_with("_se"), ~ round(.x, 3))
  )

# ---- Table 6 in Appendix ---------------------------------------

# S1 block (S5 is S1 with N=2000 here)
cate_table %>% 
  filter(scenario %in% c("S1", "S5")) %>%
  mutate(scenario = case_when(
    scenario == "S1" ~ "S1 (N=1000)",
    scenario == "S5" ~ "S1 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S2 block (S6 is S2 with N=2000 here)
cate_table %>% 
  filter(scenario %in% c("S2", "S6")) %>%
  mutate(scenario = case_when(
    scenario == "S2" ~ "S2 (N=1000)",
    scenario == "S6" ~ "S2 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S3 block (S7 is S3 with N=2000 here)
cate_table %>% 
  filter(scenario %in% c("S3", "S7")) %>%
  mutate(scenario = case_when(
    scenario == "S3" ~ "S3 (N=1000)",
    scenario == "S7" ~ "S3 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S4 block (S8 is S4 with N=2000 here)
cate_table %>% 
  filter(scenario %in% c("S4", "S8")) %>%
  mutate(scenario = case_when(
    scenario == "S4" ~ "S4 (N=1000)",
    scenario == "S8" ~ "S4 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# ---- Table 8 in Appendix ---------------------------------------

# Prepare display data: combine estimate + SE in parentheses
fmt <- function(x, se) sprintf("%.3f (%.3f)", x, se)

# S1 block (S13 (k=5), S9 (k=20) here)
df_wide <-  cate_table %>% 
  filter(scenario %in% c("S13", "S9")) %>%
  mutate(
    K_label = recode(scenario, "S13" = "K5", "S9" = "K20"),
    scenario_label = "S1",   # both map to S1
    RMSE_fmt  = fmt(RMSE,     RMSE_se),
    Bias_fmt  = fmt(AbsBias,  Bias_se),
    Cov_fmt   = fmt(Coverage, Cov_se),
    Width_fmt = fmt(Width,    Wid_se)
  ) %>%
  select(scenario_label, method, K_label, RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt) %>%
  pivot_wider(
    names_from  = K_label,
    values_from = c(RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt),
    names_glue  = "{K_label}_{.value}"
  ) %>%
  select(
    scenario_label, method,
    K5_RMSE_fmt,  K5_Bias_fmt,  K5_Cov_fmt,  K5_Width_fmt,
    K20_RMSE_fmt, K20_Bias_fmt, K20_Cov_fmt, K20_Width_fmt
  )

df_wide %>%
  kableExtra::kbl(format   = "latex",
      col.names = c("Scenario", "Method",
                    "RMSE", "|Bias|", "Cov.", "Width",
                    "RMSE", "|Bias|", "Cov.", "Width"),
      booktabs  = TRUE,
      align     = c("l", "l", rep("r", 8)),
      caption   = "Covariate size comparison for CATE."
  ) %>%
  kableExtra::kable_classic(full_width = FALSE) %>%
  kableExtra::add_header_above(c(" " = 2, "$K = 5$" = 4, "$K = 20$" = 4),
                   escape = FALSE) %>%
  kableExtra::collapse_rows(columns = 1, valign = "middle", latex_hline = "major")

# S2 block (S14 (k=5), S10 (k=20) here)
df_wide <-  cate_table %>% 
  filter(scenario %in% c("S14", "S10")) %>%
  mutate(
    K_label = recode(scenario, "S14" = "K5", "S10" = "K20"),
    scenario_label = "S2",   # both map to S1
    RMSE_fmt  = fmt(RMSE,     RMSE_se),
    Bias_fmt  = fmt(AbsBias,  Bias_se),
    Cov_fmt   = fmt(Coverage, Cov_se),
    Width_fmt = fmt(Width,    Wid_se)
  ) %>%
  select(scenario_label, method, K_label, RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt) %>%
  pivot_wider(
    names_from  = K_label,
    values_from = c(RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt),
    names_glue  = "{K_label}_{.value}"
  ) %>%
  select(
    scenario_label, method,
    K5_RMSE_fmt,  K5_Bias_fmt,  K5_Cov_fmt,  K5_Width_fmt,
    K20_RMSE_fmt, K20_Bias_fmt, K20_Cov_fmt, K20_Width_fmt
  )

df_wide %>%
  kableExtra::kbl(format   = "latex",
      col.names = c("Scenario", "Method",
                    "RMSE", "|Bias|", "Cov.", "Width",
                    "RMSE", "|Bias|", "Cov.", "Width"),
      booktabs  = TRUE,
      align     = c("l", "l", rep("r", 8)),
      caption   = "Covariate size comparison for CATE."
  ) %>%
  kableExtra::kable_classic(full_width = FALSE) %>%
  kableExtra::add_header_above(c(" " = 2, "$K = 5$" = 4, "$K = 20$" = 4),
                   escape = FALSE) %>%
  kableExtra::collapse_rows(columns = 1, valign = "middle", latex_hline = "major")

# S3 block (S15 (k=5), S11 (k=20) here)
df_wide <-  cate_table %>% 
  filter(scenario %in% c("S15", "S11")) %>%
  mutate(
    K_label = recode(scenario, "S15" = "K5", "S11" = "K20"),
    scenario_label = "S3",   # both map to S1
    RMSE_fmt  = fmt(RMSE,     RMSE_se),
    Bias_fmt  = fmt(AbsBias,  Bias_se),
    Cov_fmt   = fmt(Coverage, Cov_se),
    Width_fmt = fmt(Width,    Wid_se)
  ) %>%
  select(scenario_label, method, K_label, RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt) %>%
  pivot_wider(
    names_from  = K_label,
    values_from = c(RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt),
    names_glue  = "{K_label}_{.value}"
  ) %>%
  select(
    scenario_label, method,
    K5_RMSE_fmt,  K5_Bias_fmt,  K5_Cov_fmt,  K5_Width_fmt,
    K20_RMSE_fmt, K20_Bias_fmt, K20_Cov_fmt, K20_Width_fmt
  )

df_wide %>%
  kableExtra::kbl(format   = "latex",
      col.names = c("Scenario", "Method",
                    "RMSE", "|Bias|", "Cov.", "Width",
                    "RMSE", "|Bias|", "Cov.", "Width"),
      booktabs  = TRUE,
      align     = c("l", "l", rep("r", 8)),
      caption   = "Covariate size comparison for CATE."
  ) %>%
  kableExtra::kable_classic(full_width = FALSE) %>%
  kableExtra::add_header_above(c(" " = 2, "$K = 5$" = 4, "$K = 20$" = 4),
                   escape = FALSE) %>%
  kableExtra::collapse_rows(columns = 1, valign = "middle", latex_hline = "major")

# S4 block (S16 (k=5), S12 (k=20) here)
df_wide <-  cate_table %>% 
  filter(scenario %in% c("S16", "S12")) %>%
  mutate(
    K_label = recode(scenario, "S16" = "K5", "S12" = "K20"),
    scenario_label = "S4",   # both map to S1
    RMSE_fmt  = fmt(RMSE,     RMSE_se),
    Bias_fmt  = fmt(AbsBias,  Bias_se),
    Cov_fmt   = fmt(Coverage, Cov_se),
    Width_fmt = fmt(Width,    Wid_se)
  ) %>%
  select(scenario_label, method, K_label, RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt) %>%
  pivot_wider(
    names_from  = K_label,
    values_from = c(RMSE_fmt, Bias_fmt, Cov_fmt, Width_fmt),
    names_glue  = "{K_label}_{.value}"
  ) %>%
  select(
    scenario_label, method,
    K5_RMSE_fmt,  K5_Bias_fmt,  K5_Cov_fmt,  K5_Width_fmt,
    K20_RMSE_fmt, K20_Bias_fmt, K20_Cov_fmt, K20_Width_fmt
  )

df_wide %>%
  kableExtra::kbl(format   = "latex",
      col.names = c("Scenario", "Method",
                    "RMSE", "|Bias|", "Cov.", "Width",
                    "RMSE", "|Bias|", "Cov.", "Width"),
      booktabs  = TRUE,
      align     = c("l", "l", rep("r", 8)),
      caption   = "Covariate size comparison for CATE."
  ) %>%
  kableExtra::kable_classic(full_width = FALSE) %>%
  kableExtra::add_header_above(c(" " = 2, "$K = 5$" = 4, "$K = 20$" = 4),
                   escape = FALSE) %>%
  kableExtra::collapse_rows(columns = 1, valign = "middle", latex_hline = "major")


# ---- ATE table (forests + classicals) ---------------------------------------
ate_table <- all_results %>%
  dplyr::filter(estimand == "ATE") %>%
  dplyr::filter(method %in% method_order) %>%
  group_by(scenario, method) %>%
  summarise(
    n_reps   = n(),
    
    # RMSE: sqrt of mean squared error
    # MCSE via delta method: MCSE(sqrt(MSE)) = MCSE(MSE) / (2 * RMSE)
    mse      = mean(signed_err^2,  na.rm = TRUE),
    mse_se   = sd(signed_err^2,    na.rm = TRUE) / sqrt(sum(!is.na(signed_err))),
    RMSE     = sqrt(mse),
    RMSE_se  = mse_se / (2 * RMSE),
    
    # AbsBias: a plain mean, so MCSE = sd / sqrt(B)
    AbsBias  = mean(absbias,       na.rm = TRUE),
    Bias_se  = sd(absbias,         na.rm = TRUE) / sqrt(sum(!is.na(absbias))),
    
    # Coverage
    Coverage = mean(coverage,      na.rm = TRUE),
    Cov_se   = sd(coverage,        na.rm = TRUE) / sqrt(sum(!is.na(coverage))),
    
    # Width
    Width    = mean(width,         na.rm = TRUE),
    Wid_se   = sd(width,           na.rm = TRUE) / sqrt(sum(!is.na(width))),
    
    .groups  = "drop"
  ) %>%
  dplyr::select(-mse, -mse_se) %>%          # drop intermediate columns
  arrange(scenario, method)

method_order <- c(
  "CF (MSE)",
  "CF (LMS)",
  "CF (MAD)",
  "CF (MSD)", 
  "Diff.Mean",
  "Lin.OLS", 
  "R.Adj",
  "EIF",
  "WAQ"
)

ate_table <- ate_table %>%
  dplyr::filter(method %in% method_order) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  arrange(scenario, method) %>%
  mutate(
    across(c(RMSE, AbsBias, Coverage, Width), ~ round(.x, 2)),
    across(ends_with("_se"), ~ round(.x, 3))
  )


# ---- Table 7 in Appendix ---------------------------------------

# S1 block (S5 is S1 with N=2000 here)
ate_table %>%
  filter(scenario %in% c("S1", "S5")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S1" ~ "S1 (N=1000)",
    scenario == "S5" ~ "S1 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S2 block (S6 is S2 with N=2000 here)
ate_table %>%
  filter(scenario %in% c("S2", "S6")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S2" ~ "S2 (N=1000)",
    scenario == "S6" ~ "S2 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S3 block (S7 is S3 with N=2000 here)
ate_table %>%
  filter(scenario %in% c("S3", "S7")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S3" ~ "S3 (N=1000)",
    scenario == "S7" ~ "S3 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S4 block (S8 is S4 with N=2000 here)
ate_table %>%
  filter(scenario %in% c("S4", "S8")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S4" ~ "S4 (N=1000)",
    scenario == "S8" ~ "S4 (N=2000)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")


# ---- Table 9 in Appendix ---------------------------------------

# S1 block (S9 (K=20), S13 (K=5) here)
ate_table %>%
  filter(scenario %in% c("S13", "S9")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S13" ~ "S1 (K=5)",
    scenario == "S9" ~ "S1 (N=20)")
  ) %>%
  dplyr::select(-c(n_reps)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S2 block (S10 (K=20), S14 (N=5) here)
ate_table %>%
  filter(scenario %in% c("S14", "S10")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S14" ~ "S2 (K=5)",
    scenario == "S10" ~ "S2 (N=20)")
  ) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S3 block (S11 (K=20), S15 (N=5) here)
ate_table %>%
  filter(scenario %in% c("S11", "S15")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S15" ~ "S3 (K=5)",
    scenario == "S11" ~ "S3 (N=20)")
  ) %>%
  dplyr::select(-c(n_reps, scenario)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# S4 block (S12 (K=20), S16 (N=5) here)
ate_table %>%
  filter(scenario %in% c("S12", "S16")) %>%
  filter(!method %in% c("CF (LMS)", "CF (MAD)")) %>% # focus on MSE, MSD
  mutate(scenario = case_when(
    scenario == "S16" ~ "S4 (K=5)",
    scenario == "S12" ~ "S4 (N=20)")
  ) %>%
  filter(!method %in% c("Diff.Median", "HL", "CF (LMS)", "CF (MAD)")) %>%
  dplyr::select(-c(n_reps, scenario)) %>%
  mutate(
    RMSE = paste0(RMSE, " (", RMSE_se, ")"),
    AbsBias = paste0(AbsBias, " (", Bias_se, ")"),
    Coverage = paste0(Coverage, " (", Cov_se, ")"),
    Width = paste0(Width, " (", Wid_se, ")")
  ) %>%
  dplyr::select(-c(RMSE_se, Bias_se, Cov_se, Wid_se)) %>%
  kableExtra::kbl(format   = "latex",
      booktabs = TRUE,
      linesep  = "")

# ---- Save raw summaries -----------------------------------------------------------
print(all_results)
print(cate_table, n = Inf)
print(ate_table,  n = Inf)
saveRDS(all_results, here::here(paste0(path_out, "all_results.rds")))
saveRDS(cate_table, here::here(paste0(path_out, "summary_CATE.rds")))
saveRDS(ate_table,  here::here(paste0(path_out, "summary_ATE.rds")))

