###############################################################################
# Create plots of the simulation study; CATE and ATE summaries
###############################################################################

source(here::here("01_config.R"))


all_results <- readRDS(file=here::here("simulation/aggregate_results/all_results.rds"))

path_out <- "output/sim_plots/"

# ---- Boxplots ---------------------------------------------------------------

# Dark2 palette anchored: MSD = orange (highlighted)
method_colours <- c(
  "CF (MSD)"           = "#d95f02"   # orange — main contributed split rule (highlighted)
)


# shared theme to avoid repetition
box_theme <- theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y      = element_text(size = 11),
    axis.title       = element_text(size = 12),
    strip.text       = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 11),
    legend.title     = element_text(size = 12),
    legend.position  = "none",
    strip.background = element_rect(fill = "grey92")
  )





# --- CATE plot: Abs. Bias + Coverage + Width ---
p_cate_main <- all_results %>%
  filter(estimand == "CATE") %>%
  filter(scenario %in% c("S1", "S2", "S3", "S4")) %>%
  mutate(method = factor(method, levels = c("CF (MSE)", "CF (MSD)", "CF (MAD)", "CF (LMS)"))) %>%
  dplyr::select(scenario, method, absbias, coverage, width) %>%
  pivot_longer(c(absbias, coverage, width),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("absbias", "coverage", "width"),
                         labels = c("Abs. Bias", "Coverage", "Width"))) %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  geom_hline(data = data.frame(metric = factor("Coverage",
                                               levels = c("Coverage", "Width")),
                               yintercept = 0.95),
             aes(yintercept = yintercept),
             linetype = "dashed", colour = "black", linewidth = 0.4) +
  facet_grid(metric ~ scenario, scales = "free_y") +
  labs(x = NULL, y = NULL, fill = "Method") +
  box_theme +
  scale_fill_manual(values = method_colours)


ggsave(here::here(paste0(path_out, "boxplot_CATE_main.pdf")),
       p_cate_main, width = 10, height = 6)


# --- ATE precision: Sq. Error (-> RMSE) + Abs. Error (-> AbsBias); no coverage ---
# sq_err per replication aggregates to MSE -> RMSE
# absbias per replication aggregates to AbsBias
p_ate_main <- all_results %>%
  filter(estimand == "ATE") %>%
  filter(scenario %in% c("S1", "S2", "S3", "S4")) %>%
  filter(!method %in% c("CF (MAD)", "CF (LMS)")) %>% # focus on MSE, MSD
  mutate(method = factor(method, levels = c("CF (MSE)", "CF (MSD)", "Diff.Mean", "Lin.OLS", "R.Adj", "EIF", "WAQ"))) %>%
  mutate(sq_err = signed_err^2) %>%
  dplyr::select(scenario, method, absbias) %>%
  pivot_longer(c(absbias),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric,
                         levels = c("absbias"),
                         labels = c("Abs. Error"))) %>%
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_boxplot(outlier.size = 0.4, linewidth = 0.3) +
  facet_grid(metric ~ scenario, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  box_theme +
  scale_fill_manual(values  = method_colours)


ggsave(here::here(paste0(path_out, "boxplot_ATE_main.pdf")),
       p_ate_main, width = 10, height = 6)



