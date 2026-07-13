############################
# Create plots and tables for the empirical application based on 
# Progresa (De la O (2013), Ghosh et al. (2026))
############################

rm(list = ls())

path_in <- "application/Progresa"
path_out <- "output/application_plots/"

# load config file
source("01_config.R")

CF_estimates <- readRDS(here::here(path_in, "results/progresa_CF_estimates.rds"))

MSE.estims <- CF_estimates$CF_MSE
MSD.estims <- CF_estimates$CF_MSD


# Add method labels
MSE.estims$pred_cate <- MSE.estims$pred_cate %>% mutate(method = "CF (MSE)", ATE=MSE.estims$pred_ate[1], ATE.SE=MSE.estims$pred_ate[2])
MSD.estims$pred_cate <- MSD.estims$pred_cate %>% mutate(method = "CF (MSD)", ATE=MSD.estims$pred_ate[1], ATE.SE=MSD.estims$pred_ate[2])

# Combine and compute CIs
plot_df <- bind_rows(MSE.estims$pred_cate, MSD.estims$pred_cate) %>%
  mutate(
    ci_lower = estimate - 1.96 * std.err,
    ci_upper = estimate + 1.96 * std.err,
    panel = ifelse(ci_lower > 0 & !(ci_upper == ci_lower), "stat. significant effect", "no stat. significant effect")
  )

ate_df <- plot_df %>%
  group_by(method) %>%
  summarise(ATE = mean(ATE), 
            ATE.SE = mean(ATE.SE))  # or unique(ATE) if it's constant within method

##### Output plots #####

fig_cate <- ggplot(plot_df, aes(x = id, y = estimate, color = method)) +
  geom_linerange(
    aes(ymin = ci_lower, ymax = ci_upper),
    alpha = 0.3,
    linewidth = 0.3
  ) +
  geom_point(
    size = 0.5,
    alpha = 0.8
  )+
  scale_color_manual(values = c("CF (MSE)" = "#1b9e77", "CF (MSD)" = "#d95f02")) +
  #geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(
    data = ate_df,
    aes(yintercept = ATE),
    linetype = "solid", color = "black"
  ) +
  geom_hline(
    data = ate_df,
    aes(yintercept = ATE + 1.96*ATE.SE),
    linetype = "dashed", color = "black"
  ) +
  geom_hline(
    data = ate_df,
    aes(yintercept = ATE - 1.96*ATE.SE),
    linetype = "dashed", color = "black"
  ) +
  facet_wrap(~ method, scales = "free", ncol=1, nrow=2) +
  labs(
    x = "ID",
    y = "CATE Estimate",
    color = "Method",
    #title = "CATE Estimates with 95% CIs"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top", 
    strip.text = element_blank())

ggsave(here::here(paste0(path_out, "plot_CATE_progresa.pdf")), fig_cate,
       width = 7, height = 3.5)



# Histograms
fig_hist <- ggplot(plot_df, aes(x = estimate, fill = method)) +
  geom_histogram(
    aes(y = after_stat(density)),
    position = "identity",
    alpha = 0.45,
    bins = 40
  ) +
  scale_fill_manual(values = c("CF (MSE)" = "#1b9e77", "CF (MSD)" = "#d95f02")) +
  geom_vline(xintercept = mean(ate_df$ATE),
             linetype = "dashed",
             linewidth = 0.8
  ) +
  #scale_color_manual(values = c("CF (MSE)" = "black", "CF (MSD)" = "black")) +
  labs(
    x = "CATE estimates",
    y = "Density",
    fill = "Method",
    color = "Method",
    #title = "Distribution of CATE estimates"
  ) +
  theme_minimal() +
  theme(legend.position = "top")

ggsave(here::here(paste0(path_out, "histogram_CATE_progresa.pdf")), fig_hist,
       width = 7, height = 3.5)



##### Output tables #####

# ATE competitors
res_table <- readRDS(here::here(path_in, "results/progresa_ATE_competitors.rds")) 

# Add CF ATEs
cf_values <- data.frame(
  methods = c('CF (MSE)',
              'CF (MSD)'),
  estimate = c(MSE.estims$pred_ate[1], MSD.estims$pred_ate[1]),
  std.err = c(MSE.estims$pred_ate[2], MSD.estims$pred_ate[2]),
  ci_lo = c(MSE.estims$pred_ate[1], MSD.estims$pred_ate[1]) - 1.96 * c(MSE.estims$pred_ate[2], MSD.estims$pred_ate[2]),
  ci_hi = c(MSE.estims$pred_ate[1], MSD.estims$pred_ate[1]) + 1.96 * c(MSE.estims$pred_ate[2], MSD.estims$pred_ate[2]),
  row.names = NULL
)

res_table <-  rbind(res_table, cf_values) %>% 
  mutate(ci_width = ci_hi - ci_lo)


saveRDS(res_table, here::here(paste0("output/application_tables/ATE_table_progresa.rds")))

# latex code
res_table %>%
  xtable::xtable(digits=3, caption = "")

