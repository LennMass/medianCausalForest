############################
# Descriptive statistics for the empirical application based on 
# Progresa (De la O (2013), Ghosh et al. (2026))
############################


# load config file
source("01_config.R")

# where to save generated materials
path_out <- "output/application_plots"

set.seed(12345) 

Data <- rio::import(here::here("application/Progresa/data/PROGRESA.csv"))

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_paper <- function() {
  theme_classic(base_size = 11, base_family = "serif") +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      axis.title    = element_text(size = 9),
      axis.text     = element_text(size = 8),
      legend.position = "none",
      strip.text    = element_text(face = "bold", size = 9)
    )
}

# ── Summary stats ─────────────────────────────────────────────────────────────
stats <- Data %>%
  summarise(sk = moments::skewness(pri2000s, na.rm = TRUE))

# ── Panel A: overall density ──────────────────────────────────────────────────
pA <- ggplot(Data, aes(x = pri2000s)) +
  geom_density(fill = "grey80", color = "grey30", linewidth = 0.5, alpha = 0.8) +
  geom_rug(color = "grey50", alpha = 0.3, linewidth = 0.3) +
  labs(
    title    = "A. Full sample",
    subtitle = sprintf("Skewness = %.2f", stats$sk),
    x = "pri2000s", y = "Density"
  ) +
  theme_paper()

# ── Panel B: density by treatment ─────────────────────────────────────────────
plot_data <- Data %>%
  mutate(Group = if_else(treatment == 1, "Treatment", "Control"))

pB <- ggplot(plot_data, aes(x = pri2000s)) +
  geom_density(fill = "grey80", color = "grey30", linewidth = 0.5, alpha = 0.8) +
  geom_rug(color = "grey50", alpha = 0.3, linewidth = 0.3) +
  facet_wrap(~Group, ncol = 1) +
  labs(
    title = "B. By treatment status",
    x = "pri2000s", y = "Density"
  ) +
  theme_paper() +
  theme(strip.background = element_rect(fill = "grey92", color = NA))

# ── Combine & save ────────────────────────────────────────────────────────────
fig <- pA | pB

ggsave("outcome_distribution_progresa.pdf", fig,
       path = here::here(path_out),
       width = 7, height = 3.5)


# Codebook
# Outcome:
# pri2000s: support rates for the incumbent party as shares of the eligible voting population in the 2000 election (pri2000s) as the outcome

# Treatment:
# treatment: 1 - when the precinct includes a village assigned to early treatment; 0 - when the precinct includes a village assigned to late treatment.

# Covariates
# avgpoverty: the average poverty level in a precinct ,
# pobtot1994: the total precinct population in 1994,
# votos1994: the total number of voters who turned out in the previous election
# pri1994: total number of votes cast in the previous election for PRI (incumbent party, Institutional Revolutionary Party) 
# pan1994: total number of votes cast in the previous election for PAN (right-wing party, National Action Party)
# prd1994: total number of votes cast in the previous election for PRD (left-wing party, Party of the Democratic Revolution)
# villages: (as factors)


# ── Sample sizes ──────────────────────────────────────────────────────────────
n_total   <- nrow(Data)
n_control <- sum(Data$treatment == 0)
n_treat   <- sum(Data$treatment == 1)

# ── Variables ─────────────────────────────────────────────────────────────────
outcome    <- "pri2000s"
covariates <- c("avgpoverty", "pobtot1994", "votos1994",
                "pri1994", "pan1994", "prd1994", "villages")
all_vars   <- c(outcome, covariates)

var_labels <- c(
  pri2000s   = "PRI vote share, 2000 (\\%)",
  avgpoverty = "Avg. poverty level",
  pobtot1994 = "Total population (1994)",
  votos1994  = "Voter turnout (1994)",
  pri1994    = "PRI votes (1994)",
  pan1994    = "PAN votes (1994)",
  prd1994    = "PRD votes (1994)",
  villages   = "No. of villages"
)

# ── Build table ───────────────────────────────────────────────────────────────
tab <- Data %>%
  group_by(treatment) %>%
  summarise(across(all_of(all_vars),
                   list(mean   = ~mean(.x,   na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE),
                        sd     = ~sd(.x,     na.rm = TRUE)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(-treatment,
               names_to  = c("variable", "stat"),
               names_sep = "__",
               values_to = "value") %>%
  pivot_wider(names_from  = c(treatment, stat),
              values_from = value,
              names_glue  = "{treatment}_{stat}") %>%
  mutate(
    variable = dplyr::recode(variable, !!!var_labels),
    order    = match(variable, unname(var_labels))
  ) %>%
  arrange(order) %>%
  dplyr::select(variable,
                `0_mean`, `0_median`, `0_sd`,
                `1_mean`, `1_median`, `1_sd`)

# ── Render ────────────────────────────────────────────────────────────────────
tab %>%
  kbl("latex", booktabs = TRUE, digits = 2, escape = FALSE,
      col.names = c("Variable",
                    "Mean", "Median", "SD",
                    "Mean", "Median", "SD"),
      caption = sprintf(
        "Descriptive statistics by treatment status.
         Total sample: $N = %d$ precincts (%d control, %d treatment).",
        n_total, n_control, n_treat),
      label = "tab:desc_stats") %>%
  add_header_above(c(" "                 = 1,
                     "Control (late)"    = 3,
                     "Treatment (early)" = 3)) %>%
  pack_rows("Outcome", 1, 1,
            bold = TRUE, italic = FALSE,
            latex_gap_space = "0.3em") %>%
  pack_rows("Pre-treatment covariates", 2, nrow(tab),
            bold = TRUE, italic = FALSE,
            latex_gap_space = "0.3em") %>%
  add_footnote(sprintf(
    "Sample sizes: $N = %d$ (total), $N_0 = %d$ (control), $N_1 = %d$ (treatment).",
    n_total, n_control, n_treat),
    notation = "none") %>%
  kable_styling(latex_options = c("hold_position"))