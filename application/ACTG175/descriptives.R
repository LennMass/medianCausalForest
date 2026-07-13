############################
# Descriptive statistics for the empirical application based on 
# ACTG175 (Hammer et al. (1996), Leqi and Kennedy (2022))
############################


# load config file
source("01_config.R")

# where to save generated materials
path_out <- "output/application_plots"

set.seed(12345) 

# Load the ACTG 175 dataset from speff2trial-package
data(ACTG175)

# Binary treatment: A = 0 for zidovudine only (arms == 0), A = 1 for combination therapies (arms != 0)
ACTG175$A <- as.integer(ACTG175$arms != 0)
Y <- ACTG175$cd496

# Remove subjects with missing outcomes
complete <- !is.na(Y)
dat <- ACTG175[complete, ]

# Covariates
X <- dat[, c("cd40",       # baseline CD4 count
             "cd80",       # baseline CD8 count
             "age",        # age
             "wtkg",       # weight in kg
             "karnof",     # Karnofsky score
             "race",       # race indicator (0 = white, 1 = non-white)
             "gender",     # gender (0 = female, 1 = male)
             "hemo",       # hemophilia indicator
             "homo",       # homosexual activity indicator
             "drugs",      # drug use indicator
             "symptom",    # symptomatic indicator
             "z30",        # previous zidovudine use (> 30 days)
             "preanti")]   # previous antiretroviral use (days)


# Overwrite A (treatment) and Y (outcome) for complete cases
# Binary treatment: A = 0 for zidovudine only (arms == 0), A = 1 for combination therapies (arms != 0)
A <- dat$A
# Outcome: CD4 count at 96 +/- 5 weeks
Y <- dat$cd496

# Known propensity score (randomized trial: ~75% assigned to combination therapy)
ps <- rep(0.75, nrow(dat))

cat("Number of observations:", length(Y), "\n")
cat("Treatment distribution:\n")
print(table(A))
cat("\nCovariate dimensions:", dim(X), "\n")


#### plot outcome variable ####

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_paper <- function() {
  theme_classic(base_size = 11, base_family = "serif") +
    theme(
      plot.title      = element_text(face = "bold", size = 11),
      plot.subtitle   = element_text(size = 9, color = "grey40"),
      axis.title      = element_text(size = 9),
      axis.text       = element_text(size = 8),
      legend.position = "none",
      strip.text      = element_text(face = "bold", size = 9)
    )
}

# ── Panel A: overall density ──────────────────────────────────────────────────
stats_A <- dat %>%
  summarise(sk  = skewness(cd496, na.rm = TRUE),
            med = median(cd496, na.rm = TRUE))

pA <- ggplot(dat, aes(x = cd496)) +
  geom_density(fill = "grey80", color = "grey30", linewidth = 0.5, alpha = 0.8) +
  geom_rug(color = "grey50", alpha = 0.3, linewidth = 0.3) +
  labs(
    title    = "A. Full sample",
    subtitle = sprintf("Skewness = %.2f", stats_A$sk),
    x = "CD4 count at 96 weeks", y = "Density"
  ) +
  theme_paper()

# ── Panel B: density by treatment ─────────────────────────────────────────────
plot_data <- dat %>%
  mutate(Group = if_else(A == 1, "Combination therapy", "Zidovudine only"))

stats_B <- plot_data %>%
  group_by(Group) %>%
  summarise(sk  = skewness(cd496, na.rm = TRUE),
            med = median(cd496, na.rm = TRUE)) %>%
  mutate(label = sprintf("Skewness = %.2f", sk))

pB <- ggplot(plot_data, aes(x = cd496)) +
  geom_density(fill = "grey80", color = "grey30", linewidth = 0.5, alpha = 0.8) +
  geom_rug(color = "grey50", alpha = 0.3, linewidth = 0.3) +
  facet_wrap(~Group, ncol = 1) +
  labs(
    title = "B. By treatment status",
    x = "CD4 count at 96 weeks", y = "Density"
  ) +
  theme_paper() +
  theme(strip.background = element_rect(fill = "grey92", color = NA))

# ── Combine & save ────────────────────────────────────────────────────────────
fig <- pA | pB
ggsave(path = here::here(path_out),
       filename = "outcome_distribution_actg175.pdf",
       fig,
       width = 7, height = 3.5)


##### summary tables ####


# ── Sample sizes ──────────────────────────────────────────────────────────────
n_total   <- nrow(dat)
n_control <- sum(dat$A == 0)
n_treat   <- sum(dat$A == 1)

# ── Variables ─────────────────────────────────────────────────────────────────
outcome    <- "cd496"
covariates <- c("cd40", "cd80", "age", "wtkg", "karnof",
                "race", "gender", "hemo", "homo",
                "drugs", "symptom", "z30", "preanti")
all_vars   <- c(outcome, covariates)

var_labels <- c(
  cd496   = "CD4 count at 96 weeks",
  cd40    = "Baseline CD4 count",
  cd80    = "Baseline CD8 count",
  age     = "Age (years)",
  wtkg    = "Weight (kg)",
  karnof  = "Karnofsky score",
  race    = "Non-white",
  gender  = "Male",
  hemo    = "Hemophilia",
  homo    = "Homosexual activity",
  drugs   = "History of drug use",
  symptom = "Symptomatic",
  z30     = "Prior zidovudine use ($>$30 days)",
  preanti = "Prior antiretroviral use (days)"
)

# ── Build table ───────────────────────────────────────────────────────────────
tab <- dat %>%
  group_by(A) %>%
  summarise(across(all_of(all_vars),
                   list(mean   = ~mean(.x,   na.rm = TRUE),
                        median = ~median(.x, na.rm = TRUE),
                        sd     = ~sd(.x,     na.rm = TRUE)),
                   .names = "{.col}__{.fn}")) %>%
  pivot_longer(-A,
               names_to  = c("variable", "stat"),
               names_sep = "__",
               values_to = "value") %>%
  pivot_wider(names_from  = c(A, stat),
              values_from = value,
              names_glue  = "{A}_{stat}") %>%
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
         Total sample: $N = %d$ subjects (%d zidovudine only, %d combination therapy).",
        n_total, n_control, n_treat),
      label = "tab:desc_stats_actg175") %>%
  add_header_above(c(" "                          = 1,
                     "Control (Zidovudine only)"    = 3,
                     "Treatment (Combination therapy)" = 3)) %>%
  pack_rows("Outcome", 1, 1,
            bold = TRUE, italic = FALSE,
            latex_gap_space = "0.3em") %>%
  pack_rows("Pre-treatment covariates", 2, nrow(tab),
            bold = TRUE, italic = FALSE,
            latex_gap_space = "0.3em") %>%
  add_footnote(sprintf(
    "Sample sizes: $N = %d$ (total), $N_0 = %d$ (zidovudine only), $N_1 = %d$ (combination therapy).",
    n_total, n_control, n_treat),
    notation = "none") %>%
  kable_styling(latex_options = c("hold_position"))