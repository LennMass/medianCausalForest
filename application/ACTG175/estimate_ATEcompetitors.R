############################
# Estimation of ATE competitors for the empirical application based on 
# ACTG 175 (Hammer et al. (1996), Leqi and Kennedy (2022))
############################

rm(list = ls())

path_out <- "application/ACTG175/"

# load config file
source("01_config.R")

# specify number of cores
n.cores <- parallel::detectCores() - 2

# ------------------------------------------------
# Real data example: ACTG175
# ------------------------------------------------
# Load the ACTG 175 dataset
data(ACTG175)

# Binary treatment: A = 0 for zidovudine only (arms == 0), A = 1 for combination therapies (arms != 0)
ACTG175$A <- as.integer(ACTG175$arms != 0)

# Outcome: CD4 count at 96 +/- 5 weeks
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

Z <- dat$A
Y <- dat$cd496


set.seed(42)

B <- 1000 # bootstrap reps
alpha <- 0.05 # CI level
zc <- qnorm(1 - alpha/2)

# EIF and WAQ from Athey et al. (2023)
y0 <- Y[Z == 0]; y1 <- Y[Z == 1]
eif <- as.numeric(parTreat::eif_additive(y0, y1)$tau)
waq <- as.numeric(parTreat::waq(y0, y1)$tau)

That <- c(stat_unadj(Y, Z), stat_adj(Y, Z, X), eif, waq)

n  <- length(Z); n1 <- sum(Z == 1)
Tperm <- do.call(rbind, 
                 pbmclapply(1:B, 
                            function(b){
                              set.seed(b)
                              zpi <- integer(n); zpi[sample.int(n, n1, replace = FALSE)] <- 1
                              # EIF and WAQ from Athey et al. (2023)
                              y0.zpi <- Y[zpi == 0]; y1.zpi <- Y[zpi == 1]
                              eif.zpi <- as.numeric(parTreat::eif_additive(y0.zpi, y1.zpi)$tau)
                              waq.zpi <- as.numeric(parTreat::waq(y0.zpi, y1.zpi)$tau)
                              c(stat_unadj(Y, zpi), stat_adj(Y, zpi, X), eif.zpi, waq.zpi)
                            },
                            mc.cores = n.cores
                 )
)
dim(Tperm)

# Two-sided (1 - alpha) CI via permutation-based critical values
c_abs <- apply(abs(na.omit(Tperm)), 
               2, 
               quantile, 
               probs = 1 - alpha, 
               names = FALSE)
CI_lo <- That - c_abs; CI_hi <- That + c_abs

# output
out <- data.frame(
  methods = c('Difference-in-Means',
              'Lin estimator', 
              'Rosenbaum adj',
              'eif (Athey et al., 2023)', 
              'waq (Athey et al., 2023)'),
  estimate = as.numeric(That),
  std.err = as.numeric(CI_hi - CI_lo)/(2*zc),
  ci_lo = as.numeric(CI_lo),
  ci_hi = as.numeric(CI_hi),
  row.names = NULL
)

print(out, digits = 4)

saveRDS(out, here::here(paste0(path_out, "results/ACTG175_ATE_competitors.rds")))
