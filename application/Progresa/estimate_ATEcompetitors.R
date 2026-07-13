############################
# Estimation of ATE competitors for the empirical application based on 
# Progresa (De la O (2013), Ghosh et al. (2026))
############################

rm(list = ls())

path_out <- "application/Progresa/"

# load config file
source("01_config.R")

# specify number of cores
n.cores <- parallel::detectCores() - 2

# ------------------------------------------------
# Real data example: Progresa data
# ------------------------------------------------
set.seed(42)
Progressa = read.csv(paste0(path_out, "data/PROGRESA.csv"))
Y = Progressa$pri2000s
Z = Progressa$treatment
X = with(Progressa, data.frame(villages, pri1994, pan1994, prd1994,
                               votos1994, avgpoverty, pobtot1994))
X$villages = as.factor(X$villages)
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

saveRDS(out, paste0(path_out, "results/progresa_ATE_competitors.rds"))



