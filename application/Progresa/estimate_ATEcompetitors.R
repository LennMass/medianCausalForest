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


stat_unadj <- function(y, z, trim = 0.1, wins = 0.1) {
  y1 <- y[z==1]; y0 <- y[z==0]
  c(
    DM = mean(y1) - mean(y0)
  )
}

stat_adj <- function(y, z, x) {
  n <- length(y)
  if (length(z) != n) stop("y and z must have equal length.")
  X <- if (is.null(dim(x))) data.frame(x = x) else as.data.frame(x)
  if (nrow(X) != n) stop("nrow(x) must equal length(y).")
  
  # coerce character -> factor; leave existing factors as is
  X[] <- lapply(X, function(col) if (is.character(col)) factor(col) else col)
  
  # split controls
  is_num <- vapply(X, is.numeric, TRUE)
  num_names <- names(X)[is_num]
  fac_names <- names(X)[!is_num]
  
  # drop NAs
  df_un <- cbind.data.frame(y = y, z = z, X)
  df_un <- df_un[stats::complete.cases(df_un), , drop = FALSE]
  N <- nrow(df_un); m <- sum(df_un$z == 1L)
  
  # numeric and factor terms
  num_terms_reg <- if (length(num_names)) paste(num_names, collapse = " + ") else ""
  fac_terms_reg <- if (length(fac_names)) paste0("as.factor(", fac_names, ")", collapse = " + ") else ""

  
  # ---------- tau.interact (Lin, 2013): center ONLY numeric X, interact z with centered numerics; add factors as as.factor(...) ----------
  if (length(num_names)) {
    Xc <- X
    Xc[num_names] <- lapply(Xc[num_names], function(col) col - mean(col, na.rm = TRUE))
    cn <- paste0(num_names, "_c")
    names(Xc)[match(num_names, names(Xc))] <- cn
    df_int <- cbind.data.frame(y = y, z = z, Xc)
    
    num_terms_int <- paste(cn, collapse = " + ")
    fac_terms_int <- if (length(fac_names)) paste0(" + ", paste0("as.factor(", fac_names, ")", collapse = " + ")) else ""
    f_int <- stats::as.formula(paste0("y ~ z * (", num_terms_int, ")", fac_terms_int))
    fit2 <- stats::lm(f_int, data = df_int)
    tau.interact <- unname(coef(fit2)["z"])
  } else {
    # no numeric covariates to interact with
    rhs_int <- paste(Filter(nzchar, c("z", fac_terms_reg)), collapse = " + ")
    fit2 <- stats::lm(stats::as.formula(paste0("y ~ ", rhs_int)), data = df_un)
    tau.interact <- unname(coef(fit2)["z"])
  }
  
  # ---------- tau.Radj (Ghosh et al., 2026): root of rank balance after adjusting (y - tau z) on UNcentered X with as.factor for factors ----------
  ctrl_rhs <- paste(
    Filter(nzchar, c(num_terms_reg, fac_terms_reg)),
    collapse = " + "
  )
  aux <- function(tau) {
    y_shift <- df_un$y - tau * df_un$z
    fit <- stats::lm(stats::as.formula(
      if (nzchar(ctrl_rhs)) paste0("y_shift ~ ", ctrl_rhs) else "y_shift ~ 1"
    ), data = transform(df_un, y_shift = y_shift))
    r <- resid(fit)
    sum(rank(r)[df_un$z == 1L]) - m * (N + 1) / 2
  }
  tau.Radj <- uniroot(aux, c(0, 1), extendInt = "yes")$root
  
  c(tau.interact = tau.interact, tau.Radj = tau.Radj)
}

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