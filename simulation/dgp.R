# ---- DGP (based on Wager and Athey, 2018, Procedure 1) ------------------------------------
sigma        <- function(x) 1 + 1 / (1 + exp(-20 * (x - 1/3)))
tau_function <- function(X) sigma(X[, 1]) * sigma(X[, 2])
heavy_noise <- function(n, df = 3) rt(n, df = df)
small_region_smooth <- function(X, r = 1.2) as.integer(X[, 1]^2 + X[, 2]^2 > r^2)

simulate_data <- function(n, d,
                          noise,        # "normal" | "heavy"   (on Y0)
                          tau_effect,   # "het" | "het.sparse" (on E[Y1-Y0|X])
                          ite_noise   = c("none","symmetric","skewed"),
                          scale_treat = 1,
                          scale_U=1) {                 # lambda in Y1 = tau + lambda*Y0
  
  ite_noise <- match.arg(ite_noise)
  
  X <- matrix(runif(n * d), n, d)
  W <- rbinom(n, 1, 0.5)
  
  tau <- switch(tau_effect,
                "het"        = tau_function(X),
                "het.sparse" = tau_function(X) + 10 * small_region_smooth(X))
  
  Y0 <- switch(noise, "normal" = rnorm(n), "heavy" = heavy_noise(n))
  
  # random component of ITE, U _||_ (X, Y0)
  U <- switch(ite_noise,
              "none"      = rep(0, n),
              "symmetric" = rt(n, df = 3) / sqrt(3),                 # symmetric
              "skewed"    = rlnorm(n, 0, 1) - exp(0.5))               # mean-zero, right-skewed
  
  Y1 <- tau + scale_U*U + scale_treat * Y0
  Y  <- W * Y1 + (1 - W) * Y0
  
  data.frame(X, W = W, Y = Y, tau = tau)                              # CATE target unchanged
}