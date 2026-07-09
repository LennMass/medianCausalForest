# ---- Classical ATE estimators (estimate + SE) -----------------------------


# Rosenbaum rank-adjusted: same root-finding recipe as the application file.
rosen_root <- function(y, w, xdf) {
  df <- cbind.data.frame(y = y, w = w, xdf)
  N  <- nrow(df); m <- sum(w == 1)
  aux <- function(tau) {
    df$ys <- df$y - tau * df$w
    r <- resid(lm(ys ~ . - y - w, data = df))
    sum(rank(r)[w == 1]) - m * (N + 1) / 2
  }
  uniroot(aux, c(-10, 10), extendInt = "yes")$root
}

# Returns a data.frame: method | estimate | std.err.
# Closed-form SEs where natural; a small nonparametric bootstrap otherwise.
classical_ATE <- function(Y, W, X, R_boot = 99L) {
  
  y1 <- Y[W == 1]; y0 <- Y[W == 0]
  n1 <- length(y1); n0 <- length(y0); n <- n1 + n0
  Xdf <- as.data.frame(X)
  
  # ----- point estimates -----
  est_DM      <- mean(y1) - mean(y0)
  
  Xc      <- scale(Xdf, center = TRUE, scale = FALSE)
  fit_lin <- lm(Y ~ W * ., data = cbind.data.frame(Y = Y, W = W, Xc))
  est_Lin <- unname(coef(fit_lin)["W"])
  
  est_Radj <- rosen_root(Y, W, Xdf)
  
  eif <- parTreat::eif_additive(y0, y1)
  waq <- parTreat::waq(y0, y1)
  est_EIF <- as.numeric(eif$tau); est_WAQ <- as.numeric(waq$tau)
  
  # ----- analytic SEs -----
  se_DM  <- sqrt(var(y1)/n1 + var(y0)/n0)
  se_Lin  <- sqrt(sandwich::vcovHC(fit_lin, type = "HC2")["W", "W"])
  wt     <- suppressWarnings(wilcox.test(y1, y0, conf.int = TRUE, exact = FALSE))
  se_EIF <- tryCatch(as.numeric(eif$se), error = function(e) NA_real_)
  se_WAQ <- tryCatch(as.numeric(waq$se), error = function(e) NA_real_)
  
  # ----- bootstrap SEs (Rosenbaum) -----
  B <- matrix(NA_real_, R_boot, 1, dimnames = list(NULL, c("Radj")))
  for (b in seq_len(R_boot)) {
    idx <- sample.int(n, n, replace = TRUE)
    Yb  <- Y[idx]; Wb <- W[idx]; Xb <- Xdf[idx, , drop = FALSE]
    y1b <- Yb[Wb == 1]; y0b <- Yb[Wb == 0]
    if (length(y1b) < 2 || length(y0b) < 2) next
    B[b, "Radj"]    <- tryCatch(rosen_root(Yb, Wb, Xb),
                                error = function(e) NA_real_)
  }
  se_Radj    <- sd(B[, "Radj"],    na.rm = TRUE)
  
  data.frame(
    method   = c("DM", "Lin", "RosenbaumAdj", "EIF", "WAQ"),
    estimate = c(est_DM, est_Lin, est_Radj, est_EIF, est_WAQ),
    std.err  = c(se_DM, se_Lin, se_Radj, se_EIF, se_WAQ),
    stringsAsFactors = FALSE
  )
}