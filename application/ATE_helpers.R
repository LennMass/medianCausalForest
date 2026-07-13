############################
# ATE competitor functions based on Ghosh et al. (2026)
############################


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
