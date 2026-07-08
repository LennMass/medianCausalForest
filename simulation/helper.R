# ---- Evaluation helpers ---------------------------------------------------
rmse_eval      <- function(pred, actual) sqrt(mean((pred - actual)^2))
absbias_eval   <- function(pred, actual) abs(pred - actual)
cov_width_eval <- function(est, se, truth, z = 1.96) {
  lo <- est - z * se; hi <- est + z * se
  list(coverage = as.integer(hi >= truth & lo <= truth),
       width    = hi - lo)
}