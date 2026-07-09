###############################################################################
# Main simulation runner
# Compares CT vs MSD on CATE, and forest-based ATE vs classical ATE estimators
# across the 4 scenarios (S1 to S4) of the paper.
###############################################################################

source(here::here("01_config.R"))


# ---- Scenario grid ----------------------------------------------------------

# Adjust to the setting of interest
scenarios <- data.frame(
  id          = paste0("S", 1),                # scenario name (S1, S2, ...)
  N           = c(1000),                       # 1000, 2000 
  D           = c(5),                          # 5, 10, 20
  noise       = c("normal"),                   # normal, heavy
  tau_effect  = c("het"),                      # het, het.sparse 
  ite_noise   = c("none"),                     # none, skewed   
  scale_treat = c(1),                          # 1 
  scale_U = c(1),                              # 1, 2
  description = c(""),                         # optional description 
  stringsAsFactors = FALSE
)


# ---- Global sim settings ----------------------------------------------------
sims           <- 1                            # replications per seed
setseed_vec    <- c(6)                         # 6, 15609, 100, 5020, 42
k              <- 4                            # cross-fitting folds
ntrees         <- 2000                         # number of trees in forest
varianceMethod <- "BLB"                        # bootsrap of little bags
bag.size       <- 10                           # bag size (trees per bag) for BLB
ntrees.bag     <- 4000                         # total number of trees for BLB variance forest (more than ntrees to enhance coverage)
s_ratio        <- 0.8                          # sample ratio for BLB
minsize_blb    <- 2                            # minimum leaf size for BLB
hss            <- 2                            # half sample size for BLB: 2 is proper half sample (1 would be full sample)
methods_forest <- c("CT", "MSD", "MAD", "LMS") # CF splitting rules

path_out <- "simulation/raw_results/"
dir.create(here::here(path_out), recursive = TRUE, showWarnings = FALSE)


# ---- Main loop --------------------------------------------------------------
for (s_i in seq_len(nrow(scenarios))) {
  
  sc <- scenarios[s_i, ]
  cat(sprintf("\n=== %s  N=%d  D=%d  noise=%s  tau=%s  itenoise=%s ===\n",
              sc$id, sc$N, sc$D, sc$noise, sc$tau_effect, sc$ite_noise))
  
  for (seed in setseed_vec) {
    
    set.seed(seed)
    rows <- vector("list", sims)
    
    for (sim_i in seq_len(sims)) {
      
      # ---- 1. data ----
      d <- simulate_data(sc$N, sc$D, sc$noise, sc$tau_effect, sc$ite_noise, sc$scale_treat, sc$scale_U)
      X <- d[, !(names(d) %in% c("Y", "W", "tau"))]
      truth_ATE <- mean(d$tau)
      
      
      # ---- 2. cross-fit nuisances (Y.hat, W.hat) ----
      idx <- sample.int(k, sc$N, replace = TRUE)
      Y.hat <- numeric(sc$N)
      W.hat <- numeric(sc$N)
      for (i in seq_len(k)) {
        
        train <- idx != i
        test  <- idx == i
        
        ## --- Outcome model E[Y|X] ---
        fitY <- ranger::ranger(
          y             = d$Y[train],
          x             = X[train, , drop = FALSE],
          num.trees     = 500,
          min.node.size = 5,
          num.threads   = 1,            # avoid nested parallelism with outer loop
          seed          = 1000L + i
        )
        Y.hat[test] <- predict(fitY, data = X[test, , drop = FALSE])$predictions
        
        ## --- Propensity model P(W=1|X) ---
        fitW <- ranger::ranger(
          y             = factor(d$W[train]),     # factor -> classification mode
          x             = X[train, , drop = FALSE],
          num.trees     = 500,
          min.node.size = 10,
          probability   = TRUE,                   # return P(W=1|X)
          num.threads   = 1,
          seed          = 2000L + i
        )
        W.hat[test] <- predict(fitW, data = X[test, , drop = FALSE])$predictions[, "1"]
      }
      
      W.hat <- pmin(pmax(W.hat, 1e-3), 1 - 1e-3)
      Y.res <- d$Y - Y.hat
      W.res <- d$W - W.hat
      
      # ---- 3. forests (CT and MSD) ----
      forest_rows <- list()
      for (rule in methods_forest) {
        
        # divisor for sample.size.total; scales with sample size N
        sst <- ifelse(length(d$Y) >= 2000, 10, 5) 
        
        est <- robustCausalTree::robustCausalForest(
          Y = d$Y, W = d$W,
          Y.hat = Y.hat, W.hat = W.hat,
          Y.res = Y.res, W.res = W.res,
          X = X, 
          k=k, idx=idx, 
          ntrees=ntrees,
          varianceMethod = varianceMethod,
          s_ratio = s_ratio,
          minsize_blb = minsize_blb,
          nbag = bag.size,
          sample.size.total=floor(nrow(X)/sst),
          ntrees.bag = ntrees.bag,
          half_sample_size=hss, 
          cv.opt = rule, 
          split.R = rule
        )
        
        cate_cv <- cov_width_eval(est$pred_cate$estimate,
                                  est$pred_cate$std.err, d$tau)
        ate_cv  <- cov_width_eval(est$pred_ate["estimate"],
                                  est$pred_ate["std.err"], truth_ATE)
        
        forest_rows[[length(forest_rows) + 1]] <- data.frame(
          method     = rule, estimand = "CATE",
          rmse       = rmse_eval(est$pred_cate$estimate, d$tau),
          absbias    = mean(abs(est$pred_cate$estimate - d$tau)),
          coverage   = mean(cate_cv$coverage),
          width      = mean(cate_cv$width),
          signed_err = NA_real_,
          stringsAsFactors = FALSE
        )
        forest_rows[[length(forest_rows) + 1]] <- data.frame(
          method     = rule, estimand = "ATE",
          rmse       = NA_real_,
          absbias    = abs(est$pred_ate["estimate"] - truth_ATE),
          coverage   = ate_cv$coverage,
          width      = ate_cv$width,
          signed_err = est$pred_ate["estimate"] - truth_ATE,
          stringsAsFactors = FALSE
        )
      }
      
      # ---- 4. classical ATE estimators ----
      cat_ate <- classical_ATE(d$Y, d$W, X, R_boot = 99L)
      lo <- cat_ate$estimate - 1.96 * cat_ate$std.err
      hi <- cat_ate$estimate + 1.96 * cat_ate$std.err
      cat_rows <- data.frame(
        method     = cat_ate$method, estimand = "ATE",
        rmse       = NA_real_,
        absbias    = abs(cat_ate$estimate - truth_ATE),
        coverage   = as.integer(lo <= truth_ATE & hi >= truth_ATE),
        width      = hi - lo,
        signed_err = cat_ate$estimate - truth_ATE,
        stringsAsFactors = FALSE
      )
      
      # ---- 5. assemble this replication ----
      one_sim <- rbind(do.call(rbind, forest_rows), cat_rows)
      one_sim$sim <- sim_i
      rows[[sim_i]] <- one_sim
      
      if (sim_i %% 2 == 0)
        cat(sprintf("  seed %d  sim %d/%d   %s\n",
                    seed, sim_i, sims, format(Sys.time(), "%H:%M:%S")))
    }
    
    out <- do.call(rbind, rows)
    out$scenario <- sc$id
    out$seed     <- seed
    
    fn <- sprintf("%sresults_%s_seed%d.rds", path_out, sc$id, seed)
    saveRDS(out, here::here(fn))
  }
}