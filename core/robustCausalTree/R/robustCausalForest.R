#' Estimate Heterogeneous Treatment Effect via Causal Forest with Median Splitting Rules
#'
#' Estimate heterogeneous treatment effect based on random forests.
#' 
#' @param Y outcome variable 
#' @param W treatment variable
#' @param X covariate matrix
#' @param Y.hat cross-fitted estimates of E[Y|X]
#' @param W.hat cross-fitted estimates of E[W|X]
#' @param Y.res outcome residuals
#' @param W.res treatment residuals
#' @param k number of folds for cross-fitting
#' @param idx vector of indices for cross-fitting
#' @param ntrees number of trees in forest
#' @param minsize_estim minimum size of leaf during forest training
#' @param varianceMethod Either "IJ" for infinitesimal Jackknife or "BLB" for Bootstrap of Little Bags. 
#' @param nbag bag size for BLB for variance estimates (just for BLB)
#' @param ntrees.bag total number of trees for BLB (just for BLB)
#' @param s_ratio # subsample size ratio for blb (just for BLB)
#' @param minsize_blb # min leaf size for blb (just for BLB)
#' @param split.R splitting rule (one of "CT", "MSD", "MAD", "LMS")
#' @param cv.opt cross-valdiation option (one of "CT", "MSD", "MAD", "LMS")
#' @param ATE_estimation # should ATE be computed and given out via crossfitted AIPW. Default TRUE.
#' @returns A list with two elements. The first one is \code{pred_cate} which holds
#' the predicted CATE and the corresponding standard error for each unit.
#' The second element is \code{pred_ate} which holds the ATE and its corresponding
#' standard error, calculated based on debiased CATEs through an AIPW estimator. 
#'
#' @rdname robustCausalForest
#' @export
#' @aliases robustCausalForest
#'
robustCausalForest <- function(Y, # outcome variable
                               W, # treatment variable
                               X, # covariate matrix
                               Y.hat = NULL, 
                               W.hat = NULL, 
                               Y.res = NULL, 
                               W.res = NULL,
                               k, # number of folds for cross-fitting
                               idx, # vector of indices for cross-fitting
                               ntrees, # number of trees in a forest
                               varianceMethod="BLB", # either IJ or BLB
                               minsize_estim=5, # minimum number of obs in every leaf for training/estimation
                               nbag, # (only needed for BJB) 
                               ntrees.bag, # (only needed for BJB) 
                               s_ratio, # (only needed for BJB) 
                               minsize_blb, # (only needed for BJB)
                               half_sample_size=2, # (only needed for BJB)
                               split.alpha=0.5,
                               sample.size.total = floor(nrow(X) / 10),
                               cv.alpha=0.5, 
                               split.R = "CT", 
                               cv.opt = "CT",
                               ATE_estimation = TRUE,
                       ...){
  
  # Sample ids for cross-fitting
  D = ncol(X) # number of covariats in X
  N = nrow(X) # number of obs. in dataset
  #idx <- sample(c(1:k), N, replace = T)
  
  # Initialize storing vectors and matrices
  tau.hat <- matrix(rep(NaN, N*ntrees), ncol=ntrees, nrow=N)
  inbag.inf <- vector("list", k)
  tau.var.hat <- rep(NaN, N)

  
  
  # Prepare dataset for causalForest training ----
  
  # residualized or unresidualized outcome
  if(!is.null(Y.res)){
    # use Y-residuals to center out nuisance
    estim_dataset <- as.data.frame(cbind(Y=Y.res, X))
  } else if(is.null(Y.res)) {
    estim_dataset <- as.data.frame(cbind(Y=Y, X))
    warning("Y.res. not specified. Estimation of CATEs proceeds uncentered using Y.")
  } else (
    warning("Specify Y.res. Default is Y.res=NULL and estimation proceeds uncentered.")
  )
  
  # estimate W.hat, Y.hat, if needed (for BLB and ATE via AIPW)
  if(varianceMethod == "BLB" || ATE_estimation == TRUE){
    if(is.null(W.hat) || is.null(Y.hat) ){
      
      Y.hat <- rep(NaN, length(Y))
      W.hat <- rep(NaN, length(W))
      
      
      center_data <- data.frame(W=W, X)
      
      for (i in 1:k){
        W.model <- SimDesign::quiet(randomForest::randomForest(W~., center_data[idx != i, ], ntree=2000))
        W.hat[idx == i] <- SimDesign::quiet(predict(W.model, newdata = X[idx == i, ]))
      }
      
      center_data <- cbind(Y=Y, X)
      
      for (i in 1:k){
        Y.model <- randomForest::randomForest(Y~., center_data[idx != i, ], ntree=2000)
        Y.hat[idx == i] <- predict(Y.model, newdata = X[idx == i, ])
      }
      
      rm(center_data)
      
    }
    
  }
  
  # OOB CATE estimation ----
  
  for(i in 1:k){
    
    
    # Train honest causal forests 
    forest_honest <- robustCausalTree::causalForest(
      as.formula(paste("Y ~", paste(names(estim_dataset[, -1]), collapse = "+"), sep="")),
      data = estim_dataset[idx != i, ], 
      treatment = W[idx != i], 	# W.res not used here as it just takes W \in {0,1} for now and no continuous W. 
      double.Sample = TRUE, # honest trees
      split.Rule = split.R,
      cv.option = cv.opt,
      split.alpha=split.alpha,
      cv.alpha=cv.alpha,
      num.trees = ntrees,
      sample.size.total=sample.size.total,
      nodesize= minsize_estim,# A target for the minimum number of observations in each tree leaf.
      ncolx = D, # total number of covariates
      ncov_sample= min(round(D/2), 5) # Number of covariates randomly sampled to build each tree in the forest. At max 5 variables are considered. 
    )
    
    
    # OOB predictions
    preds <- predict(forest_honest, newdata =  X[idx == i, ], predict.all = TRUE)
    inbag.inf[[i]] <- forest_honest$inbag 
    tau.hat[idx == i, ] <- preds$individual
    
    
    if(varianceMethod=="IJ"){
      
      # Estimate variance via IJ
      tau.var.hat[idx == i] <- ranger:::rInfJack(pred=preds$individual,
                                                inbag=forest_honest$inbag,
                                                calibrate = TRUE)$var.hat
      
    } else if(varianceMethod=="BLB"){
      
      # Estimate variance  via BLB
      blb_se <- blb_cate(Y= estim_dataset$Y[idx != i],
                         W= W[idx != i],
                         W.hat = W.hat[idx != i],
                         X= X[idx != i, ],
                         Xnew = X[idx == i, ],
                         B = ntrees.bag,        # total number of trees
                         ell = nbag,       # little bag size
                         s_ratio = s_ratio, # subsample size ratio 
                         minsize = minsize_blb,  # min leaf size
                         half_sample_size=half_sample_size,
                         sample.size.total=sample.size.total,
                         split.Rule = split.R,
                         cv.option = cv.opt)
      
      tau.var.hat[idx == i] <- blb_se$tau_se
      
    } else {
      stop("varianceMethod either 'IJ' or 'BLB")
    }
    
    
    
    
    print(paste("finished run ", i, " of k=", k, "OOB folds."))
    print(Sys.time())
    
  }
  
  
  
  # prepare savings ----
  cate_est <- data.frame(estimate=rowMeans(tau.hat),
                         std.err=sqrt(tau.var.hat)) %>%
    mutate(id = row_number())
  
  
  ### ATE estimation ----
  
  
  if(ATE_estimation==TRUE){
    
    # Create forest-object to estimate doubly-robust ATE from CATEs
    # Preparation for grf:::average_treatment_effect()
    forest_obj <- list()
    forest_obj$Y.orig <- as.vector(Y)
    forest_obj$W.orig <- as.vector(W)
    forest_obj$W.hat <- W.hat
    forest_obj$Y.hat <- Y.hat
    forest_obj$tau.hat.pointwise <- cate_est$estimate
    forest_obj$tau.var.hat.pointwise <- cate_est$std.err
    #class(forest_obj) <- "causal_forest"
    
    ate_est <- avg_treat_eff(forest=forest_obj,
                             target.sample = "all", 
                             method = "AIPW")
    
  } else if(ATE_estimation==FALSE){
    
    
    ate_est <- NULL
    
  } else {
    warning("ATE_estimation needs to be TRUE/FALSE.")
  }
  
  
  
  
  
  return(list(pred_cate = cate_est, pred_ate=ate_est))
}


# Average Treatment Effect function
# Simplified version of grf::average_treatment_effect

avg_treat_eff <- function(forest,
                          target.sample = "all",
                          method = "AIPW",
                          subset = NULL,
                          debiasing.weights = NULL,
                          compliance.score = NULL,
                          num.trees.for.weights = 500) {


  clusters <- 1:NROW(forest$Y.orig)

  observation.weight <- obs_weights(forest)

  subset <- valid_subset(forest, subset)
  subset.clusters <- clusters[subset]
  subset.weights <- observation.weight[subset]

  if (length(unique(subset.clusters)) <= 1) {
    stop("The specified subset must contain units from more than one cluster.")
  }

  if (!is.null(debiasing.weights)) {
    if (length(debiasing.weights) == NROW(forest$Y.orig)) {
      debiasing.weights <- debiasing.weights[subset]
    } else if (length(debiasing.weights) != length(subset)) {
      stop("If specified, debiasing.weights must be a vector of length n or the subset length.")
    }
  }

  if (method == "AIPW" && target.sample == "all") {
    # This is the most general workflow, that shares codepaths with best linear projection
    # and other average effect estimators.

    .sigma2.hat <- function(DR.scores, tau.hat) {
      correction.clust <- Matrix::sparse.model.matrix(~ factor(subset.clusters) + 0, transpose = TRUE) %*%
        (sweep(as.matrix(DR.scores), 2, tau.hat, "-") * subset.weights)
      n.adj <- sum(rowsum(subset.weights, subset.clusters) > 0) # effective number of samples(clusters) is all with > 0 weight.

      Matrix::colSums(correction.clust^2) / sum(subset.weights)^2 *
        n.adj / (n.adj - 1)
    }


    DR.scores <- compute_DR_scores(forest,
                               subset = subset,
                               debiasing.weights = debiasing.weights,
                               compliance.score = compliance.score,
                               num.trees.for.weights = num.trees.for.weights)


    tau.hat <- weighted.mean(DR.scores, subset.weights)
    sigma2.hat <- .sigma2.hat(DR.scores, tau.hat)
    return(c(estimate = tau.hat, std.err = sqrt(sigma2.hat)))
  }
}




obs_weights <- function(forest) {
  # Case: No given sample.weights

  raw.weights <- rep(1, NROW(forest$Y.orig))

  return (raw.weights / sum(raw.weights))
}

valid_subset <- function(forest, subset) {
  if (is.null(subset)) {
    subset <- 1:NROW(forest$Y.orig)
  }
  if (is.logical(subset) && length(subset) == NROW(forest$Y.orig)) {
    subset <- which(subset)
  }
  if (!all(subset %in% 1:NROW(forest$Y.orig))) {
    stop(paste(
      "If specified, subset must be a vector contained in 1:n,",
      "or a boolean vector of length n."
    ))
  }
  subset
}


compute_DR_scores <- function(forest,
                          subset = NULL,
                          debiasing.weights = NULL,
                          num.trees.for.weights = 500,
                          ...) {
  subset <- valid_subset(forest, subset)
  W.orig <- forest$W.orig[subset]
  W.hat <- forest$W.hat[subset]
  Y.orig <- forest$Y.orig[subset]
  Y.hat <- forest$Y.hat[subset]
  tau.hat.pointwise <- forest$tau.hat.pointwise

  binary.W <- all(forest$W.orig %in% c(0, 1))

  if (is.null(debiasing.weights)) {
    if (binary.W) {
      debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
    }# else {
      # Start by learning debiasing weights if needed.
      # The goal is to estimate the variance of W given X. For binary treatments,
      # we get a good implicit estimator V.hat = e.hat (1 - e.hat), and
      # so this step is not needed. Note that if we use the present CAPE estimator
      # with a binary treatment and set V.hat = e.hat (1 - e.hat), then we recover
      # exactly the AIPW estimator of the CATE.
    #   clusters <- if (length(forest$clusters) > 0) {
    #     forest$clusters
    #   } else {
    #     1:length(forest$Y.orig)
    #   }
    #   variance_forest <- grf::regression_forest(forest$X.orig,
    #                                             (forest$W.orig - forest$W.hat)^2,
    #                                             clusters = clusters,
    #                                             sample.weights = forest$sample.weights,
    #                                             num.trees = num.trees.for.weights,
    #                                             ci.group.size = 1,
    #                                             seed = forest$seed,
    #                                             num.threads = forest$num.threads)
    #   V.hat <- predict(variance_forest)$predictions
    #   debiasing.weights.all <- (forest$W.orig - forest$W.hat) / V.hat
    #   debiasing.weights <- debiasing.weights.all[subset]
    # }
  } else if (length(debiasing.weights) == length(forest$Y.orig)) {
    debiasing.weights <- debiasing.weights[subset]
  } else if (length(debiasing.weights) != length(subset))  {
    stop("If specified, debiasing.weights must have length n or |subset|.")
  }

  # Form AIPW scores. Note: We are implicitly using the following
  # estimates for the regression surfaces E[Y|X, W=0/1]:
  # Y.hat.0 <- Y.hat - W.hat * tau.hat.pointwise
  # Y.hat.1 <- Y.hat + (1 - W.hat) * tau.hat.pointwise
  Y.residual <- Y.orig - (Y.hat + tau.hat.pointwise * (W.orig - W.hat))

  tau.hat.pointwise + debiasing.weights * Y.residual
}

# BLB-style CATE estimation using honest causal trees ----
# Based on Athey & Wager (2016, 2018) and GRF (2019)


blb_cate <- function(
    Y,
    W,
    W.hat,
    X,
    Xnew = NULL,
    B = 5000,       # total number of trees
    ell = 10,       # little bag size
    s_ratio = 0.8,
    minsize = 2,
    half_sample_size=2, # 2: proper half-sample; 1: full sample; 1.5: 3/4 sample
    sample.size.total=sample.size.total,
    split.Rule = split.R,
    cv.option = cv.opt
) {
  
  if (is.null(Xnew)) Xnew <- X
  
  if (half_sample_size > 2 || half_sample_size < 1) stop("half_sample_size should be between 1 (full sample) and 2 (proper half sample).")
  
  n      <- nrow(X)
  d      <- ncol(X)
  n_new  <- nrow(Xnew)
  n_half <- floor(n / half_sample_size)
  
  # Number of little bags (groups)
  G <- floor(B / ell)
  if (G < 2) stop("Need at least two little bags")
  
  # Storage:
  # For each group g, we store tree-level predictions at Xnew
  # Dimension: n_new × ell
  group_preds <- vector("list", G)
  
  # Step 1: grow trees inside half-sample little bags ----
  
  for (g in seq_len(G)) {
    
    # --- draw half-sample Hg (theoretical H_g) ---
    Hg <- sample.int(n, n_half, replace = FALSE)
    
    Yg <- Y[Hg]
    Wg <- W[Hg]
    Xg <- X[Hg, , drop = FALSE]
    
    # storage for tree-level predictions
    preds_g <- matrix(NA, n_new, ell)
    
    for (b in seq_len(ell)) {
      
      # --- draw subsample I_b \in H_g ---
      # smaller subsample from half-sample
      s <- floor(n_half*s_ratio)
      idx <- sample.int(n_half, s, replace = FALSE)
      
      # treatment balance check (avoid C-level crashes)
      if (length(unique(Wg[idx])) < 2){
        stop("Treatment-control unbalanced.")
      }
      
      train_df <- as.data.frame(cbind(Y=Yg[idx], Xg[idx, , drop = FALSE]))
      
      # train one tree on subsample s drawn from half-sample H_g
      forest_half <- robustCausalTree::causalForest(
        as.formula(paste("Y  ~", paste(names(train_df[, -1]), collapse = "+"), sep="")),
        data = train_df, 
        treatment = Wg[idx], 	# W.res not used here as it just takes W \in {0,1} for now and no continuous W. 
        double.Sample = TRUE, # honest trees
        split.Rule = split.Rule,
        cv.option = cv.option,
        num.trees = 1,
        sample.size.total=sample.size.total,
        nodesize= minsize,# A target for the minimum number of observations in each tree leaf.
        ncolx = d, # total number of covariates
        ncov_sample= min(round(d/2), 5) # Number of covariates randomly sampled to build each tree in the forest. At max 5 variables are considered to speed up training. 
      )
      
      
      pr <- predict(forest_half, newdata = Xnew, predict.all=TRUE)
      
      
      preds_g[, b] <- pr$individual
    }
    
    group_preds[[g]] <- preds_g
  }
  
  # Step 2: ANOVA decomposition (Section 4.1 in GRF paper) ----

  # Group means:  
  group_means <- sapply(
    group_preds,
    function(M) rowMeans(M, na.rm = TRUE)
  )
  # dimension: n_new × G
  
  # Overall mean: 
  overall_mean <- rowMeans(group_means, na.rm = TRUE)
  
  # ---- Between-group variance (half-sampling variance) ----
  between_var <- rowMeans(
    (group_means - overall_mean)^2,
    na.rm = TRUE
  )
  
  # ---- Within-group variance (Monte Carlo noise) ----
  within_var <- rowMeans(
    sapply(group_preds, function(M) {
      rowMeans((M - rowMeans(M, na.rm = TRUE))^2, na.rm = TRUE)
    }),
    na.rm = TRUE
  )
  
  # ---- H-BLB variance estimator (ANOVA formula)  ----
  H_hat_raw <- between_var - within_var / (ell - 1) # mean should be close to zero, otherwise increase B
  # Bayesian ANOVA shrinkage (Bayesian ANOVA formula)
  # shrinks more when G is small
  # converges to the classical estimator as B approaches infinity
  H_hat <- (G - 1) / G * (between_var - within_var / ell)
  H_hat <- pmax(H_hat, 0)   # Truncate at zero (posterior support)
  
  if (mean(H_hat_raw < 0) > 0.05) {
    warning(paste(mean(H_hat_raw < 0), " of H_hat values truncated at zero — consider increasing B."))
  }
  
  # Step 3: V-BLB estimator: pseudo-outcome for curvature ----
  
  # pseudo-outcome for curvature
  curv_y <- (W - W.hat)^2
  
  # honest regression forest
  rf_V <- grf::regression_forest(X, curv_y, honesty = TRUE)
  
  V_hat <- predict(rf_V, Xnew)$predictions
  V_inv <- 1/V_hat
  
  # Step 4: Final variance estimate (see Eq (16) in GRF paper) ----
  
  Sigma_hat <- V_inv^2 * H_hat
  se_hat    <- sqrt(Sigma_hat)
  
  
  list(
    tau_se         = se_hat,
    H_hat          = H_hat,
    between_var    = between_var,
    within_var     = within_var,
    V_hat          = V_hat,
    mean_H_hat_raw = mean(H_hat_raw)
  )
}





