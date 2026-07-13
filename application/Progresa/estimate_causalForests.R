############################
# Estimation of CF models for the empirical application based on 
# Progresa (De la O (2013), Ghosh et al. (2026))
############################

rm(list = ls())

path_out <- "application/Progresa/"

# load config file
source("01_config.R")

# ------------------------------------------------
# Real data example: Progresa data
# ------------------------------------------------
set.seed(42)
Progressa = read.csv(paste0(path_out, "data/PROGRESA.csv"))
Y = Progressa$pri2000s
Z = Progressa$treatment
X = Progressa %>% dplyr::select(villages, pri1994, pan1994, prd1994,
                                votos1994, avgpoverty, pobtot1994) 
X <- as.data.frame(X)

#### Estimate Causal Forests ####


# Define params ----
k <- 4 # k-fold cross fitting # 75% training, 25% test 
ntrees <- 2000 # number of trees in forest
varianceMethod <- "BLB"
bag.size <- 10 # bag size (trees per bag)
ntrees.bag <- 15000 # total number of trees for BLB (more than the default ntrees=2000 for point prediction of CATE to enhance precision)
s_ratio <- 0.8 
minsize_blb <- 2

# Sample ids for cross-fitting
D = ncol(X) # number of covariats in X
N = nrow(X) # number of obs. in dataset
idx <- sample(c(1:k), N, replace = T)

# Centering ----
center_dataset_Y <- cbind(Y=Y, X)
center_dataset_W <- cbind(W=Z, X)
Y.hat <- rep(NaN, N)
Y.res <- rep(NaN, N)
W.hat <- rep(NaN, N)
W.res <- rep(NaN, N)

for(i in 1:k){
  
  Y.model <- randomForest::randomForest(Y~., center_dataset_Y[idx != i, ], ntree=ntrees)
  Y.hat[idx == i] <- predict(Y.model, newdata = X[idx == i, ])
  Y.res[idx == i] <- center_dataset_Y$Y[idx == i] - Y.hat[idx == i]
  
  W.model <- SimDesign::quiet(randomForest::randomForest(W~., center_dataset_W[idx != i, ], ntree=ntrees))
  W.hat[idx == i] <- SimDesign::quiet(predict(W.model, newdata = X[idx == i, ]))
  W.res[idx == i] <- center_dataset_W$W[idx == i] - W.hat[idx == i]
  
}

# CF (MSE) ----
split.R = "CT"
cv.opt = "CT"
MSE.estims <- robustCausalTree::robustCausalForest(Y=Y,
                                                   W=Z,
                                                   Y.hat=Y.hat, 
                                                   W.hat=W.hat, 
                                                   Y.res=Y.res, 
                                                   W.res=W.res, 
                                                   X = X,
                                                   k=k, idx=idx, 
                                                   ntrees=ntrees,
                                                   varianceMethod = "BLB",
                                                   s_ratio = s_ratio,
                                                   minsize_blb = minsize_blb,
                                                   nbag = bag.size,
                                                   sample.size.total=floor(nrow(X)/5),
                                                   ntrees.bag = ntrees.bag,
                                                   half_sample_size=2, 
                                                   split.R = split.R, 
                                                   cv.opt = cv.opt)


# CF (MSD) ----
split.R = "MSD"
cv.opt = "MSD"
MSD.estims <- robustCausalTree::robustCausalForest(Y=Y,
                                                   W=Z,
                                                   Y.hat=Y.hat, 
                                                   W.hat=W.hat, 
                                                   Y.res=Y.res, 
                                                   W.res=W.res, 
                                                   X = X,
                                                   k=k, idx=idx, 
                                                   ntrees=ntrees,
                                                   varianceMethod = "BLB",
                                                   s_ratio = s_ratio,
                                                   minsize_blb = minsize_blb,
                                                   nbag = bag.size,
                                                   sample.size.total=floor(nrow(X)/5),
                                                   ntrees.bag = ntrees.bag,
                                                   half_sample_size=2, 
                                                   split.R = split.R, 
                                                   cv.opt = cv.opt)

out <- list(CF_MSE = MSE.estims, 
            CF_MSD = MSD.estims)


saveRDS(out, paste0(path_out, "results/progresa_CF_estimates.rds"))
