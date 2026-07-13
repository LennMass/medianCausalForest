############################
# Estimation of CF models for the empirical application based on 
# ACTG175 (Hammer et al. (1996), Leqi and Kennedy (2022))
############################

rm(list = ls())

path_out <- "application/ACTG175/"

# load config file
source("01_config.R")


set.seed(12345) 

# Load the ACTG 175 dataset from speff2trial package
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

A <- dat$A
Y <- dat$cd496



#### Estimate Causal Forest ####

# Define params ----

k <- 4 # k-fold cross fitting # 75% training, 25% test 
ntrees <- 2000 # number of trees in forest
varianceMethod <- "BLB"
bag.size <- 10 # bag size (trees per bag)
ntrees.bag <- 5000 # total number of trees for BLB
s_ratio <- 0.8 

minsize_blb <- 2

# Sample ids for cross-fitting
D = ncol(X) # number of covariats in X
N = nrow(X) # number of obs. in dataset
idx <- sample(c(1:k), N, replace = T)

# Centering ----

center_dataset_Y <- cbind(Y=Y, X)
center_dataset_W <- cbind(W=A, X)
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

# CT (MSE) ----
split.R = "CT"
cv.opt = "CT"
MSE.estims <- robustCausalTree::robustCausalForest(Y=Y,
                                                   W=A,
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


# MSD ----
split.R = "MSD"
cv.opt = "MSD"
MSD.estims <- robustCausalTree::robustCausalForest(Y=Y,
                                                   W=A,
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


saveRDS(out, here::here(paste0(path_out, "results/ACTG175_CF_estimates.rds")))



