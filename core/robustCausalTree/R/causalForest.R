#' Causal Effect Regression and Estimation Forests (Tree Ensembles)
#'
#' Build a random causal forest by fitting a user selected number of
#' \code{causalTree} models to get an ensemble of \code{rpart} objects.
#'
#' @param formula a \link{formula}, with a response and features but no
#' interaction terms.  If this a a data frome, that is taken as the model frame
#' (see \code{\link{model.frame}).}
#'
#' @param data an optional data frame that includes the variables
#'   named in the formula.
#'
#' @param weights optional case weights.
#'
#' @param treatment a vector that indicates the treatment status of
#' each observation. 1 represents treated and 0 represents control.
#' Only binary treatment supported in this version.
#'
#' @param na.action the default action deletes all observations for which
#'   \code{y} is missing, but keeps those in which one or more predictors
#'   are missing.
#'
#'
#' @param split.Rule causalTree splitting options, one of \code{"TOT"},
#' \code{"CT"}, \code{"fit"}, \code{"tstats"}, four splitting rules in
#' \code{causalTree}.  Note that the \code{"tstats"} alternative does
#' not have an associated cross-validation method \code{cv.option};
#' see Athey and Imbens (2016)
#'   for a discussion.  Note further that \code{split.Rule} and
#' \code{cv.option} can mix and match.
#'
#' @param double.Sample boolean option, \code{TRUE} or \code{FALSE},
#' if set to True, causalForest will build honest trees.
#'
#' @param split.Honest boolean option, \code{TRUE} or \code{FALSE}, used
#' to decide the splitting rule of the trees.
#'
#' @param split.Bucket boolean option, \code{TRUE} or \code{FALSE},
#' used to specify whether to apply the discrete method in splitting the tree.
#' If set as \code{TRUE}, in splitting a node, the observations in a leaf
#' will be be partitioned into buckets, with each bucket containing
#' \code{bucketNum} treated and \code{bucketNum} control units, and where
#' observations are ordered prior to partitioning. Splitting will take
#' place by bucket.
#'
#' @param bucketNum number of observations in each bucket when set
#' \code{split.Bucket} = \code{TRUE}.  However, the code will override
#' this choice in order to guarantee that there are at least \code{minsize}
#' and at most \code{bucketMax} buckets.
#'
#' @param bucketMax Option to choose maximum number of buckets to use in
#' splitting when set \code{split.Bucket} = \code{TRUE}, \code{bucketNum}
#' can change by choice of \code{bucketMax}.
#'
#' @param cv.option cross validation options, one of \code{"TOT"},
#' \code{"matching"}, \code{"CT"}, \code{"fit"}, four cross validation
#' methods in \pkg{causalTree}.  There is no \code{cv.option} for
#' the \code{split.Rule} \code{"tstats"}; see Athey and Imbens (2016)
#' for discussion.
#'
#' @param cv.Honest boolean option, \code{TRUE} or \code{FALSE}, only
#' used for \code{cv.option} as \code{"CT"} or \code{"fit"}, to specify
#' whether to apply honest risk evalation function in cross validation.
#' If set \code{TRUE}, use honest risk function, otherwise use adaptive
#' risk function in cross validation.  If set \code{FALSE}, the user
#' choice of \code{cv.alpha} will be set to 1.  If set
#' \code{TRUE}, \code{cv.alpha}
#' will default to 0.5, but the user choice of \code{cv.alpha} will be
#' respected.  Note that honest cv estimates within-leaf variances and
#' may perform better with larger leaf sizes and/or small number of
#' cross-validation sets.
#'
#' @param minsize in order to split, each leaf must have at least
#' \code{minsize} treated cases and \code{minsize} control cases.
#' The default value is set as 2.
#'
#' @param propensity propensity score used in \code{"TOT"} splitting
#' and \code{"TOT"}, honest \code{"CT"} cross validation methods.
#' The default value is the proportion of treated cases in all observations.
#' In this implementation, the propensity score is a constant for the whole
#'   dataset.  Unit-specific propensity scores are not supported; however,
#' the user may use inverse propensity scores as case weights if desired.
#'
#' @param control a list of options that control details of the
#'   \code{rpart} algorithm.  See \code{rpart.control}.
#'
#' @param split.alpha scale parameter between 0 and 1, used in splitting
#' risk evaluation function for \code{"CT"}. When \code{split.Honest = FALSE},
#' \code{split.alpha} will be set as 1.  For \code{split.Rule}=\code{"tstats"},
#' if \code{split.Honest}=\code{TRUE}, \code{split.alpha} is used in
#' calculating the risk function, which determines the order of
#' pruning in cross-validation.
#'
#' @param cv.alpha scale paramter between 0 and 1, used in cross validation
#' risk evaluation function for \code{"CT"} and \code{"fit"}.  When
#'   \code{cv.Honest = FALSE}, \code{cv.alpha} will be set as 1.
#'
#' @param cost a vector of non-negative costs, one for each variable in
#'   the model. Defaults to one for all variables. These are scalings to
#'   be applied when considering splits, so the improvement on splitting
#'   on a variable is divided by its cost in deciding which split to
#'   choose.
#'
#' @param sample.size.total Sample size used to build each tree in the
#' forest (sampled randomly with replacement).
#'
#' @param sample.size.train.frac Fraction of the sample size used for
#' building each tree (training). For eexample,  if the sample.size.total is
#' 1000 and frac =0.5 then, 500 samples will be used to build the tree and
#' the other 500 samples will be used the evaluate the tree.
#'
#' @param mtry Number of data features used to build a tree
#' (This variable is not used presently).
#'
#' @param nodesize Minimum number of observations for treated and
#' control cases in one leaf node
#'
#' @param num.trees Number of trees to be built in the causal forest
#'
#' @param ncolx Total number of covariates
#'
#' @param ncov_sample Number of covariates randomly sampled to
#' build each tree in the forest
#'
#' @param \dots arguments to \code{rpart.control} may also be
#' specified in the call to \code{causalForest}.  They are
#' checked against the
#' list of valid arguments.
#' The parameter \code{minsize} is implemented differently in
#' \code{causalTree} than in \code{rpart}; we require a minimum of \code{minsize}
#' treated observations and a minimum of \code{minsize} control
#' observations in each leaf.
#'
#' @details
#' CausalForest builds an ensemble of CausalTrees (See Athey and Imbens,
#' \emph{Recursive Partitioning for Heterogeneous Causal
#' Effects} (2016)), by repeated random sampling of the data with replacement.
#' Further, each tree is built using a randomly sampled subset of all available
#' covariates. A causal forest object is a list of trees. To predict, call R's
#' predict function with new test data and the causalForest object (estimated
#' on the training data) obtained after calling the causalForest function.
#' During the prediction phase, the average value over all tree predictions
#' is returned as the final prediction by default.
#' To return the predictions of each tree in the forest for each test
#' observation, set the flag \code{predict.all=TRUE}
#' CausalTree differs from \code{rpart} function from \pkg{rpart} package in
#' splitting rules and cross validation methods. Please check Athey
#' and Imbens, \emph{Recursive Partitioning for Heterogeneous Causal
#' Effects} (2016) and Stefan Wager and Susan Athey, \emph{Estimation and
#' Inference of Heterogeneous Treatment Effects using Random Forests
#' } for more details.
#'
#' @returns An object of class \code{rpart}.  See \code{rpart.object}.
#' @references
#' Breiman L., Friedman J. H., Olshen R. A., and Stone, C. J. (1984)
#' \emph{Classification and Regression Trees.}
#' Wadsworth.
#'
#' Athey, S and G Imbens (2016)  \emph{Recursive Partitioning for
#' Heterogeneous Causal Effects}.  http://arxiv.org/abs/1504.01132
#'
#' Wager,S and Athey, S (2015) \emph{Estimation and Inference of Heterogeneous
#' Treatment Effects using Random Forests}
#' http://arxiv.org/abs/1510.04342
#'
#' @seealso
#' \code{\link{causalTree}}
#' \code{\link{honest.causalTree}},
#' \code{rpart.control}, \code{rpart.object},
#' \code{summary.rpart}, \code{rpart.plot}
#'
#' @rdname causalForest
#' @export
#' @aliases init.causalForest
init.causalForest <- function(formula, data, treatment, weights=FALSE,
  cost=FALSE, num.trees,ncov_sample) {
  num.obs <- nrow(data)
  trees <- vector("list", num.trees)
  inbag <- matrix(0, num.obs, num.trees)
  cov_sample <- matrix(0,num.trees,ncov_sample)
  inbag.Est <- matrix(0, num.obs, num.trees)
  nameall_sample <- matrix(0,num.trees,ncov_sample+2) #2 end cols for y,w,no tau_true
  fsample<-vector("list",num.trees)
  causalForestobj <- list(trees = trees, formula=formula, data=data, treatment=treatment, weights=weights, cost=cost, ntree = num.trees, inbag = inbag,cov_sample=cov_sample, fsample=fsample,nameall_sample=nameall_sample,inbag.Est=inbag.Est)
  class(causalForestobj) <- "causalForest"
  return(causalForestobj)
}

#' @rdname causalForest
#' @export
#' @aliases predict.causalForest
#' @param object a \code{causalTree} object
#' @param newdata new data to predict
#' @param predict.all If TRUE, return predicted individual effect for
#' each observations. Otherwise, return the average effect.
#' @param type the type of returned object
#'
predict.causalForest <- function(object,newdata, predict.all = FALSE, type="vector",...) {
  if (!inherits(object, "causalForest")) stop("Not a legitimate \"causalForest\" object")
  individual <- sapply(object$trees, function(tree.fit) {
    predict(tree.fit, newdata=newdata, type="vector")
  })

  #replace sapply with a loop if needed
  # print(dim(individual))
  aggregate <- rowMeans(individual)
  if (predict.all) {
    list(aggregate = aggregate, individual = individual)
  } else {
    aggregate
  }
}

#' @rdname causalForest
#' @export
#' @aliases causalForest
#'
causalForest <- function(formula, data, treatment,
                         na.action = na.causalTree,
                         split.Rule="CT", double.Sample =TRUE, split.Honest=TRUE,
                         split.Bucket=FALSE, bucketNum = 5,
                         bucketMax = 100, cv.option="CT", cv.Honest=TRUE, minsize = 2L,
                         propensity, control, split.alpha = 0.5, cv.alpha = 0.5,
                         sample.size.total = floor(nrow(data) / 10), sample.size.train.frac = .5,
                         mtry = ceiling(ncol(data)/3), nodesize = 1, num.trees=nrow(data),
                         cost=FALSE, weights=FALSE,ncolx,ncov_sample) {
  # data.table<-..x<-..treeRange<-..estimRange<-NULL
  # do not implement subset option of causalTree, that is inherited from rpart but have not implemented it here yet
  vars <- all.vars(formula)
  y <- vars[1]
  x <- vars[2:length(vars)]
  treatmentdf <- data.frame(treatment)
  # if(class(data)[1]=="data.table"){
  if(inherits(data,"data.table",TRUE)==1){
    treatmentdt <- data.table(treatment)
    datax<-data[,..x]
    datay<-data[,y,with=FALSE]
    data <- cbind(datax,datay, treatmentdt)
  # }else if(class(data)=="data.frame"){
  }else if(inherits(data,"data.frame")){
    data <- data[, c(x, y)]
    data <- cbind(data, treatmentdf)
  }

  num.obs <-nrow(data)

  causalForest.obj <- init.causalForest(formula=formula, data=data, treatment=treatment, weights=weights, cost=cost, num.trees=num.trees,ncov_sample=ncov_sample)

  sample.size <- min(sample.size.total, num.obs)
  if (double.Sample) {
    train.size <- round(sample.size.train.frac*sample.size)
    est.size <- sample.size - train.size
  }

  # print("Building trees ...")

  for (tree.index in 1:num.trees) {

    # print(paste("Tree", as.character(tree.index)))

    full.idx <- sample.int(num.obs, sample.size, replace = FALSE)

    if(double.Sample) {
      train.idx <- full.idx[1:train.size]
      reestimation.idx <- full.idx[(train.size+1):sample.size]
    }

    #randomize over the covariates for splitting (both train and reestimation)
    cov_sample<-sample.int(ncolx)
    cov_sample<-cov_sample[1:ncov_sample]

    #modify the y=f(x) equation accordingly for this tree
    #and modify the colnames
    fsample<-""
    nextx<-""
    nameall_sample<-c()
    for (ii in 1:(ncov_sample)){
      nextxindex <- cov_sample[ii]
      nextx <- x[[nextxindex]]
      if (ii==1) {
        fsample <-nextx
        name<- nextx
      }
      if (ii>1) {
        fsample <- paste0(fsample,"+",nextx)
        name <- c(name, nextx)
      }
    }
    nameall_sample <- c( name,y, "w") #, "tau_true")

    #store this var subset for each tree (need it during testing/predict stage)
    causalForest.obj$cov_sample[tree.index,]<-cov_sample
    #also store the formula & colnames of X for each tree (need it during testing/predict stage)
    causalForest.obj$nameall_sample[tree.index,]<-nameall_sample
    causalForest.obj$fsample[[tree.index]]<-fsample

    # if(class(data)[1]=="data.table"){
    if(inherits(data,"data.table",TRUE)==1){
      if (double.Sample) {
        dataTree <- data.table(data[train.idx,])
        dataEstim <- data.table(data[reestimation.idx,])
      }else{
        dataTree <- data.table(data[full.idx,])
      }
      #pick relevant covariates for tree
      treeRange<-c(cov_sample,(ncolx+1):ncol(dataTree))
      estimRange<-c(cov_sample,(ncolx+1):ncol(dataEstim))
      dataTree <- dataTree[,..treeRange]
      if (double.Sample) {
        dataEstim <- dataEstim[,..estimRange]
      }
    # }else if(class(data)=="data.frame"){
    }else if(inherits(data,"data.frame")){
      if (double.Sample) {
          dataTree <- data.frame(data[train.idx,])
          dataEstim <- data.frame(data[reestimation.idx,])
        }else{
          dataTree <- data.frame(data[full.idx,])
        }
        #pick relevant covariates for tree
        dataTree <- dataTree[,c(cov_sample,(ncolx+1):ncol(dataTree))]
        if (double.Sample) {
          dataEstim <- dataEstim[,c(cov_sample,(ncolx+1):ncol(dataEstim))]
        }
    }


    #change colnames to reflect the sampled cols
    names(dataTree)=nameall_sample
    if(double.Sample) {
      names(dataEstim)=nameall_sample
    }

    #save rdata for debug here, if needed
    formula<-paste0(y,"~",fsample)

    if (double.Sample) {
      tree.obj <- honest.causalTree(formula, data = dataTree,
                                    treatment = treatmentdf[train.idx,],
                                    est_data=dataEstim, est_treatment=treatmentdf[reestimation.idx,],
                                    split.Rule=split.Rule, split.Honest= split.Honest, split.Bucket=split.Bucket,
                                    bucketNum = bucketNum,
                                    bucketMax = bucketMax, cv.option="CT", cv.Honest=T,
                                    minsize = nodesize,
                                    split.alpha = 0.5, cv.alpha = 0.5, xval=0,
                                    HonestSampleSize=est.size, cp=0)
    }else {
      tree.obj <- causalTree(formula, data = dataTree, treatment = treatmentdf[full.idx,],
                             na.action = na.causalTree,
                             split.Rule=split.Rule, split.Honest= split.Honest, split.Bucket=split.Bucket,
                             bucketNum = bucketNum,
                             bucketMax = bucketMax, cv.option="CT", cv.Honest=T,
                             x = FALSE, y = TRUE,
                             split.alpha = 0.5, cv.alpha = 0.5,cv.gamma=0.5,split.gamma=0.5)

    }

    causalForest.obj$trees[[tree.index]] <- tree.obj
    causalForest.obj$inbag[full.idx, tree.index] <- 1
    if (double.Sample) {causalForest.obj$inbag.Est[reestimation.idx, tree.index] <- 1}
  }
  return (causalForest.obj)
}

