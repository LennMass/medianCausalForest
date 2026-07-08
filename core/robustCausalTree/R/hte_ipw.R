#' Estimate Heterogeneous Treatment Effect via Adjusted Causal Tree
#'
#' Estimate heterogeneous treatment effect via adjusted causal tree.
#' In each leaf, the treatment effect is the difference of mean outcome
#' weighted by inverse propensity scores in treatment group and
#' control group.
#' @param ps_linear a character representing name  of  a  column  that
#' stores  linearized  propensity scores.
#' @param check if TRUE, generates 100 trees and outputs most common
#' tree structures and their frequency
#' @inheritParams hte_causalTree
#' @returns predicted treatment effect and the associated tree
#' @examples
#' library(rpart)
#' library(htetree)
#' hte_ipw(outcomevariable="outcome",
#' data=data.frame("confounder"=c(0, 1, 1, 0, 1, 1),
#' "treatment"=c(0,0,0,1,1,1), "prop_score"=c(0.4, 0.4, 0.5, 0.6, 0.6, 0.7),
#' "outcome"=c(1, 2, 2, 1, 4, 4)), treatment_indicator = "treatment",
#' ps_indicator = "prop_score", covariates = "confounder")



# Function 3: inverse propensity score weighting------
#::::::::::::::::::::::::::::::::::::
# Inverse Propensity Score Weighting#
#::::::::::::::::::::::::::::::::::::

hte_ipw <- function(outcomevariable,
                    # the name of outcome variabls we are interested in
                    minsize=20,
                    # minimum number of treated observations,
                    # control observations in a leaf
                    crossvalidation = 20,
                    # number of cross-validations to do
                    data,
                    # can be changed, and the defaul one defined here
                    # is edurose_mediation_20181126, the education dataset we are
                    # working on
                    treatment_indicator, # treatment variable
                    ps_indicator, # propensity scores
                    ps_linear=NULL,
                    covariates,
                    negative = FALSE,
                    # can be changed, specify the expected direction
                    # of the treatment effects
                    drawplot = TRUE,
                    # export the graph of tree structure if true
                    # no_indicater="",
                    varlabel=NULL,
                    maintitle="Heterogeneous Treatment Effect Estimation",
                    legend.x = 0.08,
                    legend.y = 0.25,
                    check=FALSE,
                    ...){
  # delete all missing values which is required in causal tree model
  # and use it as the train set in machine learning model
  # remotes::install_github("susanathey/causalTree",build = FALSE)
  trainset <-  data[!is.na(data[,outcomevariable]),]

  # set up the formula used for constructing causal tree
  # export the formula for causal tree model: Y~X
  if(nchar(ps_indicator)>0){
    covariates_ <- c(covariates,ps_indicator) # non-linear ps score
    formula <- as.formula(paste(outcomevariable," ~ ",
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- c(covariates,ps_linear) # linear ps score
  }else{
    formula <- as.formula(paste(outcomevariable," ~ ",
                                paste(covariates,collapse = '+'), collapse= "+"))
    covariates <- covariates
  }
  # if("rpart" %in% rownames(installed.packages()) == FALSE) {install.packages("data.tree")}
  # library("rpart")
  requireNamespace("rpart")
  # check=TRUE
  # implement the check
  if(check==TRUE){
    check_tree <- replicate(100,
                             {if(nchar(ps_indicator)>0){
                               # contruct tree
                               tree <- causalTree(formula,
                                                              # specify the model, outcome variable ~ covariates
                                                              data = trainset, # specify the dataset to be used
                                                              treatment = trainset[,treatment_indicator],
                                                              # specify the treatment variable, must be 0/1 indicator
                                                              split.Rule = "CT",
                                                              # specify split rule; for causal tree, use "CT"
                                                              # NOTE: there are four different splitting rules,
                                                              # they are different in the cross-validation criteria used
                                                              # to determine the tree structure
                                                              # 1 - TOT
                                                              # 2 - CT
                                                              # 3 - fit
                                                              # 4 - tstats
                                                              # 5 - totD
                                                              # 6 - ctD
                                                              # 7 - fitD
                                                              # 8 - tstatsD
                                                              cv.option = "CT", # specifify cross validation method
                                                              # and there are four different methods -- tot, ct, fit, tstats
                                                              # for causal tree, use "CT"
                                                              split.Honest = T, cv.Honest = T, split.Bucket = F,

                                                              xval = crossvalidation,
                                                              # number of cross-validations to do and the default number is 20
                                                              cp = 0,
                                                              propensity = trainset[,ps_indicator],
                                                              # specify the propensity score; if is not specified, it will use sum(treatment) / nobs as the propensity score
                                                              minsize = minsize # minimum number of treated observations, control observations in a leaf
                                                              # the default minimum size is 20, according to Jennie and Yu Xie's paper (Estimating Heterogeneous Treatment Effects with Observational Data, 2012)
                               )}else{
                                 tree <- causalTree(formula,
                                                                data = trainset,
                                                                treatment = trainset[,treatment_indicator],
                                                                split.Rule = "CT",
                                                                cv.option = "CT",
                                                                split.Honest = T, cv.Honest = F, split.Bucket = F,
                                                                xval = crossvalidation,
                                                                cp = 0,
                                                                minsize = minsize
                                 )
                               }

                               # prune this tree model to avoid the overfitting issues
                               # get the complexity parameter (cp) to be trimmed--the least important splits
                               opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
                               # recursively snipping off the least important tree based on the complexity parameter (cp)
                               opfit <- rpart::prune(tree, opcp)
                               # paste(opfit$frame$var,collapse=" ")
                               prty <- partykit::as.party(opfit)
                               opfit_tree <- data.tree::as.Node(prty)
                               # i <- 0
                               opfit_tree$Do(function(node){
                                 #i <<- i+1
                                 #node$name <- opfit$frame[i,"var"]
                                 ## below: tanvi edit
                                 if(node$isRoot){
                                   node$name <- "root"
                                 } else{
                                   node$name <- paste(node$parent$splitname, node$splitLevel)
                                 }
                               })
                               tree_structure <- capture.output(opfit_tree)
                               list(tree_structure[2:length(tree_structure)])
                             }
    )
  }

  if(check==TRUE){
    message(paste0("This generated tree structure appears ",max(table(sapply(check_tree,function(i){paste(i,collapse = "")})))," times in 100 iterations"))
    # # print(check_tree)
    message("Summary of tree structures:")
    check_tree1 <- lapply(check_tree,cbind)
    check <- lapply(unique(check_tree1),function(i){
      message(paste0("The following tree structure appears ",sum(sapply(check_tree1,function(j){identical(i,j)}))," times in 100 iterations:"))
      # print(i)
    })
  }

  # set up propensity score
  if(nchar(ps_indicator)>0){
    # contruct tree
    tree <- causalTree(formula,
                                   # specify the model, outcome variable ~ covariates
                                   data = trainset, # specify the dataset to be used
                                   treatment = trainset[,treatment_indicator],
                                   # specify the treatment variable, must be 0/1 indicator
                                   split.Rule = "CT",
                                   # specify split rule; for causal tree, use "CT"
                                   # NOTE: there are four different splitting rules,
                                   # they are different in the cross-validation criteria used
                                   # to determine the tree structure
                                   # 1 - TOT
                                   # 2 - CT
                                   # 3 - fit
                                   # 4 - tstats
                                   # 5 - totD
                                   # 6 - ctD
                                   # 7 - fitD
                                   # 8 - tstatsD
                                   cv.option = "CT", # specifify cross validation method
                                   # and there are four different methods -- tot, ct, fit, tstats
                                   # for causal tree, use "CT"
                                   split.Honest = T, cv.Honest = T, split.Bucket = F,

                                   xval = crossvalidation,
                                   # number of cross-validations to do and the default number is 20
                                   cp = 0,
                                   propensity = trainset[,ps_indicator],
                                   # specify the propensity score; if is not specified, it will use sum(treatment) / nobs as the propensity score
                                   minsize = minsize # minimum number of treated observations, control observations in a leaf
                                   # the default minimum size is 20, according to Jennie and Yu Xie's paper (Estimating Heterogeneous Treatment Effects with Observational Data, 2012)
    )}else{
      tree <- causalTree(formula,
                                     data = trainset,
                                     treatment = trainset[,treatment_indicator],
                                     split.Rule = "CT",
                                     cv.option = "CT",
                                     split.Honest = T, cv.Honest = F, split.Bucket = F,
                                     xval = crossvalidation,
                                     cp = 0,
                                     minsize = minsize
      )
    }

  # prune this tree model to avoid the overfitting issues
  # get the complexity parameter (cp) to be trimmed--the least important splits
  opcp <- tree$cptable[,1][which.min(tree$cptable[,4])]
  # recursively snipping off the least important tree based on the complexity parameter (cp)
  opfit <- rpart::prune(tree, opcp)

  # matchin in leaves and return the predicted heterogeneous treatment effect
  prty <- partykit::as.party(opfit)
  opfit_tree <- data.tree::as.Node(prty)

  hte_effect_setup <- list()
  # Matching Algorithms
  x_num <- 0

  opfit_tree$Do(function(node){
    x_num <<- x_num+1
    # extract data
    ipw_data <- data[as.numeric(rownames(node$data)),]
    wt <- (ipw_data[,treatment_indicator] - ipw_data[,ps_indicator])/(ipw_data[,ps_indicator] * (1-ipw_data[,ps_indicator]))
    # fit the model
    fit <- lm(as.formula(paste0(outcomevariable,'~',treatment_indicator)),
              data = ipw_data,
              weights = abs(wt))
    # print(summary(fit)$coef[2,])
    hte_effect_setup[[x_num]] <<- cbind(t(summary(fit)$coef[2,c(1,4,2)]),
                                        round(nrow(node$data)/nrow(trainset)*100,1))

    # keep all related numbers in the environment of node
    node$predicted <- summary(fit)$coef[2,1]
    node$pvalue <- summary(fit)$coef[2,4]
    node$standarderror <- summary(fit)$coef[2,2]
    node$samplesize <- nrow(node$data)
  })
  opfit_tree<<-opfit_tree
  hte_effect <- as.numeric(opfit_tree$Get("predicted"))
  opfit$frame$yval <- hte_effect

  # create a new variable indicating the estimated treatment effect for each unit
  hte_effect <- opfit$frame$yval[opfit$where]

  # statistics
  ttable <<- matrix(unlist(hte_effect_setup),ncol = 4,byrow = TRUE)
  ttable <<- as.data.frame(ttable)
  colnames(ttable) <<- c("Estimator","pvalue","se","SampleSize")

  st <- rep("",length(ttable$pvalue))
  st[ttable$pvalue<0.05] <- "*"
  st[ttable$pvalue<0.01] <- "**"
  st[ttable$pvalue<0.001] <- "***"
  ttable$star <<- st
  ttable$SampleSize <- round(ttable$SampleSize,1)

  # If makeing plots, the values from the original tree should be
  # adjusted to the value generated from matching methods
  # adj_effect <- table(hte_effect)%>%as.data.table
  # opfit$frame$yval[match(adj_effect$N,opfit$frame$n)] <- as.numeric(adj_effect$hte_effect)


  # if drawplots is TRUE, make plots and export the plots
  if(drawplot==TRUE){
    # makeplots(opfit,gph,trainset,covariates,outcomevariable)
    makeplots(negative=negative, opfit.=opfit,gph=gph,trainset,
              covariates,outcomevariable,data.=data,ttable,varlabel,
              maintitle,#no_indicater,
              legend.x,legend.y)
  }else{
    message(c('Drawplot = ', drawplot))
  }

  # export the results:
  output <- cbind(hte_effect)
  output <- as.data.frame(output)
  colnames(output) <- paste0(outcomevariable,"_predictedTE")

  return(list(predictedTE = output, tree = opfit, matching_table = ttable))
}
