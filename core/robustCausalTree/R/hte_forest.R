#' Estimate Heterogeneous Treatment Effect via Random Forest
#'
#' Estimate heterogeneous treatment effect via random forest.
#' In each leaf, the treatment effect is the difference of mean outcome
#' weighted by inverse propensity scores in treatment group and
#' control group.
#' @param gf a fitted generalized random forest object
#' @param ps_linear a character representing name  of  a  column  that
#' stores  linearized  propensity scores.
#' @inheritParams hte_causalTree
#' @returns A list with three elements. The first one is the predicted outcome
#' for each unit. The second is an \code{causalTree} object with the tree split
#' information. The third is a \code{data.frame} summarizing the prediction
#' results.





# Function: causal forest------
#::::::::::::::::::::::::::::::::::::
# Inverse Propensity Score Weighting#
#::::::::::::::::::::::::::::::::::::
hte_forest <- function(outcomevariable,
                       # the name of outcome variabls we are interested in
                       minsize=20,
                       # minimum number of treated observations,
                       # control observations in a leaf
                       crossvalidation = 20,
                       # number of cross-validations to do
                       data = edurose_mediation_20181126,
                       # can be changed, and the defaul one defined here
                       # is edurose_mediation_20181126, the education dataset we are
                       # working on
                       treatment_indicator = 'compcoll25', # treatment variable
                       ps_indicator = 'propsc_com25', # propensity scores
                       ps_linear = 'propsc_com25lin',
                       covariates = c(linear_terms,ps_indicator),
                       negative = FALSE,
                       # can be changed, specify the expected direction
                       # of the treatment effects
                       drawplot = TRUE,
                       # export the graph of tree structure if true
                       # con.num=1,
                       # the number of control variables used in matching
                       # no_indicater="",
                       legend.x = 0.08,
                       legend.y = 0.25,
                       gf,
                       ...){

  if(!exists("tau.forest")){
    stop("please first run grf algorithm and name is as tau.forest")
  }
  # delete all missing values which is required in causal tree model
  # and use it as the train set in machine learning model
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

    subgroup <- (trainset%>%rownames)%in%(node$data%>%rownames)
    # print(subgroup%>%sum)
    x <- grf::average_treatment_effect(gf,
                                  target.sample = "treated",
                                  subset = subgroup)

    z <- x[1]/x[2]
    p = exp(-0.717*z - 0.416*(z^2))

    hte_effect_setup[[x_num]] <<- cbind(x[1],p,x[2],
                                        nrow(node$data)/nrow(trainset)*100%>%round(.,1))

    # keep all related numbers in the environment of node
    node$predicted <- x[1]%>%as.numeric
    node$pvalue <- p
    node$standarderror <- x[2]%>%as.numeric
    node$samplesize <- nrow(node$data)
  })

  opfit_tree <-opfit_tree
  hte_effect <- opfit_tree$Get("predicted")%>%as.numeric
  opfit$frame$yval <- hte_effect

  # create a new variable indicating the estimated treatment effect for each unit
  hte_effect <- opfit$frame$yval[opfit$where]

  # statistics
  ttable <<- unlist(hte_effect_setup)%>%matrix(.,ncol = 4,byrow = TRUE)%>%
    as.data.frame%>%
    `colnames<-`(c("Estimator","pvalue","se","SampleSize"))
  st <- rep("",length(ttable$pvalue))
  st[ttable$pvalue<0.05] <- "*"
  st[ttable$pvalue<0.01] <- "**"
  st[ttable$pvalue<0.001] <- "***"
  ttable$star <<- st
  # ttable$star <<- ifelse(ttable$pvalue<0.1,"*","")
  ttable$SampleSize <- round(ttable$SampleSize,1)

  # If makeing plots, the values from the original tree should be
  # adjusted to the value generated from matching methods
  # adj_effect <- table(hte_effect)%>%as.data.table
  # opfit$frame$yval[match(adj_effect$N,opfit$frame$n)] <- as.numeric(adj_effect$hte_effect)


  # if drawplots is TRUE, make plots and export the plots
  if(drawplot==TRUE){
    # makeplots(opfit,gph,trainset,covariates,outcomevariable)
    makeplots(negative=negative, opfit.=opfit,trainset,
              covariates,outcomevariable,data.=data,ttable,
              # no_indicater,
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
