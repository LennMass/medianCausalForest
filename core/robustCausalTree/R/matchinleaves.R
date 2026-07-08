#' NN Matching in Leaves
#'
#' This intermediate function is used to adjust the heterogeneous treatment
#' effect estimated in each leaf with NN matching.
#'
#' @param trainset a data frame only containing the variables used in
#' the model and missings values are listwise deleted.
#' @param covariates a vector of column names of all
#' covariates (linear terms andpropensity score).
#' @param outcomevariable a character representing the column name
#' of the outcome variable.
#' @param hte_effect_setup a empty list to store the adjusted treatment
#' effect.
#' @param treatment_indicator a character representing the column
#' name of the treatment indicator.
#' @param con.num a number indicating the number of units from control
#' groups to be used in matching
#' @param ... further arguments passed to or from other methods.
#'
#' @returns A list for summarizing the results after matching.


# matching algorithms ----------------------------
# :::::::::::::::::::::::::::::::::::::::::
# matching function: 1:1 matching, full matching, 1:4 matching
# :::::::::::::::::::::::::::::::::::::::::
matchinleaves <- function(trainset=match_data,
                          covariates=covariates,
                          outcomevariable=outcomevariable,
                          hte_effect_setup = hte_effect_setup,
                          treatment_indicator,
                          con.num=1, # the numbers of controls in doing matching
                          ...){
  # variables used in the function
  match_results <- NULL
  standarderror <- NULL
  outcomevariable <- outcomevariable
  con.num <- con.num

  # make sure the matrix for predicting propensity score is of full rank
  covariates.match <- covariates[apply(trainset[,covariates],2,
                                       function(i) length(unique(i))>1)]
  # set up the matching rule
  match.formula <- as.formula(paste(treatment_indicator,'~',covariates[length(covariates)]))

  # print(match.formula)

  #::::::::
  # matching in the leaves
  #::::::::
  # optimal matching in leaves

  if(con.num==1){
    # print("con.num=1")
    trainset_imputed <- optmatch::fill.NAs(trainset)
    dist_matrix <- optmatch::match_on(match.formula,
                       # choose methods
                       # step 1: assign the varible that should be matched exactly, and here is the propensity score strata
                       # matchit can not deal with data with missing values. So I used fill.NAs functions to non-
                       # informatively fill in missing values in original data frame. And these data are prepared to
                       # to do optimal match.
                       # step 2: create treatment to control distance with match_on function
                       # step 3: use fullmatch to complete this optimal matching
                       # notes: exactMatch can be neglected if units are allowed to be matched across levels
                       method = "mahalanobis", data = trainset_imputed)
    match <- optmatch::pairmatch(dist_matrix,
                                 # controls = n (Optional, for example, set n = 1)
                                 # The number of controls used in matching process. But this constraints may not
                                 # be feasible. If so, R will automatically change the constraints and an warning message
                                 # will be issued.
                                 data = trainset_imputed,
                                 tol=0.001)
  }else if(con.num==4){
    # delete rows that contains missing values
    # MatchIt does not allow missing values
    matchdata_nomiss <-
      trainset[,as.list(c(outcomevariable,covariates[length(covariates)],treatment_indicator))]
    matchdata_nomiss <- na.omit(matchdata_nomiss)

    # the default setting is matching without replacement
    # (by setting replacement to TRUE to do matching with replacement)
    X  <- matchdata_nomiss[,covariates[length(covariates)]]
    Y  <- matchdata_nomiss[,outcomevariable]
    Tr  <- matchdata_nomiss[,treatment_indicator]

    match  <- Matching::Match(Y=Y, Tr=Tr, X=X, estimand = "ATE",
                              M=con.num,Weight = 2,replace=TRUE)

  }else{
    match <- optmatch::fullmatch(match.formula,
                                 data = fill.NAs(trainset))
  }

  match_results[[length(match_results)+1]] <- match

  # get the results
  if(con.num!=4){
    # get the new dataset after matching
    data.match <- cbind(optmatch::fill.NAs(trainset),match)
    # get the treatment effects


    average <- c(paste0('mean(',outcomevariable,')'))
    average_name <- c(outcomevariable)
    hte_effect_help <- data.match[which(!is.na(data.match$match)),
                                  c(treatment_indicator,outcomevariable,'match')]

    hte_effect_help <- dplyr::group_by(.data = hte_effect_help, dplyr::across(dplyr::all_of(c('match',treatment_indicator))))


    hte_effect_help <- dplyr::summarise(.data = hte_effect_help, groupMean = mean(hte_effect_help[[outcomevariable]], na.rm = TRUE))

    hte_effect_help <- as.data.frame(hte_effect_help)


    hte_effect_help <- reshape(hte_effect_help,idvar='match',timevar = treatment_indicator,direction = 'wide')

    # extract the control and treatment group

    treatment <- hte_effect_help[,3]
    control <- hte_effect_help[,2]

    hte_effect_help <- mean(treatment-control)
    # print(paste("WDIM:", hte_effect_help))

    star <- t.test(treatment-control,conf.level = .9)
    # pvalue
    pvalue <- star$p.value
    # standard error
    standarderror <- sqrt( var(treatment-control,na.rm = T)/length(na.omit(treatment-control)) )
  }else {
    pvalue <- 2*(1-pnorm(abs(match$est/match$se)))
    standarderror <- match$se
    hte_effect_help <- match$est
  }

  # print(c('The number of units forget to consider is',sum(is.na(hte_effect_setup))))
  # print(round(c(hte_effect_help,pvalue,standarderror),4))
  return(round(c(hte_effect_help,pvalue,standarderror),4))
}
