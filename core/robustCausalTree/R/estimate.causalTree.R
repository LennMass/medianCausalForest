#' estimate causal Tree
#' @inheritParams causalTree
#' @param object A tree-structured fit \code{rpart} object, such as one
#' generated as a \code{causalTree} fit.
#' @param data New data frame to be used for estimating effects within leaves.
#' @param treatment The treatment status of observations in the new
#' dataframe, where 1 represents treated and 0 represents control.
#'
#' @details
#' When the leaf contains only treated or control cases, the function will
#' trace back to the leaf's parent node recursively until the parent can
#' be used to compute causal effect. Please see Athey and Imbens
#' \emph{Machine Learning Methods for Estimating Heterogeneous Causal
#' Effects} (2015) for details.
#'
#' @returns Intermediate estimation results for an \code{causalTree} object
#'


estimate.causalTree <- function(object, data, weights, treatment, na.action = na.causalTree)
{
    if (!inherits(object, "rpart")) stop("Not a legitimate \"rpart\" object")
    # get the leaf of the object
    leaf <- as.numeric(row.names(object$frame)[which(object$frame$var == "<leaf>")])
    Terms <- object$terms
    data$tr <- treatment
    if (missing(weights)) {
        Terms <- object$terms
        m <- model.frame(Terms, data = data, na.action = na.action, treatment = tr,
                         xlev = attr(object, "xlevels"))
    } else {
        attr(Terms, "dataClasses")["weights"] <- "numeric"
        data$w <- weights
        m <- model.frame(Terms, data = data, na.action = na.action, treatment = tr,
                         weights = w, xlev = attr(object, "xlevels"))
    }

    if (!is.null(cl <- attr(Terms, "dataClasses")))
        .checkMFClasses(cl, m, TRUE)

    treatment <- m$`(treatment)`
    n <- nrow(m)
    Y <- model.response(m)
    X <- causalTree.matrix(m)
    if (missing(weights))
        wts <- rep(1, nrow(m))
    else
        wts <- model.weights(m)
    new_object <- data.table::copy(object)
    ans <- honest.est.causalTree(new_object, X, wts, treatment, Y)
    return(ans)
}
