#' Intermediate function for \code{causalTree}
#'
#' get model frame of causalTree, same as rpart
#'
#' @param formula a \link{formula}, with a response but no interaction terms.
#' If this is a data frame, it is taken as the model frame (see
#' \link{model.frame}).
#'
#' @inheritParams causalTree
#' @returns a model frame for \code{causalTree}.

model.frame.causalTree <- function(formula, ...)
{
    m <- formula$model
    if (!is.null(m)) return(m)
    oc <- formula$call
    if (substring(deparse(oc[[1L]]), 1L, 7L) == "predict") {
        m <- eval(oc$newdata)
        if (is.null(attr(m, "terms"))) {
            object <- eval(oc$object)
            m <- model.frame(object$terms, m, na.causalTree)
        }
        return(m)
    }
    while(!deparse(oc[[1L]]) %in%  c("causalTree", "causalTree::causalTree", "causalTree:::causalTree"))
    # while(!deparse(oc[[1L]]) %in%  c("htetree", "htetree::htetree", "htetree:::htetree"))
        oc <- eval(oc[[2L]])$call
    oc$subset <- names(formula$where)
    oc$method <- formula$method
    eval(oc)
}
