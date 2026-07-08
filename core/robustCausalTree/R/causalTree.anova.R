#' Intermediate function for \code{causalTree}
#'
#' @param y outcome variable
#' @param offset this can be used to specify an a priori known
#' component
#' to be included in the linear predictor during fitting. This should be
#' \code{NULL} or a numeric vector of length equal to the number of cases.
#' One or more \link{offset} terms can be included in the formula instead or as
#' well, and
#' if more than one is specified their sum is used. See \link{model.offset}.
#' @param wt optional weights
#' @returns No return value.

robustCausalTree.anova <- function(y, offset, wt)
{
    if (!is.null(offset)) y <- y - offset
      list(y = y, numresp = 1L, numy = 1L,
	 summary = function(yval, dev, wt, ylevel, digits) {

	   paste0("  causal effect=", formatg(yval, digits),
	          ", error=" , formatg(dev, digits))
         },
	 text = function(yval, dev, wt, ylevel, digits, n, use.n ) {
	     if (use.n) paste0(formatg(yval, digits), "\nn=", n) else
             formatg(yval, digits)
         })
}
