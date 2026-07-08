#' Intermediate function for \code{causalTree}
#'
#' @param x input training data
#' @param digits number of digits to be kept
#' @param format format of exported vector
#' @returns No return value, called for formatting the exported estimates
#'
## format a set of numbers using C's "g" format
## No longer exported.
formatg <- function(x, digits = getOption("digits"),
                    format = paste0("%.", digits, "g"))
{
    if (!is.numeric(x)) stop("'x' must be a numeric vector")

    temp <- sprintf(format, x)
    if (is.matrix(x)) matrix(temp, nrow = nrow(x)) else temp
}
