#' Include the Javascript Used in Shiny
#'
#' intermediate function used to include necessary javascript
#' to visualize tree structures and estimated treatment effect
#' in shiny
#'
#' @param \dots There is no required arguments in this function. But user
#' could manipulate to include different css files.
#' @returns No return value. It is used to pass the Javascript to Shiny.


bundScript <- function(...) {
  includeScript(system.file("js/bundle.js", package = "robustCausalTree"))
}
