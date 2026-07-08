#' Clear Temporary Files
#'
#' The files for shiny are saved in a temporary directory. The files can be
#' cleared manually using the `clearTemp()` function, or will automatically
#' be cleared when you close R
#'
#' @returns no return value, to unlink files under the temp folder
#'
#'


clearTemp <- function(){
  unlink(paste(tempdir(), "/shinyapp", sep=""), recursive=TRUE)
}
