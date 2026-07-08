#' Get the Current Working Directory
#'
#' get the current work directory and set it as the default directory to
#' save the shiny files temporarily
#'
#' @returns a temporary file path

getDefaultPath <- function(){
  paste(tempdir(), "/shinyapp", sep="")
}

