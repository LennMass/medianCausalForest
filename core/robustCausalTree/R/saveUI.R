#' Save Shiny UI Temporarily
#'
#' @inheritParams saveFiles
#' @returns No return value. It is used to save necessary files temporarily
#' to run Shiny App.
#'
saveUI <- function(filePath){
  bundstr <- "shinyServer(function(input, output) {
})
ui <- fluidPage(
  tags$script(src = 'bundle.js')
)"

  bundsplt <- strsplit(bundstr, "\n")
  bundvec <- unlist(bundsplt)

  fileConn <- file(paste(filePath,"/shinyapp/ui.R",sep=""))
  writeLines(bundvec, fileConn)
  close(fileConn)
}
