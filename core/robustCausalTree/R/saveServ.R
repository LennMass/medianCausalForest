#' Save Shiny Server Temporarily
#'
#' @inheritParams saveFiles
#' @returns No return value. It is used to save necessary files temporarily
#' to run Shiny App.


saveServ <- function(filePath){
  servstr <- "shinyServer(function(input, output) {
  output$distPlot <- renderPlot({
    dist <- rnorm(input$obs)
    hist(dist)
  })
})"

  servsplt <- strsplit(servstr, "\n")
  servvec <- unlist(servsplt)

  fileConn <- file(paste(filePath, "/shinyapp/server.R", sep=""))
  writeLines(servvec, fileConn)
  close(fileConn)


}
