#' Save Javascript Embedded in Shiny App
#'
#' @inheritParams saveFiles
#' @returns No return value. It is used to save necessary files temporarily
#' to run Shiny App.


saveBCSS <- function(filePath){
  fileConn <- file(paste(filePath, "/shinyapp/www/bundle.css", sep=""))
  writeLines(c("#svg-wrapper.svelte-1kzzeln{position:fixed;width:100vw;height:100vh;overflow:scroll;overflow-x:scroll}",
               "#steps-controller.svelte-8mmdi8{position:fixed;top:0;left:0}li.svelte-8mmdi8{cursor:pointer}.current.svelte-8mmdi8{font-weight:bold}",
               "/*# sourceMappingURL=bundle.css.map */"), fileConn)
  close(fileConn)
}
