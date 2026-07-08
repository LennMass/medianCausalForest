#' Save Necessary Files to Run Shiny App
#'
#' This function is to save files necessary to run Shiny app to
#' visualize causal tree and the estimated heterogeneous treatment
#' effects in an interactive way.
#'
#' @param filePath a character string representing the path name
#' to save the files temporarily.
#' @inheritParams runDynamic
#' @returns No return value. It is used to save necessary files temporarily
#' to run Shiny App.
#'



# change these three lines later

# Files needed:
# ui.R
# server.R
# bundle.css - function saveBCSS
# bundle.js - not created by this code
# causal_tree_raw.csv - made in saveFiles
# causal_tree_structure.json - made in saveFiles
# global.css - function saveGCSS
# index.html - made in saveFiles

saveFiles <- function(model, data, outcomevariable, treatment_indicator, propensity_score="", filePath="") {

  if (nchar(filePath) == 0){
    filePath <- tempdir()
  }
  # Create and save causal tree structure

  opfit <- model$tree
  prty <- partykit::as.party(opfit)
  opfit_tree <- data.tree::as.Node(prty)
  opfit_tree$attributesAll

  # get node labels :::::::::::::::::::::::::::::::
  nn_name <- 0
  opfit_tree$Do(function(node){
    opfit$frame$var[which(opfit$frame$var=="<leaf>")] <- NA
    nn_name <<- 1 + nn_name
    node$splitkey <- opfit$frame$var[nn_name]
  }
  )

  # get split label names :::::::::::::::::::::::
  # only if there is label to be replaced can this algorithm be used
  if( length(attr(data,"var.label"))>0 ){

    opfit$frame$var <- as.character(opfit$frame$var)
    opfit$frame$var <- attr(data,"var.label")[match(opfit$frame$var,colnames(data))]
  }else{
    opfit$frame$var <- as.character(opfit$frame$var)
  }

  nn_name <- 0
  opfit_tree$Do(function(node){
    opfit$frame$var[which(opfit$frame$var=="<leaf>")] <- NA
    nn_name <<- 1 + nn_name
    node$splitname <- opfit$frame$var[nn_name]
  }
  )

  # get split rules and split names (the same with splitkey_label)
  sp_name <- 0
  sp <- list()
  opfit_tree$Do(function(node){
    sp_name <<- 1 + sp_name
    sp[[sp_name]] <<- paste(node$splitname, # got from the previous step
                            opfit_tree$Get('splitLevel')[as.numeric(attr(node$children,'name'))])
    names(sp[[sp_name]]) <<- as.numeric(attr(node$children,'name'))
  }
  )

  sp <- unlist(sp)[order( as.numeric(names(unlist(sp))) )]
  sp <- sp[-which(is.na(names(sp)))]

  sp <- c(opfit$frame$var[1],sp)
  names(sp) <- c()
  sp[1] <- 'Root'

  # change names ::::::::::::::::
  nn_name <- 0
  opfit_tree$Do(function(node){
    nn_name <<- 1 + nn_name
    node$name <- sp[nn_name]
  }
  )

  # splitkey labels: the same as the split name :::::::::::::::::
  opfit_tree$Do(function(node){
    node$splitkey_label <- node$splitname
  })


  opfit_tree_v <- opfit_tree

  opfit_tree_v$Get(function(node){
    node$data <- c()
    node$fitted <- c()
    node$nodeinfo <- c()
    node$partyinfo <- c()
    node$terms <- c()
    node$datasetplots <- c()
    node$fill <- c()
    node$HTE <- c()
    node$SizeOfNode <- c()
    node$tooltip <- c()
    node$splitlevels <- c()
    node$split <- c()
    node$splitname <- c()
    node$treatment_effects <- c()
  })

  xxx <- jsonlite::toJSON(data.tree::ToListExplicit(opfit_tree_v,unname = TRUE,
                                   nameName = 'name'),
                    pretty=TRUE, auto_unbox = TRUE) ### THIS WORKED

  dir.create(paste(filePath, "/shinyapp", sep=""))
  dir.create(paste(filePath, "/shinyapp/www", sep=""))
  write(xxx, paste(filePath, "/shinyapp/www/causal_tree_structure.json", sep=""))


  # Save data, index.html, bundle.css, global.css, ui.R, server.R
  my_data <- data.frame(data)
  names(my_data)[names(my_data) == outcomevariable] <- "outcome"
  names(my_data)[names(my_data) == treatment_indicator] <- "treatment"
  if (nchar(propensity_score) != 0){
    names(my_data)[names(my_data) == propensity_score] <- "propscore"
  }

  write.csv(my_data, paste(filePath, "/shinyapp/www/causal_tree_raw.csv", sep=""), row.names = FALSE)
  saveInd(filePath)
  saveBCSS(filePath)
  saveGCSS(filePath)
  saveUI(filePath)
  saveServ(filePath)

  # Save bundle.js
  if (nchar(propensity_score) != 0){
    p <- system.file("shiny_ipw", "www", "bundle.js", package = "robustCausalTree")
    fileConn <- file(paste(filePath, "/shinyapp/www/bundle.js", sep=""))
    writeLines(readLines(p), fileConn)
    close(fileConn)
  } else{
    p <- system.file("shiny", "www", "bundle.js", package = "robustCausalTree")
    fileConn <- file(paste(filePath, "/shinyapp/www/bundle.js", sep=""))
    writeLines(readLines(p), fileConn)
    close(fileConn)
  }


}

