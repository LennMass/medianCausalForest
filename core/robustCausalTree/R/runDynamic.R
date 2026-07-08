#' Visualize Causal Tree and Treatment Effects via Shiny
#'
#' @param model a tree model constructed by \code{hte_causalTree, hte_matchinleaves,
#' or hte_ipw}.
#' @param data a data frame containing the variables in the model.
#' @param outcomevariable a character representing the column name
#' of the outcome variable.
#' @param treatment_indicator a character representing the column name
#' for the treatment variable in the causal setup.
#' @param propensity_score a character representing the column name of
#' the propensity score.
#'
#' @return a Shiny page.
#'


runDynamic <- function(model, data, outcomevariable, treatment_indicator, propensity_score=""){
  temppath <- tempdir()
  saveFiles(model, data, outcomevariable, treatment_indicator, propensity_score, temppath)

  shiny::runApp(appDir = paste(temppath, "/shinyapp/", sep=""),
         launch.browser = getOption("shiny.launch.browser", interactive()))
  # print(0)
  unlink(paste(tempdir(), "/shinyapp", sep=""), recursive=TRUE)

}

