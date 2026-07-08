#' Getting Distribution in Treatment and Control Groups
#'
#' Getting the density of distribution in treatment and control groups,
#' which will be displayed in the
#'
#' @param treatment A character representing the name of treatment indicator.
#' @param outcome A character representing the name of outcome variable.
#'
#' @returns vector of corresponding densities for each value of outcome vector


# get density of the distribution in different groups------
getDensities <- function(treatment, outcome){
  df <- data.frame(x=treatment, y=outcome)

  # Find the proportion of points at each possible value of the outcome
  grouped <- dplyr::group_by(df, y)

  if(dplyr::n_groups(grouped) < 20){
    grouped_data <- data.frame(dplyr::summarize(grouped, n= dplyr::n()))
    grouped_data$dens <- grouped_data$n / nrow(df)

    # Based on above calculations, create an array that has the density at each data point
    grouped_data$lowwaprop <- round(grouped_data$y, 3)
    rownames(grouped_data) <- grouped_data$y
    y_str <- as.character(round(df$y, 3))
    densities <- grouped_data[y_str,'dens']
  } else {
    d <- density(outcome, n=100*(max(outcome) - min(outcome)),from=min(outcome), to=max(outcome))
    dens_df <- data.frame(dens=d$y, row.names=round(d$x, 2))
    ch <- as.character(round(outcome,2))
    densities <- dens_df[ch, ]
  }

  return(densities)

}
