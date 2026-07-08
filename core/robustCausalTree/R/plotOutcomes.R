#' Intermediate function for \code{hte_plot_line}
#'
#' Plots the
#' different least squares models used to estimate heterogeneous
#' treatment effects(HTE) at each node. At each node, this
#' visualization aims to show how the estimated treatment effect
#' differs when using ordinary least squares and weighted least
#' squares methods. The weighted least squares method in this
#' package uses
#' inverse propensity scores as weights, in order to reduce
#' bias due to confounding variables.
#'
#' @param treatment a character representing the column name
#' for the treatment variable in the causal setup
#' @param outcome a character representing the column name
#' of the outcome variable.
#' @param propscores a character representing the column name of
#' the propensity score.
#' @param confInt a logical value indicating whether adding the 95%
#' confidence interval. The default is set as \code{TRUE}.
#' @param colbyWt a logical value indicating whether the points are
#' are colored according to inverse propensity scores. The default is
#' set as FALSE.
#' @param xlab,ylab,title Characters representing the name for x axis,
#' y axis, and main title for each node.
#' @param gamma,lambda numbers indicating the bias level used in
#' sensitivity analysis
#' @param ... further arguments passed to or from other methods.
#'
#' @returns A summary table after adjusting the estimates with
#' inverse probability weighting (ipw).


# plot outcomes for each node -----------
plotOutcomes <- function(treatment, outcome, propscores,
                         confInt=TRUE, colbyWt=FALSE,
                         # xlab='Treatment',
                         # ylab='Outcome',
                         # title='Treatment Effects: Weighted and Unweighted',
                         ylab='',
                         xlab='',
                         title='',
                         gamma=0,
                         lambda=0,
                         ...
){
  x <- treatment
  y <- outcome
  wt <- abs(treatment - propscores)/(propscores*(1-propscores))

  if (colbyWt) {
    xy <- cbind(x, y)
    df <- data.frame(xy)
    xy_wt <- df[rep(row.names(df), wt), 1:2]
    x <- xy_wt$x
    y <- xy_wt$y
  }

  df <- data.frame(x, y)
  cols <-  colorRampPalette(c("#FFD700", "#ff800e", "#ff4500", "#8B0000"))(100)

  control <- subset(df, x==0)
  dens_control <- getDensities(control$x, control$y)
  col_control <- cols[dens_control*100]

  treated <- subset(df, x==1)
  dens_treat <- getDensities(treated$x, treated$y)
  col_treat <- cols[dens_treat*100]

  col <- c(col_control, col_treat)

  ## Plot points, colored by density
  # layout(matrix(c(1, 2), 1, 2, byrow=TRUE), widths = c(4, 1), heights = c(1,1))
  # par(mar = c(3, 4, 2, 1), mgp=c(2,0.5,0))
  plot(y~x, col=col , main=list(title, font=2, cex=0.7), pch=19, xlab=xlab, ylab=ylab, xaxp=c(0,1,1), axes=FALSE)
  axis(2, at=round(seq(min(outcome), max(outcome), (max(outcome)-min(outcome))/5),2) )

  ## Create models with ordinary least squares and weighted least squares
  model_og <- lm(outcome~treatment)
  model_wt <- lm(outcome~treatment, weights=wt)

  ## Plot models and corresponding estimated treatment effects
  bias=gamma*lambda
  ## Plotting confidence interval
  if (confInt) {
    summ1 <- summary(model_og)
    SE_int1 <- coef(summ1)[1,2]
    SE_x1 <- coef(summ1)[2,2]
    summ2 <- summary(model_wt)
    SE_int2 <- coef(summ2)[1,2]
    SE_x2 <- coef(summ2)[2,2]

    polygon(c(0, 0, 1, 1), c(model_og$coefficients[1]+SE_int1,
                             model_og$coefficients[1]-SE_int1,
                             model_og$coefficients[2]-SE_x1-bias + model_og$coefficients[1]-SE_int1,
                             model_og$coefficients[2]+SE_x1-bias + model_og$coefficients[1]+SE_int1),
            col=rgb(1, 0, 0,0.25), border=NA)
    polygon(c(0, 0, 1, 1), c(model_wt$coefficients[1]+SE_int2,
                             model_wt$coefficients[1]-SE_int2,
                             model_wt$coefficients[2]-SE_x2-bias + model_wt$coefficients[1]-SE_int2,
                             model_wt$coefficients[2]+SE_x2-bias + model_wt$coefficients[1]+SE_int2),
            col=rgb(0, 0, 1,0.25), border=NA)

  }

  ## Plotting bars to show treatment effect
  rect(0.99, model_og$coefficients[2]-bias + model_og$coefficients[1],
       1.01, model_og$coefficients[1])
  text(x = 0.98,
       y = model_og$coefficients[1],
       labels = paste("RAW = ", round(model_og$coefficients[2]-bias,3)), cex=0.5, adj=1)

  rect(0.99, model_wt$coefficients[2]-bias + model_wt$coefficients[1],
       1.01, model_wt$coefficients[1], border='blue')
  text(x = 0.98,
       y = 0.5*model_wt$coefficients[2]-bias + model_wt$coefficients[1],
       labels = paste("IPW = ", round(model_wt$coefficients[2]-bias,3)), col='blue', cex=0.5, adj=1)

  ## Plotting lines to show least squares models
  lines(seq(0,1,by=0.1), (model_og$coefficients[2]-bias)*seq(0,1,by=0.1)+ model_og$coefficients[1])
  lines(seq(0,1,by=0.1), (model_wt$coefficients[2]-bias)*seq(0,1,by=0.1)+ model_wt$coefficients[1], col="blue")

  summary(model_og)
  summary(model_wt)
}
