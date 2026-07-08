#' Visualize the Estimated Results
#'
#' The function \code{hte_plot_line} takes a model created by
#' causal tree, as well as the adjusted version, and plots the
#' different least squares models used to estimate heterogeneous
#' treatment effects(HTE) at each node.  At each node, this
#' visualization aims to show how the estimated treatment effect
#' differs when using ordinary least squares and weighted least
#' squares methods. The weighted least squares method in this
#' package uses
#' inverse propensity scores as weights, in order to reduce
#' bias due to confounding variables.
#'
#' @param gamma,lambda numbers indicating the bias level used in
#' sensitivity analysis
#' @param ... further arguments passed to or from other methods.
#' @inheritParams hte_plot
#' @returns No return value, used for plotting the estimated results
#' with lines.
#'


# function 2: plot Outcome -------
hte_plot_line <- function(model,data,treatment_indicator = NULL,
                          outcomevariable,propensity_score,
                          plot.title = "Visualization of the Tree",
                          gamma=0,
                          lambda=0,
                          ...){
  # model <- xxx
  # data <- edurose_mediation_20181126
  # treatment_indicator <- "compcoll25"
  # outcomevariable <- "lowwaprop"
  # propensity_score <- "propsc_com25_rf"

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # convert the tree model to
  opfit <- model$tree
  prty <- partykit::as.party(opfit)
  opfit_tree <- data.tree::as.Node(prty)

  # add the necessary statistics
  # add estimated treatment effect ::::::::::::::::::::::::::::::
  i <- 0
  opfit_tree$Do(function(node){
    opfit$frame$yval[1] <- mean(model$predictedTE[,1],na.rm = TRUE)
    i <<- i+1
    node$HTE <- opfit$frame$yval[i]
  }
  )


  # add treatment indicator ::::::::::::::::::::::::::::::
  if(is.null(treatment_indicator)==FALSE){
    opfit_tree$Do(function(node){
      loc <- match(rownames(node$data),
                   rownames(data))
      node$data$treatment <- data[,treatment_indicator][loc]
    })}else{
      stop("Fails to get the treatment indicator")
    }

  # setup layout matrix ::::::::::::::::::::::::::::::::::
  layout.matrix <-
    matrix(as.numeric( t(data.tree::ToDataFrameTypeCol(opfit_tree)) ),
           nrow = opfit_tree$height)

  layout.matrix <- apply(layout.matrix,2,function(ccol){
    list(ccol,ccol)
  })
  layout.matrix <- matrix(unlist(layout.matrix),nrow = opfit_tree$height)
  layout.matrix[is.na(layout.matrix)] <- 0

  # for each node, keep the space for the central ones
  lapply(1:(nrow(opfit$frame)),function(i){
    loc <- which(layout.matrix==i,arr.ind = T)
    loc_i <- loc[c(nrow(loc)%/%2, (nrow(loc)%/%2+1)),]
    layout.matrix[loc_i] <<- i
    loc_0 <- loc[-c(nrow(loc)%/%2, (nrow(loc)%/%2+1)),]
    layout.matrix[loc_0] <<- 0
  })

  # set up the spaces for arrows from each node
  arr <- data.tree::ToDataFrameNetwork(opfit_tree)
  k <- nrow(opfit$frame)
  lapply(sort(as.numeric(unique(arr$from))),function(i){
    i <- as.numeric(i)
    loc_start <- which(layout.matrix==i,arr.ind = T)
    end_point <- arr$to[arr$from==i]
    loc_end <- lapply(end_point,function(x){
      which(layout.matrix==as.numeric(x),arr.ind = T)})

    k <<- k+1
    layout.matrix[unique(loc_start[,1]),
                  c( min(loc_end[[1]][,2]):(min(loc_start[,2])-1)
                  )] <<- k
    k <<- k+1
    layout.matrix[unique(loc_start[,1]),
                  c( max(loc_end[[2]][,2]):(max(loc_start[,2])+1)
                  )] <<- k
  })
  layout.matrix <- rbind(max(layout.matrix)+1,layout.matrix)
  layout.matrix <- cbind(layout.matrix,max(layout.matrix)+1,max(layout.matrix)+1)
  layout.matrix <- rbind(layout.matrix,max(layout.matrix)+1)

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

  # get split rules and split names
  sp_name <- 0
  sp <- list()
  opfit_tree$Do(function(node){
    sp_name <<- 1 + sp_name
    # node$name <-
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

  # change names
  nn_name <- 0
  opfit_tree$Do(function(node){
    nn_name <<- 1 + nn_name
    node$name <- sp[nn_name]
  }
  )


  # add split rule ::::::::::::::::::::::::::::::
  split_rule <- list()
  opfit_tree$Do(function(node){
    split_rule[[length(split_rule)+1]]<<-paste(node$splitname,node$splitlevels)
  })
  split_rule <- matrix(unlist(split_rule)[-stringr::str_which(unlist(split_rule),"NA")],byrow = T,ncol = 2)

  # makeplot()
  # plot.new()
  layout(mat = layout.matrix,
         heights = 1,
         widths = 1
  )

  i <- 0
  opfit_tree$Do(function(node){
    i <<- i+1
    par(mar = c(1,1,1,1))
    # node <- opfit_tree
    # use function hist to generate histogram in each group ------

    loc_treat <- which(node$data$treatment==1)
    loc_control <- which(node$data$treatment==0)

    wt <- (node$data$treatment - node$data[,propensity_score])/(node$data[,propensity_score] * (1-node$data[,propensity_score]))
    y <- node$data[,outcomevariable]*wt

    y[wt>0] <- (y[wt>0] / sum(wt[wt>0])) * length(y[wt>0])
    y[wt<0] <- (y[wt<0] / sum(wt[wt<0])) * length(y[wt<0])
    wt_mean <- c(mean(y[wt<0]),mean(y[wt>0]))

    control <- subset(node$data,treatment==0,select= outcomevariable)[,1]
    treat <- subset(node$data,treatment==1,select= outcomevariable)[,1]
    plotOutcomes(node$data$treatment, node$data[,outcomevariable],
                 node$data[,propensity_score],gamma=gamma,lambda=lambda)
  })


  lapply( 1:((max(layout.matrix)-2-nrow(opfit$frame))%/%2),function(k){
    plot.new()
    segments(0,0,0,0.4)
    segments(0,0.4,1,0.4)
    text(0.5,0.5,split_rule[k,1],cex = 0.85)

    plot.new()
    segments(0,0.4,1,0.4)
    segments(1,0.4,1,0)
    text(0.5,0.5,split_rule[k,2],cex = 0.85)
  } )
  cols <-  colorRampPalette(c("#FFD700", "#ff800e", "#ff4500", "#8B0000"))(100)
  plot.new()
  text(.5,.5,plot.title,font=2,cex=2)
  # plot.new()
  par(mar = c(1, 1, 2, 2), mgp=c(2,1,0))
  # s2dverification::ColorBar(seq(0,1,.01), cols)
  # image(y=2:30,z=t(2:30), col=color[2:30], axes=FALSE, main="Slope", cex.main=.8)
  image(y=seq(0,1,.01),z=t(seq(0,1,.01)), col=cols,  axes=FALSE, cex.main=.5)
  axis(4,cex.axis=0.5,mgp=c(0,.5,0))
  plot.new()
  text(.5,.5,"Red line is the difference between two groups weighted by inverse propensity score. Shadow is the confidence interval.",cex=1)
}
