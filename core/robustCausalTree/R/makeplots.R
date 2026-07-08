#' Visualize Causal Tree and the Estimated Results
#'
#' An intermediate function used for plotting
#'
#' @param negative a logical value indicating whether we expect the
#' treatment effect to be negative. The default is set as FALSE.
#' @param opfit. tree structure generated from causal tree algorithm.
#' @param trainset a data frame only containing the variables used in
#' the model and missings values are listwise deleted.
#' @param covariates a vector of column names of all
#' covariates (linear terms andpropensity score).
#' @param outcomevariable a character representing the column name
#' of the outcome variable.
#' @param data. a data frame containing the variables in the model.
#' @param hte_effect_setup a empty list to store the adjusted treatment
#' effect.
#' @param varlabel a named vector containing variable labels.
#' @param maintitle a character string indicating the main title displayed
#' when plotting the tree and results. The default is set as
#' "Heterogeneous Treatment Effect Estimation".
#' @param legend.x,legend.y x and y coordinate to position the legend.
#' The default is set as (0.08, 0.25).
#' @param ... further arguments passed to or from other methods.
#'
#' @returns A plot visualizing the tree and estimated treatment effect in each
#' node.


# Making plots function ---------------------------------------------------
# plot the tree model (heatmap)
makeplots <- function(negative,opfit.=opfit,trainset,covariates,outcomevariable,data.=data,
                      hte_effect_setup,varlabel,maintitle,legend.x=0.8,legend.y=0.25,...
){
  # opfit. <- opfit.
  # outcomevariable <- outcomevariable
  if(nrow(opfit.$frame)>1){
    # pdf(paste0(gph,outcomevariable,"_causaltree",no_indicater,".pdf")) # the name could be changed
    # define the color palette to be used in this head map
    # rbPal <- colorRampPalette(c('blue','red')) # the color palette could be changed
    # rbPal <- colorRampPalette(c('gold','royal blue')) # the color palette could be changed
    rbPal <- colorRampPalette(c('golden rod','blue')) # the color palette could be changed
    # get the color for each of the leaves
    {if(negative==FALSE){
      od <- 1:length(opfit.$frame$yval)
    }else{
      od <- length(opfit.$frame$yval):1
    }}

    cols <- rbPal(length(opfit.$frame$yval))[od]
    # change variables to its labels
    if(!is.null(varlabel)){
      opfit.$frame$var <- as.character(opfit.$frame$var)
      # *change here
      # get the label names and change the name
      opfit.$frame$var[opfit.$frame$var%in%names(varlabel)] <- na.omit(varlabel[match(opfit.$frame$var,names(varlabel))])
    }else{
      opfit.$frame$var <- as.character(opfit.$frame$var)
    }

    # make plots
    if(length(ttable)>0){
      rpart.plot::prp(opfit., # tree model estimated from causal tree
                      type = 1, # could be delted. draw a split label at each split and a node label at each leaf
                      # extra =100, # display extra information at the nodes, could be deleted
                      nn = FALSE, # could be deleted. display the node numbers and the default value is FALSE
                      box.palette = cols, # specify the colors in each node
                      suffix = paste0(ttable$star,'\n(',round(ttable$se,3),
                                      ')\n',round(ttable$SampleSize,1),"%"), # specify the label in each nodes
                      yesno = 2, # specify if "yes" and "no" are shown
                      yes.text = "y", # specify the text for "yes"
                      no.text = "n", # specify the text for "no"
                      trace = TRUE,
                      varlen = 0,
                      col = "white",
                      # branch.col = cols, # specify branch colors
                      # split.col = cols, # specify the color of the split label text
                      # nn.col = cols,# specify the color of the node numbers
                      main=maintitle # specify the title displayed at the top
      )
    }else{
      rpart.plot::prp(opfit., # tree model estimated from causal tree
                      type = 1, # could be delted. draw a split label at each split and a node label at each leaf
                      extra =100, # display extra information at the nodes, could be deleted
                      nn = FALSE, # could be deleted. display the node numbers and the default value is FALSE
                      box.palette = cols, # specify the colors in each node
                      yesno = 2, # specify if "yes" and "no" are shown
                      yes.text = "y", # specify the text for "yes"
                      no.text = "n", # specify the text for "no"
                      trace = TRUE,
                      varlen = 0,
                      col = "white",
                      # branch.col = cols, # specify branch colors
                      # split.col = cols, # specify the color of the split label text
                      # nn.col = cols,# specify the color of the node numbers
                      main=paste0("Heterogeneous Treatment Effects: ",
                                  attr(edurose_mediation_20181126,"var.labels")[as.numeric(na.omit(match(outcomevariable,colnames(edurose_mediation_20181126))) )]) # specify the title displayed at the top
      )
    }
    # legend and notes
    mtext(c("\nText in Squares: HTE & sample size (%); Color of Squares: Blue:largest treatment effects & Yellow:smallest treatment effects; Number in Parentheses: Standard Error"),
          cex=0.5)
    legend(legend.x,legend.y, legend = c("HTE","sample size (%)",
                                         "largest treatment effects",
                                         "smallest treatment effects"),
           title = "Numbers/Colors:",cex = .6,box.lty=0,#pch = c(0,0,0,0),
           fill=c("white","white","blue","golden rod"))

    # alternatively, use the following command to make the plot
    # rpart.plot(opfit.,roundint = FALSE,type = 4,col = cols,main = paste0("Heterogeneous Treatment Effects: ",outcomevariable))
    # dev.off()

    # results[[paste0(outcomevariable,"_color_",length(results)+1)]] <<- opfit.
    cols <- cols[order(order(opfit.$frame$yval))]
    # results[[paste0(outcomevariable,"_cols_",length(results)+1)]] <<- cols
  }
}
