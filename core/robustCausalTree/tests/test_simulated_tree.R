#rm(list=ls())
# save(edurose_mediation_20181126,file = "data/edurose_mediation_20181126.rda",compress = "xz")

#library("htetree")
# load("data/edurose_mediation_20181126.rda")

# construct the simulated data based on Athey's data
#install.packages("data.table",repos = NULL,
#                 type = "source")
#library("data.table")
#install.packages("causalTree",
#                 repos = "https://jiahui1902.github.io/drat/",
#                 type = "source")
# install.packages("htetree",
#                  repos = "https://jiahui1902.github.io/drat/",
#                  type = "source")
# library(causalTree)

# library(rpart)
# library(htetree)
# hte_causalTree(outcomevariable="outcome",
#                data=data.frame("confounder"=c(0, 1, 1, 0, 1, 1),
#                                "treatment"=c(0,0,0,1,1,1),
#                                "prop_score"=c(0.4, 0.4, 0.5, 0.6, 0.6, 0.7),
#                                "outcome"=c(1, 2, 2, 1, 4, 4)),
#                treatment_indicator = "treatment",
#                ps_indicator = "prop_score", covariates = "confounder")
# # causalTree(y~ x1 + x2 + x3 + x4, data = simulation.1,
# #     treatment = simulation.1$treatment,
# #     split.Rule = "CT", cv.option = "CT", split.Honest = TRUE, cv.Honest = TRUE,
# #     split.Bucket = F, xval = 5,
# #     cp = 0, minsize = 20, propensity = 0.5)
# data("simulation.1")
# # estimate the propensity score
# fit <- glm(treatment~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,
#            data=simulation.1,
#            family = "binomial")
# simulation.1$ps_score <- predict(fit,type = "response")
# linear_terms <- paste0("x",1:10)
# 
# # estimate the model with our package
# set.seed(1)
# lb <- c(paste0("var",1:10),"propensity score")
# names(lb) <- c(paste0("x",1:10),"ps_score")
# 
# fit_drawplot <- htetree::hte_ipw(outcomevariable = 'y',
#                minsize=20,crossvalidation = 40,negative = TRUE,
#                data = simulation.1,
#                ps_indicator = "ps_score",
#                covariates = c(linear_terms, "ps_score"),
#                drawplot = TRUE,treatment_indicator = "treatment",
#                # no_indicater = '_IPW_simulation',
#                legend.x = 0.1,legend.y = 0.25,varlabel = lb)
# 
# fit_noplot <- htetree::hte_ipw(outcomevariable = 'y',
#                minsize=20,crossvalidation = 40,negative = TRUE,
#                data = simulation.1,
#                ps_indicator = "ps_score",
#                covariates = c(linear_terms, "ps_score"),
#                drawplot = TRUE,treatment_indicator = "treatment",
#                # no_indicater = '_IPW_simulation',
#                legend.x = 0.1,legend.y = 0.25,varlabel = lb)

# hte_plot_line(model = xxx,data = simulation.1,
#               treatment_indicator = "treatment",
#               outcomevariable = 'y',
#               propensity_score = "ps_score",gamma = 0.5,lambda = 0.5)
# hte_plot_line(model = xxx,data = simulation.1,
#               treatment_indicator = "treatment",
#               outcomevariable = 'y',
#               propensity_score = "ps_score")

# library(htetree)
# hte_plot(model = fit_drawplot,data = simulation.1,
#          treatment_indicator = "treatment",
#          outcomevariable = 'y',
#          propensity_score = "ps_score")
#
# hte_plot_line(model = fit_drawplot,data = simulation.1,
#               treatment_indicator = "treatment",
#               outcomevariable = 'y',
#               propensity_score = "ps_score")


# ps_indicator = 'ps_score'
# covs <- c("x1", "x2", "x3",
#           "x4", "x5", "x6", "x7",
#           "x8", "x9", "x10", "ps_score")
# xxx2 <- hte_matchinginleaves(outcomevariable = 'y',
#                              data = simulation.1,
#                              drawplot = TRUE,
#                              ps_indicator = "ps_score",
#                              treatment_indicator = "treatment",
#                              covariates=covs,
#                              con.num=4)

## Dynamic visualization for IPW model (using Shiny)
# THIS IS NOT SHOWING UP (shows a blank page)
# The runDynamic function runs the visualization without saving any of the files
# runDynamic(fit_drawplot, simulation.1,
#            outcomevariable = "y", treatment_indicator = "treatment",
#            propensity_score="ps_score")

# The files for runDynamic are saved in the temporary directory
# The files can be cleared manually using the clearTemp() function, or will automatically be cleared when you close R
# clearTemp()
