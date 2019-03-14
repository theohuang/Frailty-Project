###################################################################################################
### binary outcome model evaluation functions 
### with optional inverse probability of censoring weighting (IPCW) if the binary outcome 
### is a dichotomized version of a time-to-event outcome 
###################################################################################################
# calc.oe: function for calculating calibration (observed cases / expected cases)
# calc.auc: function for calculating ROC-AUC
# calc.brier: function for calculating Brier score (mean squared error)
# km.weights: function for getting Kaplan-Meier IPC weights
# cox.weights: function for getting Cox IPC weights
# perf.boot: function for getting 95% bootstrap confidence intervals for performance measures
# format.perf.table: function for formatting matrix of performance measures from perf.boot
# multiplot: function for plotting multiple ggplot graphs on the same page
# plot.scores: function for making stratified density plots of predicted probabilities 
# plot.oe: function for plotting expected vs observed probabilities by quantiles of predicted risk
###################################################################################################

### Written by Zoe Guan
### edited by Theo for frailty work
## Last updated: December 10, 2018

library(survival)
library(ggplot2)
library(ROCR)

##### function for calculating calibration (observed cases / expected cases)
### input
# scores: vector of predicted probabilities
# outcomes: vector of outcomes (1 for case, 0 for control)
# weights: vector of IPC weights; ignore if there is no censoring 
### output
# observed cases / expected cases
calc.oe = function(scores, outcomes, weights=NULL) {
  if (is.null(weights)) {
    weights = rep(1, length(scores))
  }
  ind = which(outcomes>=0 & scores>=0)
  return(sum(outcomes[ind]*weights[ind])/sum(scores[ind]))
}

##### function for calculating discrimination
### input
# scores: vector of predicted probabilities
# outcomes: vector of outcomes (1 for case, 0 for control)
# weights: vector of IPC weights; ignore if there is no censoring 
### output
# ROC-AUC 
calc.auc = function(scores, outcomes, weights=NULL) {
  if (is.null(weights)) {
    return(unlist(performance(prediction(scores, outcomes), measure="auc")@y.values))
  }
  ind.cases = which(outcomes==1 & scores>=0)
  ind.controls = which(outcomes==0 & scores>=0)
  # for each case-control pair, check if score of case is larger than score of control
  comparisons = as.vector(outer(scores[ind.cases], scores[ind.controls], ">"))
  # check for ties
  ties = as.vector(outer(scores[ind.cases], scores[ind.controls], "=="))/2
  w = as.vector(outer(weights[ind.cases], weights[ind.controls]))
  return(sum((comparisons+ties)*w)/sum(w))
}

##### function for calculating Brier score (mean squared error)
### input
# scores: vector of predicted probabilities
# outcomes: vector of outcomes (1 for case, 0 for control)
# weights: vector of IPC weights; ignore if there is no censoring 
### output
# Brier score
calc.brier = function(scores, outcomes, weights=NULL, root = FALSE) {
  if (is.null(weights)) {
    weights = rep(1, length(scores))
  }
  ind = which(outcomes>=0 & scores>=0)
  if(root == FALSE){
    return(mean((outcomes[ind]-scores[ind])^2 * weights[ind]))
  } else if(root == TRUE){
    return(sqrt(mean((outcomes[ind]-scores[ind])^2 * weights[ind])))
  }
}

##### function for getting Kaplan-Meier IPC weights
### input
# dataset: dataframe with observed event statuses and times 
# must have columns Time (minimum of censoring time and event time) and y (event status: 1 if event observed, 0 otherwise)
# t: cutoff time used for dichotomizing the time-to-event outcome
### output
# vector of weights with size equal to the number of rows in dataset
km.weights = function(dataset, t) {
  # estimate censoring distribution using Kaplan-Meier
  km = summary(survfit(Surv(Time, y==0) ~ 1, data=dataset))
  survest = stepfun(km$time, c(1, km$surv))
  # set each person's weight to their inverse probability of not being censored by min(Time, t)
  weights = 1/survest(pmin(dataset$Time, t))
  # if censored before time t, set weight to 0
  weights[which(dataset$y==0 & dataset$Time<t | is.na(dataset$Time))] = 0
  return(weights)
}


##### function for getting Cox IPC weights
### input
# dataset: dataframe with observed event statuses and times 
# must have columns Time (minimum of censoring time and event time), y (event status: 1 if event observed, 0 otherwise), and
# columns corresponding to the names in covars
# covars: vector of names of covariates to be included in Cox model for censoring
# t: cutoff time used for dichotomizing the time-to-event outcome
### output
# vector of weights with size equal to the number of rows in dataset
cox.weights = function(dataset, covars, t) {
  # estimate censoring distribution using Cox proportional hazards model
  fit.cens.cox = coxph(as.formula(paste("Surv(Time, y==0) ~", paste(covars, collapse=" + "))), data=dataset)
  cox.test = cox.zph(fit.cens.cox)
  if (cox.test$table[nrow(cox.test$table), "p"]<0.05) {
    print("Proportional hazards assumption does not hold.")
  }
  ### set each person's weight to their inverse probability of not being censored by min(Time, t)
  #===========
  # weights = 1/sapply(1:nrow(dataset), 
  #                    function(i) summary(survfit(fit.cens.cox, newdata=dataset[i,]), times=min(t, dataset$Time[i]))$surv)
  #===========
  # first, set everyone's weight to their inverse probability of not being censored by t
  surv.t = as.vector(summary(survfit(fit.cens.cox, newdata=dataset), times=t)$surv)
  weights = 1/surv.t
  # update weights for people who had an event before t
  unique.times = unique(dataset$Time[which(dataset$y==1 & dataset$Time<t)])
  surv.dist = summary(survfit(fit.cens.cox, newdata=dataset), times=unique.times)
  for (i in 1:length(unique.times)) {
    ind.i = which(dataset$y==1 & dataset$Time==unique.times[i])
    weights[ind.i] = 1/as.vector(surv.dist$surv[which(surv.dist$time==unique.times[i]), ind.i])
  }
  # if censored before time t, set weight to 0
  weights[which(dataset$y==0 & dataset$Time<t | is.na(dataset$Time))] = 0
  return(weights)
}



##### function for getting 95% bootstrap confidence intervals for performance measures
### input
# score.matrix: matrix of predicted probabilities where each column corresponds to a different model
# can be a vector of predicted probabilities if there is only one model to evaluate
# outcomes: vector of outcomes (1 for case, 0 for control)
# cens.dist: character string specifying the type of IPC weights to use
# "none" if there is no censoring/IPCW should not be used
# "km" if Kaplan-Meier IPC weights should be used (if cens.dist="km", then dataset and t cannot be NULL)
# "cox" if Cox IPC weights should be used (if cens.dist="cox", then dataset, covars, and t cannot be NULL)
# dataset: if cens.dist="km" or "cox", dataframe with observed event statuses and times 
# must have columns Time (minimum of censoring time and event time) and y (event status: 1 if event observed, 0 otherwise)
# if cens.dist="cox", must also have columns corresponding to the names in covars
# covars: if cens.dist="cox", vector of names of covariates to be included in Cox model for censoring
# t: if cens.dist="km" or "cox", cutoff time used for dichotomizing the time-to-event outcome
# model.names: vector of model names corresponding to the columns of score.matrix 
# nboot: number of bootstrap samples
# seed: seed for random number generator
### output
# 9xm matrix where m is the number of models (number of columns in score.matrix)
# rows 1-3 correspond to calibration in observed sample, 2.5th percentile of bootstrap calibration estimates, and
# 97.5th percentile of bootstrap calibration estimates
# rows 4-6 correspond to ROC-AUC in observed sample, 2.5th percentile of bootstrap AUC estimates, and
# 97.5th percentile of bootstrap ROC-AUC estimates
# rows 7-9 correspond to Brier score in observed sample, 2.5th percentile of bootstrap Brier score estimates, and
# 97.5th percentile of bootstrap Brier score estimates
# column i corresponds to the ith model
perf.boot = function(score.matrix, outcomes, cens.dist="none", dataset=NULL, covars=NULL, t=NULL, 
                     model.names=NULL, nboot=200, seed=1) {
  set.seed(seed)
  score.matrix = as.matrix(score.matrix)
  n.models = ncol(score.matrix)
  
  # matrices for storing bootstrap results
  oe = matrix(ncol=n.models, nrow=nboot+1)
  auc = matrix(ncol=n.models, nrow=nboot+1)
  brier = matrix(ncol=n.models, nrow=nboot+1)
  
  # get weights if using IPCW
  weights = NULL
  if (cens.dist=="km") {
    if (is.null(dataset) | is.null(t)) {
      stop("dataset and t cannot be NULL.")
    }
    weights = km.weights(dataset, t)
  } else if (cens.dist=="cox") {
    if (is.null(dataset) | is.null(t) | is.null(covars)) {
      stop("dataset, covars, and t cannot be NULL.")
    }
    weights = cox.weights(dataset, covars, t)
  }
  
  # evaluate performance in observed sample
  oe[1, ] = apply(score.matrix, 2, function(x) calc.oe(x, outcomes, weights))
  auc[1, ] = apply(score.matrix, 2, function(x) calc.auc(x, outcomes, weights))
  brier[1, ] = apply(score.matrix, 2, function(x) calc.brier(x, outcomes, weights, root = TRUE))
  
  for (i in 2:(nboot+1)) {
    # generate bootstrap sample
    ind.boot = sample(1:length(outcomes), length(outcomes), replace=T)
    score.matrix.boot = as.matrix(score.matrix[ind.boot, ])
    outcomes.boot = outcomes[ind.boot]
    # re-calculate weights for bootstrap sample if using IPCW
    weights.boot = NULL
    if (cens.dist=="km") {
      weights.boot = km.weights(dataset[ind.boot,], t)
    } else if (cens.dist=="cox") {
      weights.boot = cox.weights(dataset[ind.boot,], covars, t)
    }
    # evaluate performance in bootstrap sample
    oe[i, ] = apply(score.matrix.boot, 2, function(x) calc.oe(x, outcomes.boot, weights.boot))
    auc[i, ] = apply(score.matrix.boot, 2, function(x) calc.auc(x, outcomes.boot, weights.boot))
    brier[i, ] = apply(score.matrix.boot, 2, function(x) calc.brier(x, outcomes.boot, weights.boot, root = TRUE))
  }
  
  # calculate 95% CIs for performance measures using 2.5th and 97.5th percentiles of bootstrap estimates
  oe.summary = apply(oe, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  auc.summary = apply(auc, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  brier.summary = apply(brier, 2, function(x) c(x[1], quantile(x[-1], 0.025), quantile(x[-1], 0.975)))
  perf.matrix = as.matrix(rbind(oe.summary, auc.summary, brier.summary))
  rownames(perf.matrix) = c("oe", "oe.lo", "oe.hi",
                            "auc", "auc.lo", "auc.hi",
                            "brier", "brier.lo", "brier.hi")
  colnames(perf.matrix) = model.names
  if (is.null(colnames(perf.matrix)) & !is.null(colnames(score.matrix))) {
    colnames(perf.matrix) = colnames(score.matrix)
  } else if (is.null(colnames(perf.matrix))) {
    colnames(perf.matrix) = paste("model", 1:n.models)
  }
  return(perf.matrix)
}

##### function for formatting output of perf.boot
### input
# perf.matrix: 9xm matrix formatted in the same way as the output of perf.boot, where
# m is the number of models being evaluated
# digits: vector of 3 numbers corresponding to the number of decimal places to be displayed for 
# O/E, ROC-AUC, and Brier score respectively
### output
# mx3 matrix 
# row i corresponds to the ith model
# columns 1-3 correspond to character strings for O/E, ROC-AUC, and Brier score with 95% bootstrap CIs
format.perf.table = function(perf.matrix, digits=c(2, 2, 3)) {
  oe.matrix = matrix(round(perf.matrix[1:3, ], digits[1]), ncol=ncol(perf.matrix))
  auc.matrix = matrix(round(perf.matrix[4:6, ], digits[2]), ncol=ncol(perf.matrix))
  brier.matrix = matrix(round(perf.matrix[7:9, ], digits[3]), ncol=ncol(perf.matrix))
  formatted.matrix = cbind(apply(oe.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ),
                           apply(auc.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ),
                           apply(brier.matrix, 2, function(x) paste0(x[1], " (", x[2], ", ",  x[3], ")") ))
  colnames(formatted.matrix) = c("O/E", "AUC", "Root Brier Score")
  rownames(formatted.matrix) = colnames(perf.matrix)
  return(formatted.matrix)
}


##### function for plotting multiple ggplot graphs on the same page
# code from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


##### function for making stratified density plots of predicted probabilities 
### input
# score.matrix: matrix of predicted probabilities where each column corresponds to a different model
# can be a vector of predicted probabilities if there is only one model to evaluate
# groups: vector of values for stratification variable
# model.names: vector of model names corresponding to the columns of score.matrix 
# alpha.val: ggplot alpha parameter (transparency value)
# xlimits: vector of size 2 specifying x-axis limits
# colors: vector specifying colors for the strata
# legend.title: name of stratification variable 
# legend.pos: ggplot legend.position parameter (vector of size 2 specifying legend position)
# cols: multiplot cols parameter
# layout: multiplot layout parameter
# show.plot: TRUE if the plot should be displayed
### output 
# if score.matrix contains only one model, the function returns a ggplot object corresponding to
# a density plot of the model scores stratified by the grouping variable
# if score.matrix has m>1 columns, the function returns a list of m ggplot objects corresponding to 
# density plots of the model scores stratified by the grouping variable
plot.scores = function(score.matrix, groups, model.names=NULL, alpha.val=0.3, 
                       xlimits=c(0, 1), colors=NULL, legend.title=NULL, legend.pos=c(0.8, 0.5), 
                       cols=1, layout=NULL, show.plot=T) {
  
  score.matrix = as.matrix(score.matrix)
  
  # remove rows with NAs for the stratification variable
  ind.na = which(is.na(groups))
  if (length(ind.na)>0) {
    score.matrix = score.matrix[-ind.na, ]
    groups = groups[-ind.na]
    warning(paste0("Removed ", length(ind.na), " observations with NAs in the groups vector."))
  }
  
  # set plotting parameters that are unspecified
  if (is.null(legend.title)) {
    legend.title = "group"
  }
  if (is.null(model.names)) {
    model.names = paste("model", 1:ncol(score.matrix))
  } 
  if (is.null(colors) | length(colors)!=length(unique(groups))) {
    colors = 1:length(unique(groups))
  }
  
  # make dataframe for ggplot
  df = data.frame(group=factor(groups))
  # make list for storing plots for each model
  plot.list = vector("list", length = ncol(score.matrix))
  # create plot with legend for first model 
  df$scores = score.matrix[, 1]
  plot.list[[1]] = ggplot(df, aes(x=scores)) + geom_density(aes(fill=group), alpha=alpha.val) +
    xlab(model.names[1]) + xlim(xlimits) + theme(legend.position = legend.pos) + 
    scale_fill_manual(values=colors, name=legend.title)
  if (length(plot.list)==1) {
    if (show.plot) {
      multiplot(plotlist=plot.list, cols=cols, layout=layout)
    }
    return(plot.list[[1]])
  }
  # create plots without a legend for the other models
  for (i in 2:length(plot.list)) {
    df$scores = score.matrix[, i]
    plot.list[[i]] = ggplot(df, aes(x=scores)) + geom_density(aes(fill=group), alpha=alpha.val) +
      xlab(model.names[i]) + xlim(xlimits) + 
      scale_fill_manual(values=colors) + guides(fill=FALSE)
  }
  if (show.plot) {
    multiplot(plotlist=plot.list, cols=cols, layout=layout)
  }
  return(plot.list)
  
}


##### function for plotting expected vs observed probabilities by quantiles of predicted risk
### input
# score.matrix: matrix of predicted probabilities where each column corresponds to a different model
# can be a vector of predicted probabilities if there is only one model to evaluate
# outcomes: vector of outcomes (1 for case, 0 for control)
# weights: vector of IPC weights; ignore if there is no censoring 
# n.bins: number of quantiles to use for binning
# model.names: vector of model names corresponding to the columns of score.matrix 
# xlimits: vector of size 2 specifying x-axis limits
# ylimits: vector of size 2 specifying y-axis limits
# bar.width: geom_errorbar width parameter
# threshold: minimum proportion of observed cases for the plotted bins; bins that do not meet the threshold 
# will not be plotted
# cols: multiplot cols parameter
# layout: multiplot layout parameter
# show.plot: TRUE if the plot should be displayed
### output 
# if score.matrix contains only one model, the function returns a ggplot object corresponding to
# a plot of expected vs observed probabilities for the specified quantiles with Wald confidence intervals
# if score.matrix has m>1 columns, the function returns a list of m ggplot objects corresponding to 
# plots of expected vs observed probabilities for the specified quantiles with Wald confidence intervals
plot.oe = function(score.matrix, outcomes, weights=NULL, n.bins=10, model.names=NULL,
                   xlimits=c(0, 1), ylimits=c(0, 1), bar.width=0.03, threshold=0.001,
                   cols=1, layout=NULL, show.plot=T) {
  
  score.matrix = as.matrix(score.matrix)
  
  # make list for storing plots for each model
  plot.list = vector("list", length = ncol(score.matrix))
  
  # set plotting parameteres that are unspecified
  if (is.null(weights)) {
    weights = rep(1, length(outcomes))
  }
  if (is.null(model.names)) {
    model.names = paste("model", 1:ncol(score.matrix))
  }
  
  for (j in 1:ncol(score.matrix)) {
    
    # get scores for model j
    scores = score.matrix[, j]
    ind = which(outcomes>=0 & scores>=0)
    scores = scores[ind]
    outcomes = outcomes[ind]
    weights = weights[ind]
    
    # sort scores into quantiles
    bins = cut(scores, labels = 1:n.bins, 
               breaks = quantile(scores, probs = seq(0, 1, length = n.bins+1), na.rm = TRUE), 
               include.lowest = TRUE)
    
    # compute O/E in each quantile
    df = data.frame(bin=1:n.bins, prop.obs=NA, prop.pred=NA, oe.lo=NA, oe.hi=NA)
    for (i in 1:n.bins) {
      ind.i = which(bins==i)
      size.i = length(ind.i)
      scores.i = scores[ind.i]
      outcomes.i = outcomes[ind.i]
      weights.i = weights[ind.i]
      df$prop.pred[i] = sum(scores.i)/size.i
      df$prop.obs[i] = sum(outcomes.i*weights.i)/size.i
      # get bounds for 95% Wald CI
      df$oe.lo[i] = pmax(0, df$prop.obs[i] - qnorm(0.975)*sqrt(df$prop.obs[i]*(1-df$prop.obs[i])/size.i))
      df$oe.hi[i] = pmin(1, df$prop.obs[i] + qnorm(0.975)*sqrt(df$prop.obs[i]*(1-df$prop.obs[i])/size.i))
    }
    
    df2 = df[which(df$prop.obs>=threshold), ]
    if (nrow(df2) < nrow(df)) {
      warning(paste0("Omitting ", nrow(df)-nrow(df2), " bin(s) with observed proportion falling below ", threshold, "."))
    }
    
    plot.list[[j]] = ggplot(df2, aes(x=prop.pred, y=prop.obs)) + 
      geom_point() +
      geom_errorbar(aes(ymin=oe.lo, ymax=oe.hi), width=bar.width) +
      geom_abline(linetype = 2, color = 'grey50') +
      coord_cartesian(xlim=xlimits, ylim=ylimits) + 
      labs(x="Predicted Probability", y="Observed Probability") +
      ggtitle(model.names[j])
    
  }
  
  if (show.plot) {
    multiplot(plotlist=plot.list, cols=cols, layout=layout)
  }
  
  if (length(plot.list)==1) {
    return(plot.list[[1]])
  } else {
    return(plot.list)
  }
  
}




