# helper functions for LinCDE package

#' countIndex
#'
#' This function counts the samples in each bin.
#'
#' @param yIndex vector of indices representing which bins the responoses fall into.
#' @param numberBin the number of bins for response discretization.
#'
#' @return The function returns a vector indicating the number of responses in each bin, of length \code{numberBin}.
#'
#' @export
countIndex = function(yIndex, numberBin) {
  counts = rep(0, numberBin)
  temp = table(yIndex)
  counts[as.numeric(names(temp))] = temp
  return (counts)
}

#' dfToLambda
#'
#' This function computes the regularization parameter corresponding to the given degrees of freedom.
#'
#' @param z sufficient statistics matrix, of dimension \code{numberBin} x \code{splineDf}.
#' @param counts vector indicating the number of responses in each bin, of length \code{numberBin}.
#' @param splineDf the number of sufficient statistics.
#' @param df degrees of freedom.
#' @param numberBin the number of bins for response discretization.
#' @param penalty separate penalty factors applied to each coefficient, of length \code{splineDf}.
#'
#' @return The function returns the regularization parameter.
#'
#' @export
dfToLambda = function(z, counts, splineDf, df, numberBin, penalty = NULL){
  if(var(penalty) > 0){stop("currently only allow constant penalty.factors")}
  if(splineDf <= df){return (0)}
  if(df == 0){return (1e9)}

  # estimate the weight
  # model = suppressWarnings(glm(counts ~ z, family = poisson()))
  # weight = model$fitted.values/numberBin
  #
  # Hessian = 2 * (t(z) %*% diag(weight) %*% z - t(z) %*% weight %*% weight %*% z * numberBin / sum(counts))
  # svdHessian = svd(Hessian)
  # Omega = diag(1/sqrt(svdHessian$d+1e-3)) %*% t(svdHessian$u) %*% diag(penalty) %*% svdHessian$u %*% diag(1/sqrt(svdHessian$d+1e-3))
  # eigenVal = svd(Omega)$d
  #
  # lambdaMax = sum(1/eigenVal)/df; lambdaMin = 0; lambda = (lambdaMax+lambdaMin)/2; value =  sum(1/(1+eigenVal*lambda))
  # while(abs(value - df) > 1e-3){
  #   if(value > df){lambdaMin = lambda; lambda = (lambdaMax+lambdaMin)/2; value = sum(1/(1+eigenVal*lambda))}
  #   else {lambdaMax = lambda; lambda = (lambdaMax+lambdaMin)/2; value = sum(1/(1+eigenVal*lambda))}
  # }
  #
  # lambda/2

  model = suppressWarnings(glm(counts~z, family = poisson()))
  prob = model$fitted.values/sum(counts)
  Hessian = (t(z) %*% diag(prob) %*% z - t(z) %*% prob %*% prob %*% z) * sum(counts)
  eigenVal = svd(Hessian)$d
  lambdaMax = sum(eigenVal)/df; lambdaMin = 0; lambda = (lambdaMax+lambdaMin)/2
  value = sum(eigenVal/(eigenVal+lambda))
  while(abs(value - df) > 1e-3){
    if(value > df){lambdaMin = lambda; lambda = (lambdaMax+lambdaMin)/2; value = sum(eigenVal/(eigenVal+lambda))}
    else {lambdaMax = lambda; lambda = (lambdaMax+lambdaMin)/2; value = sum(eigenVal/(eigenVal+lambda))}
  }
  lambda/numberBin
}


#' importanceScore
#'
#' This function computes the importance score for a LinCDE boosting model.
#'
#' @param tree a LinCDE boosting model.
#' @param d the number of covariates.
#' @param n.trees the number of trees in the LinCDE boosting model.
#'
#' @return The function returns a length \code{d} vector of covariate importances.
#'
#' @export
importanceScore = function(model, d, n.trees = 1, var.names = NULL){
  if(is.null(var.names)){var.names = paste("X", seq(1,d), sep = "")}
  importance = rep(0, d); names(importance) = var.names
  for(l in 1:n.trees){
    currTree = model[[l]]
    index = sort(unique(currTree$SplitVar))
    importance[index[-1]] = importance[index[-1]] + aggregate(x = currTree$ErrorReduction, by = list(currTree$SplitVar), FUN = sum)$x[-1]
  }
  return (importance)
}

#' constructSplitPoint
#'
#' This function computes the list of candidate covariate splits.
#'
#' @param X input matrix, of dimension nobs x nvars; each row represents an observation vector.
#' @param numberSplit split numbers for each covariate (length nvars). Each variable's range is divided into \code{numberSplit-1} intervals containing approximately the same number of observations.
#'
#' @return The function returns a candidate split list (length nvars). Each element is a vector corresponding to a certain variable's candidate splits (length of \code{numberSplit - 1} or 1 + the number of unique values).
#'
#' @export
constructSplitPoint = function(X, numberSplit){
  # preprocess
  d = dim(X)[2]
  if(length(numberSplit) == 1){numberSplit = rep(numberSplit,d)}

  # construct candidate splits
  splitPoint = list()
  for(i in 1:d){
    # unique values
    numberUnique = length(unique(X[,i]))
    if(numberUnique < (numberSplit[i]-1)){
      valueUnique = sort(unique(X[,i]), decreasing = FALSE)
      splitPoint[[i]] = c(valueUnique[1],(valueUnique[-numberUnique] + valueUnique[-1])/2, valueUnique[numberUnique])
    } else {
      # quantiles; remove potential duplicates
      splitPoint[[i]] = unique(quantile(X[,i], probs = seq(0,1,length.out = numberSplit[i]))[-1])
    }
    names(splitPoint[[i]]) = NULL
  }
  return (splitPoint)
}


