# helper functions for LinCDE package

#' countIndex
#'
#' This function counts samples in each bin.
#'
#' @param yIndex vector of indices representing which bins the responses fall into.
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
#' @return The function returns the regularization parameter corresponding to the desired degrees of freedom.
#'
#' @export
dfToLambda = function(z, counts, splineDf, df, numberBin, penalty = NULL){
  if(sd(penalty) > 0){stop("currently only allow constant penalty.factors")}
  if(splineDf <= df){return (0)}
  if(df == 0){return (1e9)}

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
#' This function computes the importance scores for a LinCDE boosting model.
#'
#' @param model a LinCDE boosting model.
#' @param d the number of covariates.
#' @param n.trees the number of trees in the LinCDE boosting model.
#' @param var.names vector of variable names.
#'
#' @return The function returns a vector of covariate importances.
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
  return (sqrt(importance))
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

#' constructBasis
#'
#' This function construct response basis functions.
#'
#' @param splitMidPointY vector of representative response values in each bin after discretization.
#' @param y vector of responses.
#' @param splineDf the number of response basis functions.
#' @param basis a character or a function specifying sufficient statistics, i.e., spline basis. For \code{basis = "Gaussian"}, y, y^2 are used. For \code{basis = "nsTransform"}, transformed natural cubic splines are used. If \code{basis} is a function, it should take a vector of response values and output a basis matrix: each row stands for a response value and each column stands for a basis function.
#' @param newKnots if \code{TRUE}, quantiles of \code{y} are used as the knots for splines. If \code{FALSE}, quantiles of \code{splitMidPointY} are used as the knots for splines.
#' @param eps a non-negative scalar added to the \code{penalty} for computational stability. Default is 0.001.
#'
#' @return This function returns a list of values.
#' \itemize{
#' \item z: the spline basis matrix.
#' \item penalty: a vector of penalty.factors.
#' }
#'
#' @export
constructBasis = function(splitMidPointY, y = NULL, splineDf = 2, basis, newKnots = F, eps = 0.001){
  result = list()
  if(is.character(basis)){
    if(basis == "Gaussian"){
      z = cbind(splitMidPointY, splitMidPointY^2); result$z=z
      result$penalty = c(1,1)
    } else if (basis == "nsTransform"){
      if(newKnots){
        knotsNs = quantile(y, probs = seq(1, splineDf-1)/splineDf)
        boundaryKnotsNs = range(y)
      } else{
        knotsNs = quantile(splitMidPointY, probs = (1:(splineDf+1))/(splineDf+2))
        boundaryKnotsNs = knotsNs[c(1,splineDf+1)]
        knotsNs = knotsNs[-c(1,splineDf+1)]
      }
      zOrg = splines::ns(x=splitMidPointY, knots = knotsNs, intercept = FALSE, Boundary.knots = boundaryKnotsNs)
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = splines::ns(x=xSeq, knots = attributes(zOrg)$knots, intercept = attributes(zOrg)$intercept, Boundary.knots = attributes(zOrg)$Boundary.knots)
      derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
      OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN)
      svdOmegaN = svd(OmegaN)
      penalty = svdOmegaN$d/sum(svdOmegaN$d) * splineDf + eps # + eps(0.001) for robustness
      result$zTransformMatrix = svdOmegaN$u %*% diag(sqrt(penalty[splineDf]/penalty))
      z=zOrg%*%result$zTransformMatrix; attributes(z)=attributes(zOrg); result$z=z
      result$penalty = rep(1, splineDf)
    }
  } else if(is.function(basis)){
    result$z = basis(splitMidPointY)
    result$penalty = rep(1, dim(z)[2])
  }
  return (result)
}


