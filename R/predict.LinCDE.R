#' predict.LinCDE
#'
#' This function makes predictions from a LinCDE model.
#' @method predict LinCDE
#'
#' @param object a LinCDE model.
#' @param X input matrix for prediction, of dimension nobs x nvars; each row represents an observation vector.
#' @param y response vector for prediction, of length nobs.
#' @param splitPointYTest vector of response break points. Default is NULL.
#' @param densityOnly a logical value. If \code{TRUE}, only return predicted conditional densities \code{density}.
#' @param ... other parameters.
#'
#' @return The function returns a list of values.
#' \itemize{
#' \item density: conditional density predictions at \code{y}, vector of length nobs. If \code{y = NULL}, return \code{NA}.
#' \item testLogLikelihood: average test log-likelihood at \code{y}. If \code{y = NULL}, return \code{NA}.
#' \item cellProb: cell probability prediction matrix, of dimension nobs x number of discretization bins; each row represents a cell probability vector.
#' \item yDiscretized: vector of discretized responses/response residuals, of length number of discretization bins. If \code{splitPointYTest = NULL}, discretized responses/response residuals from the LinCDE model \code{object} are used.
#' }
#'
#' @export
predict.LinCDE = function(object, ..., X = NULL, y = NULL, splitPointYTest = NULL, densityOnly = TRUE){
  # pre-process and initialization
  result = list()
  d = length(object$var.names)
  if((object$depth == 0) && is.null(X)){
    if(!is.null(y)){X = matrix(runif(length(y)), ncol = 1); colnames(X) = "pseudoX"} else {X = matrix(runif(1), ncol = 1); colnames(X) = "pseudoX"}
  }
  if("data.frame" %in% class(X)){
    if(sum(!lapply(X, class) %in% c("numeric", "integer")) > 0){stop("currently only numeric covariates are allowed")}
    X = as.matrix(X)}
  if(is.null(dim(X))){
    if(d == 1){X = matrix(X, ncol = 1)} # one covariate
    if(d > 1){X = matrix(X, nrow = 1)} # one data point
  }
  n = dim(X)[1]
  if(!is.null(colnames(X))){
    if(min(colnames(X) %in% object$var.names) == 0)
      stop("please input all covariates")
    X = X[,object$var.names, drop = FALSE]
  } else {
    if(dim(X)[2] != d)
      stop("please input all covariates in the order of the training model's var.names")
    colnames(X) = object$var.names
    # print("no input var.names and var.names from the LinCDE.boost model are used")
  }
  if(sum(!apply(X, 2, class) %in% c("numeric", "integer")) > 0){stop("currently only numeric covariates are allowed")}
  if(!is.null(y)){full = TRUE} else {full = FALSE; y = rnorm(n)}
  z = object$z; splineDf = dim(z)[2]; depth = object$depth; basis = object$basis
  centering = object$centering; centeringMethod = object$centeringMethod

  # centering
  if(centering){
    yMean = predict(object$centeringModel, newdata = as.data.frame(X))
    y = y - yMean; result$yMean = yMean
  }
  result$centering = centering

  splitMidPointY = object$splitMidPointY; h = splitMidPointY[2] - splitMidPointY[1]; numberBin = length(splitMidPointY); splitPointY = c(splitMidPointY - h/2, splitMidPointY[numberBin] + h/2)
  if(!is.null(splitPointYTest)){
    numberBin = length(splitPointYTest)-1; splitPointY = splitPointYTest
    splitMidPointY = (splitPointYTest[2:(numberBin+1)] + splitPointYTest[1:numberBin])/2; h = splitMidPointY[2] - splitMidPointY[1]
    if(is.character(basis)){
      if(basis == "Gaussian"){
        z = cbind(splitMidPointY, splitMidPointY^2)
      } else if(basis == "nsTransform"){
        z = splines::ns(x=splitMidPointY, knots = attributes(object$z)$knots, intercept = attributes(object$z)$intercept, Boundary.knots = attributes(object$z)$Boundary.knots)
        z = z %*% object$zTransformMatrix
      }
    } else if(is.function(basis)){
      z = basis(splitMidPointY)
    }
  }

  # prior
  if(!is.null(object$prior)){
    cellProb0 = t(matrix(object$prior(X = X[rep(seq(1,n), rep(numberBin, n)),,drop = FALSE], y = rep(splitMidPointY, n)), nrow = numberBin))
    if(n == 1){cellProb0 = cellProb0 / sum(cellProb0)
    } else {cellProb0 = diag(1/apply(cellProb0, 1, sum)) %*% cellProb0}
  } else{cellProb0 = matrix(1/numberBin, nrow = n, ncol = numberBin)}

  # response discretization
  if(full){yIndex = as.numeric(cut(y, c(-Inf,splitPointY[-c(1,numberBin+1)],Inf)))}

  if(object$depth == 0){
    cellProb0 = cellProb0 %*% diag(c(exp(z %*% object$trees)))
    if(n == 1){cellProb0 = cellProb0 / sum(cellProb0)
    } else {cellProb0 = diag(1/apply(cellProb0, 1, sum)) %*% cellProb0}
    object$trees = NULL
  }

  if(is.null(object$shrinkage)){shrinkage = 1} else {shrinkage = object$shrinkage}
  trees = object$trees; numTree = length(trees)
  # initialization
  beta = matrix(0, nrow = n, ncol = splineDf)

  if(numTree == 0){
    if(densityOnly){return(cellProb0[cbind(seq(1,n), yIndex)]/h)}
    result = list()
    if(full){
      result$density = cellProb0[cbind(seq(1,n), yIndex)]/h; result$testLogLikelihood = mean(log(result$density))}
    else{result$density = NA; result$testLogLikelihood = NA}
    result$cellProb = cellProb0; result$yDiscretized = splitMidPointY
    return(result)
  }

  # transfer LinCDE trees to gbm trees
  gbmTrees = list()
  for(i in 1:numTree){
    currTree = trees[[i]]
    temp = list()
    temp[[1]] = as.integer(currTree[,1] - 1) # index - 1
    temp[[2]] = currTree[,2]
    indexLeaf = which(currTree[,1] <= 0)
    temp[[2]][indexLeaf] = seq(1, length(indexLeaf)) # predict the index of the terminal nodes
    temp[[3]] = as.integer(currTree[,3] - 1) # index - 1
    temp[[4]] = as.integer(currTree[,4] - 1) # index - 1
    temp[[5]] = as.integer(currTree[,5] * 0 - 1) # missing splits, set as -1
    temp[[6]] = currTree[,5]
    temp[[7]] = currTree[,6]
    temp[[8]] = currTree[,2]
    gbmTrees[[i]] = temp
  }

  # constrcut gbm trees
  # create the basics of a gbm.object
  data = data.frame(y, X)
  if(n < 100){data = rbind(data, data[rep(1,100-n),,drop = FALSE] + matrix(rnorm((100-n)*(d+1)),ncol = (d+1)))}
  gbmObject = suppressWarnings(gbm::gbm(y~., data = data[sample(dim(data)[1], 100),,drop = FALSE], distribution = "gaussian", cv.folds = 0, n.trees = 1, interaction.depth = depth)) #2^interaction.depth
  gbmObject$trees = gbmTrees
  gbmObject$n.trees = numTree
  # compute the matrix of memberships (n by numTree)
  nodeMembership = suppressWarnings(matrix(gbm::predict.gbm(gbmObject, newdata = data[1:n, ,drop = FALSE], n.trees = seq(1:numTree), single.tree = TRUE), ncol = numTree))

  # prediction
  testLogLikelihood = numeric(0)
  for(i in 1:numTree){
    currTree = trees[[i]]
    indexLeaf = which(currTree$SplitVar == 0)
    betaTemp = currTree[indexLeaf,7:(splineDf+6)]
    beta = beta + betaTemp[nodeMembership[,i],] * shrinkage

    cellProb = exp(as.matrix(beta) %*% t(z)) * cellProb0 # homogeneous prior
    if(full){
    temp = mean(log(cellProb[cbind(seq(1,n), yIndex)]/h)) - mean(log(apply(cellProb, 1, sum)))
    testLogLikelihood = c(testLogLikelihood, temp)
    }
  }
  cellProb = exp(as.matrix(beta) %*% t(z)) * cellProb0 # homogeneous prior
  # compute normalizing constants
  a0 = -log(apply(cellProb, 1, sum)) - log(h)
  cellProb = cellProb * exp(a0) * h
  # compute conditional density
  result = list()
  if(densityOnly){return(cellProb[cbind(seq(1,n), yIndex)]/h)}
  if(full){
    result$density = cellProb[cbind(seq(1,n), yIndex)]/h; result$testLogLikelihood = mean(log(result$density))}
  else{result$density = NA; result$testLogLikelihood = NA}
  result$cellProb = cellProb; result$yDiscretized = splitMidPointY; result$testLogLikelihoodHistory = testLogLikelihood
  return(result)
}

#' LinCDETreePredict
#'
#' This function conducts prediction of a LinCDE tree.
#'
#' @param X input matrix for prediction, of dimension nobs x nvars; each row represents an observation vector.
#' @param tree a LinCDE tree.
#' @param currentRow an integer specifying the node position in \code{tree}. If the \code{currentRow}-th row of the \code{tree} is not a terminal node thus a split, the function \code{LinCDETreePredict} is applied to the left and right children of the split. If the \code{currentRow}-th row of the \code{tree} is a terminal node, predictions are computed.
#'
#' @return The function returns a list of values.
#' \itemize{
#' \item membership: length n vector: each element represents the row (terminal node) that a sample falls into.
#' }
#'
#' @export
LinCDETreePredict = function(X, tree, currentRow = 1){
  if(is.null(dim(X))){X = matrix(X, nrow = 1)}
  if(tree$SplitVar[currentRow] == 0){
    membership = rep(currentRow, dim(X)[1]); return(membership)
  }
  indexLeft = which(X[,tree$SplitVar[currentRow]] <= tree$SplitCodePred[currentRow])
  if(length(indexLeft) > 0){
    membershipLeft = LinCDETreePredict(X = X[indexLeft,,drop = FALSE], tree = tree, currentRow = tree$LeftNode[currentRow])
  } else {
    membership = LinCDETreePredict(X = X, tree = tree, currentRow = tree$RightNode[currentRow])
    return(membership)
  }
  if(length(indexLeft) < dim(X)[1]){
    membershipRight = LinCDETreePredict(X = X[-indexLeft, ,drop = FALSE], tree = tree, currentRow = tree$RightNode[currentRow])
    membership = rep(0, dim(X)[1])
    membership[indexLeft] = membershipLeft
    membership[-indexLeft] = membershipRight
    return(membership)
  } else {
    return(membershipLeft)
  }
}
