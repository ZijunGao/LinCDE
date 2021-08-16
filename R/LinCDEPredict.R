#' LinCDE prediction
#'
#' This function makes predictions from a LinCDE boosting model.
#'
#' @param X input matrix for prediction, of dimension nobs x nvars; each row represents an observation vector.
#' @param y response vector for prediction, of length nobs.
#' @param trees a LinCDE boosting model.
#' @param splitPointYTest vector of response break points. Default is NULL.
#'
#' @return The function returns a list of values.
#' \itemize{
#' \item density: conditional density predictions at \code{y}, vector of length nobs. If \code{y = NULL}, return \code{NA}.
#' \item testLogLikelihood: average test log-likelihood at \code{y}. If \code{y = NULL}, return \code{NA}.
#' \item cellProb: cell probability prediction matrix, of dimension nobs x number of discretization bins; each row represents a cell probability vector.
#' \item yDiscretized: vector of discretized responses/response residuals, of length number of discretization bins. If \code{splitPointYTest = NULL}, discretized responses/response residuals from the LinCDE boosting model \code{trees} are used.
#' }
#'
#' @export
LinCDEPredict = function(X, y = NULL, trees, splitPointYTest = NULL){
  # pre-process
  if(class(X)[1] == "data.frame"){X = as.matrix(X)}
  if(is.null(dim(X))){X = matrix(X, nrow = 1)} # predict on one observation
  n = dim(X)[1]; d = dim(X)[2];
  if(!is.null(y)){full = TRUE} else {full = FALSE; y = rnorm(n)}
  z = trees$z; k = dim(z)[2]; depth = trees$depth; type = trees$type
  augmentation = trees$augmentation; augmentationMethod = trees$augmentationMethod

  if(augmentation){
    if(augmentationMethod == "linearRegression"){
      yMean = cbind(1, X) %*% trees$augmentationModel$coefficients # order!!!
    } else if(augmentationMethod == "randomForest"){
      yMean = predict(trees$augmentationModel, newdata = X);
    }
    y = y - yMean
  }
  splitMidPointY = trees$splitMidPointY; h = splitMidPointY[2] - splitMidPointY[1]; numberBin = length(splitMidPointY); splitPointY = c(splitMidPointY - h/2, splitMidPointY[numberBin] + h/2); splitPointY = c(splitPointY, Inf);
  if(!is.null(splitPointYTest[1])){
    if(type == "Gaussian"){
      numberBin = length(splitPointYTest)-1
      splitMidPointY = (splitPointYTest[2:(numberBin+1)] + splitPointYTest[1:numberBin])/2
      z = cbind(splitMidPointY, splitMidPointY^2)
      splitPointY = c(splitPointYTest, Inf)
      penalty = c(1,1); k = 2; h = splitMidPointY[2] - splitMidPointY[1]
    } else if (type == "bsTransform"){
      numberBin = length(splitPointYTest)-1
      knotsBs = quantile(splitMidPointY, probs = (0:(k-2))/(k-2))
      splitMidPointY = (splitPointYTest[2:(numberBin+1)] + splitPointYTest[1:numberBin])/2
      z = splines::bs(x=splitMidPointY, df = k, knots = knotsBs[-c(1,k-1)], degree = 3, intercept = FALSE, Boundary.knots = knotsBs[c(1,k-1)])
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = splines::bs(x=xSeq, knots = knotsBs[-c(1,k-1)], intercept = FALSE, Boundary.knots = knotsBs[c(1,k-1)])
      derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
      OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN)
      svdOmegaN = svd(OmegaN)
      z = z %*% svdOmegaN$u
      splitPointY = c(splitPointYTest, Inf)
      h = splitMidPointY[2] - splitMidPointY[1]
    } else if (type == "nsTransform"){
      numberBin = length(splitPointYTest)-1
      knotsNs = quantile(splitMidPointY, probs = (1:(k+1))/(k+2))
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = splines::ns(x=xSeq, knots = knotsNs[-c(1,k+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,k+1)])
      derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
      OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN)
      svdOmegaN = svd(OmegaN)
      splitMidPointY = (splitPointYTest[2:(numberBin+1)] + splitPointYTest[1:numberBin])/2
      z = splines::ns(x=splitMidPointY, knots = knotsNs[-c(1,k+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,k+1)])
      z = z %*% svdOmegaN$u
      splitPointY = c(splitPointYTest, Inf)
      h = splitMidPointY[2] - splitMidPointY[1]
    }
  }

  # prior
  if(class(trees$prior) == "character"){
    if(trees$prior[1] == "uniform"){
      cellProb0 = matrix(rep(1/numberBin, numberBin), nrow = 1)
    } else if (trees$prior[1] == "Gaussian"){
      cellProb0 = matrix(dnorm(splitMidPointY, mean = as.numeric(trees$prior[2]), sd = as.numeric(trees$prior[3])), nrow = 1)
      cellProb0 = cellProb0/sum(cellProb0)
    } else if (trees$prior[1] == "LindseyMarginal"){
      betaPrior = as.numeric(trees$prior[-1])
      cellProb0 = matrix(exp(z %*% betaPrior), nrow=1)
      cellProb0 = cellProb0/sum(cellProb0)
    }
    cellProb0 = cellProb0[rep(1, n),]
  } else if(class(trees$prior) == "function"){
    cellProb0 = t(matrix(trees$prior(X = X[rep(seq(1,n), rep(numberBin, n)),], y = rep(splitMidPointY, n)), nrow = numberBin))
    cellProb0 = diag(1/apply(cellProb0, 1, sum)) %*% cellProb0
  }

  if(is.null(trees$shrinkage)){shrinkage = 1} else {shrinkage = trees$shrinkage}
  trees = trees$trees; numTree = length(trees)

  # initialization
  beta = matrix(0, nrow = n, ncol = k)

  # response discretization
  if(full){
    ySort = sort(y, decreasing = FALSE); yIndex = rep(0, n); counter = 1
    numberStop = min(n+1, which(ySort >= splitPointY[numberBin+1]))-1
    i=1
    while(i <= numberStop){
      if(ySort[i] >= splitPointY[counter+1]){
        counter = counter + 1
        if(ySort[i] >= splitPointY[counter+1]){next}
      }
      yIndex[i] = counter; i = i+1
    }
    if(numberStop < n){yIndex[(1+numberStop) : n] = numberBin}
    yIndex[order(y, decreasing = FALSE)] = yIndex
  }

  if(numTree == 0){
    result = list()
    if(full){
      # result$density = cellProb0[yIndex]/h; result$testLogLikelihood = mean(log(result$density))
      result$density = cellProb0[cbind(seq(1,n), yIndex)]/h; result$testLogLikelihood = mean(log(result$density))}
      else{result$density = NA; result$testLogLikelihood = NA}
    result$cellProb = cellProb0; result$yDiscretized = (splitPointY[-c(1, numberBin+2)] + splitPointY[-c(numberBin+1, numberBin+2)])/2
    result$augmentation = augmentation
    if(augmentation){result$yMean = yMean}
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
  if(n < 100){data = rbind(data, data[rep(1,100-n),] + matrix(rnorm((100-n)*(d+1)),ncol = (d+1)))}
  gbmObject = suppressWarnings(gbm::gbm(y~., data = data[sample(dim(data)[1], 100),], distribution = "gaussian", cv.folds = 0, n.trees = 1, interaction.depth = depth)) #2^interaction.depth
  gbmObject$trees = gbmTrees
  gbmObject$n.trees = numTree
  # compute the matrix of memberships (n by numTree)
  nodeMembership = suppressWarnings(matrix(gbm::predict.gbm(gbmObject, newdata = data[1:n, ], n.trees = seq(1:numTree), single.tree = TRUE), ncol = numTree))

  # prediction
  testLogLikelihood = numeric(0)
  for(i in 1:numTree){
    currTree = trees[[i]]
    indexLeaf = which(currTree$SplitVar == 0)
    betaTemp = currTree[indexLeaf,7:(k+6)]
    beta = beta + betaTemp[nodeMembership[,i],] * shrinkage

    cellProb = exp(as.matrix(beta) %*% t(z)) * cellProb0 # homogeneous prior -
    if(full){
    temp = mean(log(cellProb[cbind(seq(1,n), yIndex)]/h)) - mean(log(apply(cellProb, 1, sum)))
    testLogLikelihood = c(testLogLikelihood, temp)
    }
  }
  # beta = beta * shrinkage
  cellProb = exp(as.matrix(beta) %*% t(z)) * cellProb0 # homogeneous prior
  # compute normalizing constants
  a0 = -log(apply(cellProb, 1, sum)) - log(h)
  cellProb = cellProb * exp(a0) * h
  # compute conditional density
  result = list()
  if(full){
    result$density = cellProb[cbind(seq(1,n), yIndex)]/h; result$testLogLikelihood = mean(log(result$density))}
  else{result$density = NA; result$testLogLikelihood = NA}
  result$cellProb = cellProb; result$yDiscretized = (splitPointY[-c(1, numberBin+2)] + splitPointY[-c(numberBin+1, numberBin+2)])/2; result$testLogLikelihoodHistory = testLogLikelihood
  result$augmentation = augmentation
  if(augmentation){result$yMean = yMean}
  return(result)
}

#' LinCDETreePredict
#'
#' This function conducts prediction of a LinCDE tree.
#'
#' @param X input matrix for prediction, of dimension nobs x nvars; each row represents an observation vector.
#' @param tree a LinCDE tree.
#' @param currentRow the row index of the current node.
#'
#' @return The function returns the \code{membership} vector (length nobs): each element represents the row (terminal node) that a sample falls into.
#'
#' @export
LinCDETreePredict = function(X, tree, currentRow = 1){
  if(is.null(dim(X))){X = matrix(X, nrow = 1)}
  if(tree$SplitVar[currentRow] == 0){
    membership = rep(currentRow, dim(X)[1]); return(membership)
  }
  indexLeft = which(X[,tree$SplitVar[currentRow]] <= tree$SplitCodePred[currentRow])
  if(length(indexLeft) > 0){
    membershipLeft = LinCDETreePredict(X = X[indexLeft, ], tree = tree, currentRow = tree$LeftNode[currentRow])
  } else {
    membership = LinCDETreePredict(X = X, tree = tree, currentRow = tree$RightNode[currentRow])
    return(membership)
  }
  if(length(indexLeft) < dim(X)[1]){
    membershipRight = LinCDETreePredict(X = X[-indexLeft, ], tree = tree, currentRow = tree$RightNode[currentRow])
    membership = rep(0, dim(X)[1])
    membership[indexLeft] = membershipLeft
    membership[-indexLeft] = membershipRight
    return(membership)
  } else {
    return(membershipLeft)
  }
}




