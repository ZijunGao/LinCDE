#' LinCDE boosting
#'
#' This function implements the boosting algorithm of conditional density estimation with shallow LinCDE trees as base-learners, i.e., LinCDE boosting.
#'
#' @param X input matrix, of dimension nobs x nvars; each row represents an observation vector.
#' @param y response vector, of length nobs.
#' @param splitPoint a list of candidate splits, of length nvars. Each element is a vector corresponding to a variable's candidate splits (including the left and right end points). The list's elements are ordered the same as \code{X}'s columns. An alternative input is candidate split numbers, a scalar if all variables share the same number of candidate splits, a vector of length nvars if variables have different numbers of candidate splits. If candidate split numbers are given, each variable's range is divided into \code{splitPoint-1} intervals containing approximately the same number of observations. Default is 20. Note that if a variable has fewer unique values than the desired number of intervals, split intervals corresponding to each unique value are created.
#' @param numberBin the number of bins for response discretization. Default is 40. The response range is divided into \code{numberBin} equal-width bins.
#' @param z sufficient statistics, i.e., spline basis. For \code{z = "Gaussian"}, y, \eqn{y^2} are used. For \code{z = "bsTransform"}, transformed cubic B-splines are used. For \code{z = "nsTransform"}, transformed natural cubic splines are used. Default is "nsTransform".
#' @param prior the type of the initial carrying density. For \code{prior = "uniform"}, the uniform distribution over the response range is used. For \code{prior = "Gaussian"}, the Gaussian distribution with the marginal response mean and standard deviation is used. For \code{prior = "LindseyMarginal"}, the marginal response density estimated by Lindsey's method based on all responses is used. The argument \code{prior} can also be a homogeneous or heterogeneous conditional density function. The conditional density function should take a covariate matrix X, a response vector y, and output the densities at pairs {(Xi, yi)}. See the LinCDE vignette for examples. Default is "Gaussian".
#' @param splineDf the number of sufficient statistics. If \code{z = "Gaussian"}, \code{splineDf} is set to 2. Default is 10.
#' @param df the ridge Poisson regression's degrees of freedom. \code{df} is used for determining the ridge regularization hyper-parameter. If \code{z = "Gaussian"}, no penalization is implemented. If \code{df = splineDf}, there is no ridge penalization. Default is 2.
#' @param penalty vector of penalties applied to each sufficient statistics' coefficient. Default is 1 for each coefficient.
#' @param terminalSize the minimum number of observations in a terminal node. Default is 20.
#' @param depth the number of splits of each LinCDE tree. The number of terminal nodes is \code{depth + 1}. For \code{depth = 1}, an additive model is fitted. Default is 2.
#' @param n.trees the number of LinCDE trees to fit. Default is 100.
#' @param shrinkage the shrinkage parameter applied to each tree in the expansion, value in (0,1]. Default is 0.1.
#' @param augmentation  If TRUE, a conditional mean model is fitted first, and LinCDE boosting is applied to the residuals. The augmentation is recommended for responses whose conditional support varies wildly. See the LinCDE vignette for examples. Default is FALSE.
#' @param augmentationMethod  conditional mean estimator. If \code{augmentationMethod = "linearRegression"}, a regression model is fitted to the response. If \code{augmentationMethod = "randomForest"}, a random forest model is fitted. Default is randomForest. Applies only to \code{augmentation = TRUE}.
#' @param minY the user-provided left end of the response range. Default is \code{min(y)}.
#' @param maxY the user-provided right end of the response range. Default is \code{max(y)}.
#'
#' @return LinCDE boosting returns \code{trees}: a list of LinCDE trees. An \code{importanceScore} score vector measuring the contribution of each covariate to the objective is also available. For predictions based on the fitted LinCDE boosting model, please refer to the function \code{LinCDEPredict}.
#'
#' @export
LinCDEBoosting = function(y, X, splitPoint = 20, z = "nsTransform", splineDf = 10, prior = "Gaussian", numberBin = 40, df = 2, penalty = NULL, terminalSize = 20, depth = 2, n.trees = 100, shrinkage = 0.1, augmentation = FALSE, augmentationMethod = "randomForest", minY = NULL, maxY = NULL){
  alpha = 0.2; order = splineDf; type = z; z = NULL
  # pre-process
  n = length(y);
  if(class(X)[1] == "data.frame"){X = as.matrix(X)}
  if(is.null(dim(X))){X = matrix(X, ncol = 1)} # one covariate
  d = dim(X)[2]
  # candidate splits
  if(class(splitPoint) != "list"){
    numberSplit = splitPoint; splitPoint = constructSplitPoint(X, numberSplit)
  }

  # response discretization
  if(augmentation){
    if(augmentationMethod == "randomForest"){
      meanModel = randomForest::randomForest(y = y, x = X)
    } else if (augmentationMethod == "linearRegression"){
      meanModel = lm(y ~ X)
    }
    yPredict = predict(meanModel)
    y = y - yPredict
  } else {
    if(!is.null(minY)){y = c(y, minY)} else {y = c(y, min(y)-0.1)}
    if(!is.null(maxY)){y = c(y, maxY)} else {y = c(y, max(y)+0.1)}
  }
  splitPointY = seq(quantile(y, probs = 0), quantile(y, probs = 1), length.out = (numberBin + 1)); h = splitPointY[2] - splitPointY[1]
  y = y[1:n]
  splitMidPointY = (splitPointY[1 : numberBin] + splitPointY[2 : (numberBin + 1)])/2

  # compute yIndex
  ySort = sort(y, decreasing = FALSE); yIndex = rep(0, n); counter = 1
  numberStop = min(n+1,which(ySort >= splitPointY[numberBin+1]))-1
  i=1
  while(i <= numberStop){
    if(ySort[i] >= splitPointY[counter+1]){
      counter = counter+1
      if(ySort[i] >= splitPointY[counter+1]){next}
    }
    yIndex[i] = counter; i = i+1
  }
  if(numberStop < n){
    yIndex[(1+numberStop) : n] = numberBin
  }
  yIndex[order(y, decreasing = FALSE)] = yIndex

  # contruct sufficient statistics
  if(is.null(z)){
    if(type == "Gaussian"){
      z = cbind(splitMidPointY, splitMidPointY^2)
      penalty = c(1,1); order = 2
    } else if(type == "bsTransform") {
      knotsBs = quantile(splitMidPointY, probs = (0:(order-2))/(order-2))
      z = splines::bs(x=splitMidPointY, df = order, knots = knotsBs[-c(1,order-1)], degree = 3, intercept = FALSE, Boundary.knots = knotsBs[c(1,order-1)])
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = splines::bs(x=xSeq, knots = knotsBs[-c(1,order-1)], intercept = FALSE, Boundary.knots = knotsBs[c(1,order-1)])
      derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
      OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN)
      svdOmegaN = svd(OmegaN)
      z = z %*% svdOmegaN$u
      penalty = svdOmegaN$d/sum(svdOmegaN$d) * order + 0.1 # +0.1 for robustness
      penalty = penalty/sum(penalty) * order # sum(penalty) = order, as in glmnet
    } else if (type == "nsTransform"){
      knotsNs = quantile(splitMidPointY, probs = (1:(order+1))/(order+2))
      z = splines::ns(x=splitMidPointY, knots = knotsNs[-c(1,order+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,order+1)])
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = splines::ns(x=xSeq, knots = knotsNs[-c(1,order+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,order+1)])
      derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
      OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN)
      svdOmegaN = svd(OmegaN)
      z = z %*% svdOmegaN$u
      penalty = svdOmegaN$d/sum(svdOmegaN$d) * order + 0.1 # +0.1 for robustness
      penalty = penalty/sum(penalty) * order # sum(penalty) = order, as in glmnet
    }
  }

  # compute lambda
  if(is.null(df)){df = order}
  if(is.null(penalty)){penalty = rep(1,order)}
  counts = countIndex(yIndex = yIndex, numberBin = numberBin)
  if(type == "Gaussian"){lambda = 0
  } else if(type == "bsTransform" || type == "nsTransform"){lambda = dfToLambda(z = z, counts = counts, order = order, df = df, numberBin = numberBin, penalty = penalty)/n}

  # boosting
  betaName = paste0(rep("beta", order), 1:order)
  trees = list()
  if(class(prior) == "character"){
    if(prior == "uniform"){
      cellProb = matrix(1/numberBin, nrow = n, ncol = numberBin)
    } else if (prior == "Gaussian"){
      cellProb = matrix(dnorm(splitMidPointY, mean = mean(y), sd = sd(y)), nrow=1)
      cellProb = cellProb/sum(cellProb)
      cellProb = cellProb[rep(1, n),]
    } else if (prior == "LindseyMarginal"){
      counts = countIndex(yIndex = yIndex, numberBin = numberBin)
      modelPrior = suppressWarnings(glmnet::glmnet(x = z, y = counts, family = "poisson", lambda = 0, alpha = 0))
      cellProb = matrix(exp(z %*% as.vector(modelPrior$beta)), nrow=1)
      cellProb = cellProb/sum(cellProb)
      cellProb = cellProb[rep(1, n),]
    }
  } else if (class(prior) == "function"){
    cellProb = t(matrix(prior(X = X[rep(seq(1,n),rep(numberBin, n)),], y = rep(splitMidPointY, n)), nrow = numberBin))
    cellProb = diag(1/apply(cellProb, 1, sum)) %*% cellProb
  }

  # tree initialization
  if(n.trees == 0){
    result = list(); result$trees = NULL; result$splitMidPointY = splitMidPointY; result$z = z; result$depth = depth
    if(class(prior) == "character"){
      if(prior == "uniform"){
        result$prior = prior
      } else if(prior == "Gaussian"){
        result$prior = c(prior, mean(y), sd(y))
      } else if(prior == "LindseyMarginal"){
        result$prior = c(prior, as.vector(modelPrior$beta))
      }
    } else if(class(prior) == "function"){
      result$prior = prior
    }
    result$importanceScore = NULL; result$shrinkage = shrinkage; result$type = type; result$augmentation = augmentation
    if(augmentation){result$augmentationMethod = augmentationMethod; result$augmentationModel = meanModel}
    return (result)
  }
  rootNode = data.frame(matrix(rep(0, 6 + order), nrow = 1))
  indexList = list(); indexList[[1]] = seq(1, n)
  for(l in 1:n.trees){
    # recursive partitioning
    tree = LinCDEBoostingHelper(X = X, yIndex = yIndex, z = z, currDepth = 0, tree = rootNode, indexList = indexList, currNodeIndex = 1, cellProb = cellProb, improvementQueue =  collections::PriorityQueue(), splitPoint = splitPoint, numberBin = numberBin, order = order, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage)
    if(is.null(tree)){
      if(l==1) {print("no heterogeneity of conditional density is detected; the prior is returned")
        result = list(); result$trees = NULL; result$splitMidPointY = splitMidPointY; result$z = z; result$depth = depth
        if(prior == "uniform"){
          result$prior = "uniform"
        } else if(prior == "Gaussian"){
          result$prior = c("Gaussian", mean(y), sd(y))
        }
        result$importanceScore = rep(0, d); result$shrinkage = shrinkage; result$type = type; result$augmentation = augmentation
        if(augmentation){result$augmentationMethod = augmentationMethod; result$augmentationModel = meanModel}
        return(result)
      }
      print("iteration stopped after growing trees:"); print(l-1); n.trees = l-1; break
    }
    trees[[l]] = tree$tree
    colnames(trees[[l]]) = c("SplitVar", "SplitCodePred", "LeftNode", "RightNode", "ErrorReduction", "Weight", betaName)
    cellProb = tree$cellProb
  }

  # post-process
  result = list(); result$trees = trees; result$splitMidPointY = splitMidPointY; result$z = z; result$depth = depth
  if(class(prior) == "character"){
    if(prior == "uniform"){
      result$prior = prior
    } else if(prior == "Gaussian"){
      result$prior = c(prior, mean(y), sd(y))
    } else if(prior == "LindseyMarginal"){
      result$prior = c(prior, as.vector(modelPrior$beta))
    }
  } else if(class(prior) == "function"){
    result$prior = prior
  }
  result$importanceScore = importanceScore(tree = trees, d = d, n.trees = n.trees); result$shrinkage = shrinkage; result$type = type; result$augmentation = augmentation
  if(augmentation){result$augmentationMethod = augmentationMethod; result$augmentationModel = meanModel}
  return(result)
}


# This function is a helper for the function LinCDEBoosting. This helper grows one LinCDE tree with heterogeneous carrying density.
LinCDEBoostingHelper = function(X, yIndex, z = NULL, currDepth, tree, indexList, currNodeIndex, cellProb, improvementQueue, splitPoint, numberBin = 20, order = 4, lambda = 0, penalty = NULL, df = 2, alpha = 0.2, terminalSize = 20, depth = 2, shrinkage = 0.1){
  # initial split
  if(currDepth == 0){
    # preprocess
    n = length(yIndex)
    # initialization
    currNode = rep(0, order + 6)
    # estimate psi''
    counts = countIndex(yIndex = yIndex, numberBin = numberBin)
    offset = log(apply(cellProb, 2, mean))
    modelBeforeSplit = suppressWarnings(glmnet::glmnet(x = z, y = counts, offset = offset, family = "poisson", lambda = lambda * n, alpha = 0, penalty.factor = penalty))
    prob = as.vector(exp(z %*% as.vector(modelBeforeSplit$beta))) * exp(offset); prob = prob/sum(prob)
    temp = t(z) %*% prob
    covMatrix = (t(z) %*% (prob * z) - temp %*% t(temp))
    covMatrix = covMatrix + diag(lambda * numberBin * penalty)
    covMatrixInv = solve(covMatrix)
    # find the best split
    currSplit = LinCDESplit(X, yIndex-1, cellProb, z, covMatrixInv, splitPoint, numberBin)
    # insignificant
    if(currSplit[3] <= qchisq(p = (1-alpha), df = df)/2){return(NULL)}
    # update
    currNode[1] = currSplit[1]
    currNode[2] = currSplit[2]
    currNode[5] = currSplit[3]
    currNode[6] = n
    currNode[7:(order+6)] = as.vector(modelBeforeSplit$beta)
    tree[1,] = currNode
    result = LinCDEBoostingHelper(X = X, yIndex = yIndex, z = z, currDepth = 1, tree = tree, indexList = indexList, currNodeIndex = 1, cellProb = cellProb, improvementQueue = improvementQueue, splitPoint = splitPoint, numberBin = numberBin, order = order, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage)
    return(result)
  }

  # split
  # preprocess
  currIndex = indexList[[currNodeIndex]]; n = length(currIndex)
  currNode = tree[currNodeIndex,]; currNode = as.vector(t(currNode))
  currNode[3] = dim(tree)[1] + 1
  currNode[4] = currNode[3] + 1
  tree[currNodeIndex,] = currNode
  tree = rbind(tree, 0);  tree = rbind(tree, 0); # add two new nodes to the tree

  # split in the left child
  leftNode = rep(0, order + 6)
  indexLeft = currIndex[which(X[currIndex, currNode[1]] <= currNode[2])]
  nLeft = length(indexLeft)
  # estimate psi''
  counts = countIndex(yIndex = yIndex[indexLeft], numberBin = numberBin)
  offset = log(apply(cellProb[indexLeft, ], 2, mean))
  modelBeforeSplit = suppressWarnings(glmnet::glmnet(x = z, y = counts, offset = offset, family = "poisson", lambda = lambda * nLeft, alpha = 0, penalty.factor = penalty))
  prob = as.vector(exp(z %*% as.vector(modelBeforeSplit$beta))) * exp(offset); prob = prob/sum(prob)
  temp = t(z) %*% prob
  covMatrix = (t(z) %*% (prob * z) - temp %*% t(temp))
  covMatrix = covMatrix + diag(lambda * numberBin * penalty)
  covMatrixInv = solve(covMatrix)
  # find the best split
  currSplit = LinCDESplit(X[indexLeft, ], yIndex[indexLeft]-1, cellProb[indexLeft, ], z, covMatrixInv, splitPoint, numberBin)
  # insignificant
  if(currSplit[3] <= qchisq(p = (1-alpha), df = df)/2){
    indexList[[currNode[3]]] = indexLeft
    leftNode[1] = 0
    leftNode[6] = nLeft
    leftNode[7:(order+6)] = as.vector(modelBeforeSplit$beta)
    tree[currNode[3],] = leftNode
  } else {
    # update
    improvementQueue$push(currNode[3], priority = currSplit[3])
    indexList[[currNode[3]]] = indexLeft
    leftNode[1] = currSplit[1]
    leftNode[2] = currSplit[2]
    leftNode[5] = currSplit[3]
    leftNode[6] = nLeft
    leftNode[7:(order+6)] = as.vector(modelBeforeSplit$beta)
    tree[currNode[3],] = leftNode
  }

  # split in the right child
  rightNode = rep(0, order + 6)
  indexRight = currIndex[which(X[currIndex, currNode[1]] > currNode[2])]
  nRight = length(indexRight)
  # estimate psi''
  counts = countIndex(yIndex = yIndex[indexRight], numberBin = numberBin)
  offset = log(apply(cellProb[indexRight, ], 2, mean))
  modelBeforeSplit = suppressWarnings(glmnet::glmnet(x = z, y = counts, offset = offset, family = "poisson", lambda = lambda * nRight, alpha = 0, penalty.factor = penalty))
  prob = as.vector(exp(z %*% as.vector(modelBeforeSplit$beta))) * exp(offset); prob = prob/sum(prob)
  temp = t(z) %*% prob
  covMatrix = (t(z) %*% (prob * z) - temp %*% t(temp))
  covMatrix = covMatrix + diag(lambda * numberBin * penalty)
  covMatrixInv = solve(covMatrix)
  # find the best split
  currSplit = LinCDESplit(X[indexRight, ], yIndex[indexRight]-1, cellProb[indexRight, ], z, covMatrixInv, splitPoint, numberBin)
  # insignificant
  if(currSplit[3] <= qchisq(p = (1-alpha), df = df)/2){
    indexList[[currNode[4]]] = indexRight
    rightNode[1] = 0
    rightNode[6] = nRight
    rightNode[7:(order+6)] = as.vector(modelBeforeSplit$beta)
    tree[currNode[4],] = rightNode
  } else {
    # update
    improvementQueue$push(currNode[4], priority = currSplit[3])
    indexList[[currNode[4]]] = indexRight
    rightNode[1] = currSplit[1]
    rightNode[2] = currSplit[2]
    rightNode[5] = currSplit[3]
    rightNode[6] = nRight
    rightNode[7:(order+6)] = as.vector(modelBeforeSplit$beta)
    tree[currNode[4],] = rightNode
  }
  # choose the next split
  # maximum number of nodes reached
  if(currDepth == depth){
    # change the leaf node ErrorReduction to zero
    while(improvementQueue$n > 0){
      leafIndex = improvementQueue$pop()
      tree[leafIndex, c(1,2,5)] = 0
    }
    # update cell probabilities
    leafIndex = which(tree[,1] == 0)
    for(i in 1:length(leafIndex)){
      cellProbIndex = indexList[[leafIndex[i]]]
      probDelta = as.vector(z %*% as.vector(t(tree[leafIndex[i], 7:(6+order)]))) * shrinkage
      cellProb[cellProbIndex,] = t(t(cellProb[cellProbIndex,]) * exp(probDelta))
      cellProb[cellProbIndex,] = cellProb[cellProbIndex,]/apply(as.matrix(cellProb[cellProbIndex,], nrow = length(cellProbIndex)), 1, sum)
    }
    result = list(); result$tree = tree; result$cellProb = cellProb
    return(result)
  }
  # no more splits
  if(improvementQueue$n == 0){
    # update cell probabilities
    leafIndex = which(tree[,1] == 0)
    for(i in 1:length(leafIndex)){
      cellProbIndex = indexList[[leafIndex[i]]]
      probDelta = as.vector(z %*% as.vector(t(tree[leafIndex[i], 7:(6+order)]))) * shrinkage
      cellProb[cellProbIndex,] = t(t(cellProb[cellProbIndex,]) * exp(probDelta))
      cellProb[cellProbIndex,] = cellProb[cellProbIndex,]/apply(as.matrix(cellProb[cellProbIndex,], nrow = length(cellProbIndex)), 1, sum)
    }
    result = list(); result$tree = tree; result$cellProb = cellProb
    return(result)
  }
  # recursion
  nextSplit = improvementQueue$pop()
  result = LinCDEBoostingHelper(X = X, yIndex = yIndex, z = z, currDepth = currDepth + 1, tree = tree, indexList = indexList, currNodeIndex = nextSplit, cellProb = cellProb, improvementQueue = improvementQueue, splitPoint = splitPoint, numberBin = numberBin, order = order, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage)
  return (result)
}



