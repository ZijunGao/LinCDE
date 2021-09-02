#' LinCDE.boost
#'
#' This function implements LinCDE boosting: a boosting algorithm of conditional density estimation with shallow LinCDE trees as base-learners.
#'
#' @param X input matrix, of dimension nobs x nvars; each row represents an observation vector.
#' @param y response vector, of length nobs.
#' @param splitPoint a list of candidate splits of length nvars or a scalar/vector of candidate split numbers. If \code{splitPoint} is a list, each object is a vector corresponding to a variable's candidate splits (including the left and right end points). The list's objects should be ordered the same as \code{X}'s columns. An alternative input is candidate split numbers, a scalar if all variables share the same number of candidate splits, a vector of length nvars if variables have different numbers of candidate splits. If candidate split numbers are given, each variable's range is divided into \code{splitPoint-1} intervals containing approximately the same number of observations. Default is 20. Note that if a variable has fewer unique values than the desired number of intervals, split intervals corresponding to unique values are created.
#' @param numberBin the number of bins for response discretization. Default is 40. The response range is divided into \code{numberBin} equal-width bins.
#' @param basis a character or a matrix specifying sufficient statistics, i.e., spline basis. For \code{basis = "Gaussian"}, y, y^2 are used. For \code{basis = "nsTransform"}, transformed natural cubic splines are used. Default is "nsTransform".
#' @param prior a character or a function specifying initial carrier density. For \code{prior = "uniform"}, the uniform distribution over the response range is used. For \code{prior = "Gaussian"}, the Gaussian distribution with the marginal response mean and standard deviation is used. For \code{prior = "LindseyMarginal"}, the marginal response density estimated by Lindsey's method based on all responses is used. The argument \code{prior} can also be a homogeneous or heterogeneous conditional density function. The conditional density function should take a covariate matrix X, a response vector y, and output a vector of conditional densities f(yi | Xi). See the LinCDE vignette for examples. Default is "Gaussian".
#' @param splineDf the number of sufficient statistics/spline basis. If \code{z = "Gaussian"}, \code{splineDf} is set to 2. Default is 10.
#' @param df approximate degrees of freedom. \code{df} is used for determining the ridge regularization parameter. If \code{basis = "Gaussian"}, no penalization is implemented. If \code{df = splineDf}, there is no ridge penalization. Default is 6.
#' @param penalty vector of penalties applied to each sufficient statistics' coefficient.
#' @param terminalSize the minimum number of observations in a terminal node. Default is 20.
#' @param depth the number of splits of each LinCDE tree. The number of terminal nodes is \code{depth + 1}. If \code{depth = 1}, an additive model is fitted. Default is 1.
#' @param n.trees the number of trees to fit. Default is 100.
#' @param shrinkage the shrinkage parameter applied to each tree in the expansion, value in (0,1]. Default is 0.1.
#' @param subsample subsample ratio of the training samples in (0,1]. Default is 1.
#' @param centering  a logical value. If \code{TRUE}, a conditional mean model is fitted first, and LinCDE boosting is applied to the residuals. The centering is recommended for responses whose conditional support varies wildly. See the LinCDE vignette for examples. Default is \code{FALSE}.
#' @param centeringMethod a character or a function specifying the conditional mean estimator. If \code{centeringMethod = "linearRegression"}, a regression model is fitted to the response. If \code{centeringMethod = "randomForest"}, a random forest model is fitted. If \code{centeringMethod} is a function, the call \code{centeringMethod(X)} should return a conditional mean model with a predict function. Default is "randomForest". Applies only to \code{centering = TRUE}.
#' @param minY the user-provided left end of the response range. Default is NULL.
#' @param maxY the user-provided right end of the response range. Default is NULL.
#' @param alpha a hyperparameter in (0,1] to early stop the boosting. A smaller \code{alpha} is more likely to induce early stopping. If \code{alpha = 1}, no early stopping will be conducted.
#' @param verbose a logical value. If \code{TRUE}, progress and performance are printed. Default is \code{TRUE}.
#'
#' @return This function returns a LinCDE object consisting of a list of values.
#' \itemize{
#' \item trees: a list of LinCDE trees.
#' \item importanceScore: a named vector measuring the contribution of each covariate to the objective.
#' \item splitMidPointY: the vector of discretized bins' mid-points.
#' \item z: the spline basis matrix.
#' \item zTransformMatrix: the transformation matrix (of dimension splineDf x splineDf) multiplied by the standard natural cubic spline basis if \code{basis = "nsTransform"}.
#' \item prior: the prior function. The call \code{prior(X, Y)} should return a vector of prior conditional densities f(yi | Xi).
#' \item basis/depth/shrinkage/centering/centeringMethod: values inherited from the input arguments. If \code{centering} is \code{FALSE}, no \code{centeringMethod} is returned.
#' \item centeringModel: a centering model with a predict function. If \code{centering} is \code{FALSE}, no \code{centeringModel} is returned.
#' }
#'
#' @export
LinCDE.boost = function(y, X, splitPoint = 20, basis = "nsTransform", splineDf = 10, minY = NULL, maxY = NULL, numberBin = 40, df = 4, penalty = NULL, prior = "Gaussian", depth = 1, n.trees = 100, shrinkage = 0.1, terminalSize = 20, alpha = 0.2, subsample = 1, centering = FALSE, centeringMethod = "randomForest",verbose = TRUE, newKnots = FALSE){ # !!! to be deleted: newXXX
  # pre-process and initialization
  result = list(); result$depth = depth; result$shrinkage = shrinkage; result$basis = basis; result$centering = centering; result$centeringMethod = centeringMethod; result$y = y; result$X = X
  n = length(y)
  if(class(X)[1] == "data.frame"){X = as.matrix(X)}
  if(is.null(dim(X))){X = matrix(X, ncol = 1)} # a single covariate
  d = dim(X)[2]
  if(!is.null(colnames(X))){result$var.names = colnames(X)
  } else {result$var.names = paste("X", seq(1,d), sep = "")}
  if(is.null(df)){df = splineDf}
  if(is.null(penalty)){penalty = rep(1,splineDf)}
  if(class(splitPoint) != "list"){splitPoint = constructSplitPoint(X, splitPoint)}

  separate = function(x){min(abs(diff(unique(sort(x)))))} # make unique (X intervals: (]) # !!!
  resolution = min(apply(X, 2, separate)) # !!!
  X = X + runif(length(X),-resolution/10,0) # !!!

  # centering and response discretization
  if(centering){
    if(class(centeringMethod) == "character"){
      if(centeringMethod == "randomForest"){
        meanModel = randomForest::randomForest(y = y, x = X)
      } else if (centeringMethod == "linearRegression"){meanModel = lm(y ~ X)}
    }else if(class(centeringMethod) == "function"){meanModel=centeringMethod(y, X)}
    result$centeringModel = meanModel
    yPredict = predict(meanModel); y = y - yPredict
    y = c(y, 1.02*min(y)-0.02*max(y), 1.02*max(y)-0.02*min(y))
  } else {
    if(!is.null(minY)){y = c(y, minY)} else {y = c(y, 1.02*min(y)-0.02*max(y))}
    if(!is.null(maxY)){y = c(y, maxY)} else {y = c(y, 1.02*max(y)-0.02*min(y))}
  }
  splitPointY = seq(quantile(y, probs = 0), quantile(y, probs = 1), length.out = (numberBin + 1)); h = splitPointY[2] - splitPointY[1]
  splitMidPointY = (splitPointY[1:numberBin] + splitPointY[2:(numberBin+1)])/2
  result$splitMidPointY = splitMidPointY
  y = y[1:n]
  yIndex = as.numeric(cut(y, c(-Inf, splitPointY[-c(1,numberBin+1)], Inf)))

  # contruct sufficient statistics
  if(basis == "Gaussian"){
    z = cbind(splitMidPointY, splitMidPointY^2); result$z=z
    penalty = c(1,1); splineDf = 2
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
    # zTemp = splines::ns(x=xSeq, knots = knotsNs[-c(1,splineDf+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,splineDf+1)]) # !!!
    zTemp = splines::ns(x=xSeq, knots = attributes(zOrg)$knots, intercept = attributes(zOrg)$intercept, Boundary.knots = attributes(zOrg)$Boundary.knots) # !!!
    derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
    OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN)
    svdOmegaN = svd(OmegaN)
    penalty = svdOmegaN$d/sum(svdOmegaN$d) * splineDf + 0.001 # + 0.001 for robustness
    result$zTransformMatrix = svdOmegaN$u %*% diag(sqrt(penalty[splineDf]/penalty))
    z=zOrg%*%result$zTransformMatrix; attributes(z)=attributes(zOrg); result$z=z
    penalty = rep(1, splineDf)
  }

  # compute ridge regularization parameters
  counts = countIndex(yIndex = yIndex, numberBin = numberBin)
  if(basis == "Gaussian"){lambda = 0
  }else if(basis == "nsTransform"){
    lambda = dfToLambda(z = z, counts = counts, splineDf = splineDf, df = df, numberBin = numberBin, penalty = penalty)/n
    }

  # prior
  if(class(prior) == "character"){
    if(prior == "uniform"){
      cellProb = matrix(1/numberBin, nrow = n, ncol = numberBin)
      priorFunction = function(XTest, yTest){
        result=(yTest<=tail(splitPointY, n=1))*(yTest>=splitPointY[1])/numberBin/h
        return(result)
      }
    } else if (prior == "Gaussian"){
      cellProb = matrix(dnorm(splitMidPointY, mean = mean(y), sd = sd(y)), nrow=1)
      cellProb = cellProb/sum(cellProb); cellProb = cellProb[rep(1,n), ]
      priorFunction = function(XTest, yTest){
        result = dnorm(yTest, mean = mean(y), sd = sd(y)); return(result)
      }
    } else if (prior == "LindseyMarginal"){
      modelPrior = suppressWarnings(glmnet::glmnet(x = z, y = counts, family = "poisson", lambda = 0, alpha = 0, standardize = FALSE, intercept = TRUE)) # !!!
      # modelPrior = suppressWarnings(glmnet::glmnet(x = z, y = counts, family = poisson(), lambda = 0, alpha = 0, standardize = FALSE, intercept = TRUE)) # !!!
      cellProb = matrix(exp(z %*% as.vector(modelPrior$beta)), nrow=1)
      cellProb = cellProb/sum(cellProb); cellProb = cellProb[rep(1, n),]
      priorFunction = function(XTest, yTest){
        yTestIndex = as.numeric(cut(yTest, c(-Inf, splitPointY[-c(1,numberBin+1)], Inf), right = FALSE))
        result=cellProb[cbind(rep(1,length(yTest)), yTestIndex)]/h; return(result)
      }
    }
  } else if (class(prior) == "function"){
    cellProb = t(matrix(prior(X = X[rep(seq(1,n),rep(numberBin, n)),], y = rep(splitMidPointY, n)), nrow = numberBin))
    cellProb = diag(1/apply(cellProb, 1, sum)) %*% cellProb; priorFunction = prior
  }
  result$prior = priorFunction

  # boosting
  # depth = 0 gives the marginal density estimate of Lindsey's method with no prior
  if (depth == 0){
    modelBeforeSplit = glmnet::glmnet(x = z, y = counts, family = "poisson", nlambda = 100, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
    # modelBeforeSplit = glmnet::glmnet(x = z, y = counts, family = poisson(), lambda = lambda*n, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
    lambdaIndex = which.min(abs(modelBeforeSplit$lambda - n * lambda))
    result$trees = as.vector(modelPrior$beta[, lambdaIndex]); result$prior = NULL; result$importanceScore = NULL
    class(result) = "LinCDE"; return(result)
  }
  # tree initialization
  betaName = paste0(rep("beta", splineDf), 1:splineDf); trees = list()
  if(n.trees == 0){class(result) = "LinCDE"; return(result)}
  rootNode = data.frame(matrix(rep(0, 6 + splineDf), nrow = 1))
  if(subsample == 1){
    indexList = list(); indexList[[1]] = seq(1, n)
    for(l in 1:n.trees){
      # recursive partitioning
      tree = LinCDE.boost.helper(X = X, yIndex = yIndex, z = z, currDepth = 0, tree = rootNode, indexList = indexList, currNodeIndex = 1, cellProb = cellProb, improvementQueue =  collections::PriorityQueue(), splitPoint = splitPoint, numberBin = numberBin, splineDf = splineDf, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage)
      if(is.null(tree)){
        if(l==1) {print("no heterogeneity of conditional density is detected"); class(result) = "LinCDE"; return(result)}
        print("boosting stopped after trees:"); print(l-1); n.trees = l-1; break}
      trees[[l]] = tree$tree
      colnames(trees[[l]]) = c("SplitVar", "SplitCodePred", "LeftNode", "RightNode", "ErrorReduction", "Weight", betaName)
      cellProb = tree$cellProb
      if(verbose && (l %% 10 == 0)){
        currError = mean(log(cellProb[cbind(seq(1,n), yIndex)])) - log(h)
        print(paste(l, "trees ", "training error:", signif(currError, digits = 3), sep = " "))
      }
    }
  } else {
    indexList = list(); indexList[[1]] = seq(1, floor(n * subsample))
    if(floor(n*subsample)<=100){stop("insufficient number of subsamples")}
    for(l in 1:n.trees){
      indexSubsample = sample(n, floor(n*subsample), replace = FALSE)
      # recursive partitioning
      tree = LinCDE.boost.helper(X = X[indexSubsample, ], yIndex = yIndex[indexSubsample], z = z, currDepth = 0, tree = rootNode, indexList = indexList, currNodeIndex = 1, cellProb = cellProb[indexSubsample, ], improvementQueue =  collections::PriorityQueue(), splitPoint = splitPoint, numberBin = numberBin, splineDf = splineDf, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage)
      if(is.null(tree)){
        if(l==1) {print("no heterogeneity of conditional density is detected"); class(result) = "LinCDE"; return(result)
        }
        print("iteration stopped after trees:"); print(l-1); n.trees = l-1; break}
      trees[[l]] = tree$tree
      colnames(trees[[l]]) = c("SplitVar", "SplitCodePred", "LeftNode", "RightNode", "ErrorReduction", "Weight", betaName)
      cellProb[indexSubsample, ] = tree$cellProb
      membership = LinCDETreePredict(X = X[-indexSubsample,], tree = trees[[l]], currentRow = 1)
      probDelta = matrix(as.matrix(trees[[l]][,-seq(1,6)]) %*% t(z), ncol = numberBin) * shrinkage
      temp = matrix(cellProb[-indexSubsample,] * exp(probDelta)[membership,], ncol = numberBin)
      cellProb[-indexSubsample,] = temp/apply(temp, 1, sum)
    }
  }

  # post-process
  result$trees = trees; result$importanceScore = importanceScore(model = trees, d = d, n.trees = n.trees, var.names = result$var.names)
  class(result) = "LinCDE"; return(result)
}

#' LinCDE.boost.helper
#'
#' This function is a helper for the function LinCDE.boost. This helper grows one LinCDE tree with heterogeneous carrier density.
#'
#' @param X input matrix, of dimension nobs x nvars; each row represents an observation vector.
#' @param yIndex discretized response vector, of length nobs.
#' @param z sufficient statistics matrix, of dimension \code{numberBin} x \code{splineDf}. Default is NULL.
#' @param currDepth the current depth of the LinCDE tree.
#' @param tree the current LinCDE tree.
#' @param indexList the list of candidate splits.
#' @param currNodeIndex the index of the current node.
#' @param cellProb cell probability matrix, of dimension nobs x \code{numberBin}.
#' @param improvementQueue the queue recording improvements of candidate splits.
#' @param splitPoint a list of candidate splits, of length nvars.
#' @param numberBin the number of bins for response discretization. Default is 20.
#' @param splineDf the number of sufficient statistics. Default is 10.
#' @param lambda the regularization parameter for Poisson regression. Default is 0.
#' @param penalty separate penalty factors applied to each coefficient, of length \code{splineDf}. Default is NULL.
#' @param df degrees of freedom. Default is 6.
#' @param alpha parameter for the stopping rule of splitting. Default is 0.2.
#' @param terminalSize the minimum number of observations in a terminal node. Default is 20.
#' @param depth the number of splits of each LinCDE tree. The number of terminal nodes is \code{depth + 1}. For \code{depth = 1}, an additive model is fitted. Default is 1.
#' @param shrinkage the shrinkage parameter applied to each tree in the expansion, value in (0,1]. Default is 0.1.
#'
#' @return This function returns a LinCDE tree \code{trees} and the updated cell probability matrix \code{cellProb}.
#'
#' @export
LinCDE.boost.helper = function(X, yIndex, z = NULL, currDepth, tree, indexList, currNodeIndex, cellProb, improvementQueue, splitPoint, numberBin = 20, splineDf = 10, lambda = 0, penalty = NULL, df = 6, alpha = 0.2, terminalSize = 20, depth = 1, shrinkage = 0.1){
  result = list()
  # initial split
  if(currDepth == 0){
    n = length(yIndex); currNode = rep(0, splineDf + 6)
    # estimate psi''
    counts = countIndex(yIndex = yIndex, numberBin = numberBin)
    offset = log(apply(cellProb, 2, mean))
    modelBeforeSplit = glmnet::glmnet(x = z, y = counts, family = "poisson", offset = offset, nlambda = 100, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
    # modelBeforeSplit = glmnet::glmnet(x = z, y = counts, family = poisson(), offset = offset, lambda = lambda * n, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
    lambdaIndex = which.min(abs(modelBeforeSplit$lambda - n * lambda))
    prob = as.vector(exp(z %*% as.vector(modelBeforeSplit$beta[, lambdaIndex]))) * exp(offset); prob = prob/sum(prob)
    temp = t(z) %*% prob
    covMatrix = (t(z)%*%(prob*z) - temp%*%t(temp))+diag(lambda*numberBin*penalty)
    covMatrixInv = solve(covMatrix)
    # find the optimal split
    currSplit = LinCDESplit(X, yIndex-1, cellProb, z, covMatrixInv, splitPoint, numberBin)
    # early stopping
    if(currSplit[3] <= (qchisq(p = (1-alpha), df = df)/2)){return(NULL)}
    # update
    currNode[1]=currSplit[1]; currNode[2]=currSplit[2]; currNode[5]=currSplit[3]; currNode[6]=n; currNode[7:(splineDf+6)]=as.vector(modelBeforeSplit$beta[,lambdaIndex])
    tree[1,] = currNode
    result = LinCDE.boost.helper(X = X, yIndex = yIndex, z = z, currDepth = 1, tree = tree, indexList = indexList, currNodeIndex = 1, cellProb = cellProb, improvementQueue = improvementQueue, splitPoint = splitPoint, numberBin = numberBin, splineDf = splineDf, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage)
    return(result)
  }

  # split
  currIndex = indexList[[currNodeIndex]]; n = length(currIndex)
  currNode = tree[currNodeIndex,]; currNode = as.vector(t(currNode))
  currNode[3] = dim(tree)[1] + 1; currNode[4] = currNode[3] + 1
  tree[currNodeIndex,] = currNode; tree = rbind(tree, 0);  tree = rbind(tree, 0)

  # split in the left child
  leftNode = rep(0, splineDf + 6)
  indexLeft=currIndex[which(X[currIndex, currNode[1]]<=currNode[2])]
  nLeft = length(indexLeft)
  # estimate psi''
  counts = countIndex(yIndex = yIndex[indexLeft], numberBin = numberBin)
  offset = log(apply(cellProb[indexLeft, ], 2, mean))
  modelBeforeSplit = glmnet::glmnet(x = z, y = counts, offset = offset, family = "poisson", nlambda = 100, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
  # modelBeforeSplit = glmnet::glmnet(x = z, y = counts, offset = offset, family = poisson(), lambda = nLeft * lambda, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
  lambdaIndex = which.min(abs(modelBeforeSplit$lambda - nLeft * lambda)) #!!!
  prob = as.vector(exp(z %*% as.vector(modelBeforeSplit$beta[, lambdaIndex]))) * exp(offset); prob = prob/sum(prob)
  temp = t(z) %*% prob
  covMatrix = (t(z)%*%(prob*z) - temp%*%t(temp))+diag(lambda*numberBin*penalty)
  covMatrixInv = solve(covMatrix)
  # find the optimal split
  currSplit = LinCDESplit(as.matrix(X[indexLeft, ]), yIndex[indexLeft]-1, cellProb[indexLeft, ], z, covMatrixInv, splitPoint, numberBin)
  # insignificant
  if(currSplit[3] <= qchisq(p = (1-alpha), df = df)/2){
    indexList[[currNode[3]]] = indexLeft; leftNode[1] = 0; leftNode[6] = nLeft; leftNode[7:(splineDf+6)] = as.vector(modelBeforeSplit$beta[,lambdaIndex])
    tree[currNode[3],] = leftNode
  } else {
    # update
    improvementQueue$push(currNode[3], priority = currSplit[3])
    indexList[[currNode[3]]] = indexLeft; leftNode[1] = currSplit[1]; leftNode[2] = currSplit[2]; leftNode[5] = currSplit[3]; leftNode[6] = nLeft; leftNode[7:(splineDf+6)] = as.vector(modelBeforeSplit$beta[,lambdaIndex])
    tree[currNode[3],] = leftNode
  }

  # split in the right child
  rightNode = rep(0, splineDf + 6)
  indexRight = currIndex[which(X[currIndex, currNode[1]] > currNode[2])]
  nRight = length(indexRight)
  # estimate psi''
  counts = countIndex(yIndex = yIndex[indexRight], numberBin = numberBin)
  offset = log(apply(cellProb[indexRight, ], 2, mean))
  modelBeforeSplit = glmnet::glmnet(x = z, y = counts, offset = offset, family = "poisson", nlambda = 100, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
  # modelBeforeSplit = glmnet::glmnet(x = z, y = counts, offset = offset, family = poisson(), lambda = lambda * nRight, alpha = 0, penalty.factor = penalty, intercept = TRUE, standardize = FALSE) #!!!
  lambdaIndex = which.min(abs(modelBeforeSplit$lambda - nRight * lambda))
  prob = as.vector(exp(z %*% as.vector(modelBeforeSplit$beta[, lambdaIndex]))) * exp(offset); prob = prob/sum(prob)
  temp = t(z) %*% prob
  covMatrix = (t(z)%*%(prob * z) - temp%*%t(temp))+diag(lambda*numberBin*penalty)
  covMatrixInv = solve(covMatrix)
  # find the optimal split
  currSplit = LinCDESplit(as.matrix(X[indexRight, ]), yIndex[indexRight]-1, cellProb[indexRight, ], z, covMatrixInv, splitPoint, numberBin)
  # insignificant
  if(currSplit[3] <= qchisq(p = (1-alpha), df = df)/2){
    indexList[[currNode[4]]] = indexRight; rightNode[1] = 0; rightNode[6] = nRight; rightNode[7:(splineDf+6)] = as.vector(modelBeforeSplit$beta[,lambdaIndex])
    tree[currNode[4],] = rightNode
  } else {
    # update
    improvementQueue$push(currNode[4], priority = currSplit[3])
    indexList[[currNode[4]]] = indexRight; rightNode[1] = currSplit[1]; rightNode[2] = currSplit[2]; rightNode[5] = currSplit[3]; rightNode[6] = nRight; rightNode[7:(splineDf+6)] = as.vector(modelBeforeSplit$beta[,lambdaIndex])
    tree[currNode[4],] = rightNode
  }

  # maximum number of nodes reached
  if(currDepth == depth){
    # change the leaf node ErrorReduction to zero
    while(improvementQueue$n > 0){
      leafIndex = improvementQueue$pop(); tree[leafIndex, c(1,2,5)] = 0
    }
    # update cell probabilities
    leafIndex = which(tree[,1] == 0)
    for(i in 1:length(leafIndex)){
      cellProbIndex = indexList[[leafIndex[i]]]
      probDelta = as.vector(z %*% as.vector(t(tree[leafIndex[i], 7:(6+splineDf)]))) * shrinkage
      cellProb[cellProbIndex,] = t(t(cellProb[cellProbIndex,]) * exp(probDelta))
      cellProb[cellProbIndex,] = cellProb[cellProbIndex,]/apply(as.matrix(cellProb[cellProbIndex,], nrow = length(cellProbIndex)), 1, sum)
    }
    result$tree = tree; result$cellProb = cellProb; return(result)
  }
  # no more splits
  if(improvementQueue$n == 0){
    # update cell probabilities
    leafIndex = which(tree[,1] == 0)
    for(i in 1:length(leafIndex)){
      cellProbIndex = indexList[[leafIndex[i]]]
      probDelta = as.vector(z %*% as.vector(t(tree[leafIndex[i], 7:(6+splineDf)]))) * shrinkage
      cellProb[cellProbIndex,] = t(t(cellProb[cellProbIndex,]) * exp(probDelta))
      cellProb[cellProbIndex,] = cellProb[cellProbIndex,]/apply(as.matrix(cellProb[cellProbIndex,], nrow = length(cellProbIndex)), 1, sum)
    }
    result$tree = tree; result$cellProb = cellProb; return(result)
  }
  # recursion
  nextSplit = improvementQueue$pop()
  result = LinCDE.boost.helper(X = X, yIndex = yIndex, z = z, currDepth = currDepth + 1, tree = tree, indexList = indexList, currNodeIndex = nextSplit, cellProb = cellProb, improvementQueue = improvementQueue, splitPoint = splitPoint, numberBin = numberBin, splineDf = splineDf, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage)
  return (result)
}



