# LinCDE boosting
# boosting algorithm of conditional density estimation with shallow LinCDE tree as base-learners

# load helper.R 

# function: boosting algorithm of conditional density estimation with shallow LinCDE tree as base-learners
# input:
  # X: matrix of covariates (n by d)
  # y: vector of response (length n)
  # splitPoint: splits of X
  # numberBin: number of bins for the response discretization. Divide the range of y into numberBin equal width bins
  # type: type of sufficient statistics. "bs": cubic B-splines, "ns": natural cubic splines, "Gaussian": y, y^2
  # order: order of splines
  # df: degree of freedom of the sufficient statistics
  # penalty: vector of penalties applied to each coefficient. Default is 1 for all variables.
  # alpha: significance level for splitting-stopping test (primarily for LinCDE trees)
  # terminalSize: the minimum number of observations in a terminal node
  # depth: the maximum depth of the tree. No split corresponds to depth 0
  # n.trees: the total number of trees to fit (at least 2)
  # shrinkage: the shrinkage parameter applied to each tree in the expansion (value in (0,1])
  # minY, maxY: support of the responses; take min(min(Y), minY) and max(max(Y), maxY)
# output:
  # trees: a list of trees. 
  # tree: dataframe of nodes. Node:
    # SplitVar: index of which variable is used to split. 0 indicates a terminal node.
    # SplitCodePred: if the split variable is continuous then this component is the split point. If the split variable is categorical then this component contains the index of object$c.split that describes the categorical split. If the node is a terminal node then this is the prediction.
    # LeftNode: the index of the row corresponding to the left node.
    # RightNode: the index of the row corresponding to the right node.
    # ErrorReduction: the reduction in the loss function as a result of splitting this node.
    # Weight: the total weight of observations in the node. If weights are all equal to 1 then this is the number of observations in the node.
    # beta_i (1 <= i <= order): if the node is a leaf, estimated natural parameters in the leaf.
  # splitMidPointY: center points in each bin
  # z: sufficient statistics (matrix numberBin by k)
  # importance score: vector of scores measuring the contribution of each covariate to the objective
  # shrinkage: the shrinkage parameter applied to each tree in the expansion (value in (0,1])
  # time
  # prior: prior distribution, {"uniform", "Gaussian"}, default "Gaussian"

LinCDEBoosting = function(y, X, splitPoint, numberBin = 20, type = "Gaussian", order = 2, df = 2,  penalty = NULL, alpha = 1, terminalSize = 20, depth = 2, z = NULL, n.trees = 100, shrinkage = 0.1, minY = NULL, maxY = NULL, prior = "uniform"){
  # pre-process
  n = length(y); d = dim(X)[2]; splitPoint = cbind(splitPoint, max(X)); X = X + runif(length(X), -1e-6, 1e-6)
  
  # response discretization
  add = 0
  if(!is.null(minY)){y = c(y, minY)} 
  if(!is.null(maxY)){y = c(y, maxY)}
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
    } else if(type == "bs") {
      knotsBs = quantile(splitMidPointY, probs = (1:(order-3))/(order-2))
      z = bs(x=splitMidPointY, df = order, knots = knotsBs, degree = 3, intercept = FALSE)
    } else if(type == "bsTransform") {
      knotsBs = quantile(splitMidPointY, probs = (0:(order-2))/(order-2))
      z = bs(x=splitMidPointY, df = order, knots = knotsBs[-c(1,order-1)], degree = 3, intercept = FALSE, Boundary.knots = knotsBs[c(1,order-1)])
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = bs(x=xSeq, knots = knotsBs[-c(1,order-1)], intercept = FALSE, Boundary.knots = knotsBs[c(1,order-1)])
      derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
      OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN) 
      svdOmegaN = svd(OmegaN)
      z = z %*% svdOmegaN$u
      penalty = svdOmegaN$d/sum(svdOmegaN$d) * order + 0.1 # +0.1 for robustness
      penalty = penalty/sum(penalty) * order # sum(penalty) = order, as in glmnet
    } else if (type == "ns"){
      knotsNs = quantile(splitMidPointY, probs = (1:(order-1))/(order-0))
      z = ns(x=splitMidPointY, df = order, knots = knotsNs, intercept = FALSE)
    } else if (type == "nsTransform"){
      knotsNs = quantile(splitMidPointY, probs = (1:(order+1))/(order+2))
      z = ns(x=splitMidPointY, knots = knotsNs[-c(1,order+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,order+1)])
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = ns(x=xSeq, knots = knotsNs[-c(1,order+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,order+1)])
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
  if(type == "Gaussian"){
    lambda = 0
  } else if(type == "bs" || type == "ns" || type == "bsTransform" || type == "nsTransform"){lambda = dfToLambda(z = z, counts = counts, order = order, df = df, numberBin = numberBin, penalty = penalty)/n} 
  
  # boosting  
  betaName = paste0(rep("beta", order), 1:order)
  trees = list()
  if(prior == "uniform"){
    cellProb = matrix(1/numberBin, nrow = n, ncol = numberBin)
  } else if (prior == "Gaussian"){
    cellProb = matrix(dnorm(splitMidPointY, mean = mean(y), sd = sd(y)), nrow = 1)
    cellProb = cellProb/sum(cellProb)
    cellProb = cellProb[rep(1, n),]
  }
  # tree initialization
  rootNode = data.frame(matrix(rep(0, 6 + order), nrow = 1)) 
  indexList = list(); indexList[[1]] = seq(1, n)
  for(l in 1:n.trees){
    # recursive partitioning
    tree = LinCDEBoostingHelper(X = X, yIndex = yIndex, z = z, currDepth = 0, tree = rootNode, indexList = indexList, currNodeIndex = 1, cellProb = cellProb, improvementQueue =  PriorityQueue(), splitPoint = splitPoint, numberBin = numberBin, order = order, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage, time = time)
    if(is.null(tree)){
      print("iteration stopped after growing trees:"); print(l-1); n.trees = l-1; break
    }
    trees[[l]] = tree$tree
    colnames(trees[[l]]) = c("SplitVar", "SplitCodePred", "LeftNode", "RightNode", "ErrorReduction", "Weight", betaName)
    cellProb = tree$cellProb # update
  }

  # post-process
  result = list(); result$trees = trees; result$splitMidPointY = splitMidPointY; result$z = z; result$depth = depth
  if(prior == "uniform"){
    result$prior = "uniform"
  } else if(prior == "Gaussian"){
    result$prior = c("Gaussian", mean(y), sd(y)) 
  }
  result$importanceScore = importanceScore(tree = trees, d = d, n.trees = n.trees); result$shrinkage = shrinkage; result$type = type
  return(result)
}


# function: helper for function LinCDEBoosting. Grow one LinCDE tree with heterogeneous carrying density. 
# input:
  # X: matrix of covariates (n by d)
  # yIndex: vector of the bin indices that responses belong to (length n): value from 1 to numberBin
  # z: matrix of sufficient statistics evaluated at center-points of each bin (numberBin by k)
  # currDepth: current depth of the tree
  # tree: dataframe of nodes
  # indexList: list of indices of samples at each node
  # currNodeIndex: the index of current node
  # cellProb: matrix of cell probabilities (n by numberBin matrix)
  # improvementQueue: vector of improvements at each node
  # splitPoint: splits of X (matrix d by s, an artificial column with max(X) is added to the right)
  # numberBin: number of bins for the response discretization. Divide the range of y into numberBin equal width bins
  # order: number of sufficient statistics
  # lambda: regularization parameter 
  # penalty: vector of penalties applied to each coefficient. Default is 1 for all variables.
  # df: degree of freedom (used for splitting-stopping test)
  # alpha: significance level for splitting-stopping test
  # terminalSize: the minimum number of observations in a terminal node
  # depth: the maximum depth of the tree. No split corresponds to depth 0.
  # shrinkage: shrinkage parameter (value in (0,1]) 
# output:
  # tree: dataframe of nodes
  # cellProb: cell probabilities (n by numberBin matrix)
  
LinCDEBoostingHelper = function(X, yIndex, z = NULL, currDepth, tree, indexList, currNodeIndex, cellProb, improvementQueue, splitPoint, numberBin = 20, order = 2, lambda = 0, penalty, df = 2, alpha = 0.1, terminalSize = 20, depth = 5, shrinkage = 1, time){
  # initial split
  if(currDepth == 0){
    # preprocess
    n = length(yIndex)
    # initialization
    currNode = rep(0, order + 6)
    # estimate psi''
    counts = countIndex(yIndex = yIndex, numberBin = numberBin)
    offset = log(apply(cellProb, 2, mean))
    modelBeforeSplit = glmnet(x = z, y = counts, offset = offset, family = "poisson", lambda = lambda * n, alpha = 0, penalty.factor = penalty)
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
    result = LinCDEBoostingHelper(X = X, yIndex = yIndex, z = z, currDepth = 1, tree = tree, indexList = indexList, currNodeIndex = 1, cellProb = cellProb, improvementQueue = improvementQueue, splitPoint = splitPoint, numberBin = numberBin, order = order, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage, time = time)
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
  modelBeforeSplit = glmnet(x = z, y = counts, offset = offset, family = "poisson", lambda = lambda * nLeft, alpha = 0, penalty.factor = penalty)
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
  modelBeforeSplit = glmnet(x = z, y = counts, offset = offset, family = "poisson", lambda = lambda * nRight, alpha = 0, penalty.factor = penalty)
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
    result = list(); result$tree = tree; result$cellProb = cellProb; result$time = time
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
    result = list(); result$tree = tree; result$cellProb = cellProb; result$time = time
    return(result)
  }
  # recursion
  nextSplit = improvementQueue$pop()
  result = LinCDEBoostingHelper(X = X, yIndex = yIndex, z = z, currDepth = currDepth + 1, tree = tree, indexList = indexList, currNodeIndex = nextSplit, cellProb = cellProb, improvementQueue = improvementQueue, splitPoint = splitPoint, numberBin = numberBin, order = order, lambda = lambda, penalty = penalty, df = df, alpha = alpha, terminalSize = terminalSize, depth = depth, shrinkage = shrinkage, time = time)
  return (result)
}



