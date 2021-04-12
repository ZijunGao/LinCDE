# LinCDE prediction
# prediction of LinCDE trees or LinCDE boosting

# function to predict the conditional density of a sample given LinCDE trees
# input:
  # X: matrix of covariates (n by d)
  # y: vector of responses (length n)
  # trees: 
      # trees: list of LinCDE trees
      # z: matrix of sufficient statistics evaluated at center-points of each bin (numberBin by k)
      # splitMidPointY: center points in each bin
      # depth: maximum depth of the tree
# output:
  # result: a vector of conditional densities (length n)

LinCDEPredict = function(X, y = NULL, trees, boosting = FALSE, splitPointYTest = NULL){
  # pre-process
  if(is.null(dim(X))){X = matrix(X, nrow = 1)}
  n = dim(X)[1]; d = dim(X)[2];   
  if(!is.null(y)){full = TRUE} else {full = FALSE; y = rnorm(n)}
  z = trees$z; k = dim(z)[2]; depth = trees$depth; type = trees$type
  splitMidPointY = trees$splitMidPointY; h = splitMidPointY[2] - splitMidPointY[1]; numberBin = length(splitMidPointY); splitPointY = c(splitMidPointY - h/2, splitMidPointY[numberBin] + h/2); splitPointY = c(splitPointY, Inf); 
  if(!is.null(splitPointYTest[1])){
    if(type == "Gaussian"){
      numberBin = length(splitPointYTest)-1
      z = cbind(splitMidPointY, splitMidPointY^2)
      penalty = c(1,1); k = 2
    }else if(type == "bs"){
      numberBin = length(splitPointYTest)-1
      knotsBs = quantile(splitMidPointY, probs = (1:(k-3))/(k-2)) # previous knots -> previous splines
      splitMidPointY = (splitPointYTest[2:(numberBin+1)] + splitPointYTest[1:numberBin])/2
      z = bs(x=splitMidPointY, df = k, knots = knotsBs, degree = 3, intercept = FALSE)
      splitPointY = c(splitPointYTest, Inf)
      h = splitMidPointY[2] - splitMidPointY[1]
    } else if (type == "bsTransform"){
      numberBin = length(splitPointYTest)-1
      knotsBs = quantile(splitMidPointY, probs = (0:(k-2))/(k-2)) # previous knots -> previous splines
      splitMidPointY = (splitPointYTest[2:(numberBin+1)] + splitPointYTest[1:numberBin])/2
      z = bs(x=splitMidPointY, df = k, knots = knotsBs[-c(1,k-1)], degree = 3, intercept = FALSE, Boundary.knots = knotsBs[c(1,k-1)])
      # transformation
      xSeq = min(splitMidPointY) + (seq(1,1000)-0.5)/1000 * (max(splitMidPointY) - min(splitMidPointY)); hTemp = xSeq[2] - xSeq[1]
      zTemp = bs(x=xSeq, knots = knotsBs[-c(1,k-1)], intercept = FALSE, Boundary.knots = knotsBs[c(1,k-1)])
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
      zTemp = ns(x=xSeq, knots = knotsNs[-c(1,k+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,k+1)])
      derivative = apply(zTemp, 2, diff, diff = 3)/hTemp^3
      OmegaN = t(derivative) %*% derivative * hTemp; OmegaN = OmegaN/max(OmegaN) 
      svdOmegaN = svd(OmegaN)
      splitMidPointY = (splitPointYTest[2:(numberBin+1)] + splitPointYTest[1:numberBin])/2
      z = ns(x=splitMidPointY, knots = knotsNs[-c(1,k+1)], intercept = FALSE, Boundary.knots = knotsNs[c(1,k+1)])
      z = z %*% svdOmegaN$u
      splitPointY = c(splitPointYTest, Inf) 
      h = splitMidPointY[2] - splitMidPointY[1]
    }
  }
  
  # prior
  if(trees$prior[1] == "uniform"){
    cellProb0 = matrix(rep(1/numberBin, numberBin), nrow = 1)
  } else if (trees$prior[1] == "Gaussian"){
    cellProb0 = matrix(dnorm(splitMidPointY, mean = as.numeric(trees$prior[2]), sd = as.numeric(trees$prior[3])), nrow = 1)
    cellProb0 = cellProb0/sum(cellProb0)
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
  
  # transfer LinCDE tree to gbm tree
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
  # creat the basics of a gbm.object
  data = data.frame(y, X)
  if(n < 100){data = rbind(data, data[rep(1,100-n),] + matrix(rnorm((100-n)*(d+1)),ncol = (d+1)))}
  gbmObject = gbm(y~., data = data[sample(dim(data)[1], 100),], distribution = "gaussian", cv.folds = 1, n.trees = 1, interaction.depth = (depth+1-boosting)) #2^interaction.depth
  gbmObject$trees = gbmTrees
  gbmObject$n.trees = numTree
  # compute the matrix of memberships (n by numTree)
  nodeMembership = matrix(predict.gbm(gbmObject, newdata = data[1:n, ], n.trees = seq(1:numTree), single.tree = TRUE), ncol = numTree) 
  
  # prediction
  for(i in 1:numTree){
    currTree = trees[[i]]
    indexLeaf = which(currTree$SplitVar == 0)
    betaTemp = currTree[indexLeaf,7:(k+6)]
    beta = beta + betaTemp[nodeMembership[,i],]
  }
  beta = beta * shrinkage
  cellProb = exp(as.matrix(beta) %*% t(z)) * cellProb0[rep(1,n),] # homogeneous prior
  # compute normalizing constants 
  a0 = -log(apply(cellProb, 1, sum)) - log(h)
  cellProb = cellProb * exp(a0) * h
  # compute conditional density
  result = list()
  if(full){
    result$density = cellProb[cbind(seq(1,n), yIndex)]/h; result$testLogLikelihood = mean(log(result$density))} 
  else{result$density = NA; result$testLogLikelihood = NA}
  result$cellProb = cellProb; result$splitPointY = splitPointY[-(numberBin+2)]
  return(result)
}

# previous name: LinCDETree.predict 
LinCDETree.predict = LinCDEPredict

