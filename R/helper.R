# helper functions for LinCDE

# function: count the samples in each bin
# input:
  # yIndex: vector of response indices
  # numberBin: number of bins
# output:
  # counts: number of samples in each bin

countIndex = function(yIndex, numberBin) {
  counts = rep(0, numberBin)
  for(i in 1:numberBin) {counts[i] = sum(yIndex == i)}
  return (counts)
}


# function: compute the regularization parameter corresponding to the
# input:
  # z: sufficient statistics (numberBin by k matrix)
  # counts: number of samples in each bin
  # order: number of sufficient statistics
  # df: degree of freedom
  # numberBin: number of bins
  # nLambda: the number of lambda values
  # penalty: vector of penalties applied to each coefficient
# output:
  # lambda: tuning parameter

# test example:
  # dfToLambda2(z = z, counts = counts, order = 10, df = 1, numberBin = 40, penalty = c(2,rep(1,order-1)))

dfToLambda = function(z, counts, order, df, numberBin, penalty = NULL){
  if(order <= df){return (0)}
  if(df == 0){return (1e9)}

  # estimate the weight
  model = suppressWarnings(glm(counts ~ z, family = "poisson"))
  weight = model$fitted.values/numberBin

  Hessian = 2 * (t(z) %*% diag(weight) %*% z - t(z) %*% weight %*% weight %*% z * numberBin / sum(counts))
  svdHessian = svd(Hessian)
  Omega = diag(1/sqrt(svdHessian$d+1e-3)) %*% t(svdHessian$u) %*% diag(penalty) %*% svdHessian$u %*% diag(1/sqrt(svdHessian$d+1e-3))
  eigenVal = svd(Omega)$d

  lambdaMax = sum(1/eigenVal)/df; lambdaMin = 0; lambda = (lambdaMax+lambdaMin)/2; value =  sum(1/(1+eigenVal*lambda))
  while(abs(value - df) > 1e-3){
    if(value > df){lambdaMin = lambda; lambda = (lambdaMax+lambdaMin)/2; value = sum(1/(1+eigenVal*lambda))}
    else {lambdaMax = lambda; lambda = (lambdaMax+lambdaMin)/2; value = sum(1/(1+eigenVal*lambda))}
  }

  return (lambda/2)
}


# function: compute the importance score from a LinCDE tree
# input:
  # tree: a LinCDE tree
  # d: number of covariates
# output:
  # importance: vector of scores measuring the contribution of each covariate to the objective

importanceScore = function(tree, d, n.trees = 1){
  importance = rep(0, d)
  for(l in 1:n.trees){
    currTree = tree[[l]]
    index = sort(unique(currTree$SplitVar))
    importance[index[-1]] = importance[index[-1]] + aggregate(x = currTree$ErrorReduction, by = list(currTree$SplitVar), FUN = sum)$x[-1]
  }
  return (importance)
}

# 04/23/2021
# function: compute the list of candidate splits
# input:
  # X: variable matrix (nobs by nvar)
  # numberSplit: vector or scalar of the resulting split intervals' numbers (length nvar or 1)
# output:
# splitPoint: a candidate split list (length nvars). Each element is a vector corresponding to a certain variable's candidate splits (length of numberSplit or 1 + the number of unique values).
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


