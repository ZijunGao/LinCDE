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



