# example of LinCDE

# load helper.R, LinCDEBoosting.R, LinCDEPredict.R
# input:
  # setting: type of settings: "GM" or "GLM"
  # m: number of simulations
# output:
  # log-likelihood R^2
  # conditional density plots at landmarks

LinCDEExample = function(setting = "GLM", m = 2){
  library("latex2exp")
  library("lattice")
  set.seed(318)
  
  # preprocess
  if (m < 2){m = 2}
  
  # data generation parameters
  n = 1000; nTest = 1000; d = 20
  sd = 0.5; trueSplit1 = 0.2; trueSplit2 = 0.2
  beta = rep(0,d); beta0 = 0; varBeta = rep(0,d); varBase = 0.5; beta[1:5] = c(0.5,0,0,0,0); varBeta[1:5] = c(0,0.5,0,0,0)/sqrt(4) # 1
  prop = 0.5; mu = c(-0.5,0.5); sigma = c(0.25,0.25) # Gaussian
  
  # hyper-parameters
  terminalSize = 20; numberBin = 40; numberBinTest = 50; prior = "Gaussian"
  order = 10; numberSplit = 20; expand = 0.3
  if(setting == "GLM"){
    times = 300; best = 3; shrinkageEta = 0.02; alpha2 = 0.2; df = 1
  } else if (setting == "GM"){
    times = 200; best = 3; shrinkageEta = 0.02; alpha2 = 0.2; df = 2
  }
  
  # test data
  XTest = matrix(runif(nTest * d, -1, 1), ncol = d) 
  if(setting == "GLM"){
    yTest = (XTest %*% beta + beta0 + XTest[,1] * XTest[,2]) + rnorm(nTest,0,1) * (varBase + (XTest %*% varBeta))
  } else if (setting == "GM"){
    groupIndexTest = rbinom(nTest, 1, prop)
    yTest = groupIndexTest * rnorm(nTest, mu[1], (1+0.5*(XTest[,3])) * sigma[1]) + (1-groupIndexTest) * rnorm(nTest, mu[2], (1-0.5*(XTest[,3])) * sigma[2])
    yTest = XTest %*% beta/2 + beta0 + (XTest[,2]<=trueSplit1) * yTest + (XTest[,2]>trueSplit1) * rnorm(nTest, 0, sqrt(0.25 + 0.25^2))
  }
  splitPointYTest = seq(quantile(yTest, probs = 0)-expand, quantile(yTest, probs = 1)+expand, length.out = (numberBinTest + 1))
  splitMidPointYTest = (splitPointYTest[1 : numberBinTest] + splitPointYTest[2 : (numberBinTest + 1)])/2
  indexYTest = as.numeric(cut(yTest, breaks = c(-Inf, splitPointYTest[-c(1,numberBinTest+1)], Inf)))
  hTest = splitPointYTest[2] - splitPointYTest[1]
  
  
  # profile data
  if(setting == "GLM"){
    XProfile = matrix(runif(3^2 * d,-1,1), nrow = 3^2, ncol = d)
    XProfile[,1] = rep(c(-0.6,0,0.6), 3)
    XProfile[,2] = rep(c(0.6,0,-0.6), rep(3,3))
    XProfile = XProfile[rep(seq(1,9), rep(numberBinTest, 9)), ]
  } else if (setting == "GM"){
    XProfile = matrix(runif(3^2 * d,-1,1), nrow = 3^2, ncol = d)
    XProfile[1:(3^2-1),1] = rep(c(-0.5,0.5), 4)
    XProfile[1:(3^2-1),2] = rep(rep(c(-0.5,0.5), rep(2,2)),2)
    XProfile[1:(3^2-1),3] = rep(c(-0.5,0.5), c(4,4))
    XProfile[3^2,1:3] = rep(0,3)
    XProfile = XProfile[rep(seq(1,9), rep(numberBinTest, 9)), ]
  }
  yProfile = rep(splitMidPointYTest, 9)
  
  # record
  myGaussianNullLL = myCxxBoostingLL = rep(0, m)
  myCxxBoostingImportance = matrix(0, nrow = m, ncol = d)
  GaussianNullProfile = predictCxxBoostingProfile = matrix(0, nrow = m, ncol = 3^2 * numberBinTest)
  
  for(l in 1 : m){
    # data generation
    X = matrix(runif(n * d, -1, 1), ncol = d) 
    if(setting == "GLM"){
      y = (X %*% beta + beta0 + X[,1] * X[,2]) + rnorm(n,0,1) * (varBase + (X %*% varBeta))  } else if (setting == "GM"){
        groupIndex = rbinom(n, 1, prop)
        y = groupIndex * rnorm(n, mu[1], (1+0.5*(X[,3])) * sigma[1]) + (1-groupIndex) * rnorm(n, mu[2], (1-0.5*(X[,3])) * sigma[2])
        y = X %*% beta/2 + beta0 + (X[,2]<=trueSplit1) * y + (X[,2]>trueSplit1) * rnorm(n, 0, sqrt(0.25 + 0.25^2))
    }
    
    # candidate splits
    splitPoint = matrix(0, nrow = d, ncol = numberSplit + 1)
    for (i in 1:d) {splitPoint[i,] = quantile(X[,i], probs = seq(0,1,length.out = (1+numberSplit)))}
    splitPoint = splitPoint[,-c(1,dim(splitPoint)[2])]
  
    # null model
    # fit a Gaussian density on the whole data
    meanTemp = mean(y); sdTemp = sd(y)
    myGaussianNullLL[l] = mean(dnorm(yTest, mean = meanTemp, sd = sdTemp, log = TRUE))
    GaussianNullProfile[l,] = dnorm(yProfile, mean = meanTemp, sd = sdTemp)
    
    # LinCDE boosting
    myCxxBoosting = LinCDEBoosting(y = y, X = X, splitPoint = splitPoint, numberBin = numberBin, type = "nsTransform", order = order, df = df, alpha = alpha2, terminalSize = 20, depth = best, n.trees = times, shrinkage = shrinkageEta, minY = splitPointYTest[1], maxY = splitPointYTest[numberBinTest+1], prior = prior)
    predictCxxBoosting = LinCDEPredict(X = XTest, y = yTest, trees = myCxxBoosting, splitPointYTest = splitPointYTest)
    myCxxBoostingLL[l] = predictCxxBoosting$testLogLikelihood
    myCxxBoostingImportance[l,] = myCxxBoosting$importanceScore
    predictCxxBoostingProfile[l,] = LinCDEPredict(X = XProfile, y = yProfile, trees = myCxxBoosting, boosting = TRUE, splitPointYTest = splitPointYTest)$density
    
    if(l %% 1 == 0){print(l)}
  }
  
  # R^2
  if(setting == "GLM"){oracle = -0.6939}else if(setting == "GM") {oracle = -0.7725} 
  myCxxBoostingR2 = 1-(oracle - myCxxBoostingLL)/(oracle - myGaussianNullLL)
  print("LinCDE boosting R^2"); print(myCxxBoostingR2)
  
  # lattice plot of conditional densities
  if(setting == "GLM"){trueProfile = dnorm(yProfile, mean = (XProfile %*% beta + beta0 + XProfile[,1] * XProfile[,2]), sd = (varBase + (XProfile %*% varBeta)))
  } else if(setting == "GM"){
    trueProfile = rep(0, length(yProfile))
    index2 = which(XProfile[,2] < trueSplit1)
    trueProfile[index2] = (prop * dnorm(yProfile[index2], mean = XProfile[index2,] %*% beta + beta0 + mu[1] + 0.25, sd = 0.5 * (0.25 * XProfile[index2,3] + 0.5)) + (1-prop) * dnorm(yProfile[index2], mean = XProfile[index2,] %*% beta + beta0 + mu[2] + 0.25, sd = 0.5 * (-0.25 * XProfile[index2,3] + 0.5)))
    trueProfile[-index2] = (dnorm(yProfile[-index2], mean = 0.25 * XProfile[-index2,1], sd = sqrt(0.3))) 
  }
  latticeCDE = data.frame(rep(splitMidPointYTest, 9*2)); colnames(latticeCDE) = "y"
  latticeCDE$density = c(trueProfile, apply(predictCxxBoostingProfile, 2, mean))
  latticeCDE$method = factor(c(rep(c("truth", "LinCDE"), rep(numberBinTest*9,2))), levels=c("truth", "LinCDE"))
  latticeCDE$group = rep(rep(seq(1,9), rep(numberBinTest,9)), 2)
  # legend
  key=list(space="top", columns = 2, 
           lines=list(col=c("blue", "red"), lty=c(3,1), lwd=3),
           text=list(c("truth", "LinCDE")), cex = 1.5)
  if(setting == "GLM"){
    strip = strip.custom(bg="lightgrey", factor.levels = c(TeX('G_1: (-0.6,0.6)'),TeX('G_2: (0,0.6)'),TeX('G_3: (0.6,0.6)'),TeX('G_4: (-0.6,0)'),TeX('G_5: (0,0)'),TeX('G_6: (0.6,0)'),TeX('G_7: (-0.6,-0.6)'),TeX('G_8: (0,-0.6)'),TeX('G_9: (0.6,-0.6)')), par.strip.text=list(col="black", cex=1, font=1))
  } else if (setting == "GM"){
    strip = strip.custom(bg="lightgrey", factor.levels = c(TeX('G_1: (-0.5,-0.5,-0.5)'),TeX('G_2: (0.5,-0.5,-0.5)'),TeX('G_3: (-0.5,0.5,-0.5)'),TeX('G_4: (0.5,0.5,-0.5)'),TeX('G_5: (-0.5,-0.5,-0.5)'),TeX('G_6: (0.5,-0.5,0.5)')), par.strip.text=list(col="black", cex=0.8, font=1))
  }
  xyplot(density ~ y | factor(group), group = method, data = latticeCDE[which(latticeCDE$group %in% seq(1,9)),], type = "l", col = c("blue", "red"), lty = c(3,1), lwd = 3, ylab=list("conditional density", cex = 1.5), xlab = list("Y", cex = 1.5), key = key, aspect = 0.75, layout = c(3,3), strip = strip, scales = list(cex = 1))
}


