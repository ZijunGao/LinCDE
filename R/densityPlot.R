#' Lattice plot of estimated densities
#'
#' This function produces conditional density plots at input covariate values.
#' @param X covariate matrix to predict conditional densities at; each row represents an observation vector.
#' @param trees: a LinCDE boosting model. To only plot the true conditional densities, leave the \code{trees} blank. Default is NULL.
#' @param trueDensity true conditional density function. The true conditional densities at a covariate response pair (X,y) can be computed as \code{trueDensity(X,y)}. To only plot LinCDE's estimated conditional densities, leave the \code{trueDensity} blank. Default is NULL.
#' @param minY the minimal possible response value. If \code{trees} is not NULL, then the minimal value of the \code{trees}'s response grid will be used. Default is 0.
#' @param maxY the maximal possible response value. If \code{trees} is not NULL, then the maximal value of the \code{trees}'s response grid will be used. Default is 0.
#'
#' @export
densityPlot = function(X, minY = 0, maxY = 0, trees = NULL, trueDensity = NULL){
  # preprocess
  if(is.null(dim(X))){X = matrix(X, nrow = 1)}
  n = dim(X)[1]; d = dim(X)[2]
  XProfile = X[rep(seq(1,n), rep(100, n)), ]
  if(!is.null(trees)){
    splitMidPointY = trees$splitMidPointY
    h = splitMidPointY[2] - splitMidPointY[1]
    splitPointYTest = seq(min(splitMidPointY) - h, max(splitMidPointY)+h, length.out = 101)} else {splitPointYTest = seq(minY, maxY, length.out = 101)}
  splitMidPointYTest = (splitPointYTest[1:100] + splitPointYTest[2:101])/2
  yProfile = rep(splitMidPointYTest, n)
  if(!is.null(trueDensity)){
    if(!is.null(trees)){
      trueProfile = trueDensity(XProfile, yProfile)
      LinCDEBoostingProfile = matrix(LinCDEPredict(X = XProfile, y = yProfile, trees = trees, splitPointYTest = splitPointYTest)$density, nrow = 1)
      # lattice plot of the true density
      latticeCDE = data.frame(rep(splitMidPointYTest, n*2)); colnames(latticeCDE) = "y"
      latticeCDE$density = c(trueProfile, apply(LinCDEBoostingProfile, 2, mean))
      latticeCDE$method = factor(c(rep(c("truth", "LinCDE"), rep(100*n,2))), levels=c("truth", "LinCDE"))
      latticeCDE$group = rep(rep(seq(1,n), rep(100,n)), 2)
      key=list(space="top", columns = 2,
               lines=list(col=c("blue", "red"), lty=c(3,1), lwd=3),
               text=list(c("truth", "LinCDE")), cex = 1.5)

      strip = lattice::strip.custom(bg="lightgrey", factor.levels = paste("R", seq(1,n)), par.strip.text=list(col="black", cex=0.8, font=1))

      plot = lattice::xyplot(density ~ y | factor(group), group = method, data = latticeCDE, type = "l", col = c("blue", "red"), lty = c(3,1), lwd = 3, ylab=list("conditional density", cex = 1.5), xlab = list("Y", cex = 1.5), key = key, aspect = 0.75, layout = c(3,ceiling(n/3)), strip = strip, scales = list(cex = 1), main = "LinCDE's estimated conditional densities")
    } else {
      trueProfile = trueDensity(XProfile, yProfile)
      # lattice plot of the true density
      latticeCDE = data.frame(rep(splitMidPointYTest, n)); colnames(latticeCDE) = "y"
      latticeCDE$density = trueProfile
      latticeCDE$group = rep(seq(1,n), rep(100,n))

      strip = lattice::strip.custom(bg="lightgrey", factor.levels = paste("R", seq(1,n)), par.strip.text=list(col="black", cex=0.8, font=1))

      plot = lattice::xyplot(density ~ y | factor(group), data = latticeCDE, type = "l", col = c("blue"), lty = 3, lwd = 3, ylab=list("conditional density", cex = 1.5), xlab = list("Y", cex = 1.5), aspect = 0.75, layout = c(3,ceiling(n/3)), strip = strip, scales = list(cex = 1), main = "true conditional densities")
    }
  } else {
    LinCDEBoostingProfile = matrix(LinCDEPredict(X = XProfile, y = yProfile, trees = trees, splitPointYTest = splitPointYTest)$density, nrow = 1)

    latticeCDE = data.frame(rep(splitMidPointYTest, n)); colnames(latticeCDE) = "y"
    latticeCDE$density = apply(LinCDEBoostingProfile, 2, mean)
    latticeCDE$group = rep(seq(1,n), rep(100,n))

    strip = lattice::strip.custom(bg="lightgrey", factor.levels = paste("R", seq(1,n)), par.strip.text=list(col="black", cex=0.8, font=1))

    plot = lattice::xyplot(density ~ y | factor(group), data = latticeCDE, type = "l", lty = 1, lwd = 3, ylab=list("conditional density", cex = 1.5), xlab = list("Y", cex = 1.5), col = "red", aspect = 0.75, layout = c(3,ceiling(n/3)), strip = strip, scales = list(cex = 1), main = "LinCDE's estimated conditional densities")
  }
  print(plot)
}
