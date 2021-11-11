#' summary.LinCDE
#'
#' This function provides relative covariate influences/importances of a LinCDE model.
#'
#' @param object a LinCDE model.
#' @param cBars the number of bars to plot. If \code{order=TRUE} the only the variables with the \code{cBars} largest relative influence will appear in the barplot. If \code{order=FALSE} then the first \code{cBars} variables will appear in the plot. In either case, the function will return the relative influence of all of the variables.
#' @param plotit an indicator as to whether the plot is generated.
#' @param order an indicator as to whether the plotted and/or returned relative influences are sorted.
#' @param normalize if \code{TRUE} then returns the normalized influence summing to \code{100}.
#' @param ... other parameters.
#'
#' @return If the \code{object$centering} is \code{FALSE}, this function returns a data frame where the first component is the variable name and the second is the computed relative influence/importance. If the \code{object$centering} is \code{TRUE}, this function returns two data frames \code{location} and \code{beyondLocation}. In each data frame, the first component is the variable name and the second is the computed relative influence/importance. The \code{location} data frame records the relative influence/importance from \code{object$centeringModel}, and the \code{beyondLocation} data frame records the relative influence/importance from the LinCDE model fitted to the residuals of \code{object$centeringModel}.
#'
#' @export
summary.LinCDE = function(object, ..., cBars = length(object$importanceScore), plotit = TRUE, order = TRUE, normalize = TRUE){
  if(is.null(object$importanceScore)){stop("no relative covariate importance available. Probably due to a zero-depth LinCDE model.")}
  rel.inf = object$importanceScore
  if (order){index = order(-rel.inf)} else {index = seq(1,length(rel.inf))}
  if (cBars == 0){cBars = min(10, length(rel.inf))}
  if (cBars>length(rel.inf)){cBars = length(rel.inf)}
  if (normalize){rel.inf = 100 * rel.inf/max(1e-6, sum(rel.inf))}
  if (plotit) {
    if(object$centering){
      if(class(object$centeringMethod) == "character"){
        rel.inf.centering = caret::varImp(object$centeringModel)[names(rel.inf),1]
        ratio = max(0, object$errorReduction["centering"]/object$errorReduction["beyond.centering"])
        rel.inf.total = rbind(rel.inf.centering/sum(rel.inf.centering) * sqrt(ratio), rel.inf/sum(rel.inf))
        rownames(rel.inf.total) = c("location", "beyond location")
        colnames(rel.inf.total) = object$var.names
        if (order){index = order(-(apply(rel.inf.total, 2, sum)))} else {index = seq(1,dim(rel.inf.total)[2])}
        if (normalize){rel.inf.total = 100 * rel.inf.total/max(1e-6, sum(rel.inf.total))}
        barplot(rel.inf.total[1,index[cBars:1]], horiz = TRUE,
                col = rev(rainbow(cBars, start = 0/6, end = 1/6)),
                names = colnames(rel.inf.total)[index[cBars:1]],
                xlab = "Relative influence",
                xlim = c(0, max(apply(rel.inf.total, 2, sum))), las = 2)
        barplot(rel.inf.total[2,index[cBars:1]],
                offset=rel.inf.total[1,index[cBars:1]],
                add = T, axes = F, axisnames = F, horiz=T,
                col = rbind(rainbow(cBars, start = 3/6, end = 4/6)), las = 2)
        legend("bottomright", legend = c("location", "beyond location"), col = c("red", "blue"), lwd = 10)
      } else {
        barplot(rel.inf[index[cBars:1]], horiz = TRUE,
                col = rainbow(cBars, start = 3/6, end = 4/6),
                names = object$var.names[index[cBars:1]], xlab = "Relative influence", las = 2)
      }
    } else {
      barplot(rel.inf[index[cBars:1]], horiz = TRUE,
              col = rainbow(cBars, start = 3/6, end = 4/6),
              names = object$var.names[index[cBars:1]], xlab = "Relative influence (beyond location)", las = 2)
    }
  }

  if(object$centering){
    result = list()
    if(class(object$centeringMethod) == "character"){
      result$location = data.frame(var = colnames(rel.inf.total)[index], rel.inf = rel.inf.total[1,index])
      rownames(result$location) = NULL
      result$beyondLocation = data.frame(var = colnames(rel.inf.total)[index], rel.inf = rel.inf.total[2,index])
      rownames(result$beyondLocation) = NULL
    } else {
      print("no available relative covariate importance from the centering model")
      result$beyondLocation = data.frame(var = names(rel.inf)[index], rel.inf = rel.inf[index])
      rownames(result$beyondLocation) = NULL
    }
  } else{
    result = data.frame(var = names(rel.inf)[index], rel.inf = rel.inf[index])
    rownames(result) = NULL
  }
  return(result)
}
