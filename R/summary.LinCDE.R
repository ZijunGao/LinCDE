#' summary.LinCDE
#'
#' This function provides the relative covariate influence/importance of a LinCDE boosting model.
#'
#' @param object a LinCDE boosting model.
#' @param cBars the number of bars to plot. If \code{order=TRUE} the only the variables with the \code{cBars} largest relative influence will appear in the barplot. If \code{order=FALSE} then the first \code{cBars} variables will appear in the plot. In either case, the function will return the relative influence of all of the variables.
#' @param plotit an indicator as to whether the plot is generated.
#' @param order an indicator as to whether the plotted and/or returned relative influences are sorted.
#' @param normalize if TRUE then returns the normalized influence summing to \code{100}.
#'
#' @return This function returns a data frame where the first component is the variable name and the second is the computed relative influence/importance.
#'
#' @export
summary.LinCDE = function(object, cBars = length(object$importanceScore), plotit = TRUE, order = TRUE, normalize = TRUE, ...){
  rel.inf = object$importanceScore
  if (order){index = order(-rel.inf)} else {index = seq(1,length(rel.inf))}
  if (cBars == 0){cBars = min(10, length(rel.inf))}
  if (cBars>length(rel.inf)){cBars = length(rel.inf)}
  if (normalize){rel.inf = 100 * rel.inf/max(1e-6, sum(rel.inf))}
  if (is.null(names(rel.inf))){names(rel.inf) = paste("X", seq(1, length(rel.inf)), sep = "")}
  if (plotit) {
      barplot(rel.inf[index[cBars:1]], horiz = TRUE,
              col = rainbow(cBars, start = 3/6, end = 4/6),
              names = object$var.names[index[cBars:1]], xlab = "Relative influence")
  }
  result = data.frame(var = names(rel.inf)[index], rel.inf = rel.inf[index])
  rownames(result) = NULL
  return(result)
}
