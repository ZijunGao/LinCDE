#' print.LinCDE
#'
#' This function displays basic information about a LinCDE model (object).
#' @method print LinCDE
#'
#' @param x a LinCDE model (object).
#' @param ... other parameters passed to print.default.
#'
#' @return The function prints a list of values.
#' \itemize{
#' \item model: a LinCDE model. The response is denoted by \code{y}. If the covariates' names are not provided in the training stage, the covariates will be named in the format "Xj".
#' \item hyperparameters including depth, shrinkage, basis, centering (centeringMethod if centering is \code{TRUE}), and function of the prior are provided.
#' \item importanceScore: relative influences of the covariates to the response distribution.
#' }
#'
#' @export
print.LinCDE = function (x, ...){
  print.list = list()
  print.list$model = paste("y ~", paste(x$var.names, collapse = " + "))
  print.list$n.trees = length(x$trees)
  print.list = append(print.list, x)
  print.list$trees = NULL
  print.list$y = NULL; print.list$X = NULL; print.list$var.names = NULL
  print.list$z = NULL; print.list$zTransformMatrix = NULL
  print.list$penalty = NULL
  print.list$splitMidPointY = NULL
  if(!print.list$centering){print.list$centeringMethod = NULL}
  print(print.list, ...)
}
