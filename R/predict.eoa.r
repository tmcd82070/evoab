#' @title predict.eoa
#'
#' @description Prediction method for eoa objects.
#'
#' @param object an eoa object produced by function \code{eoa}
#'
#' @param newdata A data frame containing new covariate
#' values at which to do the prediction.
#'
#' @param \dots Included for compatability with other predict methods.
#'
#' @return If \code{newdata} is NULL, the predicted
#' lambda values for every row of the fitting data set.
#' If \code{newdata} is not NULL, results of applying
#' the lambda model to covarates in \code{newdata}.
#'
#' @author Trent McDonald
#'
#' @export
predict.eoa <- function(object, newdata, ...){

  X <- model.matrix(object, newdata)
  preds <- X %*% matrix(coef(object),ncol=1)
  preds <- exp(preds)
  preds

}
