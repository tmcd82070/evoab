#' @title predict.eoa
#'
#' @description Prediction method for eoa objects.
#'
#' @param object an eoa object produced by function \code{eoa}
#'
#' @param newdata A data frame containing new covariate
#' values at which to do the prediction.
#'
#' @param type Text string specifying the type of prediction.
#' \itemize{
#'   \item \code{"response"} returns the predicted value of the eoa rate parameter,
#' lambda, evaluated using the median coefficient estimates
#' in \code{object}. \code{"response"}
#' evaluates the model using \code{object$estimates[object$coef.labels]}
#' (or \code{coef(object)}) as coefficients, and returns the result after
#' taking the anti-log (exp()).
#'   \item \code{"quantile"} returns a quantile from the distribution of
#'   the eoa rate parameter, lambda, evaluated using all the posterior samples
#'   of the coefficients.  \code{"quantile"} evaluates the model using
#'   all posterior samples of the coefficients housed in \code{object$out},
#'   and returns the \code{p}-th quantile(s) from the distribution of
#'   these predicted values (in anti-log space). This provides a way to
#'   compute predictive intervals for the response.  Setting
#'   \code{type="quantile"} and \code{p=c(0.025, 0.975)} produces a
#'   95% predictive interval for the predicted values.
#' }
#'
#' @param p A vector specifing which quantiles to return when
#' \code{type="quantile"}.
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
predict.eoa <- function(object, newdata, type="response", p=c(0.025,0.5,0.975), ...){

  X <- model.matrix(object, newdata)
  if(type=="response"){
    preds <- X %*% matrix(coef(object),ncol=1)
    preds <- exp(preds)
  } else {
    betaMat <- t(as.matrix(object$out[,object$coef.labels]))
    preds <- X %*% betaMat
    preds <- exp(preds)
    preds <- t(apply( preds, 1, quantile, p=p))
  }
  preds

}
