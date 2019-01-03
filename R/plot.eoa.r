#' @export
#'
#' @title Plot rate estimates from an \code{eoa} model
#'
#' @description Plot method for \code{eoa} objects.
#'
#' @param obj An object of class "eoa".  Usually output by the
#' routine \code{eoa}.
#'
#' @param xvar The x-axis variable to plot.  Points in the plot
#' are carcasses/g/numTurbines ~ xvar and line in plot is lambda ~ xvar.
#' xvar must be in the model you want and in the data frame contained
#' in the model object. All other covariates in the model are held
#' constant at their mean value.
#'
#' @param main  The main title of the plot
#'
#' @author Trent McDonald
#'
#' @seealso \code{\link{plot.Mest}}
#'
#'
plot.eoa <- function(obj, xvar=NULL, main=NULL){

  if(is.null(xvar)){
    # xaxis is first variable after intercept in design mat
    xvar <- attr(obj$terms, "term.labels")[1]
  }

  df <- obj$data
  df$g <- df$gFac.a / (df$gFac.a+df$gFac.b)

  xmat <- model.matrix(obj)[,-1, drop=FALSE]
  xmeans <- colMeans(xmat)
  xmeans <- xmeans[names(xmeans)!=xvar]
  nms <- names(xmeans)

  xvar.seq <- seq(min(df[,xvar]), max(df[,xvar]), length=70)
  xmeans <- data.frame(matrix(xmeans, length(xvar.seq),
                              length(xmeans),byrow=T),
                       xvar.seq )
  names(xmeans) <- c(nms, xvar)

  preds <- predict(obj, newdata=xmeans)

  plot(df[,xvar],df$carcasses/df$g/exp(df$logOffset),
       pch=16, col="black", xlab=xvar, ylab="Lambda/turbine/year",
       main=main)
  lines(xvar.seq, preds, col="red",lwd=2)

  data.frame(xmeans, predicted=preds)

}
