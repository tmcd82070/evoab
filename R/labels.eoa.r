#' @export
#'
#' @title labels.eoa - Extract parameter labels from \code{eoa} objects.
#'
#' @param obj An \code{eoa} model object.  See \code{\link{eoa}}.
#'
#' @param type The type of parameter label requred or a regular expression.
#' Parameter type possibilities are "coef" or
#' "derived".  Regular expressions are used to match parameter labels using \code{grep}.
#' If \code{type} is not "coef" and not "derived" and the regular expression
#' fails to match anything, all parameter labels.
#'
#' @details Coefficient labels are for variables in the log-linear model for
#' lambda.  Derived parameter labels are non-coefficient parameters link
#' M, Mtot, and lambda.
#'
#' @author Trent McDonald
#'
#' @seealso \code{\link{coef}}, \code{\link{eoa}}
#'
#' @examples
#' # A 3 year study of 7 sites. 21 "cells". lambda change = 20/year
#' set.seed(9430834) # fixes Y and g of this example, but not the RNG's used in chains
#' ns <- 3
#'
#' ny <- 7
#' g <- data.frame(
#'  alpha = rnorm(ns*ny,70,2),
#'  beta = rnorm(ns*ny,700,25)
#' )
#' Y <- rbinom(ns*ny, c(rep(20,ny), rep(40,ny), rep(60,ny)), g$alpha/(g$alpha+g$beta))
#'
#' df <- data.frame(year=factor(c(rep("2015",ny),rep("2016",ny),rep("2017",ny))),
#'    Year=c(rep(1,ny),rep(2,ny),rep(3,ny)))
#'
#' # Uninformed eoa (use low number of iterations because it's and example)
#' eoa.1 <- eoa(Y~year, g, df, nburn = 1000, niters= 50*10, nthins = 10 )
#'
#' labels(eoa.1)
#' labels(eoa.1,"derived")
#' labels(eoa.1,"^M")  # all M parameters
#' labels(eoa.1,"\\[3\\]$") # M[3] and lambda[3]
#' labels(eoa.1,".") # all parameter labels
#'
#' plot(eoa.1$out[,labels(eoa.1)])  # trace plot of coefficients.
#'
labels.eoa <- function(obj, type="coef"){

  if(type=="coef"){
    ans <- obj$coef.labels
  } else {
    ans <- dimnames(obj$estimates)[[1]]
    if(type=="derived"){
      ans <- ans[!(ans %in% obj$coef.labels)]
    } else {
      ind <- grep(type, ans)
      if(length(ind)>0){
        ans <- ans[ind]
      }
    }
  }
  ans
}
