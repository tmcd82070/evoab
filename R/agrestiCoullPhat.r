#' @export
#'
#' @title Agresti-Coull adjustment to binomial proportion
#'
#' @description Applies the Agresti and Coull (1998) adjustment
#' to an observed number of 'successes' out of a number of 'trials'.
#'
#' @param x The number of 'successes' (an integer vector)
#'
#' @param n The number of 'trials' (an integer vector)
#'
#' @param conf Confidence interval level (a scalar between 0.5 and 1)
#'
#' @details
#' The Agresti-Coull adjusted point estimate of a binomial proportion
#' is,
#' \deqn{phat = (x+(z^2)/2) / (n+z^2)}
#' where \code{z} is the \code{(1-(1-conf)/2)} quantile of a standard
#' normal distribution.  This point estimator was actually proposed by
#' Brown et al. (2001) but named for Agresti and Coull because the latter
#' discussed performance of the confidence interval under the
#' special case of \code{conf = 0.95} (i.e., \code{z = 1.96}).
#'
#' The estimated standard error of the Agresti-Coull adjusted proportion
#' is,
#' \deqn{se = sqrt( phat*(1-phat) / (n+z^2)).}
#'
#' The Agresti-Coull confidence interval is,
#' \deqn{phat +- z*se.}
#'
#' The Agresti-Coull confidence interval has
#' excellent coverage for ratios between
#' approximately 2% and 98%. For ratios between 0% and 2% (and 98% and 100%)
#' interval coverage
#' is too high (i.e., coverage exceeds \code{conf} substantially).
#' The Agresti-Coull point estimator is biased high except when true p = 0.5.
#' Simulations by the author of this routine suggest bias
#' in the point estimate is <10% for true p's between approximately
#' 0.15 and 0.85. For ratios between 0 and 0.02 (and 0.98 and 1), bias
#' in the point estimate can exceed 100%, but these are small (and large)
#' ratios so percent bias is sensitive to small changes. The author of
#' this routine is not aware of any studies of the variance estimator.
#'
#'
#'
#'
#' @return A data frame containing the Agresti-Coull adjusted proportion,
#' standard error, and confidence limits. The data frame has the
#' following columns:
#' \enumerate{
#'   \item \code{phat} : the Agresti-Coull adjusted point estimates.
#'   \item \code{se.phat} : the Agresti-Coull estimated standard
#'   error of the point estimates \code{phat}.
#'   \item \code{ll.phat} : the lower limit of a \code{100*(1-(1-conf)/2)}%
#'   confidence interval for \code{phat}.
#'   \item \code{ul.phat} : the upper limit of a \code{100*(1-(1-conf)/2)}%
#'   confidence interval for \code{phat}.
#' }
#'
#' @author Trent McDonald
#'
#' @references
#' Agresti, A. and B. A. Coull. 1998. Approximate is Better than
#' "Exact" for Interval Estimation of Binomial Proportions.
#' The American Statistician 52: 119:126.
#'
#' Brown, L. D., T. T. Cai, and A. DasGupta. 2001. Interval Estimation
#' for a Binomial Proportion. Statistical Science 16:101-133.
#'
#' @seealso \code{\link{bayesPhat}}
#'
#' @examples
#' agrestiCoullPhat(0:5, 100)
#'
#' # Simulation: point est bias and ci coverage
#' trueP <- 0.01
#' n <- 1000
#' x <- rbinom( 1000, n, trueP)
#' agPhat <- agrestiCoullPhat( x, n )
#' muAG <- mean(agPhat$phat)
#' covAG <- mean(agPhat$ll.phat <= trueP & trueP <= agPhat$ul.phat)
#' agStats <- c(mean=muAG,
#'              relBias = abs(muAG-trueP)/trueP),
#'              coverage = covAG)
#'
agrestiCoullPhat <- function(x,n,conf=.9){
  z <- qnorm(1-(1-conf)/2,0,1)
  x <- x + z^2/2
  n <- n + z^2
  phat <- x / n

  se <- sqrt(phat*(1-phat)/n)

  ll <- phat - z*se
  ul <- phat + z*se

  data.frame(phat=phat, se.phat=se , ll.phat=ll, ul.phat=ul)
}
