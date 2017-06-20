#' @export
#'
#' @title Baysian estimation to binomial proportion using a conjugate
#' prior.
#'
#' @description This routine assumes a beta(a,b) prior for a binomial
#' proportion and computes posterior point and inteval estimates
#' for the proportion given
#' an observed number of 'successes' and number of 'trials'.
#'
#' @param x The number of 'successes' (an integer vector)
#'
#' @param n The number of 'trials' (an integer vector)
#'
#' @param a A vector of shape parameters to specify the beta prior
#'
#' @param b A vector of scale parameters to specify the beta prior
#'
#' @param conf A scalar specifying the desired confidence level
#' of the interval estimate.  This scalar should be between 0.5 and 1.
#'
#' @param ci A character string specifying the type of posterior
#' confidence interval to compute.  A value of 'eq' computes
#' posterior intervals with equal area in both tails.  A value
#' of 'hdi' computes highest density interval estimates.  See Details.
#'
#' @details
#'
#' Computations are elementary because the beta(a,b) prior is
#' conjugate for the binomial. Nearly every text on Bayesian
#' estimation shows that given values for \code{x}, \code{n},
#' \code{a}, and \code{b}, the posterior distribution of p
#' is,
#' \deqn{beta(x+a, n-x+b).}
#' Hence, the Bayes point estimator is,
#' \deqn{phat = (x+a) / (n+a+b),}
#' which is the mean of the posterior. Standard error of the
#' posterior is,
#' \deqn{se.phat=sqrt(a*b/((a+b)^2*(a+b+1))).}
#' If \code{a>1} and \code{b>1}, mode of the posterior for p is,
#' \deqn{(a-1)/(a+b-2).}
#'
#'
#' Confidence intervals can be computed two ways. The default,
#' \code{ci='eq'}, puts \code{(1-conf)/2} probability in the
#' lower tail and \code{(1-conf)/2} in the upper tail.
#' That is, the lower limit is
#' \code{qbeta((1-conf)/2,x+a,n-x+b)} and the upper limit is
#' \code{qbeta(1-(1-conf)/2,x+a,n-x+b)}.  If \code{ci='hdi'},
#' this routine finds the interval [l,u] such that
#' \code{pbeta(u,x+a,n-x+b) - pbeta(l,x+a,n-x+b)} and
#' \code{dbeta(u,x+a,n-x+b) = dbeta(l,x+a,n-x+b)}.
#'
#' The default values for \code{a} and \code{b}
#' imply use of the Jeffery's prior.  The Jeffery's prior is
#' proportional to the root of Fisher's information and
#' is equal to beta(0.5,0.5). A flat prior is beta(1,1).
#'
#' @return A data frame with number of rows equal to
#' \code{max(length(x), length(n), length(a), length(b)}
#' containing the Baysian point and interval estimates.
#' The data frame has the
#' following columns:
#' \enumerate{
#'   \item \code{phat} : the Bayes point estimates equal to the
#'   mean of the posterior distribution.
#'   \item \code{phat.mode} : if \code{a>1} and \code{b>1}, this
#'   column contains the mode of the posterior. If either
#'   \code{a<=1} or \code{b<=1}, the posterior is either multi-modal
#'   or its mode is infinity.  In this case, \code{phat.mode = NA}.
#'   \item \code{se.phat} : the standard
#'   error of the posterior distribution.
#'   \item \code{ll.phat} : the lower limit of a \code{100*(1-(1-conf)/2)}%
#'   posterior confidence interval for \code{phat}. See Details for
#'   interpretation under the allowed values of \code{ci}.
#'   \item \code{ul.phat} : the upper limit of a \code{100*(1-(1-conf)/2)}%
#'   posterior confidence interval for \code{phat}. See Details for
#'   interpretation under the allowed values of \code{ci}.
#' }
#'
#' @author Trent McDonald
#'
#'
#' @seealso \code{\link{agrestiCoullPhat}}
#'
#' @examples
#' bayesPhat(0:5, 100, 0.5, 0.5)  # Jeffery's prior
#' bayesPhat(0:5, 100, 1, 1) # flat prior
#'
#' # When prior data is available:
#' x.prior <- 5
#' n.prior <- 100
#' bayesPhat(0:5, 100, x.prior+0.5, n.prior-x.prior+0.5)
#'
#' # Simulation: point est bias and ci coverage
#' trueP <- 0.01
#' n <- 1000
#' x <- rbinom( 1000, n, trueP)
#' baPhat <- bayesPhat( x, n )
#' muBA <- mean(baPhat$phat)
#' covBA <- mean(baPhat$ll.phat <= trueP & trueP <= baPhat$ul.phat)
#' baStats <- c(mean=muBA,
#'              relBias = abs(muBA-trueP)/trueP),
#'              coverage = covBA)
#'
#'
#'
#'
bayesPhat <- function(x, n, a=0.5, b=0.5, conf=0.9, ci="eq"){
  aPost <- x + a
  bPost <- n - x + b
  phat <- aPost / (aPost + bPost)


  se <- sqrt(aPost*bPost / ((aPost+bPost)^2*(aPost+bPost+1)))

  if((aPost>1) & (bPost>1)){
    phat.mode <- (aPost-1) / (aPost + bPost - 2)
  } else {
    phat.mod <- NA
  }

  ll <- qbeta((1-conf)/2,aPost,bPost)
  ul <- qbeta(1-(1-conf)/2,aPost,bPost)
  if(ci=="hdi"){
      if((aPost>1) & (bPost>1)){
        fl <- function(q, a, b, p){
          (dbeta(q,a,b) - p)^2
        }
        hdi.obj <- function(p, a, b, conf, phat.mode){
          q1 <- optim(par=phat.mode/2,fn=fl,a=a, b=b, p=p,
                      method="Brent",lower=0, upper = phat.mode)$par
          q2 <- optim(par=(1+phat.mode)/2,fn=fl,a=a, b=b, p=p,
                      method="Brent",lower=phat.mode, upper = 1)$par
          hdi.prob <- pbeta(q2,a,b) - pbeta(q1,a,b)
          (hdi.prob - conf)^2
        }
        fp <- optim(par=dbeta(phat.mode,aPost,bPost)/2, fn=hdi.obj,
                    a=aPost, b=bPost, conf=conf, phat.mode=phat.mode,
                    method="Brent", lower=1e-10,
                    upper=dbeta(phat.mode,aPost,bPost))$par
        ll <- optim(par=ll,fn=fl,a=aPost, b=bPost, p=fp,
                    method="Brent",lower=0, upper = phat.mode)$par
        ul <- optim(par=ul,fn=fl,a=aPost, b=bPost, p=fp,
                    method="Brent",lower=phat.mode, upper = 1)$par
    } else {
      ll <- ul <- NA
    }
  }

  data.frame(phat=phat, se.phat=se, phat.mode=phat.mode, ll.phat=ll, ul.phat=ul)
}
