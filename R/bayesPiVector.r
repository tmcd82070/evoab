#' @export
#'
#' @title Baysian estimation of a multinomial
#' proportion vector assuming a conjugate
#' prior.
#'
#' @description This routine assumes you have an observation
#' from a multinomial distribution with k classes (k >= 2, parameter
#' \code{x}) and have assumed that the multinomial distribution's
#' proportion vector ("pi vector") follows a Dirichelet distribution.
#' If so, this routine estimates the proportion vector's
#' posterior distribution mean, variance, and mode.
#'
#' @param x An integer vector
#' containing the number of observed 'successes' in
#' each catagory of the multinomial. Total number of trials is \code{sum(x)}.
#' The number of catagories is K = \code{length(x)}.
#'
#' @param pseudoCounts A vector of real-valued "pseudo counts"
#' for the K catagories in the problem. This is sometimes called
#' the "concentration" parameter.
#'
#'
#' @details
#'
#' Computations are elementary because the Dirichlet(a1, a2, ..., aK) prior is
#' conjugate for the multinomial. Nearly every text on Bayesian
#' estimation shows that given values for \code{x} and \code{pseudoCounts},
#' the posterior distribution of the mulitinomial's p vector
#' is,
#' \deqn{Dirichlet(x1+a1, x2+a2, ..., xk+ak).}
#' Hence, the Bayes point estimator of the multinomial's proportions is,
#' \deqn{phat_i = (xi+ai) / sum(xi + ai),}
#' which is the mean of the posterior. Standard error of the
#' posterior is,
#' \deqn{se.phat_i=sqrt((xi+ai)*(A-ai)/(A^2*(A+1))).}
#' where A = sum(xi + ai).  If \code{(xi+ai)>1} for all i, mode of the
#' posterior for the proportion vector is,
#' \deqn{(xi+ai-1)/(A-K).}
#'
#'
#' The default value for \code{pseudoCounts}
#' corresponds to the Jeffery's prior.  The Jeffery's prior is
#' proportional to the root of Fisher's information and
#' is equal to Dirichlet(1/K,1/K, ..., 1/K).
#'
#' @return A data frame with number of rows equal to
#' \code{length(x)}
#' containing the Baysian point estimates for the proportion
#' in each catagory.
#' The data frame has the
#' following columns:
#' \enumerate{
#'   \item \code{phat} : the Bayes point estimates equal to the
#'   mean vector of the posterior distribution. This column sums to 1.0
#'   \item \code{phat.mode} : if \code{xi+ai} > 1 for all i, this
#'   column contains the mode vector of the posterior. Mode vector
#'   is the most vector of proportions with maximum likelihood.
#'   If any \code{xi+ai} < 1, \code{phat.mode = NA}.
#'   \item \code{se.phat} : the standard
#'   error vector of the posterior distribution.
#'   \item \code{psuedoCounts} : the vector of pseudoCounts
#'   associated with the Dirichlet posterior.  This vector
#'   can be used to accumulate counts over muliple calls.
#' }
#'
#' @author Trent McDonald
#'
#'
#' @seealso \code{\link{agrestiCoullPhat}}
#'
#' @examples
#' bayesPiVector(c(1,5), c(.5,.5))  # Jeffery's prior
#' bayesPiVector(c(1,5), c(1, 1)) # flat prior
#'
#' # When prior data is available:
#' x.prior <- 5
#' n.prior <- 100
#' bayesPiVector(c(1,5), c(x.prior+0.5, n.prior-x.prior+0.5))
#'
#' # Simulation: point est bias and ci coverage
#' trueP <- c(0.01, 0.04, 0.95)
#' n <- 20
#' x <- rbinom( 1000, n, trueP)
#' baPhat <- apply(x, 1, bayesPiVector, pseudoCounts=rep(1,3)/3 )
#' muBA <- mean(baPhat$phat)
#'
bayesPiVector <- function(x, pseudoCounts=rep(1,length(x))/length(x)){
  aPost <- x + pseudoCounts
  A <- sum(aPost)
  phat <- aPost / A

  se <- sqrt(aPost*(A-aPost) / (A^2*(A+1)))

  if( all(aPost>1) ){
    phat.mode <- (aPost-1) / (A - length(aPost))
  } else {
    phat.mode <- NA
  }

  data.frame(phat=phat, se.phat=se, phat.mode=phat.mode, pseudoCounts=aPost)
}
