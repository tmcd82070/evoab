#' @export
#'
#' @title compVagueSd - Compute sD's of vague priors
#'
#' @description  Compute the standard deviation of a normal
#' distribution that is big enough to be considered a 'vague'
#' prior.  This is not straight forward when there are covariates, as here.
#'
#' @param Y Vector of number of carcasses found, one element per year. If
#' multiple sites are involved, elements of Y are the total (summed) number
#' of targets found per season.
#'
#' @param beta.params A data frame containing, at a minimum, two columns named \code{$alpha}
#' and \code{$beta}.
#' These are the annual alpha and beta parameters that determine the overall annual g.
#' Length of \code{$alpha} and \code{$beta} is either one, in which case g is constant across years,
#' or one element per year.
#'
#' @param Lprior Prior for lambda.  "uniform" uses a uniform[0,LMax] prior for lambda.
#' "normal" uses a normal(Lprior.mean,Lprior.sd) prior for lambda.
#'
#' @param Lprior.mean Mean of lambda prior when Lprior == "normal".
#'
#' @param Lprior.sd Standard deviation of normal when Lprior == "normal".
#'
#' @param LMax Maximum lambda when Lprior = "uniform".
#'
#' @param conf.level Confidence level for the confidence intervals on lambda.
#'
#' @return List containing two components.  \code{$lamda.ests} is a data
#' frame containing the lambda estiamtes (point est and confidence interval).
#' In this component \code{$Mtot} is the estimated number of targets over
#' all seasons.
#' \code{$out} is the full MCMC chain object.  Use this to check convergence, etc.
#'
#'
#' @examples
#'
compVagueSd<- function(Y,alpha.vec,beta.vec,X, range.multiplier=100){

  g <- alpha.vec / (alpha.vec+beta.vec)

  if( min(Y) == 0 ){
    m.c <- 0.5
  } else {
    m.c <- 0
  }
  M <- log(Y+m.c) - log(g + min(g)*m.c)

  lm.fit <- lm( M ~ -1+X )
  lm.fit <- summary(lm.fit)$coefficients

  coef.range <- abs(lm.fit[,"Estimate"]) + 2*lm.fit[,"Std. Error"]
  coef.range <- coef.range * range.multiplier

  coef.range
}


#tmp <- compVagueSd(Y, g$alpha, g$beta, eoa$design.mat)

