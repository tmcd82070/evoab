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
#' @param alpha.vec vector of the alpha parameters of the Beta distributions
#'
#' @param beta.vec vector of the beta parameters of the Beta distributions
#'
#' @param X a design matrix upon which an approximation of inflated
#' numbers of targets is regressed.  Usually, this is meant to be the
#' design matrix from the \code{eoa} function.
#'
#' @param range.multiplier a multipiler for the range of coefficient
#' estimates to make the output standard deviations sufficiently vague.
#' Increasing this number increases vagueness.
#'
#' @return List containing two components.  \code{$vagueSd} is a
#' vector, one per parameter, of standard deviations that should
#' be large enough to call vague when used in a normal prior.
#' \code{$startA} is a vector of potential starting values for
#' the coefficients in the model.
#'
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
  if( sum(!is.na(lm.fit$coefficients)) != ncol(X) ){
    stop("Multicolinearity in the model.")
  }

  lm.fit <- summary(lm.fit)$coefficients

  if( is.na(lm.fit[,"Std. Error"]) ){
    # in absence of coef var estimate, seems reasonable to
    # assume a CV of 100%
    lm.fit[,"Std. Error"] <- lm.fit[,"Estimate"]
  }
  coef.range <- 2*lm.fit[,"Std. Error"]
  coef.range <- coef.range * range.multiplier

  # names of coef.range must match coefficient names because
  # we use them to subset outside this routine.
  names(coef.range) <- dimnames(X)[[2]]

  list(vagueSd=coef.range, startA=lm.fit)
}


#tmp <- compVagueSd(Y, g$alpha, g$beta, eoa$design.mat)

