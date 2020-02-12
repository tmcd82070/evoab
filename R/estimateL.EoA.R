#' @title estimateL.EoA - Estimate rate parameter Lambda for
#' a single-site using the EoA method
#'
#' @description This routine computes lambda, the mean number of search targets
#' out there per season,
#' using information from the number of found targets and the g-value (=probatility of
#' discovery).  The method is Bayesian and allows either an uniform prior for lambda
#' or an informed prior.
#' Estimation is direct in the sense that this routine uses numerical
#' integration to compute the posterior of lambda.
#'
#' @param X Total number of search targets found at all searched sites during the
#' entire search season.
#'
#' @param beta.params A list containing, at a minimum, components named $alpha and $beta.
#' These are the all-site alpha and beta parameters for g. In many cases, these parameters
#' are computed using function \code{\link{getFleetG}}.
#'
#' @param Lprior A string naming the prior distribution to use for lambda.
#' The following priors are implimented:
#' \enumerate{
#'    \item "normal" : uses a normal(\code{Lprior.mean},\code{Lprior.sd}) prior for lambda.
#'    \item "gamma" : uses a gamma(alpha,beta) prior for lambda, where alpha =
#' \code{Lprior.mean^2}/\code{Lprior.sd^2} and beta = \code{Lprior.mean}/
#' \code{Lprior.sd^2}. That is, alpha and beta are the method of moment estimates
#' for the shape and rate parameter of a gamma distribution.
#'    \item "jefferys" : uses a beta(L + 0.5, 0.5) function as the prior for lambda.  This prior
#'    is improper and really close to the actual Jeffery's prior for a poisson random variable.
#'    This is the Jeffery's prior implemented by Dalthorp's eoa package.
#' }
#'
#' @param Lprior.mean Mean of lambda prior when Lprior == "normal" or "gamma".
#'
#' @param Lprior.sd Standard deviation of normal when Lprior == "normal" or "gamma".
#'
#'
#' @param conf.level Confidence level for the confidence intervals on lambda.
#'
#'
#' @return List containing two components:
#' \itemize{
#'   \item  \code{$L.ests} is a data
#' frame containing the lambda estiamtes (point estimate and confidence interval).
#'   \item \code{$L.posterior} is a data frame containing the
#'   posterior, posterior cdf, prior, and likelihood. This is returned in case
#'   you want to plot them.
#' }
#'
#' @examples
#'
#' syr <- data.frame(species=c("LBBA","LBBA","LBBA"),
#'    facility=c("f1","f2","f2"),
#'    gFac.a = c( 69.9299, 63.5035,  84.6997),
#'    gFac.b = c(  736.4795,  318.3179, 759.9333 ),
#'    year = c(2015,2015,2016))
#' g <- getFleetG(syr, "LBBA"))
#'
#' eoa <- estimateL.EoA( 1, g )  # Un-informed EoA
#'
#' ieoa <- estimateL.EoA( 1, g, Lprior="normal", Lprior.mean=20, Lprior.sd=4) # Informed EoA
#'
#' # interesting plot showing movement of posterior
#' plot(ieoa$L.posterior$L, ieoa$L.posterior$pdf, type="l")
#' lines(ieoa$L.posterior$L, ieoa$L.posterior$like.pdf, col="red")
#' lines(ieoa$L.posterior$L, ieoa$L.posterior$prior.pdf, col="blue")
#' legend("topright", legend=c("prior","likelihood","posterior"), col=c("blue","red","black"), lty=1)
#'
#' @export
estimateL.EoA <- function(X,
                          beta.params,
                          Lprior="jeffreys",
                          Lprior.mean=NULL,
                          Lprior.sd=NULL,
                          conf.level=0.9){

  quants <- c((1-conf.level)/2, 0.5, 1-(1-conf.level)/2)
  zero <- 1e-6
  support.n <- 800  # MUST BE EVEN!

  ## ---- LSupport ----
  if( Lprior == "normal"){
    # Use mean and sd passed in
    support.L <- qnorm(c(zero,1-zero), Lprior.mean, Lprior.sd)
    Lmax <- support.L[2]
    Lmin <- max(zero,support.L[1])
    L.x <- seq(Lmin, Lmax, length=support.n)
    L.fx <- dnorm(L.x, Lprior.mean, Lprior.sd)
  } else if(Lprior == "gamma"){
    shape <- Lprior.mean^2 / Lprior.sd^2
    scale <- Lprior.sd^2 / Lprior.mean
    support.L <- qgamma(c(2*zero,1-2*zero), shape=shape, scale=scale)
    Lmax <- support.L[2]
    Lmin <- max(zero,support.L[1])
    L.x <- seq(Lmin, Lmax, length=support.n)
    L.fx <- dgamma(L.x, shape = shape, scale = scale )
  } else {
    #	this is the Jeffery's prior
    Lmin <- zero
    Lmax <- fmmax.ab( X, beta.params$alpha, beta.params$beta)
    L.x <- seq(Lmin, Lmax, length=support.n)
    L.fx <- beta(L.x + 0.5, 0.5)/sqrt(pi)  # not sure we need the sqrt(pi)
  }

  ## ---- MSupport ----
  Mmin <- X
  Mmax <- qpois(0.999, Lmax)
  M.x <- Mmin:Mmax

  ## ---- Prep for integration ----
  h <- L.x[2]-L.x[1]
  simp.coef <- (h/3)*c(1,rep(c(4,2),(length(L.x)/2)-1),1) # length(L.x) must be even. i.e., support.n must be even

  ## ---- Compute the posterior and check if we're close enough ----
  repeat{

    # ----- Main computations: likelihood, posterior ----
    like.X <- VGAM::dbetabinom.ab(X, size = M.x, shape1 = beta.params$alpha, shape2 = beta.params$beta)

    L.M.grid <- outer(M.x, L.x, FUN = dpois)
    L.M.grid <- L.M.grid * matrix(like.X, length(M.x), length(L.x))
    L.like <- colSums(L.M.grid) # colSums integrate out M; colSums(L.M.grid) is the likelihood; save for later
    L.post <- L.like * L.fx  # L.fx is the prior for L
    intgral <- matrix(simp.coef,1,support.n) %*% L.post
    L.post <- L.post / c(intgral)

    # ---- Check that both ends of L posterior are near zero ----

    if( L.post[which.max(L.x)] <= zero ) {
      if( L.post[which.min(L.x)] <= zero ) {
        break
      }
    }

    if( L.post[which.min(L.x)] > zero ) {
        # decrease Lmin
        Lmin.prev <- Lmin
        slp <- mean(diff(L.post[1:5]))
        if( slp < 0 ){
          Lmin <- 0.5*Lmin
        } else {
          Lmin <- (slp*Lmin - L.post[1]) / slp
        }
        cat(paste0("Pr(L=Lmin)=", round(L.post[1],6),
                   ". Expanding Lmin from ", Lmin.prev, " to ", Lmin, "\n"))
    }


    if( L.post[which.max(L.x)] > zero ) {
      # increase Lmax
      Lmax.prev <- Lmax
      slp <- mean(diff(L.post[(length(L.x)-5):length(L.x)]))
      if( slp > 0 ){
        Lmax <- 2*Lmax
      } else {
        Lmax <- (slp*Lmax - L.post[length(L.x)]) / slp
      }
      cat(paste0("Pr(L=Lmax)=", round(L.post[length(L.x)],6),
                 ". Expanding Lmax from ", Lmax.prev, " to ", Lmax, "\n"))
    }

    # ---- Recompute prior using expanded range. ----
    if( Lprior == "normal"){
      L.x <- seq(Lmin, Lmax, length=support.n)
      L.fx <- dnorm(L.x, Lprior.mean, Lprior.sd)
    } else if(Lprior == "gamma"){
      L.x <- seq(Lmin, Lmax, length=support.n)
      L.fx <- dgamma(L.x, shape = shape, scale = scale )
    } else {
      #	this is the Jeffery's prior
      L.x <- seq(Lmin, Lmax, length=support.n)
      L.fx <- beta(L.x + 0.5, 0.5)/sqrt(pi)  # not sure we need the sqrt(pi)
    }

    ## ---- New MSupport ----
    Mmin <- X
    Mmax <- qpois(0.999, Lmax)
    M.x <- Mmin:Mmax


  }

  ## ---- Summarize posterior ----
  L.cdf <- cumsum(L.post) * h

  # plots if you need to check
  #plot(L.x, L.fx)
  #plot(L.x, L.post)
  #lines(L.x, L.cdf)

  # Mean
  mu.L <- sum(L.x * L.post)*h

  # SD
  sd.L <- sqrt(sum((L.x - mu.L)^2 * L.post)*h)
  v.L <- sd.L * sd.L

  # Quantiles
  med.L <- approx(L.cdf, L.x, xout=quants, method="constant", f=1, rule=2)$y

  # Save likelihood and prior for plotting later
  intgral <- matrix(simp.coef,1,support.n) %*% L.like
  L.like <- L.like / c(intgral)
  intgral <- matrix(simp.coef,1,support.n) %*% L.fx
  L.fx <- L.fx / c(intgral)

  ans <- list(L.est=data.frame(L=med.L[2], L.mu=mu.L, L.sd=sd.L, L.lo=med.L[1],
                               L.hi=med.L[3], ci.level=conf.level),
              L.posterior = data.frame(L=L.x, pdf=L.post, cdf=L.cdf,
                                  prior.pdf=L.fx, like.pdf=L.like))

  class(ans) <- c("Lest","evoab")

  ans


}
