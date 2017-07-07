#' @export
#'
#' @title estimateM.EoA
#'
#' @description  Estimate single-site or multiple-class M (=mortalities)
#' parameter of Evidence of Absence (EoA)
#' using objective or informed priors.
#'
#' @param X Total number of carcasses found.
#'
#' @param beta.params A list or data frame containing at a minimum
#' components named \code{$alpha} and \code{$beta}.
#' These are the alpha and beta parameters of a Beta distribution which
#' is used for g=Pr(discovery).
#'
#' @param Mprior Character string specifying the prior distribution
#' for M.
#'  \itemize{
#'    \item "objective" uses an objective prior very close to the
#'    Jeffery's prior for a Poisson, i.e., sqrt(m+1)-sqrt(m).
#'    \item "normal" uses a truncated and descretized normal(Mprior.mean,Mprior.sd).
#'    \item "gamma" uses a descretized gamma with mean Mprior.mean and
#'    standard deviation Mprior.sd.  Note, this always assigns zero prior probability
#'    to M=0.
#'  }
#'
#' @param Mprior.mean Mean of M prior when Mprior == "normal" or "gamma".
#'
#' @param Mprior.sd Standard deviation of M prior when Mprior == "normal" or "gamma".
#'
#'
#' @param conf.level Confidence level for the confidence intervals on
#' posterior estimates of M and g.
#'
#'
#' @details This routine replicates the M estimates of the 'Single Year' and
#' 'Multiple Classes' modules in package \code{eoa}.  To repeat either case,
#' input the composite g parameter's "a" and "b" parameters here, along
#' with the number of carcasses "X", and specify the "objective" prior. See
#' Examples.
#'
#'
#' @return List containing the following components.
#' \itemize{
#'   \item \code{M.est} : A data
#'    frame containing the following:
#'    \itemize{
#'      \item \code{M} = usual point estimate of M = median of M posterior
#'      distribution.
#'      \item \code{M.mu} = mean of M posterior distribution
#'      \item \code{M.sd} = standard deviation of M posterior distribution
#'      \item \code{M.lo} = lower endpoint of a 100(conf.level)% posterior credible
#'      interval for M.
#'      \item \code{M.hi} = upper endpoint of a 100(conf.level)% posterior credible
#'      interval for M.
#'      \item \code{g} = mean of g posterior distribution
#'      \item \code{g.lo} = lower endpoint of a 100(conf.level)% posterior credible
#'      interval for g.  Note, this does not agree with the analogous number
#'      from package \code{eoa} because this comes from the posterior for g.
#'      \code{eoa} reports the quantile from the prior for g.
#'      \item \code{g.hi} = upper endpoint of a 100(conf.level)% posterior credible
#'      interval for g. Note, this does not agree with the analogous number
#'      from package \code{eoa} because this comes from the posterior for g.
#'      \code{eoa} reports the quantile from the prior for g.
#'      \item \code{ci.level} = confidence level of the credible intervals (same
#'      as input \code{conf.level}).
#'    }
#'
#'
#'   \item \code{M.margin} : the full posterior marginal distribution for M.
#'   This is a data frame with the following columns
#'   \itemize{
#'
#'     \item \code{M} : value of M in its support.
#'     \item \code{pdf} : posterior probability mass function for M.  Pr(M=m)
#'     \item \code{cdf} : posterior cummulative probability mass function for M.
#'     Pr(M<=m).
#'     \item \code{prior.pdf} : prior probability mass function for M.
#'     \item \code{like.pdf} : likelihood probability mass function for M.
#'   }
#'   Note, all three of the pdf columns sum to 1.0, even in the
#'   case of an improper prior. These columns can be plotted together
#'   using the plot method for \code{Mest} objects.
#'
#'   \item \code{g.margin} : the full posterior marginal distribution for g.
#'   This is a data frame with columns \code{$g}, \code{$pdf.g}, and \code{$cdf.g}
#'   corresponding to g, probability of g, and probability of being less than or
#'   equal to g, respectively. Note, \code{sum(result$g.margin$pdf.g)}
#'   \code{*diff(result$g.margin$g)[1] == 1.0}.
#'  }
#'
#' @author Trent McDonald
#'
#' @seealso \code{\link{estimateL.EoA}}, \code{\link{plot.Mest}}
#'
#' @examples
#' g.params <- list(alpha=600, beta=1200)
#' X <- 5
#' m.ests <- estimateM.EoA(X,g.params)
#' print(m.ests$M.est)
#'
#' m.ests <- estimateM.EoA(X, g.params, Mprior = "normal", Mprior.mean = 50, Mprior.sd = 30)
#' print(m.ests$M.est)
#'
#' m.ests <- estimateM.EoA(X, g.params, Mprior = "gamma", Mprior.mean = 50, Mprior.sd = 30)
#' print(m.ests$M.est)


estimateM.EoA <- function(X, beta.params, Mprior="objective", Mprior.mean, Mprior.sd,
                           conf.level=0.9){

  quants <- c((1-conf.level)/2, 0.5, 1-(1-conf.level)/2)

  ## ---- gSupport -----
  zero <- 1e-6  # this really effects accuracy of results, more than length of g.x or M.x
  support.g <- qbeta(c(zero,1-zero), beta.params$alpha, beta.params$beta)
  g.x <- seq(support.g[1], support.g[2], length=200)
  g.fx <- dbeta(g.x, beta.params$alpha, beta.params$beta)

  ## ---- MSupport ----
  if( Mprior == "normal"){
    # Use mean and sd passed in
    support.M <- qnorm(c(zero,1-zero), Mprior.mean, Mprior.sd)
    Mmax <- ceiling(support.M[2])
    Mmin <- max(0,floor(support.M[1]))
    M.x <- seq(Mmin, Mmax, by=1)
    M.fx <- dnorm(M.x, Mprior.mean, Mprior.sd)
  } else if(Mprior == "gamma"){
    shape <- Mprior.mean^2 / Mprior.sd^2
    scale <- Mprior.sd^2 / Mprior.mean
    support.M <- qgamma(c(zero,1-zero), shape=shape, scale=scale)
    Mmax <- ceiling(support.M[2])
    Mmin <- floor(support.M[1])
    M.x <- seq(Mmin, Mmax, by=1)
    M.fx <- dgamma(M.x, shape = shape, scale = scale )
  } else {
    #	The following line is the "objective" prior from eoa
    mu.g <- beta.params$alpha / (beta.params$alpha + beta.params$beta)
    Mmax <- ceiling(3*(max(X,1) / mu.g))
    Mmin <- 0
    M.x  <- Mmin:Mmax
    M.fx <- sqrt(M.x + 1) - sqrt(M.x)   # This matches numbers scraped from eoa. See plotEoaPriors.r
  }

  ## ---- thePrior -----
  prior <- c(outer(M.fx, g.fx, FUN="*"))

  ## ---- Loop until we get all of the distribution
  # require posterior pdf at (Mmax,gMax) to be less than "zero"
  support <- expand.grid(M=M.x, g=g.x)

  # The following sets up use of Simpson's rule to integrate and scale the joint
  # dist'n of M and g.  This is not completely necessary. Setting h = 1 (below), then
  # scaling post by post/sum(post) (i.e., forget the Simpson rule stuff)
  # will get all the estimates and quantiles
  # correct.  Only issue is that true integral of marginals will not be 1, and
  # this is important because I use a probability cutoff to decide on Mmax.

  h <- g.x[2]-g.x[1] # g values MUST be equally spaced. They are in gSupport section above
  simp.coef <- (h/3)*c(1,rep(c(4,2),(length(g.x)/2)-1),1) # length(g.x) must be even

  repeat{

    like <- dbinom(X, support$M, support$g)

    post <- prior * like
    post <- matrix(post, length(M.x), length(g.x))

    simp.coef.mat <- matrix(simp.coef, length(M.x), length(g.x), byrow=TRUE)

    intgral <- sum(simp.coef.mat * post)
    post <- post/intgral

    M.margin <- rowSums(post)
    g.margin <- colSums(post)

    # Once we go to marginals, M is discrete and g is continuous. Rescale.
    # Keep in mind, g must have an interval (h) with it.  M does not.
    M.margin <- M.margin / sum(M.margin)

    if( M.margin[length(M.x)] <= zero){
      break
    } else {
      Mmax.prev <- Mmax
      slp <- mean(diff(M.margin[(length(M.x)-5):length(M.x)]))
      if( slp > 0 ){
        Mmax <- 2*Mmax
      } else {
        Mmax <- ceiling((slp*Mmax - M.margin[length(M.x)]) / slp)
      }
      cat(paste0("Pr(M=Mmax)=", round(M.margin[length(M.x)],6),
                 ". Expanding Mmax from ", Mmax.prev, " to ", Mmax, "\n"))

      if( Mprior == "normal"){
        M.x <- seq(Mmin, Mmax, by=1)
        M.fx <- dnorm(M.x, Mprior.mean, Mprior.sd)
      } else if(Mprior == "gamma"){
        M.x <- seq(Mmin, Mmax, by=1)
        M.fx <- dgamma(M.x, shape = shape, scale = scale )
      } else {
        M.x  <- Mmin:Mmax
        M.fx <- sqrt(M.x + 1) - sqrt(M.x)   # This matches numbers scraped from eoa. See plotEoaPriors.r
      }
      # must recompute prior to accomodate expanded M
      prior <- c(outer(M.fx, g.fx, FUN="*"))
      support <- expand.grid(M=M.x, g=g.x)

    }

  }

  M.cdf <- cumsum(M.margin)
  g.cdf <- cumsum(g.margin)*h  # note inclusion of h here, and in moment computations below.

  #plot(M.x, M.fx)
  #plot(M.x, like)
  #lines(M.x, M.margin)

  # Mean
  mu.M <- sum(M.x * M.margin)
  mu.g <- sum(g.x * g.margin)*h

  # SD
  sd.M <- sqrt(sum((M.x - mu.M)^2 * M.margin))
  sd.g <- sqrt(sum((g.x - mu.g)^2 * g.margin)*h)
  v.g <- sd.g*sd.g

  pBa <- mu.g*(mu.g*(1-mu.g)/v.g -1)
  pBb <- (1-mu.g)*(mu.g*(1-mu.g)/v.g -1)
  gq <- qbeta(quants,shape1=pBa,shape2=pBb)


  # Quantiles
  med.M <- approx(M.cdf, M.x, xout=quants, method="constant", f=1, rule=2)$y

  # Save likelihood and prior for plotting later
  like <- matrix(like, length(M.x), length(g.x))
  like.M <- rowSums(like)
  like.M <- like.M / sum(like.M)
  M.fx <- M.fx / sum(M.fx)

  ans <- list(M.est=data.frame(M=med.M[2], M.mu=mu.M, M.sd=sd.M, M.lo=med.M[1],
                               M.hi=med.M[3], g=mu.g, g.lo=gq[1],
                               g.hi=gq[3], ci.level=conf.level),
              M.margin=data.frame(M=M.x, pdf=M.margin, cdf=M.cdf,
                                  prior.pdf=M.fx, like.pdf=like.M),
              g.margin=data.frame(g=g.x, pdf.g=g.margin, cdf.g=g.cdf))

  class(ans) <- c("Mest","evoab")

  ans

}

# -------------------------------------------

# Testing code
# beta.params <- data.frame(alpha = 1306.232, beta = 14095.39)
# beta.params <- data.frame(alpha = 2306.232, beta = 14095.39)
# beta.params <- data.frame(alpha = 5306.232, beta = 14095.39)
# beta.params <- data.frame(alpha = 991.5, beta = 11070)
# beta.params <- beta.params/1
#
# #print(beta.params$alpha/(beta.params$alpha+beta.params$beta))
#
#
# post <- estimateM.EoA(7, beta.params, conf.level = .95)
# print(post$M.est)
#
# beta.params <- beta.params*4
# post2 <- estimateM3.EoA(7, beta.params, conf.level = .95)
# print(post2$M.est)
#
# beta.params <- beta.params/100
# post3 <- estimateM3.EoA(1, beta.params, conf.level = .95)
# print(post3$M.est)
#
# post <- estimateM3.EoA(0, beta.params, Mprior = "normal", Mprior.mean = 50, Mprior.sd = 30, conf.level = .95)
# print(post$M.est)
#
# post <- estimateM3.EoA(0, beta.params, Mprior = "gamma", Mprior.mean = 50, Mprior.sd = 30, conf.level = .95)
# print(post$M.est)
