#' @export
#'
#' @title update.eoa - Update method for EoA objects
#'
#' @description Update (run more iterations of) an EoA model object output by \code{eoa()}.
#'
#' @param obj An object of class \code{eoa}.  See function \code{\link{eoa}}.
#'
#' @param nburns Number of additional iterations to burn.
#'
#' @param niters Number of additional iterations to do.
#'
#' @param nthins Thinning interval for additional iterations
#'
#' @param add If TRUE, the additional \code{niters/nthin} kept iterations
#' will be appended to the previous iterations. If FALSE, the additional
#' \code{niters/nthin} kept iterations will overwrite the previous iterations,
#' effectively increasing the burn-in period of the MCMC model.
#'
#' @return An \code{eoa} model object.  See \code{\link{eoa}}
#'
#' @author Trent McDonald
#'
#' @seealso \code{\link{eoa}}, \code{\link{coef.eoa}}.
#'
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
#' # Could do this if eoa.1 did not converge
#' #  1) Consider everything in eoa.1 burn-in. (nadapt+nburn+niters iterations)
#' #  2) Burn another 1000 iterations, then sample 1000 taking every 10th
#' eoa.2 <- update(eoa.1, nburns=1000, niters=1000, nthin=10)
#'
#' # Could do the following if eoa.1 converged, but autocorrelation was high.
#' #  1) Consider everything in eoa.1 burn-in. (nadapt+nburn+niters iterations)
#' #  2) Do another 1000 iterations with higher thinning
#' eoa.2 <- update(eoa.1, nburns=0, niters=1000, nthin=50)
#'
#' # Could do the following if eoa.1 converged and autocorrelation were low,
#' # but simulation error was high (i.e., need more iterations)
#' #  1) Consider the n.iters in eoa.1 "good".
#' #  2) Add another 1000 iterations with same (or different) thinning
#' eoa.2 <- update(eoa.1, nburns=0, niters=1000, nthin=50, add=TRUE)
#'
update.eoa <- function(obj, niters, nburns=0, nthins=1, add=FALSE, quiet=FALSE){

  stop("This update routine does not work yet. Needs to be converted to runjags.\n")

  # Establish parameters to monitor.
  vnames <- obj$coef.labels  # coefficients
  tmp <- varnames(obj$out)
  tmp <- gsub('\\[.+\\]','',tmp)
  tmp <- unique(tmp)
  tmp <- tmp[!(tmp %in% vnames)]
  params <- c(tmp,"a")  # all coefficients in jags are "a"

  # Burn in some more if called for
  t1<-Sys.time()
  jags <- obj$jags.model

  cat("Initial iteration:")
  print(jags$iter())
  cat("\n")

  #jags$recompile()

  cat("After recompile iteration:")
  print(jags$iter())
  cat("\n")

  #writeLines(jags$model(), "model.txt")

  if(nburns > 0){
    if(!quiet) cat(paste0("Burnin...(",nburns," iterations)\n"))
    out = update(jags, n.iter=nburns, progress.bar=ifelse(quiet, "none","text"))
    cat("After burnin iteration:")
    print(jags$iter())
    cat("\n")

  }

  # Additional sampling
  if(!quiet) cat(paste0("MCMC Sampling...(",niters," iterations)\n"))
  out = coda.samples(jags,
                     variable.names=params,
                     n.iter=niters,
                     thin=nthins,
                     progress.bar=ifelse(quiet, "none","text"))

  cat("After sampling iteration:")
  print(jags$iter())
  cat("\n")

  t2=Sys.time()
  t3 <- t2-t1
  if(!quiet) cat(paste("Execution time:", round(t3,2), attr(t3,"units"), "\n\n"))


  # Add to previous, if called for
  if(add){

  }

  # extract the coefficients and label em
  # Careful.  You can't have any other parameters that start with 'a'.
  out.coefs <- out[,grep("^a",varnames(out))]
  varnames(out)[grep("^a",varnames(out))] <- vnames
  varnames(out.coefs) <- vnames

  # Check convergence
  conv <- checkIsConverged(out.coefs,1.1,quiet)

  # Check autocorrelation
  auto <- checkIsAutocorrelated(out.coefs, 0.4, 2, quiet)

  ## ---- quantiles ----
  alpha <- c(0,0.5,1) + c(1,0,-1)*(1-obj$conf.level)/2 # median plus two extreme-er quantiles
  out.sumry <- summary(out, quantiles=alpha)
  out.sd <- out.sumry$statistics[,c("SD")]
  out.sumry <- out.sumry$quantiles
  out.sd <- cbind(out.sumry[,"50%"], out.sd)
  dimnames(out.sd)[[2]] <- c("Estimate","SD")

  # Put changed things back into obj
  obj$estimates <- out.sd
  obj$intervals <- out.sumry[,-2]
  obj$out <- out
  obj$jags.model <- jags
  obj$converged <- conv$converged
  obj$Rhats <- conv$Rhats
  obj$autoCorrelated <- auto$autoCorrelated
  obj$autoCorrs <- auto$autoCorrs

  obj
}

