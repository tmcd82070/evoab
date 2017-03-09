#' @export
#'
#' @title checkIsAutocorrelated - Check whether MCMC chains are autocorrelated.
#'
#' @description Run autocorrelation on every parameter in an mcmc.list and make
#' a decision whether they are significantly auto correlated.  Ideally, chains
#' are not autocorrelated when converged.
#'
#' @param obj An mcmc.list object from the \code{rjags} package.  These are output
#' by the \code{coda.samples} function.
#'
#' @param criterion scalar specifying the level of acceptable autocorrelation at
#' the assessed \code{lag}.  If all average lag2 autocorrelations are less than this number, the
#' whole mcmc object is deems 'random' or not autocorrelated.
#'
#' @param lag The lag at which to assess autocorrelation.  This should generally be a
#' low number like 1, 2, 5, or 10.
#'
#' @param quiet Whether to print any output on screen.
#'
#' @details This routine calculates auto correlation separately at the specified lag
#' for every chain in the mcmc.list, then averages autocorrelations across
#' chains.
#'
#' @return A list with the following components
#' \list{
#'   \item autoCorrelated: If the average autocorrelations at the specified
#'   lag are all less than
#'   \code{criterion}, the chain is declared un-autocorrelated and is set FALSE.
#'   Otherwise, this component of the returns is TRUE.
#'   \item autoCorrs: the computed autocorrelations at the specified lag.
#'   This value is the average across chains in \code{obj}.
#'  }
#'
#' @author Trent McDonald
#'
#'
#' @seealso \link{\code{checkIsConverged}}
#'
#' @examples
#' \dontrun{
#' jags = jags.model(file="someJAGSFile.txt",
#'    data=someJAGS.data, inits=someJAGS.inits,
#'    n.chains=3,n.adapt=100)
#' update(jags, n.iter=1000)
#' out = coda.samples(jags,
#'    variable.names=c("someParameter"),
#'    n.iter=1000,
#'    thin=2)
#'
#' auto <- chechIsAutocorrelated(out)
#' }


checkIsAutocorrelated <- function(obj, criterion=0.4, lag=1, quiet=FALSE){

  if(!quiet) cat(paste0("\n", paste(rep("-",4),collapse="")," Checking autocorrelation:\n"))
  vn <- varnames(obj)
  ncovars <- nvar(obj)
  nchains <- nchain(obj)
  acflag2 <- rep(NA,ncovars)
  names(acflag2) <- vn
  for( i in 1:ncovars ){
    tmp <- rep(NA,nchains)
    for(j in 1:nchains){
      tmp.2 <- acf(as.matrix(obj[[j]][,i]), lag.max = lag, plot=FALSE)
      tmp[j] <- tmp.2$acf[tmp.2$lag ==lag,,]
    }
    acflag2[i] <- mean(tmp)
    if(!quiet) cat(paste( vn[i], ": mean(acf[",lag, "]) =", round(acflag2[i],5), "\n") )
  }
  if(any(acflag2>criterion, na.rm=TRUE)){
    if(!quiet) cat(paste0("\nSome autocorrelations exceed ",criterion, ".\nInspect acfplot(mcmc.list) and increase thinning.\n"))
    ans <- TRUE
  } else {
    ans <- FALSE
  }

  list(autoCorrelated = ans, autoCorrs=acflag2)

}
