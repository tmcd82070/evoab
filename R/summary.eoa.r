#' @export
#'
#' @title summary.eoa - summary method for EoA objects
#'
#' @description Summarizes an EoA model object output by \code{eoa()}.
#'
#' @param obj An object of class \code{eoa}.  See function \code{\link{eoa}}.
#'
#' @return NULL.  This function prints a summary of the model on screen.
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
#' summary(eoa.1)
#'
summary.eoa <- function(obj){

  cat("Call:\n")
  print(obj$call)
  cat("\n")

  cat("Coefficients:\n")
  cis <- obj$intervals[labels(obj),,drop=FALSE]
  dimnames(cis)[[2]] <- paste0("CI.",dimnames(cis)[[2]])
  ests <- cbind(obj$estimates[labels(obj),,drop=FALSE], cis)
  print(ests)
  cat("\n")

  cat("Total estimated targets:\n")
  ests <- c(obj$estimates[labels(obj,"Mtot"),], obj$intervals[labels(obj,"Mtot"),])
  ests <- matrix(ests,1)
  dimnames(ests) <- list("Mtot", c("Estimate","SD",paste0("CI.",dimnames(obj$intervals)[[2]])))
  print(ests)
  cat("\n")


  cat("MCMC Information:\n")
  cat(paste(nchain(obj$out), "chains,", start(obj$out)-thin(obj$out), "burn-in steps,",
            end(obj$out)-start(obj$out), "iterations thinned by", thin(obj$out), "to", niter(obj$out), "total.\n"))

  # the following lines will print Rhats and autocorrelations
  #tmp <- cbind(obj$Rhats, obj$autoCorrs)
  #dimnames(tmp)[[2]] <- c("R-hat","Corr.Lag(2)")
  #print(tmp)

  cat("The chains likely", ifelse(obj$converged, "converged.", "DID NOT CONVERGE (try additional burn-in)."), "\n")
  if(obj$autoCorrelated){
    ncovars <- length(labels(obj))
    nchains <- nchain(obj$out)
    acf0mat <- matrix(NA, ncovars, nchains)
    for( i in 1:ncovars ){
      for(j in 1:nchains){
        tmp.2 <- acf(as.matrix(obj$out[[j]][,i]), lag.max = niter(obj$out)/2, plot=FALSE)
        tmp.2 <- data.frame(acf=tmp.2$acf,lag=tmp.2$lag)
        if( sum(tmp.2$acf>0) >0){
          tmp.2 <- tmp.2[tmp.2$acf>0,]
          acf0mat[i,j] <- min(which(diff(tmp.2$lag)>1))  # first lag acf below 0
        }
      }
    }
    cat("Within chain auto-correlation is HIGH.")
    cat(paste0(" (try nthins = ", round(mean(acf0mat,na.rm = TRUE)*thin(obj$out)), ")\n"))
  } else {
    cat("Within chain auto-correlation is acceptable.\n")
  }

}
