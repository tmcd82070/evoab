#' @export
#'
#' @name mixtureBeta
#'
#' @title Calculate parameters of a composite beta distribution arising
#' from a mixture of beta distributions.
#'
#' @description  Assuming n beta distributions, each with parameters a and b, this
#' routine computes the A and B of an overall beta distribution when the n distributions
#' are mixed.  Mixing parameters (fractions) can be specified as weights.
#'
#' @param Ba A vector containing the alpha parameter for all beta distributions.
#' @param Bb A vector containing the beta parameter for all beta distributions.
#' @param w A vector indicating the relative weight for each beta distribution
#' in the mixture distribution. The default is \code{NULL}, in which case weights are all equal.
#'
#'
#' @details
#' The assumption is that each beta distribution is independent of every
#' other, and that they are combined into a composite mixture beta distribution.
#' This routine uses mixture distribution theory to calculate a mean and a
#' variance. The mean and variance are then used in a method of moments
#' approach to calculate the parameters for the composite beta distribution.
#'
#'
#' @return A data frame (containing 1 row) with the beta mixture
#' parameter estimates and summary statistics.
#' The compoents (columns) of the returned data frame are:
#' \itemize{
#'   \item alpha = alpha parameter of the composite beta distribution,
#'   \item beta = beta parameter of the composite beta distribution,
#'   \item mean = mean of the composite beta distribution (alpha/(alpha+beta)),
#'   \item variance = variance of the composite beta (alpha*beta)/((alpha+beta)^2*(alpha+beta+1)),
#'   \item lower2.5 = 2.5-th quantile of the composite beta distribution
#'   \item upper97.5 = 97.5-th quantile of the coomposite beta distribution.
#' }
#'
#' @examples
#'
#' ## alpha parameter for the beta distribution
#' (alpha <- c(100.5,100.5,100.5,235.4,235.4))
#' ## beta parameter for the beta distribution
#' (beta <- c(234.5,234.5,234.5,2708,2708))
#' ## alpha and beta are assumed to be in the same order
#'
#' mixtureBeta(Ba=alpha,Bb=beta) ## equal weights
#'
#' weight <- 1:5
#'
#' mixtureBeta(Ba=alpha,Bb=beta,w=weight) ## unequal weights



mixtureBeta <- function(Ba,Bb,w=NULL){

    if(length(Ba)!=length(Bb)){
        stop('The length of Ba and Bb must be the same.')
    }
    if(is.null(w)){
        w <- rep(1,length(Ba))
    }

    if(length(w)!= length(Ba)){
        stop('The length of Ba, Bb and w must all be the same.')
    }

    if(any(Ba<=0)|any(Bb<=0)|any(w<=0)){
        stop('All values in Ba, Bb and w must be greater than zero.')
    }

    ## scale weighted to sum to one
    w <- w/sum(w)
    (mu <- Ba/(Ba+Bb)) ## beta means
    (sig2 <- Ba*Bb/((Ba+Bb)^2*(Ba+Bb+1))) ## beta variances

    (Eg<-sum(mu*w)) # average, expected G is the weighted average of the observed g's
    ## calculate the variance based on a mixture distribution
    (Vg <- sum(w^2*sig2))

    if(!Vg<Eg*(1-Eg)){ ## this check is need for the method of moments
        message(paste0('The mixture beta mean is ',Eg,' and the variance is ',Vg))
        stop('The variance is to big! This will produce negative values from the method of moments estimation.')
    }
    ## method of moments
    (pBa <- Eg*(Eg*(1-Eg)/Vg -1))
    (pBb <- (1-Eg)*(Eg*(1-Eg)/Vg -1))

    (gmin<-qbeta(.025,shape1=pBa,shape2=pBb))
    (gmax<-qbeta(.975,shape1=pBa,shape2=pBb))

    out <- data.frame(alpha=pBa,beta=pBb,mean=pBa/(pBa+pBb),variance=pBa*pBb/((pBa+pBb)^2*(pBa+pBb+1)),lower2.5=gmin,upper97.5=gmax)

    return(out)

} # end weightedG
