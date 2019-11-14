#' @export
#'
#' @title Compute average carcass persistence
#'
#' @description This computes average carcass persistence
#' during a set of searches, given a cummulative hazard distribution
#'
#' @param searches A vector containing julian dates of all
#' carcass searches. Only the number of searches and difference
#' between min and max matter.  Average search interval is
#' computed and average persistence probability is computed
#' between 0 and average interval.  Hence, placement of searches
#' within \code{min(searches)} and \code{max(searches)} does
#' not matter.
#'
#' @param dist Character string naming the survival distribution
#' to use. Valid values are "exponential" and "weibull".
#'
#' @param shape Shape parameter of the survival distribution.
#' If \code{dist} == "exponential", this parameter is ignored.
#'
#' @param scale Scale parameter of the survival distribution.
#' If \code{dist} == "exponential", this parameter is the mean,
#' which is 1/rate, where rate is the parameter of \code{pexp} in R.
#'
#' @return A single number.  The average probability of a carcass
#' surviving from 0 to the average search interval. i.e.,
#' average survival under the survival curve from 0 to
#' \code{diff(range(searches))/2}.
#'
#' @author Trent McDonald
#'
#' @seealso \code{\link{eoa}}
#'
#' @examples
#' avgPersist(c(0,15,30),"exponential",scale=8.6251)
#' avgPersist(c(0,15,30),"weibull",0.57792,7.1683)
#'
#' @export
#'
avgPersist <- function(searches, dist="exponential", shape, scale){

 Ir <- diff(range(searches))/2

 tmp <- seq(0,Ir,length=100)
 tmp2 <- tmp[-length(tmp)] + diff(tmp)[1]/2

 if( dist=="exponential"){
    ans <- mean(1-pexp(tmp2,1/scale))
 } else if( dist == "weibull"){
    ans <- mean(1-pweibull(tmp2,shape,scale))
 } else {
    stop("Unknown carcass persistence distribution.")
 }

 ans


}
