#' @export
#'
#' @title calcGwithK - Calculate eoa's g value
#'
#' @description Compute EoA g = Pr(discovery) given searcher efficiency,
#' carcass removal distribution, area correction, "k", and search dates.
#'
#' @param searchDates A vector of search dates.  This can either be a character
#' vector containing representations of real dates (e.g.,
#' c('2014-03-01','2014-04-01', etc.)), a vector of \code{POSIXct}
#' objects representing real dates, or a vector of integers specifying days
#' after start of the study (e.g., c(0,30,60, etc.)).  If \code{searchDates}
#' is a vector of integers, the first element must be 0.  If \code{searchDates}
#' is a vector of characters, format of the dates must be "%Y-%m-%d"
#' (or YYYY-MM-DD).
#'
#' @param crShape Shape parameter of the persistence distribution for
#' objects. If \code{crDist} == "exponential", this parameter is ignored.
#'
#' @param crScale Scale parameter of the persistence distribution for
#' objects. If \code{crDist} == "exponential", this parameter is the mean,
#' which is 1/rate, where rate is the parameter of \code{pexp} in R.
#'
#' @param crDist Character string naming the persistence distribution
#' to use. Valid values are "exponential", "weibull", "lognormal", and "loglogistic".
#' Anything else throws an error.
#'
#' @param seef A scalar specifying probability of detection. This is probability
#' of detecting an object given it is present.
#'
#' @param k A scalar specifying the proportion of objects remaining available
#' for detection
#' following the first search they became available.  E.g., \code{k} of the
#' objects missed on a search are available during the next search.
#'
#' @param a A scalar specifying the proportion of the object's spatial distribution
#' searched, sometimes called area correction proportion.
#'
#'
#' @return A scalar, the probability of object discovery, which is the product
#' of probability of an object being available and probability of detection
#' given presence.
#'
#' @details This routine assumes the arrival function of objects on
#' the study area is uniform between \code{min(searchDates)} and
#' \code{max(searchDates)}.
#'
#' @author Trent McDonald, with help from Jared Studyvin,
#' Kristen Nasmen, and Paul Rabie.  Based largely on code in package \code{eoa}
#' by Dan Dalthorp.
#'
#' @examples
#' days <- c(0,7,14,21,28,58,88)
#' calcGwithK(days, NA, 21, "exponential", 0.7, 0.2, 0.5)
#'
#' days <- c("2014-05-15","2014-06-15","2014-08-15")
#' calcGwithK(days, NA, 21, "exponential", 0.7, 0.2, 0.5)
#'
#' startDay <- c("2014-05-15")
#' endDay <- c("2014-09-15")
#' searchInterval <- 7 # days
#' nSearches <- floor(as.numeric(as.POSIXct(endDay) - as.POSIXct(startDay))/searchInterval)
#' days <- as.POSIXct(startDay) + (0:nSearches)*searchInterval*(60*60*24)
#' calcGwithK(days, 0.57792, 7.16, "weibull", 0.7, 0.2, 0.5)
#'
#' @export


calcGwithK <- function(searchDates, crShape, crScale, crDist, seef, k, a){

  ## ---- Fix up searchDates ----
  if( is.character(searchDates) ){
    searchDates <- as.POSIXct(searchDates)
    if( any(is.na(searchDates))) stop("Invalid date format in search dates vector")
  }

  if( class(searchDates)[1] == "POSIXct"){
    # This will make searchDates # days after start of study
    searchDates <- difftime( searchDates, min(searchDates), units="days" )
    searchDates <- as.numeric(searchDates)
  }

  days <- searchDates

  ## ---- Fix up CR parameters ----
  if(grepl('loglogistic',crDist,ignore.case=TRUE)){
    pda <- 1/crScale
    pdb <- crShape
  } else if(grepl('weibull',crDist,ignore.case=TRUE)){
      pda <- 1/crScale
      pdb <- crShape
  } else if(grepl('lognormal',crDist,ignore.case=TRUE)){
      pda <- crScale^2
      pdb <- log(crShape)
  } else if(grepl('exponential',crDist,ignore.case=TRUE)){
      pda <- 1/crScale
      pdb <- NA
  }

  ## ---- Some of the original functions from Dalthorp

  fn.arrivals = function( duration) {
     ans <- dunif(1:length(duration), 0, length(duration))
     ans
  }

  # These probability of persisting from day z to search on
  # day days[tt+1] = the next search, given a persistance distribution
  int_Exp_Psi_0 <-  function(z,tt, days, pda){
    return(pexp(days[tt+1]-z, pda, lower.tail=FALSE)* fn.arrivals(z))
  }
  int_Weib_Psi_0 <- function(z,tt, days, pda, pdb){
     return(pweibull(days[tt+1]-z, shape=pda, scale = pdb, lower.tail=FALSE)* fn.arrivals(z))
  }
  int_LNorm_Psi_0 <- function(z,tt, days, pda, pdb){
    return(plnorm(days[tt+1]-z, meanlog=pdb, sdlog=sqrt(pda), lower.tail=FALSE)  * fn.arrivals(z))
  }
  int_LLogis_Psi_0 <- function(z,tt, days, pda, pdb){
    return(pllogis(days[tt+1]-z,shape=pda,scale=pdb,lower.tail=FALSE) * fn.arrivals(z))
  }

  # This integrates under the persistence distribution between
  # the j-th and (j+1)-th search days, assuming the final
  # search occurs on days[tt+1], and given a distribution.
  # Hell, I know what this does but I don't know why.  Seems
  # to me that when tt=4 but we integrate from surveys 1 to 2,
  # that this somehow messes up the interpretation of k.  Why
  # do they do this search by search?  May not matter.
  fn.Psi.1 <- function(tt,j,days,crDist, pda, pdb){
    if (crDist=="exponential"){
      ans <-integrate(int_Exp_Psi_0, days[j], days[j+1],tt=tt,
        days=days,pda=pda,subdivisions=1000)$val
    }
    if (crDist=="weibull"){
      ans <-integrate( int_Weib_Psi_0,days[j],days[j+1],tt=tt,
        days=days, pda = pda, pdb = pdb,subdivisions=1000)$val
    }
    if (crDist=="loglogistic"){
      ans <- integrate( int_LLogis_Psi_0,days[j],days[j+1], tt = tt,
        days=days, pda= pda, pdb = pdb,subdivisions=1000)$val
    }
    if (crDist=="lognormal"){
      ans <- integrate( int_LNorm_Psi_0,days[j],days[j+1], tt = tt,
        days=days, pda = pda, pdb = pdb,subdivisions=1000)$val
    }
    return(ans)
  }

  ## ---- Calculate the probability of available and detected ----

  # KA: not sure why this has - 1
  nt <- length(days) - 1
  scale.init.1 <- numeric(nt)
  hold.tj <- matrix(0, nt, nt)

  # Calculate probability carcass persists until sampling given
  # sampling date
   for(tt in 1:nt){
     for(j in 1:tt){
       r <- fn.Psi.1(tt = tt,j = j,days=days, crDist, pda, pdb = pdb)
       hold.tj[tt, j] <- prod(1 - seef*k^(0:(tt - j - 1)))^(j < tt)*seef*k^(tt - j)*r
     }
     scale.init.1[tt] <- integrate(fn.arrivals, days[tt], days[tt + 1])$val
   }
   scale.init <- sum(scale.init.1)
   pr.t <- numeric(nt)

   # Trent does not understand the need for all this scaling,
   # ...in both directions of hold.tj
   for(tt in 1:nt){
      pr.t[tt] <- sum(hold.tj[tt,] / scale.init.1[tt])
   }

   prob_obs = a*sum(pr.t * scale.init.1 / scale.init)

   return(prob_obs)
}

