##function to calculate g 

##gdat is a list of data including (according to the EoA V2 user manual p.101)
## a:  spatial coverage (as a proportion)
## v:  temporal coverage (probably as a proportion)--ignored if there is not an 
##      object called arrfun available to the function
## p:  SEEF as a probability
## k:  k (probably as a proportion)
## samtype: type of search schedule:  "Formula" or "Custom"
## Isam, nsearch:  search interval and number of searches (if samtype == 'Formula')
## days: vector of search dates, starting with 0 (if samtype == 'Custom'
##persistence_distn: name of persistence distribution: "Exponential", "Weibull", "Log-Logistic", or "Lognnormal"
##pda, pdb: alpha and beta parameters of persistence distribution (see appendix H of user guide for non-standard parameterizations!)

# calcg.fixed(gdat = list(a = .1, v = .01, p = 0.5, k = 0, samtype = 'Formula', 
    # Isam = 1, nsearch = 365, persistence_distn = 'Exponential', 
    # pda = NULL, pdb = 100))

calcg.fixed <- function (gdat, arrdat = NULL) 
{
  clipProb <- 0.001
  for (nm in names(gdat)) assign(nm, gdat[[nm]])
  if (samtype == "Formula") {
    days <- c(0, 1:nsearch) * Isam
  }
  else if (samtype == "Custom") {
    Isam <- round(max(days)/(length(days) - 1), 1)
    nsearch <- length(days) - 1
  }
  else {
    warning("error in search schedule. aborting calculation.")
    return(F)
  }
  ind1 <- rep(1:nsearch, times = nsearch:1)
  ind2 <- ind1 + 1
  ind3 <- unlist(lapply(1:nsearch, function(x) x:nsearch)) + 
    1
  schedule.index <- cbind(ind1, ind2, ind3)
  schedule <- cbind(days[ind1], days[ind2], days[ind3])
  nmiss <- schedule.index[, 3] - schedule.index[, 2]
  maxmiss <- max(nmiss)
  pobs <- numeric(dim(schedule)[1])
  powk <- cumprod(c(1, rep(k, maxmiss)))
  notfind <- cumprod(1 - p * powk[-length(powk)])
  nvec <- c(1, notfind) * p
  pfind.si <- nvec * powk
  intxsearch <- unique(cbind(schedule[, 2] - schedule[, 1], 
                             schedule[, 3] - schedule[, 2]), MAR = 1)
  if (persistence_distn == "Exponential") 
    pda <- 1/pdb
  ppersu <- eoa::ppersist(persistence_distn, t_arrive0 = 0, 
                          t_arrive1 = intxsearch[, 1], t_search = intxsearch[, 
                                                                             1] + intxsearch[, 2], pda = pda, pdb = pdb)
  if (is.null(arrdat)) {
    arrvec <- (schedule[, 2] - schedule[, 1])/max(days)
    arrmiss0 <- 0
    arrmissf <- 0
    v <- 1
    arrSimplify <- T
  }
  else if (is.numeric(arrdat)) {
    if (length(arrfun) != nsearch) {
      warning(paste0("Vector of arrival probabilities is ", 
                     ifelse(length(arrfun) > nsearch, "longer ", "shorter "), 
                     "than number of searches. Aborting calculation."))
      return(F)
    }
    if (sum(is.na(arrfun)) > 0) {
      warning("Missing values not allowed in vector of arrival probabilities (arrfun).")
      return(F)
    }
    arrvec <- rep(arrfun, nsearch:1)
    arrmiss0 <- 0
    arrmissf <- 0
    v <- 1
    arrSimplify <- T
  }
  else if (is.list(arrdat) || inherits(arrdat, "R6")) {
    if (is.function(arrdat$arrfun)) {
      if (arrdat$duration < arrdat$s0 + max(days)) {
        warning("Error in arrdat: duration must be at least s0 + max(days). Period of inference shorter than monitored period. Aborting calculation.")
        return(F)
      }
      arrSimplify <- ifelse(is.null(arrdat$arrSimplify), 
                            F, arrdat$arrSimplify)
      deno <- integrate(arrfun, lower = 0, upper = arrdat$duration)$val
      arrmiss0 <- integrate(arrfun, lower = 0, upper = arrdat$s0)$val/deno
      arrmissf <- integrate(arrfun, lower = arrdat$s0 + 
                              max(days), upper = arrdat$duration)$val/deno
      v <- 1 - arrmiss0 - arrmissf
      if (arrSimplify) {
        arrvec0 <- numeric(nsearch)
        for (i in 1:nsearch) {
          arrvec0[i] <- integrate(arrfun, lower = days[i], 
                                  upper = days[i + 1])$val
        }
        arrvec0 <- arrvec0/sum(arrvec0)
        arrvec <- rep(arrvec0, nsearch:1)
      }
    }
    else if (is.numeric(arrdat$arrfun) || arrdat$arrfun == 
             "Uniform") {
      arrmiss0 <- arrdat$arrmiss0
      arrmissf <- arrdat$arrmissf
      v <- 1 - arrmiss0 - arrmissf
      arrSimplify <- T
      if (arrdat$arrfun == "Uniform") {
        arrvec <- (schedule[, 2] - schedule[, 1])/max(days)
      }
      else {
        arrvec <- rep(arrfun, nsearch:1)
      }
    }
  }
  else {
    warning("error in arrdat specification. Aborting calculation")
    return(F)
  }
  if (arrSimplify) {
    for (i in 1:length(pobs)) {
      pobs[i] <- pfind.si[nmiss[i] + 1] * ppersu[which(abs(intxsearch[, 
                                                                      1] - (schedule[i, 2] - schedule[i, 1])) < clipProb & 
                                                         abs(intxsearch[, 2] - (schedule[i, 3] - schedule[i, 
                                                                                                          2])) < clipProb), ] * arrvec[i]
    }
    pobs <- sum(pobs)
  }
  else {
    arrdeno <- integrate(arrfun, lower = s0, upper = sf)$val
    ipart <- numeric(dim(schedule)[1])
    for (i in 1:dim(schedule)[1]) {
      ipart[i] <- switch(persistence_distn, Exponential = integrate(function(tt) (1 - 
                                                                                    pexp(schedule[i, 3] - tt, rate = 1/pdb)) * arrfun(tt + 
                                                                                                                                        s0)/arrdeno, lower = schedule[i, 1], upper = schedule[i, 
                                                                                                                                                                                              2])$val, Weibull = integrate(function(tt) (1 - 
                                                                                                                                                                                                                                           pweibull(schedule[i, 3] - tt, shape = pda, scale = pdb)) * 
                                                                                                                                                                                                                             arrfun(tt + s0)/arrdeno, lower = schedule[i, 
                                                                                                                                                                                                                                                                       1], upper = schedule[i, 2])$val, Lognormal = integrate(function(tt) (1 - 
                                                                                                                                                                                                                                                                                                                                              plnorm(schedule[i, 3] - tt, meanlog = pdb, sdlog = sqrt(pda))) * 
                                                                                                                                                                                                                                                                                                                                arrfun(tt + s0)/arrdeno, lower = schedule[i, 
                                                                                                                                                                                                                                                                                                                                                                          1], upper = schedule[i, 2])$val, `Log-Logistic` = integrate(function(tt) (1 - 
                                                                                                                                                                                                                                                                                                                                                                                                                                                      actuar::pllogis(schedule[i, 3] - tt, shape = pda, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      scale = pdb)) * arrfun(tt + s0)/arrdeno, lower = schedule[i, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                1], upper = schedule[i, 2])$val, NA)
    }
    pfind <- pfind.si[nmiss + 1]
    pobs <- sum(pfind * ipart)
  }
  return(list(`Full site, full year` = pobs * a * v, `Full site, monitored period` = pobs * 
                a, `Searched area, monitored period` = pobs))
}