#' @export
#'
#' @title estimateL.EoA.MultiYear - Estimate multi-season
#' EoA and informed EoA Lambda parameter via JAGS
#'
#' @description  This routine computes lambda, the mean number of search targets per season,
#' from information about the number of targets found and the g-values (=probatility of
#' discovery) at all
#' searched sites.  The method is Bayesian and allows either an uniform prior for lambda
#' or an informed prior.
#' Estimation is performed using JAGS.
#'
#' @param Y Vector of number of carcasses found, one element per year. If
#' multiple sites are involved, elements of Y are the total (summed) number
#' of targets found per season.
#'
#' @param beta.params A data frame containing, at a minimum, two columns named \code{$alpha}
#' and \code{$beta}.
#' These are the annual alpha and beta parameters that determine the overall annual g.
#' Length of \code{$alpha} and \code{$beta} is either one, in which case g is constant across years,
#' or one element per year.
#'
#' @param Lprior Prior for lambda.  "uniform" uses a uniform[0,LMax] prior for lambda.
#' "normal" uses a normal(Lprior.mean,Lprior.sd) prior for lambda.
#'
#' @param Lprior.mean Mean of lambda prior when Lprior == "normal".
#'
#' @param Lprior.sd Standard deviation of normal when Lprior == "normal".
#'
#' @param LMax Maximum lambda when Lprior = "uniform".
#'
#' @param conf.level Confidence level for the confidence intervals on lambda.
#'
#' @return List containing two components.  \code{$lamda.ests} is a data
#' frame containing the lambda estiamtes (point est and confidence interval).
#' In this component \code{$Mtot} is the estimated number of targets over
#' all seasons.
#' \code{$out} is the full MCMC chain object.  Use this to check convergence, etc.
#'
#'
#' @examples
#' # Three year study
#' g <- data.frame(
#'    alpha = c( 69.9299, 63.5035,  84.6997),
#'    beta = c(  736.4795,  318.3179, 759.9333 )
#'    )
#' Y <- c( 0, 1, 3)
#'
#'
#' eoa <- estimateL.EoA.MultiYear( Y, g, LMax=500 )  # Un-informed EoA
#'
#' ieoa <- estimateL.EoA.MultiYear( Y, g, Lprior="normal", Lprior.mean=20, Lprior.sd=4) # Informed EoA
#'
#'
#' # Check convergence and autocorrelation
#'
#' gelman.diag(ieoa$out) # gelmanStats
#' gelman.plot(ieoa$out) # gelmanPlot
#'
#' # Traceplots
#' library(lattice)
#' plot(ieoa$out) # tracePlot, all parameters
#' xyplot(ieoa$out[,grep("^a\\[",varnames(ieoa$out))]) # nicer trace of coefficients
#'
#' # Autocorrelation functions
#' acfplot(ieoa$out[,grep("^a\\[",varnames(ieoa$out))], ylim=c(-.2,1), lag.max=300)
#' acfplot(ieoa$out[,grep("^a\\[",varnames(ieoa$out))], ylim=c(-.2,1), thin=10)
#'
#' # Density plots
#' densityplot(ieoa$out[,grep("^a\\[",varnames(ieoa$out))])
#'
#' # Correlation among coefficients
#' levelplot(ieoa$out[,grep("^a\\[",varnames(ieoa$out))][[1]])
#'
#' # QQ plots
#' qqmath(ieoa$out[,grep("^a\\[",varnames(ieoa$out))])
#'
estimateL.EoA.MultiYear <- function(lambda, beta.params, data, priors=NULL,
                                    conf.level=0.9, nburns = 500000, niters = 20000,
                                    nthins = 10, nchains = 3, nadapt = 3000,
                                    quiet=FALSE, seeds=NULL ){
	library(rjags)

  ## ---- lambdaModel ----
  # Resolve formula for lambda
  if (missing(data))
    data <- environment(lambda)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("lambda", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  names(mf)[names(mf)=="lambda"] <- "formula"
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf,"any")
  lambda.covars <- if (!is.empty.model(mt)){
    model.matrix(mt, mf, contrasts)
  }
  ncovars <- ncol(lambda.covars)
  vnames<-dimnames(lambda.covars)[[2]]


  ## ---- initialize ----
  # Make sure one beta dist'n parameter per row
	nyrs <- length(Y)
	if( length(beta.params$alpha) == 1 ){
		alpha.vec <- rep(beta.params$alpha, nyrs)
		beta.vec <- rep(beta.params$beta, nyrs)
	} else {
		alpha.vec <- beta.params$alpha
		beta.vec <- beta.params$beta
	}

	if( length(alpha.vec) != nyrs ) stop("Length of alpha and beta vectors must be 1 or equal length of Y")
	if( length(alpha.vec) != length(beta.vec) ) stop("Lengths of alpha and beta inputs must be equal") # Actually, this can't happen if beta.params is a data frame. Oh well, leave it.


  ## ---- resolvePriors ----
  # Use vague normal priors for coefficients (huge SE's) by default
  sd.n.start <- compVagueSd(Y,alpha.vec,beta.vec,lambda.covars, range.multiplier = 100)
  coefTaus <- sd.n.start$vagueSd
  coefMus <- rep(0,ncovars)
  names(coefMus)<-vnames
  if(is.vector(priors)){
    # Use informed prior for intercept parameter
    if(!all(c("mean","sd") %in% names(priors))){
      stop("Priors must contain names 'mean' and 'sd'.")
    } else {
      coefMus["(Intercept)"] <- priors["mean"]
      coefTaus["(Intercept)"] <- priors["sd"]
    }
  } else if(!is.null(priors)){
    # Use informed priors for all coefficience, priors must be a data.frame
    if( !is.data.frame(priors)){
      stop("Priors must be NULL, a vector, or a data.frame.")
    } else if(!all(c("mean","sd") %in% names(priors))){
      stop("Priors data.frame must contain 'mean' and 'sd'.")
    } else if(!any(row.names(priors) %in% vnames)){
      warning(paste("No row.names(priors) matched parameter names. Vague priors used."))
    } else {
      fnd <- names(coefMus) %in% rownames(prior)
      coefMus[fnd] <- priors[names(coefMus)[fnd], "mean"]
      coefTaus[fnd] <- priors[names(coefTaus)[fnd], "sd"]
    }
  }

  # Recall: tau of dnorm in JAGS is 1/variance
  #print(coefTaus)
  coefTaus <- 1/coefTaus^2


	## ---- bayesModelCode ----
	jagsModel <- "model{

		# Priors
    for(i in 1:ncovars){
      a[i] ~ dnorm( coefMus[i], coefTaus[i] )
    }

    # functional relations
    for(i in 1:nx){
      for(j in 1:ncovars){
        logl[i,j] <- a[j]*lambda.covars[i,j]
      }
      lambda[i] <- exp(sum(logl[i,]))
    }

		# Likelihood
		for( i in 1:nx ){
			g[i] ~ dbeta(alpha[i], beta[i])
			M[i] ~ dpois(lambda[i])
			Y[i] ~ dbin( g[i], M[i] )
		}

		Mtot <- sum(M[])

		}
	"

	JAGS.data.0 <- list ( Y = Y,
													nx = nyrs,
													ncovars = ncovars,
													coefTaus = coefTaus,
													coefMus = coefMus,
													alpha = alpha.vec,
													beta = beta.vec,
													lambda.covars = lambda.covars)


  writeLines(jagsModel, "model.txt")
  #cat(jagsModel)


  ## ---- initialValues ----
  # MCMC sample size settings:


  if( is.null(seeds) ){
    tmp.mult <- 10^(8)
    seeds <- round(runif(nchains,0,tmp.mult))
  } else if (length(seeds) != nchains) {
    stop(paste("MCMC random number seeds should be NULL or length", nchains ))
  }

  Inits <- function(x,strt,seed){
  		gg <-rbeta(x$nx, x$alpha, x$beta)
  		M <- ceiling(x$Y / gg) + 1
  		#a <- strt$startA[,"Estimate"] + rnorm(nrow(strt$startA),0,strt$startA[,"Std. Error"])
  		a <- rep(0,x$ncovars)
  		a[1] <- log(mean(M))
  		list ( a = a,
  					 M = M,
  					 g=gg,
  					 .RNG.name="base::Mersenne-Twister",
  					 .RNG.seed=seed )
  }

  inits<-vector("list",nchains)
  for(i in 1:nchains){
    inits[[i]]<-Inits(JAGS.data.0, sd.n.start, seeds[i])
  }



  # Parameters to be monitored by WinBUGS
  params <- c("a", "M", "lambda", "Mtot")


  ## ---- jagsRun ----

  # Initialize the chains and adapt
  library(rjags)
  (t1=Sys.time())
  jags = jags.model(file="model.txt",
  									data=JAGS.data.0,
  									inits=inits,
  									n.chains=nchains,
  									n.adapt=nadapt,
  									quiet = quiet)

  #   Run out in the chain a ways
  if(!quiet) cat("Burnin...\n")
  out = update(jags, n.iter=nburns, progress.bar=ifelse(quiet, "none","text"))

  # Run the MCMC chains
  if(!quiet) cat("MCMC Sampling...\n")
  out = coda.samples(jags,
  									 variable.names=params,
  									 n.iter=niters,
  									 thin=nthins,
  									 progress.bar=ifelse(quiet, "none","text"))

  (t2=Sys.time())
  t3 <- t2-t1
  if(!quiet) cat(paste("Execution time:", round(t3,2), attr(t3,"units"), "\n\n"))

  ## ---- mcmcChecking ----
  #print(summary(out))

  # extract the coefficients and label em
  out.coefs <- out[,grep("^a\\[",varnames(out))]
  varnames(out)[grep("^a\\[",varnames(out))] <- vnames
  varnames(out.coefs) <- vnames

  # Check convergence
  conv <- checkIsConverged(out.coefs,1.1,quiet)

  # Check autocorrelation
  auto <- checkIsAutocorrelated(out.coefs, 0.4, 2, quiet)



  ## ---- quantiles ----
  alpha <- c(0,0.5,1) + c(1,0,-1)*(1-conf.level)/2
  out.sumry <- summary(out, quantiles=alpha)
  out.sd <- out.sumry$statistics[,c("SD")]
  out.sumry <- out.sumry$quantiles
  out.sd <- cbind(out.sumry[,"50%"], out.sd)
  dimnames(out.sd)[[2]] <- c("Estimate","SD")
  out.sumry <- list(estimates=out.sd,
                    intervals=out.sumry[,-2]
                    )

  ## ---- Done ----
  priors.df <- data.frame(mean=coefMus, sd=1/sqrt(coefTaus))
  c(out.sumry,
    list(
    out=out,
    priors = priors.df,
    design.mat=lambda.covars),
    conv,
    auto,
    conf.level=conf.level,
    seeds=seeds
  )


}


# ==================================================================================================

# ---- Testing

# A simple example with three years, four observations
# g <- data.frame(
#   alpha = c( 69.9299, 63.5035,  84.6997, 84.6997),
#   beta = c(  736.4795,  318.3179, 759.9333, 759.9333 )
#   )
# Y <- c( 0, 1, 3, 10 )
#
# tmp.df <- data.frame(year=factor(c("2015","2016","2017","2017")),
#                      Year=c(1,2,3,3))


# A more complicated example, with 3 years, 100 observations, lambda slope = 20
set.seed(9430834)
n <- 21
n3 <- round(n/3)

g <- data.frame(
  alpha = rnorm(n,70,2),
  beta = rnorm(n,700,25)
)
Y <- rbinom(n, c(rep(20,n3), rep(40,n3), rep(60,n-2*n3)), g$alpha/(g$alpha+g$beta))

tmp.df <- data.frame(year=factor(c(rep("2015",n3),rep("2016",n3),rep("2017",n-2*n3))),
                     Year=c(rep(1,n3),rep(2,n3),rep(3,n-2*n3)))


prior <- data.frame(mean=c(log(3),10,5,5), sd=c(.5,6,7,8), x=10:13)
row.names(prior) <- c("(Intercept)","year1","year2", "bob")

eoa <- estimateL.EoA.MultiYear(Y~year, g, data=tmp.df, nburn = 700000, niters= 2000*100, nthins = 100, quiet=FALSE )  # Un-informed EoA
#tmp.1 <-apply(as.matrix(eoa$out),2,median)

#print(tmp.1)

#eoa <- estimateL.EoA.MultiYear(Y~year, g, data=tmp.df, quiet=FALSE, priors=NULL )  # Un-informed EoA
