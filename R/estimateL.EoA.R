#' @export
#'
#' @title estimateL.EoA - Estimate single-site EoA and Informed EoA
#' lambda parameters via JAGS
#'
#' @description This routine computes lambda, the mean number of search targets per season,
#' from information about the number of found targets and the g-value (=probatility of
#' discovery).  The method is Bayesian and allows either an uniform prior for lambda
#' or an informed prior.
#' Estimation is performed using JAGS.
#'
#' @param X Total number of search targets found at all searched sites during the
#' entire search season.
#'
#' @param beta.params A list containing, at a minimum, components named $alpha and $beta.
#' These are the all-site alpha and beta parameters for g. In many cases, these parameters
#' are computed using function \code{\link{getFleetG}}.
#'
#' @param Lprior A string naming the prior distribution to use for lambda.
#' "uniform" uses a uniform prior on [0,\code{LMax}].
#' "normal" uses a normal(\code{Lprior.mean},\code{Lprior.sd}) prior for lambda.
#'
#' @param Lprior.mean Mean of lambda prior when Lprior == "normal".
#'
#' @param Lprior.sd Standard deviation of normal when Lprior == "normal".
#'
#' @param LMax Maximum lambda when Lprior = "uniform".  CAREFUL: set this to a number
#' that is well beyond the upper limit of the confidence interval on lambda.  i.e.,
#' you may need to run this a couple times to make this big enough.
#'
#' @param conf.level Confidence level for the confidence intervals on lambda.
#'
#' @param seeds A vector of length \code{nchains} containing random number
#' seeds for the MCMC sampler.  If NULL, \code{nchains} random numbers are
#' generated from R's base random number generator which is controled outside
#' this routine using \code{set.seed}.  Regardless of the random number
#' sequence set outside this routine, specifying \code{seeds} will set
#' the MCMC seeds in JAGS so that exact chains can be reproduced. Note
#' that \code{set.seed} has no effect on the random number sequences used
#' in JAGS.  The seeds, whether chosen by this routine or specified, are
#' stored in the output object.
#'
#' @param nburns The number of burn-in steps to take during MCMC sampling.
#'
#' @param niters The number of sampling steps to take during MCMC sampling.
#'
#' @param nthins The amount of MCMC chain thinning to do.  Every (\code{nthins})-th
#' sampling iteration (of which there are \code{niters}) is saved, while the
#' rest are discarded. In the end, a total of \code{niters/nthins} samples from
#' the posterior are available to the user.
#'
#' @param nchains The number of MCMC sampling chains. Specify 2 or more to
#' check convergence.
#'
#' @param nadapt The number of adapting iterations to perform before
#' burn-in.  During adaptin, JAGS is trying to optimize it's proposal distribution
#' and stepsize to increase convergence speed.
#'
#' @return List containing two components:
#' \itemize{
#'   \item  \code{$lamda.ests} is a data
#' frame containing the lambda estiamtes (point estimate and confidence interval).
#'   \item \code{$out} is the full MCMC chain object.  Use this to check convergence, etc.
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
#' eoa <- estimateL.EoA( 1, g, LMax=500 )  # Un-informed EoA
#'
#' ieoa <- estimateL.EoA( 1, g, Lprior="normal", Lprior.mean=20, Lprior.sd=4) # Informed EoA
#'
#'
#' # To check convergence, run traceplot and Gelman stats
#' plot(ieoa$out) # tracePlot
#' gelman.diag(ieoa$out) # gelmanStats
#' gelman.plot(ieoa$out) # gelmanPlot


estimateL.EoA <- function(X, beta.params, Lprior="uniform",
                          Lprior.mean=NULL, Lprior.sd=NULL,
                          LMax=62, conf.level=0.9,
                          seeds=NULL,
                          niters = 10000,
                          nthins = 10,
                          nburn  = 50000,
                          nadapt = 2000,
                          nchains = 3){

	## ---- bayesModelCode ----
	if( Lprior == "uniform" ){
		jagsModel <- "model{

		# Priors
		lambda ~ dunif(0,lambdaMax)
		g ~ dbeta(alpha, beta)

		M ~ dpois(lambda)

		# Likelihood
		for( i in 1:nx ){
		X[i] ~ dbin( g, M )
		}

		}
		"

		JAGS.data.0 <- list ( X = X,
													nx = length(X),
													alpha = beta.params$alpha,
													beta = beta.params$beta,
													lambdaMax = LMax)

} else if (Lprior == "normal"){
	jagsModel <- "model{

	# Priors
	lambda ~ dnorm(lmean,ltau)
	g ~ dbeta(alpha, beta)

	M ~ dpois(lambda)

	# Likelihood
	for( i in 1:nx ){
	X[i] ~ dbin( g, M )
	}

	}
	"

	JAGS.data.0 <- list ( X = X,
												nx = length(X),
												alpha = beta.params$alpha,
												beta = beta.params$beta,
												lmean = Lprior.mean,
												ltau = 1/Lprior.sd^2)

}

writeLines(jagsModel, "model.txt")
#cat(jagsModel)


## ---- initialValues ----

Inits <- function(x, seed){
	if( !is.null(x$lmean) ){
		list ( lambda=abs(rnorm(1,x$lmean,1/sqrt(x$ltau))), g=rbeta(1, x$alpha, x$beta),
		       .RNG.name="base::Mersenne-Twister",
		       .RNG.seed=seed
		        )
	} else {
		list ( lambda=runif(1,0,x$lambdaMax), g=rbeta(1, x$alpha, x$beta),
		       .RNG.name="base::Mersenne-Twister",
		       .RNG.seed=seed
		       )
	}
}


if( is.null(seeds) ){
  tmp.mult <- 10^(8)
  seeds <- round(runif(nchains,0,tmp.mult))
} else if (length(seeds) != nchains) {
  stop(paste("MCMC random number seeds should be NULL or length", nchains ))
}

inits<-vector("list",nchains)
for(i in 1:nchains){
  inits[[i]]<-Inits(JAGS.data.0, seeds[i])
}



# Parameters to be monitored by WinBUGS
params <- c("M", "g", "lambda")


## ---- jagsRun ----

# Initialize the chains and adapt
library(rjags)
(t1=Sys.time())
jags = jags.model(file="model.txt",
									data=JAGS.data.0,
									inits=inits,
									n.chains=nchains,
									n.adapt=nadapt)
#   Run out in the chain a ways
cat("Burning...\n")
out = update(jags, n.iter=nburn)

# Run the MCMC chains
cat("MCMC Sampling...\n")
out = coda.samples(jags,
									 variable.names=params,
									 n.iter=niters,
									 thin=nthins)

(t2=Sys.time())
print(t2-t1)

## ---- mcmcSummary ----
print(summary(out))


## ---- quantiles ----
tmp <- as.matrix(out)
l.50 <- quantile(tmp[,"lambda"], 0.5)
alpha <- c(0,1) + c(1,-1)*(1-conf.level)/2
l.CI <- quantile(tmp[,"lambda"], alpha)

## ---- STD ----
# We may never use this
l.sd <- sd(tmp[,"lambda"])



## ---- tracePlot ----

#plot(out)


## ---- gelmanStats ----

#print(gelman.diag(out))


## ---- gelmanPlot ----

#gelman.plot(out)


list(lambda.est=data.frame(lambda=l.50, lambda.sd=l.sd, lambda.lo=l.CI[1], lambda.hi=l.CI[2], ci.level=conf.level),
     out=out,
     seeds=seeds)
}
