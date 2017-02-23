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
#' @param X Vector of number of carcasses found, one element per year. If
#' multiple sites are involved, elements of X are the total (summed) number
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
#' X <- c( 0, 1, 3)
#'
#'
#' eoa <- estimateL.EoA( X, g, LMax=500 )  # Un-informed EoA
#'
#' ieoa <- estimateL.EoA( X, g, Lprior="normal", Lprior.mean=20, Lprior.sd=4) # Informed EoA
#'
#'
#' # To check convergence, run traceplot and Gelman stats
#' plot(ieoa$out) # tracePlot
#' gelman.diag(ieoa$out) # gelmanStats
#' gelman.plot(ieoa$out) # gelmanPlot

#'
estimateL.EoA.MultiYear <- function(X, beta.params, Lprior="uniform", Lprior.mean=NULL, Lprior.sd=NULL, LMax=62, conf.level=0.9 ){
	library(rjags)

	nyrs <- length(X)
	if( length(beta.params$alpha) == 1 ){
		alpha.vec <- rep(beta.params$alpha, nyrs)
		beta.vec <- rep(beta.params$beta, nyrs)
	} else {
		alpha.vec <- beta.params$alpha
		beta.vec <- beta.params$beta
	}

	if( length(alpha.vec) != nyrs ) stop("Length of alpha and beta vectors must be 1 or equal length of X")
	if( length(alpha.vec) != length(beta.vec) ) stop("Lengths of alpha and beta inputs must be equal") # Actually, this can't happen if beta.params is a data frame. Oh well, leave it.

	## ---- bayesModelCode ----
	if( Lprior == "uniform" ){
		jagsModel <- "model{

		# Priors
		lambda ~ dunif(0.1,lambdaMax)

		# Likelihood
		for( i in 1:nx ){
			g[i] ~ dbeta(alpha[i], beta[i])
			M[i] ~ dpois(lambda)
			X[i] ~ dbin( g[i], M[i] )
		}

		Mtot <- sum(M[])

		}
		"

		JAGS.data.0 <- list ( X = X,
													nx = nyrs,
													alpha = alpha.vec,
													beta = beta.vec,
													lambdaMax = LMax)

} else if (Lprior == "normal"){
	jagsModel <- "model{

	# Priors
	lambda ~ dnorm(lmean,ltau)

	# Likelihood
	for( i in 1:nx ){
		g[i] ~ dbeta(alpha[i], beta[i])
		M[i] ~ dpois(lambda)
		X[i] ~ dbin( g[i], M[i] )
	}

	Mtot <- sum(M[])

	}
	"

	JAGS.data.0 <- list ( X = X,
												nx = nyrs,
												alpha = alpha.vec,
												beta = beta.vec,
												lmean = Lprior.mean,
												ltau = 1/Lprior.sd^2)
	#print(JAGS.data.0)

}

writeLines(jagsModel, "model.txt")
#cat(jagsModel)


## ---- initialValues ----

Inits <- function(x){
	if( !is.null(x$lmean) ){
		gg <-rbeta(x$nx, x$alpha, x$beta)
		M <- ceiling(x$X / gg) + 1
		l <- mean(M)
		list ( lambda = l,
					 M = M,
					 g=gg,
					 .RNG.name="base::Wichmann-Hill",
					 .RNG.seed=93229 )
	} else {
		gg <-rbeta(x$nx, x$alpha, x$beta)
		M <- ceiling(x$X / gg) + 1
		l <- mean(M)
		list ( lambda = l,
					 M = M,
					 g=gg,
					 .RNG.name="base::Wichmann-Hill",
					 .RNG.seed=62234 )
	}
}

inits = list(  Inits(JAGS.data.0),
							 Inits(JAGS.data.0),
							 Inits(JAGS.data.0))


# MCMC sample size settings:
niters <- 10000
nthins <- 10
nburn  <- 75000
nchains <- 3
nadapt <- 3000


# Parameters to be monitored by WinBUGS
params <- c("M", "g", "lambda", "Mtot")


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

M.50 <- quantile(tmp[,"Mtot"], 0.5)
M.CI <- quantile(tmp[,"Mtot"], alpha)

## ---- STD ----
# We may never use this
l.sd <- sd(tmp[,"lambda"])
M.sd <- sd(tmp[,"Mtot"])



list(lambda.est=data.frame(lambda=l.50, lambda.sd=l.sd, lambda.lo=l.CI[1], lambda.hi=l.CI[2],
													 Mtot=M.50, Mtot.sd=M.sd, Mtot.lo=M.CI[1], Mtot.hi=M.CI[2],
													 ci.level=conf.level), out=out)
}
