#' @export
#'
#' @title estimateM.EoA
#'
#' @description  Estimate single-site or multiple-class M (=mortalities)
#' parameter of Evidence of Absence (EoA)
#' using either the objective or an
#' informed prior.
#'
#' @param X Total number of carcasses found.
#'
#' @param beta.params A list containing at a minimum components named $alpha and $beta.
#' These are the alpha and beta parameters of a Beta distribution which
#' is used for g=Pr(discovery).
#'
#' @param Mprior Character string specifying the prior distribution
#' for M.  "objective"
#' uses the objective prior of Dalthorp
#' i.e., sqrt(m+1)-sqrt(m).  "normal" uses
#' a normal(Mprior.mean,Mprior.sd) prior for M, discretized to
#' the integers 0:Mmax.
#'
#' @param Mprior.mean Mean of M prior when Mprior == "normal".
#'
#' @param Mprior.sd Standard deviation of normal when Mprior == "normal".
#'
#' @param Mmax Maximum M.  Set this to be larger than 99-th
#' percentile of posterior M distribution. Note: run time depends strongly
#' on this parameter.  Larger Mmax = longer run-time.
#'
#' @param conf.level Confidence level for the confidence intervals on M.
#'
#' @param nburns The number of burn-in steps per chain to take during MCMC sampling.
#'
#' @param niters The number of sampling steps per chain to take during MCMC sampling.
#'
#' @param nthins The amount of MCMC chain thinning to do.  Every (\code{nthins})-th
#' sampling iteration in each chain is saved, while the
#' rest are discarded. Each chain has \code{niters} iterations. In
#' the end, a total of \code{nchains*niters/nthins} samples from
#' the posterior are available to the user.
#'
#' @param nchains The number of MCMC sampling chains. Must specify 2 or more to
#' check convergence.
#'
#' @param nadapts The number of adapting iterations to perform before
#' burn-in.  During adapting, JAGS tries to optimize it's proposal
#' distribution and stepsize to increase convergence speed.
#'
#' @param quiet Logical indicating whether to print any output on the
#' screen.
#'
#' @param seeds A vector of length \code{nchains} containing random number
#' seeds for the MCMC sampler.  If NULL, \code{nchains} random numbers are
#' generated from R's base random number generator which is controled outside
#' this routine using \code{set.seed}.  Note
#' that \code{set.seed} has no effect on the random number sequences used
#' in JAGS because JAGS is a separate package.
#' The seeds, whether chosen by this routine or specified, are
#' stored in the output object.
#'
#'
#' @details This routine replicates the M estimates of the 'Single Year' and
#' 'Multiple Classes' modules in package \code{eoa}.  To repeat either case,
#' input the composite g parameter's "a" and "b" parameters here, along
#' with the number of carcasses "X", and specify the "objective" prior. See
#' Examples.
#'
#' Due to discretizing the distribution for M, length of time this
#' routine takes to run depends on \code{Mmax}.  Larger \code{Mmax}
#' require much longer run times.
#'
#' There are two ways
#' to obtain the exact same results across multiple
#' runs.  One method is to specify
#' \code{seeds} here. This will set
#' the MCMC seeds in JAGS so that exact chains are reproduced.
#' For example, if \code{run1} is the
#' result of a previous call, \code{estimateM.EoA(X,b,seeds=run1$seeds)}
#' will replicate \code{run1} exactly.
#' The second method is to use R's default \code{set.seed}
#' just before calling this routine.
#'
#' @return List containing the following components.
#' \itemize{
#'   \item \code{M.ests} : A data
#'    frame containing the M estimates (point est and confidence interval).
#'   \item \code{out} : The full MCMC chain object.  Use this to
#'    check convergence, etc.
#'   \item \code{seeds} : The initial seeds used by JAGS.  Re-use this vector
#'   to replicated results exactly.
#'  }
#'
#' @examples
#' g.params <- list(alpha=600, beta=1200)
#' X <- 5
#' m.ests <- estimateM.EoA(X,g.params,Mmax=75)
#'
#' # Nice plots for checking convergence, correlation, density shape, etc.
#' library(lattice)
#' xyplot(m.ests$out)
#' acfplot(m.ests$out, ylim=c(-.2,1), lag.max=100)
#' densityplot(m.ests$out)
#'

estimateM.EoA <- function(X, beta.params, Mprior="objective", Mprior.mean, Mprior.sd,
                          Mmax=1000, conf.level=0.9,
                          niters = 10000,
                          nthins = 10,
                          nburns  = 20000,
                          nchains = 3,
                          nadapts = 2000,
                          quiet=FALSE,
                          seeds=NULL){

	## ---- bayesModelCode ----

	jagsModel <- "model{

	# Priors
	m ~ dcat(M.fx)
	g ~ dbeta(alpha, beta)

	M <- m-1

	# Likelihood
	X ~ dbin(g, M)
}
"

# for( i in 1:nx ){
# 	X[i] ~ dbin( g, M )
# }


writeLines(jagsModel, "model.txt")
#cat(jagsModel)

## ---- Mprior ----
M.x  <- 0:Mmax

if( Mprior == "normal"){
	# Use mean and sd passed in
	M.fx <- dnorm(M.x, Mprior.mean, Mprior.sd)
} else if(Mprior == "gamma"){
  shape <- Mprior.mean^2 / Mprior.sd^2
  scale <- Mprior.sd^2 / Mprior.mean
  M.fx <- dgamma(M.x, shape = shape, scale = scale )
} else {
	#	The following line is the "objective" prior from eoa
	M.fx <- sqrt(M.x + 1) - sqrt(M.x)   # This matches numbers scraped from eoa. See plotEoaPriors.r
}

if( is.null(seeds) ){
  tmp.mult <- 10^(8)
  seeds <- round(runif(nchains,0,tmp.mult))
} else if (length(seeds) != nchains) {
  stop(paste("MCMC random number seeds should be NULL or length", nchains ))
}


## ---- jagsObject -----



JAGS.data.0 <- list ( X = X,
											M.fx = M.fx,
											alpha = beta.params$alpha,
											beta = beta.params$beta
											 )



## ---- initialValues ----

Inits <- function(a,b,x,m.max,seed){
  g <- rbeta(1, a, b)
  m <- (x+1)/g
  m <- m + rnorm(1,0,0.15*m)
  m <- ifelse(m<0, 1, ifelse(m>m.max, m.max-1, round(m)))
  list ( m=m,
         g=g,
	       .RNG.name="base::Mersenne-Twister",
	       .RNG.seed=seed)
}

inits<-vector("list",nchains)
for(i in 1:nchains){
  inits[[i]]<-Inits(JAGS.data.0$alpha, JAGS.data.0$beta, X, Mmax, seeds[i])
}

# Parameters to be monitored by WinBUGS
params <- c("M", "g")


## ---- jagsRun ----

# Initialize the chains and adapt
library(rjags)
(t1=Sys.time())
jags = jags.model(file="model.txt",
									data=JAGS.data.0,
									inits=inits,
									n.chains=nchains,
									n.adapt=nadapts,
									quiet = quiet)

#   Run out in the chain a ways
if(!quiet) cat("Burning...\n")
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

# Check convergence
conv <- checkIsConverged(out,1.1,quiet)

# Check autocorrelation
auto <- checkIsAutocorrelated(out, 0.4, 2, quiet)

## ---- mcmcSummary ----
if(!quiet) print(summary(out))


## ---- quantiles ----
tmp <- as.matrix(out)
M.50 <- quantile(tmp[,"M"], 0.5)
M.CI <- quantile(tmp[,"M"], c((1-conf.level)/2, 1-(1-conf.level)/2))

## ---- STD ----
# We may never use this
M.sd <- sd(tmp[,"M"])



list(M.est=data.frame(M=M.50, M.sd=M.sd, M.lo=M.CI[1],
                      M.hi=M.CI[2], ci.level=conf.level),
     out=out,
     seeds=seeds)
}

