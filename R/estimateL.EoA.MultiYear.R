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
estimateL.EoA.MultiYear <- function(lambda, beta.params, data, Lprior="vague",
                                    Lprior.mean=NULL, Lprior.sd=NULL,
                                    conf.level=0.9, nburns = 500000, niters = 20000,
                                    nthins = 10 ){
	library(rjags)

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


	## ---- bayesModelCode ----
	if( Lprior == "vague" ){
	  # Recall: tau of dnorm in JAGS is 1/variance
	  coefTaus <- compVagueSd(Y,alpha.vec,beta.vec,lambda.covars, range.multiplier = 100)
	  print(coefTaus)
	  coefTaus <- 1/coefTaus^2

		jagsModel <- "model{

		# Priors
    for(i in 1:ncovars){
      a[i] ~ dnorm( 0, coefTaus[i] )
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
													ncovars = ncol(lambda.covars),
													coefTaus = coefTaus,
													alpha = alpha.vec,
													beta = beta.vec,
													lambda.covars = lambda.covars)

} else if (Lprior == "normal"){
	jagsModel <- "model{

	# Priors
  for(i in 1:ncovars){
	  a[i] ~ dnorm( 0, 4)
  }

	# functional relations
	for(j in 1:nx){
	  logl[j] <- a[1]*lambda.covars[j,1]  # model must have at least one covar/intercept
	  for(i in 2:ncovars){
	    logl[j] <- logl[i-1] + a[i]*lambda.covars[j,i]
	  }
	  lambda[j] <- exp(logl[j])
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
												ncovars = ncol(lambda.covars),
												alpha = alpha.vec,
												beta = beta.vec,
												lmean = Lprior.mean,
												ltau = 1/Lprior.sd^2,
												lambda.covars = lambda.covars)
	#print(JAGS.data.0)

}

writeLines(jagsModel, "model.txt")
#cat(jagsModel)


## ---- initialValues ----

Inits <- function(x){
	if( !is.null(x$lmean) ){
		gg <-rbeta(x$nx, x$alpha, x$beta)
		M <- ceiling(x$Y / gg) + 1
		l <- mean(M)
		list ( lambda = l,
					 M = M,
					 g=gg,
					 .RNG.name="base::Wichmann-Hill",
					 .RNG.seed=93229 )
	} else {
		gg <-rbeta(x$nx, x$alpha, x$beta)
		M <- ceiling(x$Y / gg) + 1
		a <- rep(0,ncol(x$lambda.covars))
		a[1] <- log(mean(M))
		list ( a = a,
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
nchains <- 3
nadapt <- 3000


# Parameters to be monitored by WinBUGS
params <- c("a", "M", "g", "lambda", "Mtot","logl")


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
out = update(jags, n.iter=nburns)

# Run the MCMC chains
cat("MCMC Sampling...\n")
out = coda.samples(jags,
									 variable.names=params,
									 n.iter=niters,
									 thin=nthins)

(t2=Sys.time())
t3 <- t2-t1
cat(paste("Execution time:", round(t3,2), attr(t3,"units"), "\n\n"))

## ---- mcmcChecking ----
#print(summary(out))
out.coefs <- out[,grep("^a\\[",varnames(out))]

# Check convergence
cat(paste0(paste(rep("-",4),collapse="")," Coefficient convergence:\n"))
Rhats <- gelman.diag(out.coefs)
print(Rhats)
Rhats <- Rhats$psrf[,"Upper C.I."]
if( any(Rhats > 1.1) ){
  cat("The following coefficients probably did not converge (Rhat > 1.1):\n")
  print(names(Rhats)[Rhats > 1.1])
  cat("Try increasing nburns and/or nthins or informing the coefficient's priors.\n")
} else {
  cat("The model converged. All Rhats < 1.1\n")
}

# Check autocorrelation
cat(paste0("\n", paste(rep("-",4),collapse="")," Coefficient autocorrelation:\n"))
vn <- varnames(out.coefs)
ncovars <- ncol(lambda.covars)
acflag2 <- rep(NA,ncovars)
for( i in 1:ncovars ){
  tmp <- rep(NA,nchains)
  for(j in 1:nchains){
    tmp.2 <- acf(as.matrix(out.coefs[[j]][,i]), lag.max = 2, plot=FALSE)
    tmp[j] <- tmp.2$acf[tmp.2$lag ==2,,]
  }
  acflag2[i] <- mean(tmp)
  cat(paste("Coef.", vn[i], ": mean(acf[lag 2]) =", round(acflag2[i],5), "\n") )
}
if(any(acflag2>0.4)){
  cat("\nSome autocorrelations exceed 0.4. Inspect acfplot() and increase thinning.\n")
}



## ---- quantiles ----
# tmp <- as.matrix(out)
# l.50 <- quantile(tmp[,"lambda"], 0.5)
# alpha <- c(0,1) + c(1,-1)*(1-conf.level)/2
# l.CI <- quantile(tmp[,"lambda"], alpha)
#
# M.50 <- quantile(tmp[,"Mtot"], 0.5)
# M.CI <- quantile(tmp[,"Mtot"], alpha)
#
# ## ---- STD ----
# # We may never use this
# l.sd <- sd(tmp[,"lambda"])
# M.sd <- sd(tmp[,"Mtot"])
#


list(
  #lambda.est=data.frame(lambda=l.50, lambda.sd=l.sd, lambda.lo=l.CI[1], lambda.hi=l.CI[2],
	#												 Mtot=M.50, Mtot.sd=M.sd, Mtot.lo=M.CI[1], Mtot.hi=M.CI[2],
	#												 ci.level=conf.level),
  out=out,
  priorCoefSd = 1/sqrt(coefTaus),
  design.mat=lambda.covars)
}



# ---- Testing

# A simple example with three years, four observations
g <- data.frame(
  alpha = c( 69.9299, 63.5035,  84.6997, 84.6997),
  beta = c(  736.4795,  318.3179, 759.9333, 759.9333 )
  )
Y <- c( 0, 1, 3, 10 )

tmp.df <- data.frame(year=factor(c("2015","2016","2017","2017")),
                     Year=c(1,2,3,3))


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




eoa <- estimateL.EoA.MultiYear(Y~Year, g, data=tmp.df, niters= 2000*10, nthins = 10 )  # Un-informed EoA
tmp.1 <-apply(as.matrix(eoa$out),2,median)

#print(tmp.1)
