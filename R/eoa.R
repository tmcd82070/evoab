#' @export
#'
#' @title eoa - Evidence of Absence model fitting.
#'
#' @description  This routine estimates an EoA model.  An EoA model
#' consists of a log-linear for lambda, the mean number of search targets per "cell", where
#' "cell" could be a turbine, season, year, etc.  Inputs include
#' information about the number of targets found (\code{Y}) and the
#' g-values (=probatility of
#' discovery) at all
#' searched sites.  The method is Bayesian and allows either an uniform prior for lambda
#' or an informed prior.
#' Estimation is performed using JAGS.
#'
#' @param lambda A model formula for the lambda parameters of EoA.
#' This formulat has the form \code{Y ~ X1 + X2} + etc. (exactly like \link{\code{lm}}).
#' Here,  \code{Y} is a Vector of number of carcasses found, one element per "cell". For example, if
#' multiple sites are searched in a single season, elements of Y are the number
#' of targets found at each site.  If multiple sites are searched during multiple seasons,
#' \code{Y} contains number of targets found each site X season combination (length
#' of \code{Y} is number of sites times number of seasons, assuming no missing).  Covariate
#' vectors (or matricies) \code{X1} etc. must have the same length as \code{Y} and
#' can be anything that \code{formula} handles (e.g., factors, interactions, \code{I(x^2)}, etc.)
#'
#' @param beta.params A data frame containing, at a minimum, two columns named \code{$alpha}
#' and \code{$beta}.
#' These are the alpha and beta parameters of a Beta distribution that determine the g-values
#' in each "cell".
#' Length of \code{$alpha} and \code{$beta} is either one, in which case g is assumed constant
#' over "cells",
#' or one element per "cell".
#'
#' @param data An optional data frame, list or environment
#' (or object coercible by \code{as.data.frame} to a data frame)
#' containing the variables in the model. If variables are not found in data,
#' the variables are taken from environment(formula), typically the
#' environment from which \code{eoa} is called.
#'
#' @param priors An optional data frame specifying vague priors or
#' containing information
#' about (informed) prior distributions of coefficients in the lambda model.
#'
#' If \code{priors} = NULL, vague normal priors are assumed for
#' all coefficients.  Coefficients are on a log scale.  The
#' vague normal priors chosen here have mean 0 and very large standard deviations.
#' Check the output object to ensure the choosen standard deviations
#' are sufficiently large enough to be considered vague.
#'
#' If \code{priors} is not NULL, it should be a data frame containing
#' at a minimum \code{$mean} and \code{$sd} specifying the
#' mean and standard deviation of a normal prior for coefficients in the
#' lambda model.  The \code{row.names} of \code{priors} are matched to
#' coefficients, and this provides a way to use informed priors for some
#' coefficients and vague priors for others.  For example, if \code{X1} and
#' \code{X2} are both in the model and \code{priors}
#' contains only one row with \code{row.name} == "X1",  \code{$mean} and \code{$sd}
#' values for that row are used to inform the coefficient of \code{X1},
#' but vague priors are used for all other coefficients. See examples.
#'
#' @param conf.level The confidence level for all posterior confidence intervals.
#' Part of the value returned by this routine is \code{$intervals} containing
#' posterior confidence intervals for all parameters.  \code{conf.level} specifies
#' the two-tailed confidence level for all these confidence intervals.
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
#' @param quiet Logical indicating whether to print output during estimation.
#' \code{quiet==FALSE} shows text based progress bars during estimation.
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
#' @details
#' The EoA model implemented here is as follows:
#' \list{
#'   \item Target count Y[i] is assumed to be binomial(M[i],g[i]).
#'   \item Binomial index M[i] is assumed to be poisson(lambda[i]).
#'   \item Binomial probability g[i] is assumed to be Beta(alpha[i],beta[i]).
#'   \item Poisson rate lambda[i] is constrained to be a log-linear function
#' of covariaes. That is, log(lambda[i]) = A[0] + A[1]x[1i] + A[2]x[2i] + etc.
#'   \item Beta hyper-parameters alpha[i] and beta[i] are constants.
#' }
#'
#'
#' @return List containing the following components:
#' \list{
#'   \item \code{estimates} : a matrix containing parameter estimates and
#'   standard errors.  This matrix contains  one row per parameter and
#'   two columns.  \code{rownames(x$estimates)} lists the parameters.
#'   Columns are \code{Estimate}, containing the median value of
#'   of that parameter's posterior marginal, and \code{SD}, containing
#'   the standard deviation of the parameter's posterior marginal.  The
#'   following parameters are included:
#'   \list{
#'     \item \code{M} : Mortality estimates in individual "cells". For example,
#'     \code{M[4]} is the estimate (median of posterior) of mortalities at
#'     the site represented in cell 4 of the data set.
#'
#'     \item \code{Mtot} : A derived parameter, the median of the posterior
#'     distribution of the sum of all \code{M}'s.
#'
#'     \item \code{lambda} : Mortality rates, potentially averaged over
#'     multiple cells in the problem.  Lambda's average over cells according
#'     to the model, M's do not.  M's apply to single cells.
#'
#'     \item coefficients : One or more coefficients in the log model for
#'     lambda.
#'   }
#'
#'   \item \code{intervals} : a matrix containing two-tailed posterior confidence
#'   intervals.  Same dimension as \code{estimates} but columns are lower
#'   and upper endpoints of the intervals. Confidence level of all intervals
#'   is stored in \code{conf.level}.
#'
#'   \item \code{out} : a \code{mcmc.list} object containing the output
#'   of the JAGS runs. This is the raw output from the \code{coda.samples} function
#'   of the \code{coda} package and can be used to check convergence, etc. This
#'   object can be subsetted to chains for one or more parameter by thinking of
#'   the list as a 2-D matrix and using named dimensions.  For example,
#'   \code{x$out[,c("(Intercept)","year")]} produces an mcmc.list object containing only
#'   chains for parameters named "(Intercept)" and "year".  See examples for some
#'   useful plots.
#'
#'   \item \code{priors} : the mean and standard deviation of the normal
#'   prior distributions for all coefficients. This differs from the input
#'   \code{priors} parameter because matching of terms and row names has been done.
#'   There is exactly one row in the output \code{priors} for every parameter,
#'   while that is not necessarily true for the input.
#'
#'   \item \code{design.mat} : the lambda model's design matrix containing
#'   covariates values.
#'
#'   \item \code{seeds} : the vector of MCMC random number seeds actually used
#'   in the analysis.  Specify these as input seeds to repeat, exactly, the MCMC
#'   sampling.
#'
#'   \item \code{converged} : Logical value for whether this routine thinks the
#'   MCMC chain has converged.  The MCMC sampling is deemed to have converged if all Gelman
#'   R-hats are less than some criterion. This function simply calls
#'   \link{\code{checkIsConverged}} with an mcmc.list object containing coefficients only
#'   and a criterion of 1.1.  Convergence checking is only done on lambda model coefficients.
#'
#'   \item \code{Rhats} : Gelman's R-hats for every parameter.  This is one of the output
#'   components of \link{\code{checkIsConverged}}.
#'
#'   \item \code{autoCorrelated} :  Logical value indicating whether this routine
#'   thinks the MCMC is autocorrelated. Normally, one does not want autocorrelation
#'   in the final chains. This function simply calls \link{\code{checkIsAutocorrelated}}
#'   with the mcmc.list and \code{criterion=0.4} and \code{lag=2}.
#'
#'   \item \code{autoCorrs} : a vector of autocorrelations used the judge whether
#'   the chains are autocorrelated.  This is a vector with same length as number of
#'   coefficients, and is the second component of the output of \link{\code{checkIsAutocorrelated}}.
#'   Only model coefficients are checked for autocorrelation.
#'
#'   \item \code{conf.level} : the two-tailed confidence level for all confidence intervals
#'   in \code{intervals}.  Because the entire mcmc.list is returned, one can change
#'   confidence levels without re-running by executing \code{apply(as.matrix(x$out),
#'   2, quantile, p=c(<new quantiles>))}
#'
#' }
#'
#' @author Trent McDonald
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
#' # 3 year study of 7 sites. 21 "cells". lambda change = 20/year
#' set.seed(9430834) # fixes Y and g of this example, but not the actual chains
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
#' # Repeat and get exact same answers
#' eoa.1 <- eoa(Y~year, g, df, nburn = 1000, niters= 50*10, nthins = 10, seeds=eoa.1$seeds )
#'
#' # Run Informed EoA.  Assume prior annual lambda estimates
#' # are 10, 15, and 20, all with sd=.5.
#' prior <- data.frame(mean=c(log(10),log(15/10),log(20/10)),
#'   sd=c(0.5,0.5,0.5))
#' row.names(prior) <- c("(Intercept)","year2016","year2017")
#'
#' ieoa <- eoa(Y~year, g, df, priors=prior, nburn = 1000, niters= 50*10, nthins = 10 )
#'
#' # The above chains do not converge due to the low number of iterations used here.
#' # To check convergence and autocorrelation
#'
#' gelman.diag(ieoa$out) # gelmanStats
#' gelman.plot(ieoa$out) # gelmanPlot
#'
#' # Nice traceplots
#' library(lattice)
#' plot(ieoa$out) # tracePlot, all parameters
#' xyplot(ieoa$out[,c("(Intercept)","year2016","year2017")]) # nicer trace of coefficients
#'
#' # Autocorrelation functions
#' acfplot(ieoa$out[,c("(Intercept)","year2016","year2017")], ylim=c(-.2,1), lag.max=300)
#' acfplot(ieoa$out[,c("(Intercept)","year2016","year2017")], ylim=c(-.2,1), thin=2)
#'
#' # Density plots
#' densityplot(ieoa$out[,c("(Intercept)","year2016","year2017")])
#' densityplot(ieoa$out[,c("lambda[1]","lambda[8]","lambda[15]")]) # Mean lambda each year
#'
#' # Correlation among coefficients
#' levelplot(ieoa$out[,c("(Intercept)","year2016","year2017")][[1]])
#'
#' # QQ plots
#' qqmath(ieoa$out[,c("(Intercept)","year2016","year2017")])
#'
#'
eoa <- function(lambda, beta.params, data, priors=NULL,
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
  alpha <- c(0,0.5,1) + c(1,0,-1)*(1-conf.level)/2 # median plus two extreme-er quantiles
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
    design.mat=lambda.covars,
    seeds=seeds),
    conv,
    auto,
    conf.level=conf.level
  )


}



