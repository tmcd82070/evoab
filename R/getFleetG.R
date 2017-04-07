#' @export
#'
#' @title getFleetG - Compute all-site (fleetwide) Evidence of Absence g-value.
#'
#' @description Compute a single all-site g value from site-season specifice g values.
#'
#' @param siteYearResults A data frame containing site-season results.
#' One usually obtains this by running the site-year module of the \code{eoa}
#' package on
#' every combination of year and site.  To assist with this, we at WEST have a python app
#' which copies and pastes the input parameters from CSV files into the
#' interface of
#' \code{eoa}.
#' Once site-year g-values are computed (using the auto-copy-paste routine),
#' we have a "scraper" function in R which
#' pulls parameters from all the text files that \code{eoa} produces.
#' This parameter is a data frame containing the output CSV results
#' from the "scraper".
#' At a minimum, \code{siteYearResults} must contain the following
#' columns:  \code{$species}, \code{$facility}, \code{$gFac.a} = the
#' alpha paramters of this facility's beta distribution, \code{$gFac.b} =
#' the beta parameter of this facility's beta distribution, and \code{$year}.
#'
#'
#' @param species Species abbreviation to run.
#'
#' @param weights A scalar or data frame containing weights for each
#' facility in the fleet. One row per
#' facility.  This is the "DWP" field of \code{eoa}. If \code{weights} = NULL,
#' all facilities receive equal weight.  Otherwise, \code{weights} must be
#' a data frame containing \code{$facility} and \code{$weight} columns.  This is merged
#' with \code{siteYearResults} using \code{$facility} as the key. These weights are
#' re-scaled to sum to 1.0 in the \code{mixtureBeta} function.
#'
#' @param yearWeights When more than one year at a facility is present, this
#' parameter allows different weights to be applied across years.  The length
#' of this vector must be equal to the number of years of data for all facilities that have
#' multiple years of data.  That is, an individual facility can have either 1 or length(yearWeights)
#' rows in siteYearResults.  If yearWeights == NULL, equal weights are use. These
#' weights are re-scaled to sum to 1.0 inside the \code{mixtureBeta} function.
#'
#' @return The fleet-wide g-value for a particular species.  This
#' value can then go into \code{estimateL.eoa}.
#'
#' @details
#'
#'
#' Seasons are summarized first.  That is, if a site appears more than
#' once for a species, we assume each line is a separate season
#' for that species.  g is averaged over those lines first.  At the end,
#' this average goes into the multi-site g average.
#'
#' @examples
#' syr <- data.frame(species=c("LBBA","LBBA","LBBA"),
#'    facility=c("f1","f2","f2"),
#'    gFac.a = c( 69.9299, 63.5035,  84.6997),
#'    gFac.b = c(  736.4795,  318.3179, 759.9333 ),
#'    year = c(2015,2015,2016))
#' getFleetG(syr, "LBBA")

getFleetG <- function(siteYearResults, species="LBBA", weights=NULL, yearWeights=NULL){

	df <- siteYearResults[siteYearResults$species == species, ]

	facilities <- unique(df$facility)

	cat(paste("---- Species =", species, "----\n"))


	#	compute or extract facility G's.  Must do this loop because
	# some facilities were sampled in multiple years.
	found.multiyrs <- FALSE
	fac.G <- NULL
	for( i in 1:length(facilities)){

		df.fac <- df[df$facility == facilities[i], ]

		if( nrow(df.fac) > 1){
			if(!found.multiyrs){
				cat("Found multiple years for the following facilities:\n")
				found.multiyrs <- TRUE
			}

			# We have multiple years for this facility. Average across years.
			fac.g <- mixtureBeta(df.fac$gFac.a, df.fac$gFac.b, w=yearWeights)

			#	Output some results for checking
			cat( paste0(facilities[i], ": years="))
			cat(df.fac$year)
			cat("\n")
		} else {
			mn <- df.fac$gFac.a / (df.fac$gFac.a + df.fac$gFac.b)
			fac.g <- list(alpha=df.fac$gFac.a, beta=df.fac$gFac.b, mean=mn)
		}

		fac.G <- rbind(fac.G, data.frame(facility=facilities[i], g=fac.g$mean, alpha=fac.g$alpha, beta=fac.g$beta))
	}

	# Tack on the weight vector to fac.G
	if( is.null(weights)){
		weights <- data.frame(facility=facilities, weight=1)
	} else {
		weights <- weights[,c("facility", "weight")]  # drop other columns, not really necessary
	}
	fac.G <- merge(fac.G, weights, by="facility", all.x=TRUE)

	#	Some print out to check
	cat("\n")
	cat(paste("Site-specific g's:", "\n"))
	print(fac.G)

	# Finally, compute the across-facility fleet-wide g
	ans <- mixtureBeta(fac.G$alpha, fac.G$beta, fac.G$weight)

	ans
}
