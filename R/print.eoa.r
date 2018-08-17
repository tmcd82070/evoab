#' @name print.eoa
#'
#' @title Print a Evidence of Absence model
#'
#' @description Print method for EoA models produced by \code{eoa},
#' which are of class \code{eoa}.
#'
#' @param x An estimated eoa object from \code{eoa}.
#'
#' @param \dots Included for compatibility with other print methods.  Ignored here.
#'
#' @return The input value of \code{obj} is invisibly returned.
#' @author Trent McDonald, WEST Inc. \email{tmcdonald@west-inc.com}
#'
#' @seealso \code{\link{summary.eoa}}
#' @examples
#'
#' @keywords models
#' @export

print.eoa <- function( x, ... ){

  summary(x)

}
