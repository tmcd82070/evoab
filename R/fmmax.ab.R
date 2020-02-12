#' @title fmmax.ab
#'
#' @description Compute the maximum M (mortality) we need to worry about
#' during EoA and EoAR calculations.
#'
#' @param x The number of carcasses (number of successes).
#'
#' @param pBa The alpha parameter of g's Beta distribution
#'
#' @param pBb The beta parameter of g's Beta distribution
#'
#' @return The M with less than 0.0001 probability of occurring.
#'
#' @author Dan Dalthorp
#'
#' @export
#'
#' @examples
#' beta.params <- data.frame(alpha=231.59, beta=2239.925)
#' x <- 40
#' fmmax.ab(x, beta.params$alpha, beta.params$beta)
#'
fmmax.ab <- function (x, pBa, pBb)
{
  zero <- 1e-4
  maxMmax <- 100000


  if (VGAM::pbetabinom.ab(x, maxMmax, pBa, pBb) > zero) {
    g <- pBa / (pBa + pBb)
    warning(paste0("P(X <= ", x, " | g = ", g, ", m = ", maxMmax, ") = ",
                   signif(pbinom(x, maxMmax, g), 6), ". Taking mmax = ", maxMmax, "..."))
    return(maxMmax)
  }
  mmax <- x
  while (1) {
    m <- mmax:(mmax + 100)
    mmax <- mmax + 100
    if (VGAM::pbetabinom.ab(x, size = mmax, shape1 = pBa,
                            shape2 = pBb) < 1e-04) {
      mmax <- m[min(which(VGAM::pbetabinom.ab(x, size = m,
                                              shape1 = pBa, shape2 = pBb) < 1e-04))]
      break
    }
  }
  mmax
}
