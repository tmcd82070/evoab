#' @export
#'
#' @title Plot M posterior distribution
#'
#' @description Plot method for the M estimates
#' (from \code{evoab}) posterior distribution.
#'
#' @param obj An object of class "Mest".  Usually output by the
#' routine \code{estimateM.EoA}.
#'
#' @param plot.like Logical for whether to plot the likelihood for M as
#' well as the posterior.
#'
#' @param plot.prior Logical for whether to plot the prior for M as well
#' as the posterior.
#'
#' @details  \code{plot.like} and \code{plot.prior} are additive. Set both
#' to TRUE and the prior, likelihood, and posterior will all three be plotted.
#' If neither \code{plot.like} nor \code{plot.prior} are TRUE, the posterior
#' distribution is plotted.
#'
#' @author Trent McDonald
#'
#' @seealso \code{\link{estimateM.EoA}}
#'
#'
plot.Mest <- function(obj, plot.like=FALSE, plot.prior=FALSE){

  x <- obj$M.margin$M
  fx <- obj$M.margin$pdf
  max.fx <- max(fx)
  ml <- obj$M.est$M.lo
  mh <- obj$M.est$M.hi
  M <- obj$M.est$M
  conf.level <- obj$M.est$ci.level*100

  old.par <- par()
  par(mar=c(5.1,4.1,1,1))

  rng.x <- range(x)

  if(plot.like & plot.prior){
    rng.fx <- range(fx, obj$M.margin$like.pdf, obj$M.margin$prior.pdf)
  } else if( !plot.like & plot.prior){
    rng.fx <- range(fx, obj$M.margin$prior.pdf)
  } else if( plot.like & !plot.prior){
    rng.fx <- range(fx, obj$M.margin$like.pdf)
  } else if( !plot.like & !plot.prior){
    rng.fx <- range(fx)
  }
  max.fx <- rng.fx[2]


  plot(rng.x, rng.fx, type="n", lwd=2,
       xlab="Mortalities (M)", ylab="",
       yaxt="n", ylim=c(0,max.fx*1.7), bty="n")
  axis(2, at=pretty(c(0,max.fx)))
  mtext(side=2, text="Prob. Density (pdf)", at=max.fx/2, line=2.5)

  ci.poly.x <- c(x,rev(x))
  ci.poly.y <- c(fx,rep(0,length(fx)))
  midInterval <- (ml <= ci.poly.x) & (ci.poly.x <= mh)
  ci.poly.x <- ci.poly.x[midInterval]
  ci.poly.y <- ci.poly.y[midInterval]

  doTheLegend <- FALSE
  legEntries <- c("Posterior")
  legCols <- c("black")
  if(plot.like){
    lines(x, obj$M.margin$like.pdf, lwd=2, col="blue")
    doTheLegend <- TRUE
    legEntries <- c("Likelihood",legEntries)
    legCols <- c("blue", legCols)
  }
  if(plot.prior){
    lines(x, obj$M.margin$prior.pdf, lwd=2, col="red")
    doTheLegend <- TRUE
    legEntries <- c("Prior",legEntries)
    legCols <- c("red", legCols)
  }
  if(!plot.prior & !plot.like){
    polygon(ci.poly.x, ci.poly.y, col="cornsilk", border=NA)
  }
  lines(x, fx, lwd=2 )


  oneLineHgt <- par("cxy")[2]
  segments(ml,max.fx+4*oneLineHgt,mh,max.fx+4*oneLineHgt,col='grey60',lwd=13)

  points(ml,max.fx+4*oneLineHgt,bg="red",pch=21,cex=1.8,xpd=T)
  points(mh,max.fx+4*oneLineHgt,bg="red",pch=21,cex=1.8,xpd=T)
  points(M,max.fx+4*oneLineHgt,bg='turquoise',pch=21,cex=2.2,xpd=T)

  segments(ml,max.fx,ml,0,lty=3,lwd=2,col='darkred')
  segments(mh,max.fx,mh,0,lty=3,lwd=2,col='darkred')

  text(M,max.fx+5.5*oneLineHgt,paste("Estimate:\n median=",format(M,big.mark=",")),
       cex=.8,font=2,family="sans",col='darkslategray')

  text(ml,max.fx+2*oneLineHgt,paste0("LOWER\n",conf.level,"% Conf\n",format(ml,big.mark=",")),
       cex=.7,col="darkred",font=2,family="sans",xpd=NA)

  text(mh,max.fx+2*oneLineHgt,paste0("UPPER\n",conf.level,"% Conf\n",format(mh,big.mark=",")),
       cex=.7,col="darkred",font=2,family="sans",xpd=NA)

  if(doTheLegend){
    legend("topright",legend=legEntries, lty=1, col=legCols, lwd=2,
           inset=c(0,1-(max.fx/(max.fx+6*oneLineHgt))), cex=0.75)
  }


}
