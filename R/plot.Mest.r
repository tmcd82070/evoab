#' @export
#'
#' @title Plot M posterior distribution
#'
#' @description Plot the M (from \code{evoab}) posterior distribution.
#'
#' @param obj An object of class "Mest".  Usually output by the
#' routine \code{estimateM.EoA}.
#'
#'
plot.Mest <- function(obj){

  x <- obj$M.margin$M
  fx <- obj$M.margin$pdf.M
  max.fx <- max(fx)
  ml <- obj$M.est$M.lo
  mh <- obj$M.est$M.hi
  M <- obj$M.est$M
  conf.level <- obj$M.est$ci.level*100

  old.par <- par()
  par(mar=c(5.1,4.1,1,1))

  plot(x, fx, type="l", lwd=2,
       xlab="Mortalities (M)", ylab="",
       yaxt="n", ylim=c(0,max.fx*1.7), bty="n")
  axis(2, at=pretty(c(0,max.fx)))
  mtext(side=2, text="Prob. Density (pdf)", at=max.fx/2, line=2.5)

  ci.poly.x <- c(x,rev(x))
  ci.poly.y <- c(fx,rep(0,length(fx)))
  midInterval <- (ml <= ci.poly.x) & (ci.poly.x <= mh)
  ci.poly.x <- ci.poly.x[midInterval]
  ci.poly.y <- ci.poly.y[midInterval]

  polygon(ci.poly.x, ci.poly.y, col="cornsilk", border=NA)

  oneLineHgt <- par("cxy")[2]
  segments(ml,max.fx+4*oneLineHgt,mh,max.fx+4*oneLineHgt,col='grey60',lwd=13)

  points(ml,max.fx+4*oneLineHgt,bg="red",pch=21,cex=1.8,xpd=T)
  points(mh,max.fx+4*oneLineHgt,bg="red",pch=21,cex=1.8,xpd=T)
  points(M,max.fx+4*oneLineHgt,bg='turquoise',pch=21,cex=2.2,xpd=T)

  segments(ml,max.fx,ml,0,lty=3,lwd=2,col='darkred')
  segments(mh,max.fx,mh,0,lty=3,lwd=2,col='darkred')

  text(M,max.fx+5.5*oneLineHgt,paste("Estimate:\n median=",format(M,big.mark=",")),
       cex=.8,font=2,family="sans",col='darkslategray')
  text(ml,max.fx+2*oneLineHgt,paste0("LOWER\n",conf.level,"% Conf\n",format(ml,big.mark=",")),cex=.7,col="darkred",font=2,family="sans")

  text(mh,max.fx+2*oneLineHgt,paste0("UPPER\n",conf.level,"% Conf\n",format(mh,big.mark=",")),cex=.7,col="darkred",font=2,family="sans")

  #main.label<- verboseModels
  #assign('main.label',main.label,envir = serverProtected)

  # segments(0,0,ml,0)

}
