mtTest <- function(lambda, data, offset){

  ## ---- lambdaModel ----
  # Resolve formula for lambda
  if (missing(data))
    data <- environment(lambda)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("lambda", "data", "offset"), names(mf), 0L)
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
  offset <- as.vector(model.offset(mf))
  ncovars <- ncol(lambda.covars)
  vnames<-dimnames(lambda.covars)[[2]]


ans <- c(
  list(k=1, ll=4:5),
  list(
    design.mat=lambda.covars,
    offset=offset,
    coef.labels=vnames,
    call=cl,
    data=data
  ),
  terms = mt
)

class(ans) <- "eoa"
ans

}



