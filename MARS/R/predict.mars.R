predict.mars <- function(object,newdata) {
  if(missing(newdata) || is.null(newdata)) {
    B <- as.matrix(object$B)
  }
  else {
    tt <- terms(object$formula,data=newdata)
    tt <- delete.response(tt)
    mf <- model.frame(tt,newdata)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)[,-1] # remove intercept
    B <- make_B(X,object$Bfuncs)
  }
  beta <- object$coefficients
  drop(B %*% beta)
}

make_B <- function(X, Bfuncs){
  B <- matrix(1, nrow = nrow(X),ncol = length(Bfuncs))
  for(m in 2:length(Bfuncs)) {
    for(i in 1:nrow(Bfuncs[[m]])){
      B[,m] = B[,m]*h(X[,Bfuncs[[m]][i,"v"]],Bfuncs[[m]][i,"s"],Bfuncs[[m]][i,"t"])
    }
  }
  return(B)
}


