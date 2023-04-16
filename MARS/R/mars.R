#' Multivariate Adaptive Regression Splines (MARS)
#'
#' Fit Friedman's Multivariate Adaptive Regression Splines (MARS) model.
#'
#' @param formula an R formula
#' @param data a data frame containing the data
# ....
# .....

mars <- function(formula,data,control=mars.control()) {
  cc <- match.call() # save the call
  mf <- model.frame(formula,data)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)[,-1,drop=FALSE]
  x_names <- colnames(x)
  control <- validate_mars.control(control)
  fwd <- fwd_stepwise(y,x,control)
  bwd <- bwd_stepwise(fwd,control)
  fit <- lm(y~.-1,data=data.frame(y=y,bwd$B)) # notice -1 added
  out <- c(list(call=cc,formula=formula,y=y,B=bwd$B,Bfuncs=bwd$Bfuncs,
                x_names=x_names),fit)
  class(out) <- c("mars",class(fit))
  out
}

fwd_stepwise <- function(y,x,control=mars.control()){
  # Initialize:
  N <- length(y) # sample size
  n <- ncol(x) # number of predictors = number of X
  # B: a data frame with optimal basis function as columns
  B <- init_B(N,control$Mmax)
  # Bfuncs(): a list records the optimal (m, v, t) for basis functions (regions)
  Bfuncs <- vector(mode="list", length = control$Mmax+1)
  #---------------------------------------------------
  # Looping for forward selection:
  for(i in 1:(control$Mmax/2)) { # contrast to indexing 2...Mmax in Friedman
    M <- 2*i-1
    if(control$trace) cat("M", M, "\n")
    lof_best <- Inf
    for(m in 1:M) { # choose a basis function to split
      svars <- setdiff(1:n, Bfuncs[[m]][,"v"]) # vars not in B_m >> return 1:n if Bfuncs[[m]] null
      if(control$trace) cat("M, m, svars",M,m,svars,"\n")
      for(v in svars){ # select a variable to split on
        tt <- split_points(x[,v],B[,m])
        for(t in tt) {
          Bnew <- data.frame(B[,1:M],
                             # replace parent B[,m] with Btem1,Btem2
                             Btem1=B[,m]*h(x[,v], +1, t),
                             Btem2=B[,m]*h(x[,v], -1, t))
          gdat <- data.frame(y=y,Bnew)
          lof <- LOF(y~.-1,gdat,control)
          if(lof < lof_best) {
            lof_best <- lof
            best_split <- c(m=m,v=v,t=t)
          } # end if
        } # end loop over splits
      } # end loop over variables
    } # end loop over basis functions to split
    # save optimal (m, v, t) and update basis functions
    m <- best_split["m"]; v <- best_split["v"]; t <- best_split["t"]
    Bfuncs[[M+1]] <- rbind(Bfuncs[[m]], c(s=-1, v, t))
    Bfuncs[[M+2]] <- rbind(Bfuncs[[m]], c(s=1, v, t))
    B[,M+1:2] <- cbind(B[,m]*h(x[,v], -1, t), B[,m]*h(x[,v], 1, t))
  } # end loop over i
  colnames(B) <- paste0("B",(0:(ncol(B)-1)))
  return(list(y=y,B=B, Bfuncs=Bfuncs))
}


init_B <- function(N,Mmax){
  # Input: N- # of rows; Mmax: # of basis funcs
  # output: a N by (Mmax+1) dataframe
  B <- data.frame( matrix(NA,nrow=N,ncol=(Mmax+1)) )
  B[,1] <- 1 # first column for intercept: B0
  names(B) <- c("B0",paste0("B",1:Mmax))
  return(B)
}

bwd_stepwise <- function(fwd,control=mars.control()) {
  # fwd is a list with elements y, B, and Bfuncs
  Mmax <- ncol(fwd$B)-1
  Jstar <- 2:(Mmax+1)
  Kstar <- Jstar
  dat <- data.frame(y=fwd$y, fwd$B)
  lofstar <- LOF(y~.-1, dat, control)
  for(M in (Mmax+1):2){
    b <- Inf
    L <- Kstar
    if(control$trace) cat("L: ", L, "\n")
    for(m in L) {
      K <- setdiff(L,m)
      dat <- data.frame(y=fwd$y, fwd$B[,K])
      lof <- LOF(y~., dat, control)
      if(control$trace) cat("M:K:lof", M, ":", K, ":", lof, "\n")
      if(lof < b) {
        b <- lof
        Kstar <- K
      }
      if(lof < lofstar) {
        lofstar <- lof
        Jstar <- K
      }
    }
    if(control$trace) cat("M:Jstar:lofstar", M, ":", Jstar, ":", lofstar, "\n")
  }
  Jstar <- c(1,Jstar)
  return(list(y=fwd$y,B=fwd$B[,Jstar],Bfuncs=fwd$Bfuncs[Jstar]))
}

LOF <- function(form,data,control) {
  # update this LOF to GCV
  ff <- lm(form,data)
  RSS <- sum(residuals(ff)^2)
  N <- nrow(data)
  M <- length(coef(ff))-1
  Ctilda <- sum(diag(hatvalues(ff)))+control$d*M
  return(RSS * N/(N-Ctilda)^2)
}

h <- function(x,s,t) {
  # if x>t, s=+1, this return max(0,x-t)
  # if x<t, s=-1, this return max(0,t-x)
  return(pmax(0,s*(x-t)))
}

split_points <- function(xv,Bm){
  ## Remove Points where Bm(xv) <= 0
  for(i in seq_along(xv)){
    if(Bm[[i]] <= 0){
      xv <- xv[-i]
    }
  }
  ## Order unique values
  xv <- sort(unique(xv))
  ## Remove the highest value
  xv <- xv[-length(xv)]
  return(xv)
}

#------------------------------------------------------------------------
# constructor, validator and helper for class mars.control
#------------------------------------------------------------------------
#
new_mars.control <- function(control) {
  structure(control,class="mars.control")}

validate_mars.control <- function(control) {
  stopifnot(is.integer(control$Mmax),is.numeric(control$d),
            is.logical(control$trace))
  if(control$Mmax < 2) {
    warning("Mmax must be >= 2; Reset it to 2")
    control$Mmax <- 2}
  if(control$Mmax %% 2 > 0) {
    control$Mmax <- 2*ceiling(control$Mmax/2)
    warning("Mmax should be an even integer. Reset it to ",control$Mmax)}
  control
}


#' Constructor for `mars.control` objects
#'
#' This function constructs a `mars.control` object that specifies
#' parameters used in the model fitting procedure.
#'
#' @param Mmax Maximum number of basis functions. Should be an even integer. Default value is 2.
# .....
# ...

mars.control <- function(Mmax=2,d=3,trace=FALSE) {
  Mmax <- as.integer(Mmax)
  control <- list(Mmax=Mmax,d=d,trace=trace)
  control <- validate_mars.control(control)
  new_mars.control(control)
}
