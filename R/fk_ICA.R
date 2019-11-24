f_ica <- function(w, X, h, betas, nbin = NULL){
  x <- c(X%*%w/sqrt(sum(w^2)))
  n <- length(x)
  if(is.null(nbin)){
    o <- order(x)
    xo <- x[o]
    hf <- ksum(xo, rep(1,n), xo, h, betas, 1:n)
    hf[which(hf < 1e-300)] <- 1e-300
    sum(log(hf))
  }
  else{
    m <- min(x)-1e-5
    M <- max(x)+1e-5
    ynew <- bin_wts(x, rep(1,n), nbin, m, M)
    xs <- seq(m, M, length = nbin)
    hf <- ksum(xs, ynew, xs, h, betas, 1:nbin)
    hf[which(hf < 1e-300)] <- 1e-300
    sum(log(hf)*ynew)
  }
}

#NumericVector ksum_bin(NumericVector x, NumericVector y, NumericVector x_eval, NumericVector counts, double h, int n, int n_eval, int ord, NumericVector betas){

df_ica <- function(w, X, h, betas, nbin = NULL){
  nw <- sqrt(sum(w^2))
  x <- c(X%*%w/nw)
  n <- length(x)
  if(is.null(nbin)){
    o <- order(x)
    xo <- x[o]
    hf <- ksum(xo, rep(1,n), xo, h, betas, 1:n)
    hf[which(hf < 1e-300)] <- 1e-300
    dhf <- -dksum(xo, rep(1,n), xo, h, betas, 1:n)/h
    d2 <- -dksum(xo, 1/hf, xo, h, betas, 1:n)/h
    dp <- (d2+dhf/hf)
    dp%*%(X[o,]/nw - xo%*%t(w)/nw^2)
  }
  else{
    m <- min(x)-1e-5
    M <- max(x)+1e-5
    ynew <- bin_wts(x, rep(1,n), nbin, m, M)
    xs <- seq(m, M, length = nbin)
    alloc <- cbin_alloc(x, nbin, m, M)
    hf <- ksum(xs, ynew, xs, h, betas, 1:nbin)
    hf[which(hf < 1e-300)] <- 1e-300
    dhf <- -dksum(xs, ynew, xs, h, betas, 1:nbin)/h
    ynew <- bin_wts(xs, 1/hf*ynew, nbin, m, M)
    d2 <- -dksum(xs, ynew, xs, h, betas, 1:nbin)/h
    dp <- (d2+dhf/hf)
    dp[alloc]%*%(X/nw - x%*%t(w)/nw^2)
  }
}

fk_ICA <- function(X, ncomp = 1, beta = c(.25,.25), hmult = 1.2, it = 20, nbin = NULL){
  dim <- ncomp

  ### define the objective function to be optimised, and its gradient. Functions LL() and dLL
  ### are in the c++ code. They compute the log likelihood and gradient for the projected
  ### points.
  W <- whiten(X, dim)
  Y <- W$Y
  V <- matrix(0, dim, dim)
  n <- nrow(X)
  ### the folowing is a very simple gradient descent implementation. You could also use
  ### optim() here. We can see which works better.
  for(i in 1:(dim-1)){
    v0 <- numeric(dim)
    v0[i] <- 1
    if(i>1){
      for(j in 1:(i-1)){
        v0 <- v0-(v0%*%V[,j])[1]*V[,j]
        v0 <- v0/sqrt(sum(v0^2))
      }
    }
    h <- hmult*(norm_K(beta)^2/var_K(beta)^2*8*sqrt(pi)/3/nrow(X))^(1/5)
    for(iter in 1:it){
      fval <- f_ica(v0, Y, h, beta, nbin)
      dir <- c(df_ica(v0, Y, h, beta, nbin))
      ndir <- sqrt(sum(dir^2))
      dir <- dir/ndir
      stp <- .5
      repeat{
        fnew <- f_ica(v0+stp*dir, Y, h, beta, nbin)
        if(fnew>(fval/.99999)) break
        else stp <- stp*.5
        if(stp<1e-9) break
      }
      v0 <- v0+stp*dir
      v0 <- v0/sqrt(sum(v0^2))
      if(stp<1e-9) break
    }
    v0 <- v0/sqrt(sum(v0^2))
    if(i>1){
      for(j in 1:(i-1)){
        v0 <- v0-(v0%*%V[,j])[1]*V[,j]
        v0 <- v0/sqrt(sum(v0^2))
      }
    }
    Y <- Y-Y%*%v0%*%t(v0)
    V[,i] <- v0
  }
  V[,dim] <- Null(V[,1:(dim-1)])
  list(X = X, K = W$E$vectors[,1:dim]%*%diag(1/sqrt(W$E$values[1:dim])), W = V, S = W$Y%*%V)
}

whiten <- function(X, ncomp){
  if(ncol(X)>300) E <- eigs_sym(cov(X),ncomp)
  else E <- eigen(cov(X))
  E$values[which(E$values<1e-10)] <- Inf
  Y <- t(t(X)-colMeans(X))%*%E$vectors[,1:ncomp]%*%diag(1/sqrt(E$values[1:ncomp]))
  list(E = E, Y = Y)
}

