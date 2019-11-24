### X must be centred at zero!

fk_fmdh <- function(w, X, h, beta, alpha, C, n_eval = 500){
  x <- X%*%w/sqrt(sum(w^2))
  s <- sd(x)
  n <- length(x)
  fk_md(x,rep(1,n),1.2*seq(-alpha*s, alpha*s, length = n_eval),h,beta,alpha,C)/n/h
}

fk_dfmdh <- function(w, X, h, beta, alpha, C, n_eval = 500){
  nw <- sqrt(sum(w^2))
  x <- X%*%w/nw
  s <- sd(x)
  n <- length(x)
  dp <- fk_md_dp(x,rep(1,n),1.2*seq(-alpha*s, alpha*s, length = n_eval),h,beta,alpha,C)/n/h^2
  c(dp%*%(X/nw-x%*%t(w)/nw^2))
}


fk_mdh <- function(X, v0 = NULL, hmult = 1, beta = c(.25,.25), alphamax = 1){
  n <- nrow(X)
  mn <- colMeans(X)
  X <- sweep(X,2,mn,'-')
  if(is.null(v0)) v0 <- eigs_sym(cov(X), 1)$vectors[,1]
  else v0 <- v0/sqrt(sum(v0^2))
  h <- hmult*(8*sqrt(pi)*roughness_K(beta)/var_K(beta)^2/3/n)^.2*sd(X%*%v0)
  alpha_opt <- 0
  vopt <- v0
  alpha <- 0
  b <- 0
  C <- 100000
  while(alpha < alphamax){
    v0 <- optim(v0, fk_fmdh, fk_dfmdh, X, h, beta, alpha, C, method = 'BFGS', control = list(maxit=50))$par
    v0 <- v0/sqrt(sum(v0^2))
    x <- X%*%v0
    s <- sd(x)
    if(fk_is_minim_md(x,rep(1,n),2*seq(-alpha*s, alpha*s, length = 1000),h,beta,alpha,C)){
      vopt <- v0
      alphaopt <- alpha
      b <- fk_md_b(x,rep(1,n),1.1*seq(-alpha*s, alpha*s, length = 1000),h,beta,alpha,C)+mn%*%v0
    }
    alpha <- alpha + .1
  }
  list(v = v0, b = b)
}
