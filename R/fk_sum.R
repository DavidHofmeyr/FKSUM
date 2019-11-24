fk_sum <- function(x, omega, h, x_eval = NULL, beta = c(.25,.25), nbin = NULL, type = "ksum"){
  n <- length(x)
  mn <- mean(x)
  x <- x-mn
  if(is.null(nbin)){
    o <- order(x)
    xo <- x[o]
    omo <- omega[o]
    if(is.null(x_eval)){
      if(type=="ksum") ksum(xo,omo,x,h,beta,match(x,xo))
      else if(type=="dksum") dksum(xo,omo,x,h,beta,match(x,xo))
      else if(type=="both") kndksum(xo,omo,x,h,beta,match(x,xo))
    }
    else{
      x_eval <- x_eval - mn
      o_eval <- order(x_eval)
      x_eo <- x_eval[o_eval]
      if(type=="ksum") ksum(xo,omo,x_eo,h,beta)[match(x_eval,x_eo)]
      else if(type=="dksum") dksum(xo,omo,x_eo,h,beta)[match(x_eval,x_eo)]
      else if(type=="both") kndksum(xo,omo,x_eo,h,beta)[match(x_eval,x_eo),]
    }
  }
  else{
    xo <- seq(min(x)-1e-5, max(x)+1e-5, length = nbin)
    omo <- sm_bin_wts(x, omega, nbin, xo[1], xo[nbin])
    if(is.null(x_eval)){
      if(type=="ksum") ksum(xo,omo,x,h,beta,cbin_alloc(x,nbin,xo[1],xo[nbin]))
      else if(type=="dksum") dksum(xo,omo,x,h,beta,cbin_alloc(x,nbin,xo[1],xo[nbin]))
      else if(type=="both") kndksum(xo,omo,x,h,beta,cbin_alloc(x,nbin,xo[1],xo[nbin]))
    }
    else{
      x_eval <- x_eval - mn
      if(type=="ksum") ksum(xo,omo,x_eval,h,beta,cbin_alloc(x_eval,nbin,xo[1],xo[nbin]))
      else if(type=="dksum") dksum(xo,omo,x_eval,h,beta,cbin_alloc(x_eval,nbin,xo[1],xo[nbin]))
      else if(type=="both") kndksum(xo,omo,x_eval,h,beta,cbin_alloc(x_eval,nbin,xo[1],xo[nbin]))
    }
  }
}

norm_const_K <- function(beta){
  2*sum(sapply(1:length(beta), function(k) beta[k]*factorial(k-1)))
}

roughness_K <- function(beta){
  beta_use <- beta/norm_const_K(beta)
  betakj <- beta_use%*%t(beta_use)
  kpj <- log((2^(0:(length(beta)-1)))%*%t(2^(0:(length(beta)-1))), base = 2)
  sum(betakj/2^kpj*factorial(kpj))
}

var_K <- function(beta){
  beta_use <- beta/norm_const_K(beta)
  2*sum(sapply(1:length(beta), function(k) beta_use[k]*factorial(k+1)))
}

h_Gauss_to_K <- function(h, beta){
  h*(roughness_K(beta)/var_K(beta)^2*8*sqrt(pi)/4)^(1/5)
}

plot_kernel <- function(beta, type = 'l', ...){
  sig <- sqrt(var_K(beta))
  nc <- norm_const_K(beta)
  xs <- seq(-4, 4, length = 500)
  ys <- rowSums(sweep(sapply(0:(length(beta)-1), function(k) abs(xs*sig)^k), 2, beta, '*')*exp(-abs(xs*sig)))
  plot(xs, ys/nc*sig, xlab = 'x', ylab = 'K(x)', type = type, ...)
}

norm_K <- function(beta) sqrt(roughness_K(beta))

h_K_to_Gauss <- function(h, beta){
  h/(roughness_K(beta)/var_K(beta)^2*8*sqrt(pi)/4)^(1/5)
}
