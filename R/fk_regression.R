fk_regression <- function(x, y, h = 'amise', beta = NULL, from = NULL, to = NULL, ngrid = 1000, nbin = NULL, type = 'loc-lin'){
  n <- length(x)
  if(is.null(nbin) & is.unsorted(x)){
    o <- order(x)
    xo <-  x[o]
    yo <- y[o]
  }
  else if(!is.null(nbin)){
    xo <- x
    yo <- y
  }
  if(is.null(beta)) beta <- c(.25,.25)
  if(h=='amise'){
    h0 <- sd(x)/n^.2
    if(!is.null(nbin)){
      xs <- seq(min(x), max(x), length = nbin)
      wts1 <- sm_bin_wts(x, rep(1,n), nbin, xs[1], xs[nbin])
      wtsy <- sm_bin_wts(x, y, nbin, xs[1], xs[nbin])
      sK <- ksum(xs, wts1, xs, h0, beta, 1:nbin)
      sKy <- ksum(xs, wtsy, xs, h0, beta, 1:nbin)
      dsK <- dksum(xs, wts1, xs, h0, beta, 1:nbin)
      dsKy <- dksum(xs, wtsy, xs, h0, beta, 1:nbin)
      df <- (dsK*sKy-sK*dsKy)/sK^2/h0
      sK2 <- ksum(xs, wts1, xs, h0/1000, beta, 1:nbin)
      sKy2 <- ksum(xs, wtsy, xs, h0/1000, beta, 1:nbin)
      sig2 <- mean((sKy2/sK2-sKy/sK)^2)
      h <- (sig2*roughness_K(beta)/var_K(beta)^2/sum(sK[-1]/(n-1)/h0*diff(df)^2/(xs[2]-xs[1]))/n/2)^.2
    }
    else{
      sK <- ksum(xo, rep(1,n), xo, h0, beta, 1:n)
      sKy <- ksum(xo, yo, xo, h0, beta, 1:n)
      dsK <- dksum(xo, rep(1,n), xo, h0, beta, 1:n)
      dsKy <- dksum(xo, yo, xo, h0, beta, 1:n)
      df <- (dsK*sKy-sK*dsKy)/sK^2/h0
      sig2 <- mean((yo-sKy/sK)^2)
      h <- (sig2*roughness_K(beta)/var_K(beta)^2/mean((diff(df)/diff(xo))^2)/n/2)^.2
    }
  }
  else if(h=='cv'){
    if(!is.null(nbin)){
      message('Cross validation bandwidth estimation not implemented for binned estimator.
          Switching to AMISE optimal estimate \n')
      return(fk_regression(x, y, h = 'amise', beta, from, to, ngrid, nbin, type))
    }
    else{
      if(type=='loc-lin'){
        loo_cv <- function(h){
          sum((fk_loc_lin(xo, yo, h, beta, loo=1)-yo)^2)
        }
      }
      else{
        loo_cv <- function(h){
          sum((fk_NW(xo, yo, h, beta, loo=1)-yo)^2)
        }
      }
      h <- optimise(loo_cv, sd(x)/n^.2*c(.05, 2))$minimum
    }
  }
  else if(!is.numeric(h) || h < 0) stop('"h" must be either positive numeric or one of "cv" and "amise"')
  if(!is.null(nbin)){
    if(is.null(from)) from <- min(x)
    if(is.null(to)) to <- max(x)
    if(type=='loc-lin'){
      list(x = seq(from, to, length = nbin), y = fk_loc_lin(xo, yo, h, beta, nbin = nbin, from = from, to = to), h = h)
    }
    else{
      list(x = seq(from, to, length = nbin), y = fk_NW(xo, yo, h, beta, nbin = nbin, from = from, to = to), h = h)
    }
  }
  else if(!is.null(ngrid)){
    if(is.null(from)) from <- xo[1]
    if(is.null(to)) to <- xo[n]
    if(type=='loc-lin'){
      list(x = seq(from, to, length = ngrid), y = fk_loc_lin(xo, yo, h, beta, ngrid = ngrid, from = from, to = to), h = h)
    }
    else{
      list(x = seq(from, to, length = ngrid), y = fk_NW(xo, yo, h, beta, ngrid = ngrid, from = from, to = to), h = h)
    }
  }
  else{
    if(type=='loc-lin'){
      list(x = xo, y = fk_loc_lin(xo, yo, h, beta), h = h)
    }
    else{
      list(x = xo, y = fk_NW(xo, yo, h, beta), h = h)
    }
  }
}



fk_loc_lin <- function(x, y, h, beta, nbin = NULL, ngrid = NULL, from = NULL, to = NULL, loo = 0){
  n <- length(x)
  if(!is.null(nbin)){
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = nbin)
    wts1 <- sm_bin_wts(x, rep(1,n), nbin, xs[1], xs[nbin])
    wtsx <- sm_bin_wts(x, x, nbin, xs[1], xs[nbin])
    wtsy <- sm_bin_wts(x, y, nbin, xs[1], xs[nbin])
    wtsx2 <- sm_bin_wts(x, x^2, nbin, xs[1], xs[nbin])
    wtsxy <- sm_bin_wts(x, x*y, nbin, xs[1], xs[nbin])
    sK <- ksum(xs, wts1, xs, h, beta, 1:nbin) - loo*beta[1]
    sKx <- ksum(xs, wtsx, xs, h, beta, 1:nbin) - loo*beta[1]*wtsx/wts1
    sKy <- ksum(xs, wtsy, xs, h, beta, 1:nbin) - loo*beta[1]*wtsy/wts1
    sKx2 <- ksum(xs, wtsx2, xs, h, beta, 1:nbin) - loo*beta[1]*wtsx2/wts1
    sKxy <- ksum(xs, wtsxy, xs, h, beta, 1:nbin) - loo*beta[1]*wtsxy/wts1
    (((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*xs)/(sK*sKx2-sKx^2))
  }
  else if(!is.null(ngrid)){
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = ngrid)
    sK <- ksum(x, rep(1,n), xs, h, beta)
    sKx <- ksum(x, x, xs, h, beta)
    sKy <- ksum(x, y, xs, h, beta)
    sKx2 <- ksum(x, x^2, xs, h, beta)
    sKxy <- ksum(x, x*y, xs, h, beta)
    (((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*xs)/(sK*sKx2-sKx^2))
  }
  else{
    sK <- ksum(x, numeric(n)+1, x, h, beta, 1:n) - loo*beta[1]
    sKx <- ksum(x, x, x, h, beta, 1:n) - loo*beta[1]*x
    sKy <- ksum(x, y, x, h, beta, 1:n) - loo*beta[1]*y
    sKx2 <- ksum(x, x^2, x, h, beta, 1:n) - loo*beta[1]*x^2
    sKxy <- ksum(x, x*y, x, h, beta, 1:n) - loo*beta[1]*x*y
    (((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*x)/(sK*sKx2-sKx^2))
  }
}


fk_NW <- function(x, y, h, beta, nbin = NULL, ngrid = NULL, from = NULL, to = NULL, loo = 0){
  n <- length(x)
  if(!is.null(nbin)){
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = nbin)
    wts1 <- sm_bin_wts(x, rep(1,n), nbin, xs[1], xs[nbin])
    wtsy <- sm_bin_wts(x, y, nbin, xs[1], xs[nbin])
    sK <- ksum(xs, wts1, xs, h, beta, 1:nbin) - loo*beta[1]
    sKy <- ksum(xs, wtsy, xs, h, beta, 1:nbin) - loo*beta[1]*wtsy/wts1
    sKy/sK
  }
  else if(!is.null(ngrid)){
    if(is.null(from)) from <- x[1]
    if(is.null(to)) to <- x[n]
    xs <- seq(from, to, length = ngrid)
    sK <- ksum(x, numeric(n)+1, xs, h, beta)
    sKy <- ksum(x, y, xs, h, beta)
    sKy/sK
  }
  else{
    sK <- ksum(x, numeric(n)+1, x, h, beta, 1:n) - loo*beta[1]
    sKy <- ksum(x, y, x, h, beta, 1:n) - loo*beta[1]*y
    sKy/sK
  }
}
