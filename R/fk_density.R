fk_density <- function(x, h = 'Silverman', h_adjust = 1, beta = NULL, from = NULL, to = NULL, ngrid = 1000, nbin = NULL, x_eval = NULL){

  n <- length(x)

  # computation is slightly more stable for smaller magnitude sample, so centralise at zero
  mn <- mean(x)
  x <- x - mn
  if(!is.null(x_eval)) x_eval <- x_eval - mn
  if(!is.null(from)) from <- from - mn
  if(!is.null(to)) to <- to - mn

  # sample points must be sorted unless binning is used
  if(is.null(nbin) && is.unsorted(x)){
    o <- order(x)
    xo <-  x[o]
  }
  else if(is.null(nbin)) xo <- x

  # set kernel coefficients, or normalise if provided
  if(is.null(beta)) beta <- c(.25,.25)
  else beta <- beta/norm_const_K(beta)

  # determine bandwidth value if actual value not provided
  if(!is.numeric(h)){
    if(h=='Silverman') h <- (8*sqrt(pi)/3*roughness_K(beta)/var_K(beta)^2/length(x))^.2*sd(x)
    else if(h=='mlcv'){

      # mlcv bandwidth handled differently if binning is used
      if(!is.null(nbin)){
        frm <- min(x)
        tot <- max(x)
        wts <- sm_bin_wts(x,rep(1,n),nbin,frm,tot)
        xs <- seq(frm,tot,length=nbin)
        loo_cv <- function(h){
          hf <- ksum(xs,wts,xs,h,beta,1:nbin)-beta[1]
          hf[hf<1e-300] <- 1e-300
          sum(log(hf)*wts)-n*log(h)
        }
        h <- suppressWarnings(optimise(loo_cv,sd(x)/n^.2*c(.05,5),maximum = TRUE)$maximum)
      }
      else{
        loo_cv <- function(h){
          hf <- ksum(xo,rep(1,n),xo,h,beta,1:n)-beta[1]
          hf[hf<1e-300] <- 1e-300
          sum(log(hf))-n*log(h)
        }
        h <- suppressWarnings(optimise(loo_cv,sd(x)/n^.2*c(.05,5),maximum = TRUE)$maximum)
      }
    }
    else stop('"h" must be positive numeric or one of "Silverman" or "mlcv"')
  }
  h <- h*h_adjust

  # evaluate density
  # three possibilities are binned, grid evaluation or neither
  if(!is.null(nbin)){ # binning

    # establish range for grid of bin locations
    if(is.null(from)) from <- min(x)-6*h
    if(is.null(to)) to <- max(x)+6*h

    # compute bin values
    wts <- sm_bin_wts(x, rep(1,n), nbin, from, to)

    # set evaluation grid
    xs <- seq(from, to, length = nbin)

    # evaluation of density is handled differently if evaluation at the bin locations or specified set
    # of evaluation points
    if(is.null(x_eval)){

      # evaluate the density at bin locations
      y <- ksum(xs, wts, xs, h, beta, 1:nbin)

      # ensure normalisation is correct despite binning
      y <- y/sum(y)/(xs[2]-xs[1])
    }
    else{

      # determine location of evaluation points in bins
      alloc <- cbin_alloc(x_eval, nbin, from, to)

      # evaluate density at evaluation points
      y <- ksum(xs, wts, x_eval, h, beta, alloc)/n/h
      xs <- x_eval
    }
  }
  else if(!is.null(ngrid)){ # exact evaluation on a grid

    # establish range for grid
    if(is.null(from)) from <- xo[1]-6*h
    if(is.null(to)) to <- xo[n]+6*h

    # set evaluation points to grid locations
    xs <- seq(from, to, length = ngrid)

    # evaluate density with sorted sample at grid locations
    y <- ksum(xo, rep(1,n), xs, h, beta)/n/h
  }
  else{ # exact evaluation at non-grid locations
    if(is.null(x_eval)){ # default is evaluation at the sample points
      xs <- x
      xe <- xo
    }
    else{ # otherwise ensure evaluation points are sorted for fast kernel summing
      xs <- x_eval
      if(is.unsorted(x_eval)) xe <- sort(x_eval)
      else xe <- x_eval
    }

    # evaluate density at evaluation points
    y <- ksum(xo,rep(1,n),xe,h,beta)[match(xs,xe)]/n/h
  }

  # return evaluation points (ensure to add the mean which was subtracted initially) and corresponding density values, plus bandwidth
  list(x=xs+mn, y=y, h=h)
}
