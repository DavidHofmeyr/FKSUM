f_ppr <- function(v, X, y, h, betas, loss = NULL, dloss = NULL, type = 'loc-lin'){
  p <- X%*%v/sqrt(sum(v^2))
  n <- length(p)
  o <- order(p)
  po <- p[o]
  yo <- y[o]
  if(type=='loc-lin'){
    sK <- ksum(po, rep(1, n), po, h, betas, 1:n)-betas[1]
    sKx <- ksum(po, po, po, h, betas, 1:n)-po*betas[1]
    sKy <- ksum(po, yo, po, h, betas, 1:n)-yo*betas[1]
    sKx2 <- ksum(po, po^2, po, h, betas, 1:n)-po^2*betas[1]
    sKxy <- ksum(po, po*yo, po, h, betas, 1:n)-po*yo*betas[1]
    hy <- ((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*po)/(sK*sKx2-sKx^2)
  }
  else{
    sK <- ksum(po, rep(1, n), po, h, betas, 1:n)-betas[1]
    sK[sK < 1e-20] <- 1e-20
    sKy <- ksum(po, yo, po, h, betas, 1:n)-yo*betas[1]
    hy <- sKy / sK
  }
  if(is.null(loss)) sum((yo-hy)^2)
  else sum(loss(yo, hy))
}


kLLreg <- function(x, y, h, betas, type){
  n <- length(x)
  o <- order(x)
  xo <- x[o]
  yo <- y[o]
  if(type == 'loc-lin'){
    sK <- ksum(xo, numeric(n)+1, xo, h, betas)
    sKx <- ksum(xo, xo, xo, h, betas)
    sKy <- ksum(xo, yo, xo, h, betas)
    sKx2 <- ksum(xo, xo^2, xo, h, betas)
    sKxy <- ksum(xo, xo*yo, xo, h, betas)
    (((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*xo)/(sK*sKx2-sKx^2))[rank(x)]
  }
  else{
    sK <- ksum(xo, numeric(n)+1, xo, h, betas)
    sKy <- ksum(xo, yo, xo, h, betas)
    (sKy / sK)[rank(x)]
  }
}


df_ppr <- function(v, X, y, h, betas, loss = NULL, dloss = NULL, type = 'loc-lin'){
  nv <- sqrt(sum(v^2))
  p <- X%*%v/nv
  n <- length(p)
  o <- order(p)
  po <- p[o]
  yo <- y[o]

  if(type == 'loc-lin'){
    S1 <- kndksum(po, numeric(n)+1, po, h, betas, 1:n)-rep(c(betas[1], 0), each = n)
    Sx <- kndksum(po, po, po, h, betas, 1:n)-cbind(po*betas[1], 0)
    Sy <- kndksum(po, yo, po, h, betas, 1:n)-cbind(yo*betas[1], 0)
    Sxx <- kndksum(po, po^2, po, h, betas, 1:n)-cbind(po^2*betas[1], 0)
    Sxy <- kndksum(po, po*yo, po, h, betas, 1:n)-cbind(po*yo*betas[1], 0)
    N <- ((Sxx[,1]*Sy[,1]-Sx[,1]*Sxy[,1])+(S1[,1]*Sxy[,1]-Sx[,1]*Sy[,1])*po)
    D <- (S1[,1]*Sxx[,1]-Sx[,1]^2)
    hy <- N/D
    if(is.null(loss)){
      dL <- 2*(hy-yo)
    }
    else dL <- dloss(yo, hy)
    dLD <- dL/D
    Sy.p_Spy_Sp.hy <- (Sy[,1]*po+Sxy[,1]-2*Sx[,1]*hy)*dLD
    S1.hy_Sy <- (S1[,1]*hy-Sy[,1])*dLD
    Sp_S1.p <- (Sx[,1]-S1[,1]*po)*dLD
    Sp.p_Spp <- (Sx[,1]*po-Sxx[,1])*dLD
    Spp.hy_Spy.p <- (Sxx[,1]*hy-Sxy[,1]*po)*dLD

    S_Sy.p_Spy_Sp.hy <- kndksum(po, Sy.p_Spy_Sp.hy, po, h, betas, 1:n) - betas[1]*cbind(Sy.p_Spy_Sp.hy, 0)
    S_S1.hy_Sy <- kndksum(po, S1.hy_Sy, po, h, betas, 1:n) - betas[1]*cbind(S1.hy_Sy, 0)
    S_Sp_S1.p <- kndksum(po, Sp_S1.p, po, h, betas, 1:n) - betas[1]*cbind(Sp_S1.p, 0)
    S_Sp.p_Spp <- dksum(po, Sp.p_Spp, po, h, betas, 1:n)
    S_Spp.hy_Spy.p <- dksum(po, Spp.hy_Spy.p, po, h, betas, 1:n)

    dp <- po/h*S_Sy.p_Spy_Sp.hy[,2] + po^2/h*S_S1.hy_Sy[,2] + po*yo/h*S_Sp_S1.p[,2] + yo/h*S_Sp.p_Spp
    dp <- dp + S_Spp.hy_Spy.p/h - 2*po*S_S1.hy_Sy[,1] - S_Sy.p_Spy_Sp.hy[,1] - yo*S_Sp_S1.p[,1]
    dp <- dp + dLD*(Sy[,1]*(po/h*Sx[,2]-1/h*Sxx[,2]-Sx[,1]) - 1/h*Sy[,2]*(Sxx[,1]-po*Sx[,1])-1/h*Sxy[,2]*(po*S1[,1]-Sx[,1]))
    dp <- dp + dLD*(Sxy[,1]*(Sx[,2]/h-po/h*S1[,2]+S1[,1]) + hy*(Sxx[,1]/h*S1[,2]+Sxx[,2]/h*S1[,1]-2/h*Sx[,1]*Sx[,2]))
  }
  else{
    S1 <- kndksum(po, numeric(n)+1, po, h, betas, 1:n)-rep(c(betas[1], 0), each = n)
    S1[S1[, 1] < 1e-20, 1] <- 1e-20
    Sy <- kndksum(po, yo, po, h, betas, 1:n)-cbind(yo*betas[1], 0)
    hy <- Sy[, 1] / S1[, 1]
    if(is.null(loss)){
      dL <- 2*(hy-yo)
    }
    else dL <- dloss(yo, hy)

    dp <- (dksum(po, hy * dL / S1[, 1], po, h, betas, 1:n) - yo * dksum(po, dL / S1[, 1], po, h, betas, 1:n) +
             dL / S1[, 1] * (hy * S1[, 2] - Sy[, 2])) / h
  }
  c(c(dp)%*%(X[o,]/nv-(po)%*%t(v)/nv^2))
}




fk_ppr <- function(X, y, nterms = 1, hmult = 1, betas = NULL, loss = NULL, dloss = NULL, initialisation = 'lm', type = 'loc-lin'){
  n <- nrow(X)
  d <- ncol(X)
  mu <- mean(y)
  mu_X <- colMeans(X)
  r <- y-mu
  X <- sweep(X, 2, mu_X, '-')
  hs <- numeric(nterms)
  vs <- matrix(0,nterms,d)
  fitted <- matrix(0,nterms,n)
  if(is.null(betas)) betas = c(0.25,0.25)
  for(tm in 1:nterms){
    if(initialisation=='lm') v0 <- solve(t(X)%*%X+.01*diag(ncol(X)))%*%t(X)%*%r
    else if(initialisation=='random') v0 <- rnorm(d)
    else if(is.function(initialisation)) v0 <- initialisation(X, r)
    else stop('argument "initialisation" must be a function of X and y or one of "lm" and "random".')
    v0 <- v0/sqrt(sum(v0^2))
    if(d < 100) h <- sqrt(eigen(cov(X))$values[1]) / n^.2 * hmult
    else h <- sqrt(eigs_sym(cov(X), 1)$values[1]) / n^.2 * hmult
    vs[tm,] <- optim(v0, f_ppr, df_ppr, X, r, h, betas, loss, dloss, type, method = 'L-BFGS-B', control = list(maxit = 50))$par
    vs[tm,] <- vs[tm,]/sqrt(sum(vs[tm,]^2))
    p <- X%*%vs[tm,]
    cv_fun <- function(h) f_ppr(matrix(1, 1, 1), p, r, h, betas, loss, dloss, type)
    h <- optimise(cv_fun, c(.05*h, h))$minimum
    fitted[tm,] <- kLLreg(X%*%vs[tm,], r, h, betas, type)
    r <- r - fitted[tm,]
    hs[tm] <- h
  }
  sol <- list(mu = mu, mu_X = mu_X, y = y, X = X, hs = hs, vs = vs, fitted = fitted, betas = betas, type = type)
  class(sol) <- "fk_ppr"
  sol
}


predict.fk_ppr <- function(object, Xtest = NULL, ...){
  if(is.null(Xtest)) Xtest <- object$X
  else Xtest <- sweep(Xtest, 2, object$mu_X, '-')
  betas <- object$betas
  n <- nrow(object$X)
  ntest <- nrow(Xtest)
  yhat <- object$mu
  for(tm in 1:length(object$hs)){
    h <- object$hs[tm]
    p <- object$X%*%object$vs[tm,]
    o <- order(p)
    ptest <- Xtest%*%object$vs[tm,]
    otest <- order(ptest)
    if(object$type == 'loc-lin'){
      sK <- ksum(p[o], numeric(n)+1, ptest[otest], h, betas)
      sKx <- ksum(p[o], p[o], ptest[otest], h, betas)
      sKy <- ksum(p[o], object$fitted[tm,o], ptest[otest], h, betas)
      sKx2 <- ksum(p[o], p[o]^2, ptest[otest], h, betas)
      sKxy <- ksum(p[o], p[o]*object$fitted[tm,o], ptest[otest], h, betas)
      yhat <- yhat + (((sKx2*sKy-sKx*sKxy)+(sK*sKxy-sKx*sKy)*ptest[otest])/(sK*sKx2-sKx^2))[rank(ptest)]
    }
    else{
      sK <- ksum(p[o], numeric(n)+1, ptest[otest], h, betas)
      sKy <- ksum(p[o], object$fitted[tm,o], ptest[otest], h, betas)
      yhat <- yhat + (sKy / sK)[rank(ptest)]
    }
  }
  yhat
}
