#
# Copyright (C) weiya
#
#
# 2016-12-08
# Given Auto-Covariance Function to determine the parameter of ARMA model
#

myacf2arma <- function(x, p, q)
{
  #if(!missing(order))
  #  if(!is.numeric(order) || length(order) != 2L || any(order < 0))
  #    stop("'order' must be a non-negative numeric vector of length 3")
  YW <- matrix(0,p,p)
  for (i in 1:p)
    YW[i,] <- x[seq(q+i,length.out = p,by=-1)]
  b <- x[seq(q+2,length.out = p,by=1)]
  a <- solve(YW,matrix(b))
  a <- rbind(-1, a)
  # gamma_y
  gamma_y = c(0:q)
  for (k in 0:q)
  {
    Yk <- matrix(0, p+1, p+1)
    for (i in 0:p)
    {
      idx <- seq(k-i, length.out = p+1,by=1)
      # rho(-k)=rho(k)
      idx <- abs(idx)
      Yk[i+1,] <- x[idx+1]
    }
    gamma_y[k+1] <- t(a) %*% Yk %*% a
  }
  # the coefficient of ma
  A <- diag(q-1)
  A <- cbind(0, A)
  A <- rbind(A, 0)
  C <- c(1, rep(0, q-1))
  C <- as.matrix(C)
  K = 50
  GammaK <- matrix(0, K, K)
  qy <- length(gamma_y) - 1
  for (i in (qy+1):(K-qy-1))
    GammaK[i,(i-qy):(i+qy)] <- c(rev(gamma_y[-1]),gamma_y[1],gamma_y[-1])
  for (i in 1:qy)
    GammaK[i,1:(i+qy)] <- c(rev(gamma_y[-1][0:(i-1)]),gamma_y[1],gamma_y[-1])
  for (i in (K-qy):K)
    GammaK[i,(i-qy):K] <- c(rev(gamma_y[-1]),gamma_y[1],gamma_y[-1][0:(K-i)])
  # Omega k
  Omegak <- matrix(0, qy, K)
  for (i in 1:qy)
    Omegak[i,1:qy-i+1] <- gamma_y[-1][i:qy]

  # PI
  PI <- Omegak %*% solve(GammaK) %*% t(Omegak)
  sigma2 <- gamma_y[1] - t(C) %*% PI %*% C
  bq <- (as.matrix(gamma_y[-1]) - A %*% PI %*% C)/as.numeric(sigma2)
}
