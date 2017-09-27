loglikelihood <- function(X, Y, Y.c, ziMatrix, para) {
  # get dims
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  q <- dim(ziMatrix)[2]
  
  # reformat the parameters
  theta <- matrix(para[1:((p+q+1)*m)], p+q+1, m)
  
  B.current <- theta[1:p, ]
  XB.mat <- X %*% B.current
  p.current <- exp(XB.mat) / (1 + exp(XB.mat))   # n-by-m
  idx.large <- XB.mat >= 15
  idx.small <- XB.mat <= -15
  p.current[idx.large] <- exp(15) / (1 + exp(15))
  p.current[idx.small] <- 1 - exp(15) / (1 + exp(15))
  
  # psi <- theta[p+1, ]
  x <- apply(X %*% as.matrix(B.current), 2, mean)
  x2 <- x^2
  x3 <- x^3
  Coef <- para[((p+q+1)*m + 1):length(para)]
  psi <- Coef[1] + x * Coef[2] + x2 * Coef[3] + x3 * Coef[4]
  tmp <- exp(-psi)
  idx.large <- psi > 20
  idx.small <- psi < -20
  tmp[idx.large] <- exp(-20)
  tmp[idx.small] <- exp(20)
  s1 <- t(apply(p.current, 1, '*', tmp))   # n-by-m
  s2 <- t(apply((1 - p.current), 1, '*', tmp))   # n-by-m
  
  eta.current <- theta[(p+2):(p+q+1), ]
  zeta.mat <- ziMatrix %*% eta.current
  patzero.current <- exp(zeta.mat) / (1 + exp(zeta.mat))   # n-by-m
  idx.large <- zeta.mat >= 15
  idx.small <- zeta.mat <= -15
  patzero.current[idx.large] <- exp(15) / (1 + exp(15))
  patzero.current[idx.small] <- 1 - exp(15) / (1 + exp(15))
  
  loglike <- matrix(0, n, m)
  zeroInd <- Y == 0
  Y.c.aug <- matrix(rep(Y.c, m), n, m)
  
  # log-likelihood
  loglike[zeroInd] <- log(patzero.current[zeroInd] + 
                            (1 - patzero.current[zeroInd]) * 
                            beta(s1[zeroInd], Y.c.aug[zeroInd] + s2[zeroInd]) / 
                            beta(s1[zeroInd], s2[zeroInd]))
  
  loglike[!zeroInd] <- log(1 - patzero.current[!zeroInd]) + 
    lchoose(Y.c.aug[!zeroInd], Y[!zeroInd]) + 
    lbeta(Y[!zeroInd] + s1[!zeroInd], 
          Y.c.aug[!zeroInd] - Y[!zeroInd] + s2[!zeroInd]) - 
    lbeta(s1[!zeroInd], s2[!zeroInd])
  
  return(rowSums(loglike, na.rm = TRUE))
}