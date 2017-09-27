gammaFunc <- function(para, theta.init, X, Y, K, w, rho) {
  # get dims
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  
  # reformat the parameters
  theta <- matrix(theta.init, p+q+1, m)
  
  # beta
  B.current <- theta[1:p, ]
  XB.mat <- X %*% B.current
  # E(\mu): n-by-m
  p.current <- exp(XB.mat) / (1 + exp(XB.mat))
  idx.large <- XB.mat >= 15
  idx.small <- XB.mat <= -15
  p.current[idx.large] <- exp(15) / (1 + exp(15))
  p.current[idx.small] <- 1 - exp(15) / (1 + exp(15))
  
  #*Re-parameterize over-dispersion
  psi <- theta[p+1, ]
  
  # gamma
  gamma.k <- para
  meanXB <- apply(XB.mat, 2, mean)   # m-vector
  tmp.vec <- rep(0, m)
  for (i in 1:(K+1)) {
    tmp.vec <- tmp.vec + gamma.k[i]*(meanXB^(i-1))
  }
  
  # output
  tmp.val <- psi - tmp.vec
  obj.val <- sum(w * tmp.val) + rho * sum(tmp.val^2) / 2
  return(obj.val)
}