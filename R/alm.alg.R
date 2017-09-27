alm.alg <- function(para.init, s, lambda, grpIdx, grpSize, X, Y, Y.c, 
                    ziMatrix, K, rho, tol.fista, imax.fista, 
                    tol.alm, imax.alm) {
  # get dims
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  q <- dim(ziMatrix)[2]
  
  # initial Lagrangian multiplier w's
  theta <- matrix(para.init[1:((p+q+1)*m)], p+q+1, m)
  B.current <- theta[1:p, ]
  #*Re-parameterize over-dispersion
  psi <- theta[p+1, ]
  # phi <- theta[p+1, ]
  gamma.k <- para.init[((p+q+1)*m + 1):((p+q+1)*m + K + 1)]
  x <- apply(X %*% B.current, 2, mean)   # m-vector
  tmp <- rep(gamma.k[1], m)
  for (i in 2:(K+1)) {
    tmp <- tmp + gamma.k[i]*(x^(i-1))
  }
  # tmp <- log(phi / (1-phi)) - tmp
  tmp <- psi - tmp
  w.out <- rho * tmp
  
  # w.out <- rep(0.1, m)
  
  para.out <- para.init
  para.old <- para.init
  
  for (i in 1:imax.alm) {
    # find the unconstrained minimum
    # para.out <- fista.alg(para.old, s, lambda, grpIdx, grpSize, X, Y, Y.c, 
    # 											ziMatrix, p, q, K, w.out, rho, tol.fista, imax.fista)
    para.out <- cd.fista(para.old, s, lambda, grpIdx, grpSize, X, Y, Y.c, 
                         ziMatrix, p, q, K, w.out, rho, tol.fista, imax.fista)
    
    # update the multiplier w's
    theta <- matrix(para.out[1:((p+q+1)*m)], p+q+1, m)
    B.current <- theta[1:p, ]
    #*Re-parameterize over-dispersion
    psi <- theta[p+1, ]
    # phi <- theta[p+1, ]
    gamma.k <- para.out[((p+q+1)*m + 1):((p+q+1)*m + K + 1)]
    x <- apply(X %*% B.current, 2, mean)   # m-vector
    tmp <- rep(gamma.k[1], m)
    for (i in 2:(K+1)) {
      tmp <- tmp + gamma.k[i]*(x^(i-1))
    }
    # tmp <- log(phi / (1-phi)) - tmp
    tmp <- psi - tmp
    w.out <- w.out + rho * tmp
    
    # check for convergence
    dpara <- (para.out - para.old) / para.old
    bmax <- max(abs(dpara), na.rm = TRUE)
    # idx.max <- order(abs(dpara), decreasing = TRUE)[1:5]
    # print(paste("ALM change is", bmax))
    # print(paste("Corresponding parameter index is", idx.max))
    if (bmax < tol.alm) {break}
    
    # store previous updates
    para.old <- para.out
  }
  
  return(para.out)
}
