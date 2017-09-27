cd.fista <- function(para.init, s, lambda, grpIdx, grpSize, X, Y, Y.c, 
                     ziMatrix, p, q, K, w, rho, tol.fista, imax.fista) {
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  theta.init <- para.init[1:((p+q+1)*m)]
  theta1 <- theta.init
  theta2 <- theta.init
  theta.out <- theta.init
  gamma.out <- para.init[((p+q+1)*m + 1):((p+q+1)*m + K + 1)]
  gamma1 <- gamma.out
  G <- length(grpSize)
  N <- length(theta.init)
  
  for (i in 3:imax.fista) {
    ## update gamma, by optim function
    temp <- try(optim(par = gamma.out, fn = gammaFunc, theta.init = theta.out,
                      X = X, Y = Y, K = K, w = w, rho = rho, hessian = TRUE,
                      control = list(reltol = 1e-200, maxit = 100)), TRUE)
    if (class(temp) != "try-error") {
      gamma.tmp <- as.vector(temp$par)
      idx.na <- is.na(gamma.tmp)
      idx.inf <- is.infinite(gamma.tmp)
      gamma.tmp[idx.na] <- gamma.out[idx.na]
      gamma.tmp[idx.inf] <- gamma.out[idx.inf]
      gamma.out <- gamma.tmp
    }
    
    ## update other parameters, by FISTA
    # extrapolation
    v <- theta1 + (i-2) / (i+1) * (theta1 - theta2)
    
    # proximate gradient descent
    # tmp.v <- v - s * grad(augLagrangian, v, X = X, Y = Y, Y.c = Y.c, 
    # 											ziMatrix = ziMatrix, K = K, w = w, rho = rho)
    # tmp.v <- v - s * gradAL(v, X = X, Y = Y, Y.c = Y.c, ziMatrix = ziMatrix, 
    # 												K = K, w = w, rho = rho)
    tmp.v <- v - s * gradAL.cd(c(v, gamma.out), X = X, Y = Y, Y.c = Y.c, 
                               ziMatrix = ziMatrix, K = K, w = w, rho = rho)
    
    theta.out <- tmp.v
    
    for (g in 1:G) {
      # construct the index for parameters that belong to group g
      idx <- rep(0, 0)
      groupings <- grpIdx[[g]]
      for (j in 1:length(groupings)) {
        # tmp.idx <- (groupings[j] - 1) * (p+1+q) + c(1:(p+1+q))
        tmp.idx <- (groupings[j] - 1) * (p+1+q) + 2
        idx <- c(idx, tmp.idx)
      }
      # prox operation
      theta.out[idx] <- proxGL(tmp.v[idx], s, lambda, grpSize[g])
    }
    
    # check for convergence
    dtheta <- (c(theta.out, gamma.out) - c(theta1, gamma1)) / c(theta1, gamma1)
    bmax <- max(abs(dtheta), na.rm = TRUE)
    # idx.max <- order(abs(dtheta), decreasing = TRUE)[1:5]
    # print(paste("CD change is", bmax))
    # print(paste("Corresponding parameter index is", idx.max))
    if (bmax < tol.fista) {break}
    
    # store previous updates
    gamma1 <- gamma.out
    theta2 <- theta1
    theta1 <- theta.out
  }
  
  return(c(theta.out, gamma.out))
}