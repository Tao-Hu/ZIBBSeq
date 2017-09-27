gradAL.cd <- function(para, X, Y, Y.c, ziMatrix, K, w, rho) {
  # get dims
  n <- dim(Y)[1]
  m <- dim(Y)[2]
  p <- dim(X)[2]
  q <- dim(ziMatrix)[2]
  
  # reformat the parameters
  theta <- matrix(para[1:((p+q+1)*m)], p+q+1, m)
  
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
  tmp <- exp(-psi)
  idx.large <- psi > 20
  idx.small <- psi < -20
  tmp[idx.large] <- exp(-20)
  tmp[idx.small] <- exp(20)
  # if(sum(tmp == 0) > 0) {
  #   print(paste("In function gradAL, over-dispersion parameter too large!",
  #               sum(tmp == 0)))
  # }
  # if(sum(is.infinite(tmp)) > 0) {
  #   print(paste("In function gradAL, over-dispersion parameter too small!",
  #               sum(is.infinite(tmp))))
  # }
  # if(sum(is.na(tmp)) > 0) {
  #   print(paste("In function gradAL, over-dispersion parameter NA found!",
  #               sum(is.na(tmp))))
  # }
  s1 <- t(apply(p.current, 1, '*', tmp))   # n-by-m
  s2 <- t(apply((1 - p.current), 1, '*', tmp))   # n-by-m
  
  eta.current <- theta[(p+2):(p+q+1), ]
  zeta.mat <- ziMatrix %*% eta.current
  # \pi_{ij}: n-by-m
  patzero.current <- exp(zeta.mat) / (1 + exp(zeta.mat))
  idx.large <- zeta.mat >= 15
  idx.small <- zeta.mat <= -15
  patzero.current[idx.large] <- exp(15) / (1 + exp(15))
  patzero.current[idx.small] <- 1 - exp(15) / (1 + exp(15))
  
  gamma.k <- para[((p+q+1)*m + 1):((p+q+1)*m + K + 1)]
  meanXB <- apply(XB.mat, 2, mean)   # m-vector
  
  #*Re-parameterize over-dispersion
  Y.c.aug <- matrix(rep(Y.c, m), n, m)
  # phi.aug <- matrix(rep(phi, n), n, m, byrow = TRUE)
  psi.aug <- matrix(rep(psi, n), n, m, byrow = TRUE)
  
  # mat1 = choose(s, y) * B1() / B2(), n-by-m
  mat1.tmp <- lchoose(Y.c.aug, Y) + lbeta(Y+s1, Y.c.aug-Y+s2) - lbeta(s1, s2)
  mat1 <- exp(mat1.tmp)
  idx.large <- mat1.tmp > 20
  idx.small <- mat1.tmp < -20
  mat1[idx.large] <- exp(20)
  mat1[idx.small] <- exp(-20)
  # mat2 = (1 - pi) * mat1, n-by-m
  mat2 <- (1 - patzero.current) * mat1
  # mat3 = pi * indicator(y == 0) + mat2, density for the zibb model
  mat3 <- patzero.current * (Y == 0) + mat2
  # mat4 = indicator(y == 0) - mat1
  mat4 <- (Y == 0) - mat1
  # mat5 = mat2 / mat3
  mat5 <- mat2 / mat3
  # mat6 = w + rho * (logit(phi) - \sum (gamma_k * meanXB^k)), m-vector
  mat6 <- rep(gamma.k[1], m)
  tmp.vec <- rep(0, m)
  for (i in 2:(K+1)) {
    mat6 <- mat6 + gamma.k[i]*(meanXB^(i-1))
    tmp.vec <- tmp.vec + gamma.k[i]*(i-1)*(meanXB^(i-2))
  }
  # mat6 <- w + rho * (log(phi / (1-phi)) - mat6)
  mat6 <- w + rho * (psi - mat6)
  
  # Derivative with respect to beta
  # grad.beta0, grad.beta1: m-vector
  grad.beta.tmp <- matrix(0, n, m)
  # grad.beta.tmp <- -mat5 * p.current * (1 - p.current) * (1 - phi.aug) * 
  # 	(digamma(Y+s1)-digamma(Y.c.aug-Y+s2)-digamma(s1)+digamma(s2)) / phi.aug
  tmp.val <- exp(-psi.aug)
  idx.large <- psi.aug > 20
  idx.small <- psi.aug < -20
  tmp.val[idx.large] <- exp(-20)
  tmp.val[idx.small] <- exp(20)
  grad.beta.tmp <- -mat5 * p.current * (1 - p.current) * tmp.val * 
    (digamma(Y+s1)-digamma(Y.c.aug-Y+s2)-digamma(s1)+digamma(s2))
  x0.aug <- matrix(rep(X[,1], m), n, m)
  x1.aug <- matrix(rep(X[,2], m), n, m)
  grad.beta0 <- colSums(grad.beta.tmp * x0.aug) - mat6 * tmp.vec * mean(X[,1])
  grad.beta1 <- colSums(grad.beta.tmp * x1.aug) - mat6 * tmp.vec * mean(X[,2])
  
  # # Derivative with respect to phi
  # # grad.phi: m-vector
  # grad.phi.tmp <- matrix(0, n, m)
  # grad.phi.tmp <- -mat5 * 
  # 	(digamma(s2) - digamma(Y.c.aug-Y+s2) - 
  # 	 	(digamma(Y+s1)-digamma(Y.c.aug-Y+s2)-digamma(s1)+digamma(s2))*p.current) / 
  # 	(phi.aug)^2
  # grad.phi <- colSums(grad.phi.tmp) + mat6 / (phi * (1 - phi))
  
  # Derivative with respect to psi
  # grad.psi: m-vector
  grad.psi.tmp <- matrix(0, n, m)
  grad.psi.tmp <- -mat5 * 
    (digamma(s2) - digamma(Y.c.aug-Y+s2) - 
       (digamma(Y+s1)-digamma(Y.c.aug-Y+s2)-digamma(s1)+digamma(s2))*p.current) * 
    tmp.val
  grad.psi <- colSums(grad.psi.tmp) + mat6
  
  # Derivative with respect to eta
  # grad.eta0, grad.eta1: m-vector
  grad.eta.tmp <- -(mat4 / mat3) * patzero.current * (1 - patzero.current)
  z0.aug <- matrix(rep(ziMatrix[,1], m), n, m)
  z1.aug <- matrix(rep(ziMatrix[,2], m), n, m)
  grad.eta0 <- colSums(grad.eta.tmp * z0.aug)
  grad.eta1 <- colSums(grad.eta.tmp * z1.aug)
  
  # # Derivative with respect to gamma
  # # grad.gamma: (K+1)-vector
  # tmp1 <- matrix(0, K+1, m)
  # for (i in 1:(K+1)) {
  # 	tmp1[i,] <- meanXB^(i-1)
  # }
  # grad.gamma <- rowSums(-tmp1 * matrix(rep(mat6, K+1), K+1, m, byrow = TRUE))
  
  # output
  # stack.grad.mat <- rbind(matrix(grad.beta0, 1, m), matrix(grad.beta1, 1, m), 
  # 												matrix(grad.phi, 1, m), matrix(grad.eta0, 1, m), 
  # 												matrix(grad.eta1, 1, m))
  stack.grad.mat <- rbind(matrix(grad.beta0, 1, m), matrix(grad.beta1, 1, m), 
                          matrix(grad.psi, 1, m), matrix(grad.eta0, 1, m), 
                          matrix(grad.eta1, 1, m))
  # grad.out <- c(as.vector(stack.grad.mat), grad.gamma)
  grad.out <- as.vector(stack.grad.mat)
  return(grad.out)
}